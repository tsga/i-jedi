#!/usr/bin/env python3
"""Check MOM6 C-grid <-> h-point velocity conversion in i-jedi.

Usage: check_mom6_velocity.py <data_generated/mom6 dir>

Expects, in that directory:
  MOM.res.nc               analytic C-grid u, v (from gen_mom6_testdata.py)
  ocean_hgrid.nc           supergrid
  mom6_velocity_tracer.nc  i-jedi output, velocity grid: tracer
  mom6_velocity_native.nc  i-jedi output, velocity grid: native

Checks
  1. Strict: both outputs match a numpy re-implementation of the MOM6
     conventions (2-point face averages, tripolar fold, rotation angle from
     MOM6 initialize_grid_rotation_angle) to round-off at every wet point.
  2. tracer output (eastward/northward on h-points) matches the analytic
     flow at wet cells whose four faces are wet, south of POLAR_CAP: MOM6's
     finite-difference rotation angle is inaccurate within a few degrees of
     the geographic pole on this 5-degree grid.  The same faces averaged
     WITHOUT rotation must be clearly worse, so the rotation is exercised.
  3. native output (grid-relative on faces) matches MOM.res.nc u, v on faces
     whose neighbouring cells are both fully wet (h -> face is a 2-point
     smoother, so this is tolerance-based), with land faces exactly 0 and
     variables on lonq / latq.
"""

import os
import sys

import netCDF4
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_mom6_testdata import VELOCITY_SCALE, velocity_east_north  # noqa: E402

TRACER_TOL = 0.02   # max |error| / VELOCITY_SCALE at h-points
NATIVE_TOL = 0.02   # max |error| / VELOCITY_SCALE on faces
STRICT_TOL = 1e-10  # max |error| / VELOCITY_SCALE vs the numpy reference
POLAR_CAP = 75.0    # degrees north; analytic check skipped beyond
MINIMUM_DEPTH = 0.5
MINIMUM_THICKNESS = 1e-6


def mom6_rotation(x, y):
    """cos/sin of the grid rotation at h-points, MOM6
    initialize_grid_rotation_angle (non-bugs branch), from supergrid corners."""
    lonc = x[1::2, 1::2]

    def corner(js, is_):
        lon = x[js::2, is_::2][:lonc.shape[0], :lonc.shape[1]]
        lat = y[js::2, is_::2][:lonc.shape[0], :lonc.shape[1]]
        return lon - 360.0 * np.floor((lon - lonc) / 360.0 + 0.5), lat

    lonSW, latSW = corner(0, 0)
    lonSE, latSE = corner(0, 2)
    lonNW, latNW = corner(2, 0)
    lonNE, latNE = corner(2, 2)
    lon_scale = np.cos(np.pi / 720.0 * ((latSW + latNE) + (latSE + latNW)))
    angle = np.arctan2(lon_scale * ((lonNW - lonSE) + (lonNE - lonSW)),
                       (latNW - latSE) + (latNE - latSW))
    return np.cos(angle), np.sin(angle)


def geometry_mask(d, h):
    """mask3d as built by GeometryMOM6: wet column, sealed below the first
    layer thinner than MINIMUM_THICKNESS."""
    with netCDF4.Dataset(os.path.join(d, "ocean_topog.nc")) as ds:
        wet = ds["depth"][:].data > MINIMUM_DEPTH
    return wet[None] & (np.cumsum(h < MINIMUM_THICKNESS, axis=0) == 0)


def read(path, name):
    with netCDF4.Dataset(path) as ds:
        var = ds[name]
        return var[0].data if var.dimensions[0] == "Time" else var[:].data, var.dimensions


def main():
    d = sys.argv[1]
    res = os.path.join(d, "MOM.res.nc")
    u, _ = read(res, "u")
    v, _ = read(res, "v")
    h, _ = read(res, "h")
    layer, _ = read(res, "Layer")
    nz, nj, ni = h.shape

    with netCDF4.Dataset(os.path.join(d, "ocean_hgrid.nc")) as ds:
        x = ds["x"][:].data
        y = ds["y"][:].data
    lonT = x[1::2, 1::2]
    latT = y[1::2, 1::2]
    cos, sin = mom6_rotation(x, y)
    decay = np.exp(-layer / 200.0)[:, None, None]
    uE, vN = velocity_east_north(lonT, latT)
    uE, vN = decay * uE, decay * vN

    tmask = h > 1e-5
    u_wet = u != 0.0
    v_wet = v != 0.0
    v_south = np.zeros_like(v_wet)
    v_south[:, 1:, :] = v_wet[:, :-1, :]
    interior = tmask & u_wet & np.roll(u_wet, 1, axis=2) & v_wet & v_south

    ok = True

    trc = os.path.join(d, "mom6_velocity_tracer.nc")
    uT, dims_u = read(trc, "u")
    vT, _ = read(trc, "v")
    assert dims_u[-1] == "lonh", f"tracer u should be on lonh, got {dims_u}"
    nat = os.path.join(d, "mom6_velocity_native.nc")
    uF, dims_u = read(nat, "u")
    vF, dims_v = read(nat, "v")
    assert dims_u[-2:] == ("lath", "lonq"), f"native u dims {dims_u}"
    assert dims_v[-2:] == ("latq", "lonh"), f"native v dims {dims_v}"

    # --- 1. strict: numpy reference of the MOM6 conventions --------------------
    mask = geometry_mask(d, h)
    u_avg = 0.5 * (u + np.roll(u, 1, axis=2))
    v_avg = 0.5 * (v + np.concatenate([np.zeros_like(v[:, :1]), v[:, :-1]], axis=1))
    uE_ref = cos * u_avg + sin * v_avg
    vN_ref = cos * v_avg - sin * u_avg
    err = max(np.abs(uT - uE_ref)[mask].max(), np.abs(vT - vN_ref)[mask].max())
    err /= VELOCITY_SCALE
    print(f"strict tracer: max rel error {err:.2e} (tol {STRICT_TOL})")
    if err > STRICT_TOL:
        print("FAIL: tracer output differs from the MOM6-convention reference")
        ok = False

    ug = np.where(mask, cos * uT - sin * vT, 0.0)
    vg = np.where(mask, cos * vT + sin * uT, 0.0)
    u_ref = mask * np.roll(mask, -1, axis=2) * 0.5 * (ug + np.roll(ug, -1, axis=2))
    v_ref = np.zeros_like(vg)
    v_ref[:, :-1] = mask[:, :-1] * mask[:, 1:] * 0.5 * (vg[:, :-1] + vg[:, 1:])
    v_ref[:, -1] = mask[:, -1] * mask[:, -1, ::-1] * 0.5 * (vg[:, -1] - vg[:, -1, ::-1])
    err = max(np.abs(uF - u_ref).max(), np.abs(vF - v_ref).max()) / VELOCITY_SCALE
    print(f"strict native: max rel error {err:.2e} (tol {STRICT_TOL})")
    if err > STRICT_TOL:
        print("FAIL: native output differs from the MOM6-convention reference")
        ok = False

    # --- 2. tracer output vs the analytic flow ---------------------------------
    check = interior & (latT < POLAR_CAP)[None]
    err = max(np.abs(uT - uE)[check].max(), np.abs(vT - vN)[check].max())
    err /= VELOCITY_SCALE
    err_norot = max(np.abs(u_avg - uE)[check].max(),
                    np.abs(v_avg - vN)[check].max()) / VELOCITY_SCALE
    print(f"tracer: max rel error {err:.4f} (tol {TRACER_TOL}); "
          f"without rotation {err_norot:.4f}; {check.sum()} points")
    if err > TRACER_TOL:
        print("FAIL: tracer-point velocity does not match the analytic flow")
        ok = False
    if err_norot < 5.0 * TRACER_TOL:
        print("FAIL: rotation not exercised (unrotated error too small)")
        ok = False

    # --- 3. native output vs MOM.res.nc -----------------------------------------
    u_cmp = interior & np.roll(interior, -1, axis=2)
    v_cmp = np.zeros_like(interior)
    v_cmp[:, :-1, :] = interior[:, :-1, :] & interior[:, 1:, :]
    err_u = np.abs(uF - u)[u_cmp].max() / VELOCITY_SCALE
    err_v = np.abs(vF - v)[v_cmp].max() / VELOCITY_SCALE
    land_u = np.abs(uF[~u_wet]).max(initial=0.0)
    land_v = np.abs(vF[~v_wet]).max(initial=0.0)
    print(f"native: max rel error u {err_u:.4f}, v {err_v:.4f} (tol {NATIVE_TOL}); "
          f"max |land face| u {land_u:.3g}, v {land_v:.3g}")
    if max(err_u, err_v) > NATIVE_TOL:
        print("FAIL: C-grid velocity does not match MOM.res.nc")
        ok = False
    if max(land_u, land_v) != 0.0:
        print("FAIL: masked faces are not zero")
        ok = False

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
