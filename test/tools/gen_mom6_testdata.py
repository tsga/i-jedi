#!/usr/bin/env python3
"""
Generate MOM6 geometry test data from CDL files.

Recreates ocean_hgrid.nc and ocean_topog.nc from CDL dumps of the real
72x35x25 GFDL low-resolution tripolar grid.

CDL files are committed to the repo (no git-lfs required). ncgen recreates
the netCDF files at test time.  Grid parameters (NIGLOBAL, NJGLOBAL, etc.)
are specified directly in the geometry YAML config — no MOM_input file is
generated or needed.

Usage: python3 gen_mom6_testdata.py <input_dir> <output_dir>
"""

import os
import sys
import subprocess

import numpy as np
import netCDF4


def ncgen(cdl_path, nc_path):
    result = subprocess.run(
        ["ncgen", "-o", nc_path, cdl_path],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ncgen failed for {cdl_path}:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)


def ice_concentration(lat2d):
    """Smooth sea-ice area fraction derived from 2-D latitude.

    Arctic:    linear ramp from 0 at 70 N to 1 at 80 N.
    Antarctic: linear ramp from 0 at 60 S to 1 at 70 S.
    """
    arctic    = np.clip((lat2d - 70.0) / 10.0, 0.0, 1.0)
    antarctic = np.clip((-lat2d - 60.0) / 10.0, 0.0, 1.0)
    return np.maximum(arctic, antarctic)


VELOCITY_SCALE = 0.1  # m/s


def velocity_east_north(lon, lat):
    """Smooth analytic surface flow: solid-body rotations about the polar axis
    and about the (lon=0, lat=0) axis.  Well defined at the poles and has both
    components, so a sign or ordering error in the rotation shows up."""
    lam = np.deg2rad(lon)
    phi = np.deg2rad(lat)
    uE = VELOCITY_SCALE * (np.cos(phi) - np.sin(phi) * np.cos(lam))
    vN = VELOCITY_SCALE * np.sin(lam)
    return uE, vN


def _xyz(lon, lat):
    lam, phi = np.deg2rad(lon), np.deg2rad(lat)
    return np.stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)], -1)


def _velocity_xyz(lon, lat):
    """velocity_east_north as a 3-D vector."""
    lam, phi = np.deg2rad(lon), np.deg2rad(lat)
    east = np.stack([-np.sin(lam), np.cos(lam), np.zeros_like(lam)], -1)
    north = np.stack([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)], -1)
    uE, vN = velocity_east_north(lon, lat)
    return uE[..., None] * east + vN[..., None] * north


def _tangent(p, a, b):
    """Unit vector along chord a->b, projected on the tangent plane at p."""
    t = b - a
    t = t - np.sum(t * p, -1, keepdims=True) * p
    return t / np.linalg.norm(t, axis=-1, keepdims=True)


def cgrid_velocity(hgrid_nc, ni, nj):
    """Grid-relative analytic velocity on non-symmetric MOM6 C-grid faces.

    u on the east face of cell i (supergrid (2j+1, 2i+2)) is the component
    along grid east, v on the north face of cell j (supergrid (2j+2, 2i+1))
    the component along grid north.  Grid directions come from the faces' own
    supergrid edges in 3-D, so they stay exact near the geographic pole.
    Returns u, v of shape (nj, ni).
    """
    with netCDF4.Dataset(hgrid_nc) as ds:
        x = ds["x"][:].data
        y = ds["y"][:].data
    xyz = _xyz(x, y)
    jT = 2 * np.arange(nj) + 1
    iT = 2 * np.arange(ni) + 1

    # u face: grid north along the face edge, grid east = north x up
    J, I = np.meshgrid(jT, iT + 1, indexing="ij")
    p = xyz[J, I]
    gnorth = _tangent(p, xyz[J - 1, I], xyz[J + 1, I])
    geast = np.cross(gnorth, p)
    u = np.sum(_velocity_xyz(x[J, I], y[J, I]) * geast, -1)

    # v face: grid east along the face edge, grid north = up x east
    J, I = np.meshgrid(jT + 1, iT, indexing="ij")
    p = xyz[J, I]
    geast = _tangent(p, xyz[J, I - 1], xyz[J, I + 1])
    gnorth = np.cross(p, geast)
    v = np.sum(_velocity_xyz(x[J, I], y[J, I]) * gnorth, -1)

    return u, v

def gen_mom6_restart(outdir):
    """Generate MOM.res.nc with analytically defined fields.

    Mimics a MOM6 z-coordinate restart for the 72x35x25 test grid.
    Only the variables needed for JEDI state I/O are written:
    h, Temp, Salt, u, v, ave_ssh.

    Coordinates and mask are read from soca_gridspec.72x35x25.nc (tripolar).
    """
    # ---- grid dimensions -----------------------------------------------------
    nz = 25
    nj = 35   # lath / latq size
    ni = 72   # lonh / lonq size

    # ---- read tripolar grid from soca_gridspec -------------------------------
    gridspec_nc = os.path.join(outdir, "soca_gridspec.72x35x25.nc")
    with netCDF4.Dataset(gridspec_nc) as gs:
        lonh  = gs["lonh"][0].data          # (ni,)  nominal T-cell lon axis
        lath  = gs["lath"][0].data          # (nj,)  nominal T-cell lat axis
        lonq  = gs["lonq"][0].data          # (ni,)  nominal B-grid lon axis
        latq  = gs["latq"][0].data          # (nj,)  nominal B-grid lat axis
        lat2d = gs["lat"][0].data           # (nj, ni) actual 2-D latitude
        mask  = gs["mask2d"][0].data        # (nj, ni) 1=ocean 0=land

    # Fill value for land / massless cells.  Deliberately a ridiculous number
    # (not 0.0) so that any code that fails to mask land blows up loudly
    # instead of silently consuming a plausible-looking value.
    LAND_FILL = -1.0e38

    # ---- z-level midpoint depths (72x35x25 MOM6 configuration) -------------
    layer_z = np.array([
        2.5, 7.5, 12.505, 17.545, 22.705, 28.170, 34.285, 41.610, 50.990,
        63.630, 81.165, 105.755, 140.170, 187.880, 253.155, 341.170, 458.115,
        611.290, 809.230, 1061.82, 1380.41, 1777.955, 2269.125, 2870.44,
        4429.32233,
    ])

    # Derive interface depths from layer midpoints: zi[0]=0, zi[k+1]=2*z[k]-zi[k]
    zi = np.zeros(nz + 1)
    for k in range(nz):
        zi[k + 1] = 2.0 * layer_z[k] - zi[k]
    dz = np.diff(zi)   # target layer thicknesses, shape (nz,)

    # ---- depth from ocean_topog.nc (consistent with geometry setup) ----------
    topog_nc = os.path.join(outdir, "ocean_topog.nc")
    with netCDF4.Dataset(topog_nc) as ds_topo:
        depth = ds_topo["depth"][:].data   # (nj, ni), land = 0

    # ---- h: z-coordinate layer thicknesses (vectorised) ---------------------
    # Broadcast shapes: dz (nz,1,1), depth (1,nj,ni), zi (nz+1,) → slices
    zi_lo = zi[:-1, np.newaxis, np.newaxis]   # (nz, 1, 1) upper interface
    zi_hi = zi[1:,  np.newaxis, np.newaxis]   # (nz, 1, 1) lower interface
    dz_b  = dz[:,   np.newaxis, np.newaxis]   # (nz, 1, 1)
    dep   = depth[np.newaxis, :, :]           # ( 1, nj, ni)

    h = np.where(
        dep <= 0.0,             # land: massless
        1e-10,
        np.where(
            zi_lo >= dep,       # layer below ocean bottom: massless
            1e-10,
            np.where(
                zi_hi <= dep,   # layer fully inside ocean: full target thickness
                dz_b,
                np.maximum(dep - zi_lo, 1e-10),  # partial bottom layer
            ),
        ),
    )   # shape: (nz, nj, ni)

    # ---- Temp: warm tropics, cold poles; exponential decay with depth --------
    # Use actual 2-D latitude (tripolar grid distorts near the poles).
    # Force SST to the freezing point (-1.8 C) wherever sea ice is present.
    t_surf = np.clip(28.0 * np.cos(np.deg2rad(lat2d)) - 1.8, -1.8, 30.0)
    t_surf = np.where(ice_concentration(lat2d) > 0.15, -1.8, t_surf)
    Temp = np.where(
        h > 1e-5,
        (t_surf[np.newaxis, :, :]
         * np.exp(-layer_z[:, np.newaxis, np.newaxis] / 500.0)
         + 2.0),
        LAND_FILL,
    )   # shape: (nz, nj, ni)

    # ---- Salt: subtropical maximum at surface, nearly uniform with depth -----
    s_surf = 35.0 - 2.0 * np.sin(np.deg2rad(lat2d)) ** 2
    Salt = np.where(
        h > 1e-5,
        (s_surf[np.newaxis, :, :]
         * (1.0 - 0.02 * np.exp(-layer_z[:, np.newaxis, np.newaxis] / 2000.0))),
        LAND_FILL,
    )   # shape: (nz, nj, ni)

    # ---- u, v: analytic geographic flow on the MOM6 C-grid ------------------
    # (uE, vN) is evaluated at each face position and projected on that face's
    # grid direction, so i-jedi's face->h averaging and grid->geographic
    # rotation can be checked against velocity_east_north().
    # Masked faces are 0, as in a real MOM6 restart.
    tmask = h > 1e-5                                          # (nz, nj, ni)
    u_surf, v_surf = cgrid_velocity(os.path.join(outdir, "ocean_hgrid.nc"), ni, nj)
    decay = np.exp(-layer_z / 200.0)[:, np.newaxis, np.newaxis]

    u_mask = tmask & np.roll(tmask, -1, axis=2)               # periodic in i
    u = np.where(u_mask, decay * u_surf, 0.0)                 # (Layer, lath, lonq)

    v_mask = np.zeros_like(tmask)
    v_mask[:, :-1, :] = tmask[:, :-1, :] & tmask[:, 1:, :]
    v_mask[:, -1, :] = tmask[:, -1, :] & tmask[:, -1, ::-1]   # tripolar fold
    v = np.where(v_mask, decay * v_surf, 0.0)                 # (Layer, latq, lonh)

    # ---- ave_ssh: small sinusoidal SSH with zero on land --------------------
    ave_ssh = np.where(
        mask > 0.5,
        0.1 * np.sin(np.deg2rad(lat2d)),
        LAND_FILL,
    )   # shape: (nj, ni)

    # ---- write NetCDF --------------------------------------------------------
    out_nc = os.path.join(outdir, "MOM.res.nc")
    with netCDF4.Dataset(out_nc, "w", format="NETCDF4") as ds:
        ds.createDimension("Time",      1)
        ds.createDimension("Layer",     nz)
        ds.createDimension("Interface", nz + 1)
        ds.createDimension("lath",      nj)
        ds.createDimension("lonh",      ni)
        ds.createDimension("latq",      nj)
        ds.createDimension("lonq",      ni)

        def defvar(name, dims, units, data):
            v = ds.createVariable(name, "f8", dims)
            v.units = units
            v[:] = data

        defvar("Time",      ("Time",),                        "days",          [0.0])
        defvar("Layer",     ("Layer",),                       "meter",         layer_z)
        defvar("Interface", ("Interface",),                   "meter",         zi)
        defvar("lonh",      ("lonh",),                        "degrees_east",  lonh)
        defvar("lath",      ("lath",),                        "degrees_north", lath)
        defvar("lonq",      ("lonq",),                        "degrees_east",  lonq)
        defvar("latq",      ("latq",),                        "degrees_north", latq)
        defvar("h",       ("Time", "Layer", "lath", "lonh"), "m",     h[np.newaxis])
        defvar("Temp",    ("Time", "Layer", "lath", "lonh"), "degC",  Temp[np.newaxis])
        defvar("Salt",    ("Time", "Layer", "lath", "lonh"), "PPT",   Salt[np.newaxis])
        defvar("u",       ("Time", "Layer", "lath", "lonq"), "m s-1", u[np.newaxis])
        defvar("v",       ("Time", "Layer", "latq", "lonh"), "m s-1", v[np.newaxis])
        defvar("ave_ssh", ("Time", "lath", "lonh"),          "meter", ave_ssh[np.newaxis])

    print(f"  {out_nc}")


def gen_cice6_history(outdir):
    """Generate cice6.nc with CICE6-format sea-ice fields on the 72x35x25 grid.

    Produces realistic Arctic (> 70 N) and Antarctic (< 60 S) sea-ice
    distributions: aice_h, hi_h, hs_h, sice_h.  Arctic ice is thicker
    (2 m in-ice) than Antarctic (1 m).
    """
    nj, ni = 35, 72
    FILL_LARGE = np.float32(1.e30)

    gridspec_nc = os.path.join(outdir, "soca_gridspec.72x35x25.nc")
    with netCDF4.Dataset(gridspec_nc) as gs:
        lonh  = gs["lonh"][0].data          # (ni,)
        lat2d = gs["lat"][0].data           # (nj, ni)
        mask  = gs["mask2d"][0].data        # (nj, ni) 1=ocean 0=land

    lon2d = np.tile(lonh[np.newaxis, :], (nj, 1))  # (nj, ni)
    aice  = (ice_concentration(lat2d) * mask).astype(np.float32)

    in_ice_hi = np.where(lat2d > 0.0, 2.0, 1.0)   # Arctic: 2 m, Antarctic: 1 m
    in_ice_hs = np.where(lat2d > 0.0, 0.3, 0.1)   # Arctic: 30 cm, Antarctic: 10 cm

    # Grid-cell means (= in-ice value * aice); 0 where no ice
    hi_h   = np.where(aice > 0.0, aice * in_ice_hi, 0.0).astype(np.float32)
    hs_h   = np.where(aice > 0.0, aice * in_ice_hs, 0.0).astype(np.float32)
    sice_h = np.where(aice > 0.0, np.float32(5.0), FILL_LARGE)

    out_nc = os.path.join(outdir, "cice6.nc")
    with netCDF4.Dataset(out_nc, "w", format="NETCDF4") as ds:
        ds.title       = "sea ice model output for CICE (synthetic test data)"
        ds.source      = "gen_mom6_testdata.py"
        ds.conventions = "CF-1.0"

        ds.createDimension("time", 1)
        ds.createDimension("nj",   nj)
        ds.createDimension("ni",   ni)

        vt           = ds.createVariable("time", "f4", ("time",))
        vt.units     = "days since 2021-07-01"
        vt.long_name = "model time"
        vt[:]        = np.float32(0.0)

        def coord_var(name, long_name, units, data):
            v               = ds.createVariable(name, "f4", ("nj", "ni"),
                                                fill_value=FILL_LARGE)
            v.long_name     = long_name
            v.units         = units
            v.missing_value = FILL_LARGE
            v[:]            = data.astype(np.float32)

        coord_var("TLAT", "T grid center latitude",  "degrees_north", lat2d)
        coord_var("TLON", "T grid center longitude", "degrees_east",  lon2d)

        def data_var(name, long_name, units, data, fill, coords, cell_measures=None):
            v             = ds.createVariable(name, "f4", ("time", "nj", "ni"),
                                              fill_value=fill)
            v.long_name   = long_name
            v.coordinates = coords
            v.time_rep    = "instantaneous"
            if units is not None:
                v.units = units
            if not np.isnan(float(fill)):
                v.missing_value = fill
            if cell_measures is not None:
                v.cell_measures = cell_measures
            v[0] = data

        data_var("aice_h", "ice area  (aggregate)",        "1",
                 aice,   FILL_LARGE, "TLON TLAT time", "area: tarea")
        data_var("hi_h",   "grid cell mean ice thickness", None,
                 hi_h,   np.float32(0.0), "TLAT TLON")
        data_var("hs_h",   "grid cell mean snow thickness", None,
                 hs_h,   np.float32(0.0), "TLAT TLON")
        data_var("sice_h", "bulk ice salinity",            "ppt",
                 sice_h, FILL_LARGE, "TLON TLAT time", "area: tarea")

    print(f"  {out_nc}")


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_dir> <output_dir>")
        sys.exit(1)

    indir = sys.argv[1]
    outdir = sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    # -------------------------------------------------------------------------
    # ocean_hgrid.nc and ocean_topog.nc from CDL
    # -------------------------------------------------------------------------
    for name in ("ocean_hgrid", "ocean_topog"):
        cdl = os.path.join(indir, f"{name}.cdl")
        nc  = os.path.join(outdir, f"{name}.nc")
        ncgen(cdl, nc)
        print(f"  {nc}")

    # -------------------------------------------------------------------------
    # SOCA reference gridspec (used by comparison test)
    # -------------------------------------------------------------------------
    soca_cdl = os.path.join(indir, "soca_gridspec.72x35x25.cdl")
    soca_nc  = os.path.join(outdir, "soca_gridspec.72x35x25.nc")
    ncgen(soca_cdl, soca_nc)
    print(f"  {soca_nc}")

    # -------------------------------------------------------------------------
    # MOM.res.nc — analytical restart (h, Temp, Salt, u, v, ave_ssh)
    # -------------------------------------------------------------------------
    gen_mom6_restart(outdir)

    # -------------------------------------------------------------------------
    # cice6.nc — CICE6-format sea-ice history (aice_h, hi_h, hs_h, sice_h)
    # -------------------------------------------------------------------------
    gen_cice6_history(outdir)

    print(f"MOM6 test data written to: {outdir}")


if __name__ == "__main__":
    main()
