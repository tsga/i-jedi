# Standalone unit tests (no YAML or test data needed)
ecbuild_add_test( TARGET  ijedi_unit_swap_io_member
                  SOURCES SwapIoMember.cc
                  LIBS    ijedi eckit )

# Generic in-place NetCDF update (Io/utils/UpdateNetcdf) on an MPAS-like restart
ecbuild_add_test( TARGET  ijedi_unit_update_netcdf
                  SOURCES UpdateNetcdf.cc
                  MPI     4
                  LIBS    ijedi eckit atlas NetCDF::NetCDF_C )
