# Standalone unit tests (no YAML or test data needed)
ecbuild_add_test( TARGET  ijedi_unit_swap_io_member
                  SOURCES SwapIoMember.cc
                  LIBS    ijedi eckit )
