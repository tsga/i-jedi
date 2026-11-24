add_subdirectory( ijedi )
add_subdirectory( mains )

# find cpplint.py provided by oops
find_file( CPPLINT_PY oops_cpplint.py HINTS ${oops_BINDIR} )

ecbuild_add_test( TARGET ijedi_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CPPLINT_PY}
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
