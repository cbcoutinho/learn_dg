add_subdirectory(legendre)

add_lib_sources(
                mod_assembly.f90
                mod_integration.f90
                mod_io.f90
                mod_linalg.f90
                mod_misc.f90
                fortranlib/src/lib_array.f90
                )

# Propagate lib_srcs list to parent directory
add_lib_sources()


# Add main executable and propagate to parent
add_exe_sources(
                driverA.f90
                driverB.f90
                )
