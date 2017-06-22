add_subdirectory(legendre)

add_lib_sources(
                assembly_mod.f90
                integration_mod.f90
                io_mod.f90
                legendre_mod.f90
                linalg_mod.f90
                misc_mod.f90
                fortranlib/src/lib_array.f90
                )

# Propagate lib_srcs list to parent directory
add_lib_sources()


# Add main executable and propagate to parent
add_exe_sources(
                driverA.f90
                # driverB.f90
                )
