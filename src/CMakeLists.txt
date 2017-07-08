add_subdirectory(legendre)
add_subdirectory(assembly)
add_subdirectory(integration)
add_subdirectory(io)
add_subdirectory(linalg)
add_subdirectory(misc)

add_lib_sources(
                fortranlib/src/lib_array.f90
                )

# Propagate lib_srcs list to parent directory
add_lib_sources()


# Add main executable and propagate to parent
add_exe_sources(
                driver1D.f90
                driver2D.f90
                )
