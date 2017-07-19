add_lib_sources(
                mod_legendre_c.f90
                mod_legendre.f90
                smod_legendre_1D.f90
                smod_legendre_2D.f90
                )

# Propagate lib_srcs list to parent directory
add_lib_sources()
