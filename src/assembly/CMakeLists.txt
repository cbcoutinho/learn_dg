add_lib_sources(
                mod_assembly.f90
                mod_assembly_c.f90
                smod_assemble_1D.f90
                smod_assemble_2D.f90
                )

# Propagate lib_srcs list to parent directory
add_lib_sources()
