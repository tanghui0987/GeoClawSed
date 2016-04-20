subroutine setprob()

    use geoclaw_module !geoclaw_module.f90
    use topo_module !topo_module.f90
    use qinit_module !qinit_module.f90
    use fixedgrids_module !fixedgrids_module.f90
    use refinement_module !refinement_module.f90
    use sediment_module !sediment_module.f90

    implicit none

    call set_geo()                    !# sets basic parameters g and coord system
    call set_refinement()             !# sets refinement control parameters
    call read_dtopo_settings()        !# specifies file with dtopo from earthquake
    call read_topo_settings()         !# specifies topography (bathymetry) files
    call set_qinit()                  !# specifies file with dh if this used instead
    call set_fixed_grids()            !# Fixed grid settings
    call set_sediment()               !#sets basic parameters for sediment transport
    call read_sed_settings()          !#specifies sediment (thickness and grainsize dsitribution) files

end subroutine setprob
