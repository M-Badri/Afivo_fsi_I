#include "cpp_macros.h"
program main
    use m_config
    use m_incomp_flow_field
    use m_af_all
    use m_navier_stokes_solver

    implicit none
    !@todo t_config should be changed to config_t to be consistent with the others
    type (t_config) :: config
    type(af_t)       :: tree
    type(ref_info_t) :: refine_info
    character(len=:), allocatable :: input_file

    call OMP_SET_NUM_THREADS (1)

    input_file = trim(adjustl("./inputs/my_input_file.inp"))
    call config%build (filename = input_file)
    call build_incomp_flow_field (config, tree, refine_info)
    call build_navier_stokes_solver (config, tree, refine_info)

    call p_solve (tree, refine_info)
end program main
