#include "cpp_macros.h"
module m_incomp_flow_field
    use m_config
    use m_solver_parameters
    use m_initial_condition
    use m_boundary_condition

    implicit none

    procedure(generic_init_ker), pointer :: p_init_ker
    procedure(generic_bndr_ker), pointer :: p_bndr_ker


contains


    subroutine build_incomp_flow_field (config, tree, refine_info )
        implicit none
        class(t_config),   intent(inout) :: config
        class(af_t)      , intent(inout) :: tree
        class(ref_info_t), intent(inout) :: refine_info
        integer(kind=I4) :: box_size
        integer(kind=I4) :: sol_type
        integer(kind=I4) :: coord_type
        integer(kind=I4) :: nbox(NDIM)
        integer(kind=I4) :: nobj

        call config%get(section_name='field', option_name='solution_type', val=sol_type)
        call config%get(section_name='field', option_name='coord_type',val=coord_type)
        call config%get(section_name='field', option_name='box_size',val=box_size)
        call config%get(section_name='field', option_name='domain_len', val=domain_len)
        call config%get(section_name='field', option_name='nbox',  val=nbox)
        call config%get(section_name='immersed_boundary', option_name='nobjects', val=nobj)

        !> "max_refinement" is a global variable which is available via
        !  "m_solver_parameters" module.
        call config%get(section_name='refinement',  &
        option_name='max_refinement', val=max_refinement)

        call config%get(section_name='problem', option_name='test_case', val=test_case)
        test_case = trim(adjustl(test_case))

        call af_add_cc_variable(tree, "p",      ix=i_p)
        call af_add_cc_variable(tree, "uo",     ix=i_uold)
        call af_add_cc_variable(tree, "un",     ix=i_unew)
        call af_add_cc_variable(tree, "vo",     ix=i_vold)
        call af_add_cc_variable(tree, "vn",     ix=i_vnew)
        call af_add_cc_variable(tree, "phi",    ix=i_phi)
        call af_add_cc_variable(tree, "prhs",   ix=i_prhs)
        call af_add_cc_variable(tree, "hrhsx",  ix=i_hrhsx)
        call af_add_cc_variable(tree, "hrhsy",  ix=i_hrhsy)
        call af_add_cc_variable(tree, "hphix",  ix=i_hphix)
        call af_add_cc_variable(tree, "hphiy",  ix=i_hphiy)
        call af_add_cc_variable(tree, "tmp",    ix=i_tmp)
        call af_add_cc_variable(tree, "tmp2", ix=i_tmp2)
        call af_add_cc_variable(tree, "tmp3", ix=i_tmp3)

!        call af_add_cc_variable(tree, "utmp", ix=i_utmp)
!        call af_add_cc_variable(tree, "vtmp", ix=i_vtmp)
        call af_add_cc_variable(tree, "vortz", ix=i_vortz)
!        call af_add_cc_variable(tree, "div_1", ix=i_div_1)
!        call af_add_cc_variable(tree, "div_2", ix=i_div_2)

        call af_add_cc_variable(tree, "uolder", ix=i_uolder)
        call af_add_cc_variable(tree, "volder", ix=i_volder)
        call af_add_cc_variable(tree, "ibfx",   ix=i_ibfx)
        call af_add_cc_variable(tree, "ibfy",   ix=i_ibfy)
        call af_add_cc_variable(tree, "resx",   ix=i_resx)
        call af_add_cc_variable(tree, "resy",   ix=i_resy)

        call af_add_fc_variable(tree, "f",      ix=i_f)
        call af_add_fc_variable(tree, "folder", ix=i_folder)

        !> Keeps different distance function for each object of interest.
        call af_add_cc_variable(tree, "df",     ix=i_df, n_copies=nobj)

#if NDIM == 3
        call af_add_cc_variable(tree, "wo",     ix=i_wold)
        call af_add_cc_variable(tree, "wn",     ix=i_wnew)
        call af_add_fc_variable(tree, "f",      ix=i_f)
#endif
        call af_init(tree, & ! Tree to initialize
        box_size, & ! A box contains box_size**DIM cells
        domain_len, &
        nbox * box_size,    &
        periodic=[DTIMES(.false.)], &
        coord=coord_type)

        select case (test_case)
            case("Box")
                p_init_ker => box_ker
                call set_initial (tree)
                p_bndr_ker => box_boundary_condition

            case("Channel")
                p_init_ker => channel_ker
                call set_initial (tree)
                p_bndr_ker => channel_boundary_condition

            case("Cavity")
                p_init_ker => cavity_ker
                call set_initial (tree)
                p_bndr_ker => cavity_boundary_condition

!            case("Taylor")
!                fld%p_init_ker => taylor_ker
!                call fld%set_initial ()
!                fld%p_bndr_ker => taylor_boundary_condition
!
            case default
                write(*,*) "Error: The <<", " ", test_case,">> ",  &
                "problem is not acceptable."
                stop
        end select

       !> Sets proper boundary condition regardless of the problem is being solved.
       !> Also sets the required prolongation method.
        call set_boundary (tree)

    end subroutine  build_incomp_flow_field



    subroutine set_initial(tree)
        class(af_t), intent(inout) :: tree

        call af_loop_box(tree, p_init_ker)
    end subroutine set_initial



    subroutine set_boundary(tree)
        use m_af_prolong, only: af_prolong_linear, af_prolong_quadratic
        implicit none
        class(af_t), intent(inout) :: tree

        !> According to problem determined in the configuration, set boundary condition
        call af_set_cc_methods(tree, i_uold, p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_unew, p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_vold, p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_vnew, p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_p,    p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_df,   p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_ibfx, p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_ibfy, p_bndr_ker, prolong=af_prolong_quadratic)
        call af_set_cc_methods(tree, i_vortz, p_bndr_ker, prolong=af_prolong_quadratic)
!        call af_set_cc_methods(tree, i_div_1, p_bndr_ker)
!        call af_set_cc_methods(tree, i_div_2, p_bndr_ker)
    end subroutine set_boundary



end module m_incomp_flow_field
