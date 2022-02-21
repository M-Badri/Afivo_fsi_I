#include "cpp_macros.h"
module m_helmholtz
    use m_solver_parameters
    use m_config
    use m_af_all
    use m_incomp_flow_field
    use m_initial_condition
    use m_boundary_condition, only: channel_boundary_condition
    use m_boundary_condition, only: cavity_boundary_condition
    use m_boundary_condition, only: taylor_boundary_condition

    implicit none
    private
    public :: buid_helmholtz, get_helmholtz

    type(mg_t) :: helmholtz_x
    type(mg_t) :: helmholtz_y
    integer(kind=I4) :: max_itr_helmholtz
    procedure(generic_sides_bc),      private, pointer :: pt_sides_bc
    procedure(generic_get_helmholtz), public,  pointer :: pt_get_helmholtz

    abstract interface

        subroutine generic_get_helmholtz(tree)
            import
            class(af_t), intent(inout) :: tree
        end subroutine generic_get_helmholtz


        subroutine generic_sides_bc(box, nb, iv, coords, bc_val, bc_type)
            import
            type(box_t), intent(in) :: box
            integer(kind=I4), intent(in) :: nb
            integer(kind=I4), intent(in) :: iv
            real(kind=R8), intent(in)  :: coords(NDIM, box%n_cell**(NDIM-1))
            real(kind=R8), intent(out) :: bc_val(box%n_cell**(NDIM-1))
            integer(kind=I4),  intent(out) :: bc_type
        end subroutine generic_sides_bc

    end interface



contains



    subroutine buid_helmholtz (config, tree)
        class(af_t),     intent(inout) :: tree
        class(t_config), intent(inout) :: config
        character(len=:), allocatable  :: time_method, time_scheme

        test_case = trim(adjustl(test_case))
        select case (test_case)
            case ("Cavity")
                pt_sides_bc => cavity_boundary_condition

            case ("Channel")
                pt_sides_bc => channel_boundary_condition

            case ("Box")
                pt_sides_bc => box_boundary_condition

            case ("Taylor")
                pt_sides_bc => taylor_boundary_condition
        end select
!
        helmholtz_x%i_phi    = i_hphix     ! Solution variable
        helmholtz_x%i_rhs    = i_hrhsx    ! Right-hand side variable
        helmholtz_x%i_tmp    = i_tmp2     ! Variable for temporary space
        helmholtz_x%sides_bc => pt_sides_bc
        helmholtz_x%helmholtz_lambda = (2.0_R8 * rey * (1.0_R8 / nu_f)) / dt

        helmholtz_y%i_phi    = i_hphiy     ! Solution variable
        helmholtz_y%i_rhs    = i_hrhsy    ! Right-hand side variable
        helmholtz_y%i_tmp    = i_tmp3     ! Variable for temporary space
        helmholtz_y%sides_bc => pt_sides_bc
        helmholtz_y%helmholtz_lambda = (2.0_R8 * rey * (1.0_R8 / nu_f)) / dt ! (2.0_dp * rey) / dt

        !> Initializes the Helmholtz multigrid solver using the given field (and its tree)
        call mg_init (tree, helmholtz_x)
        call mg_init (tree, helmholtz_y)

        allocate(character(50) :: time_method)
        allocate(character(50) :: time_scheme)
        call config%get (section_name='solver_time', option_name='method', &
        &   val=time_method)
        call config%get (section_name='solver_time', option_name='scheme', &
        &   val=time_scheme)
        time_method = trim (adjustl (time_method))
        time_scheme = trim (adjustl (time_scheme))

        call config%get (section_name='helmholtz', option_name='max_iteratoin', &
        &   val=max_itr_helmholtz)

        select case(time_method)

            case ("Adams_Bashforth_simp")
                select case (time_scheme)
                    case("smac")
                        pt_get_helmholtz => get_helmholtz

                    case ("mac")
                        pt_get_helmholtz => get_helmholtz

                    case default
                        write(*,*) "Error: The <<", " ", time_scheme,">> ", &
                        "time integration scheme is not acceptable for Helmholtz solver."
                        stop
                end select

            case default
                write(*,*) "Error: The <<", " ", time_method,">> ", &
                "time integration method is not acceptable for Helmholtz solver."
                stop

        end select
    end subroutine buid_helmholtz



    subroutine get_helmholtz (tree)
        implicit none
        class(af_t), intent(inout) :: tree
        integer(kind=I4) :: mg_iter

        do mg_iter = 1, max_itr_helmholtz
            call mg_fas_fmg(tree, helmholtz_x, set_residual=.true., have_guess=(mg_iter>1))
        end do



        do mg_iter = 1, max_itr_helmholtz
            call mg_fas_fmg(tree, helmholtz_y, set_residual=.true., have_guess=(mg_iter>1))
        end do

        call af_tree_copy_cc(tree, i_hphix, i_unew)
        call af_tree_copy_cc(tree, i_hphiy, i_vnew)
    end subroutine get_helmholtz

end module m_helmholtz
