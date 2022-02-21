#include "cpp_macros.h"
module m_projection
    use m_config
    use m_af_all
    use m_incomp_flow_field
    use m_solver_parameters
    use m_initial_condition
    use m_boundary_condition, only: channel_boundary_condition
    use m_boundary_condition, only: cavity_boundary_condition
    use m_boundary_condition, only: taylor_boundary_condition

    implicit none
    private
    public :: buid_projector

    type(mg_t) :: mg
    real(kind = dp) :: mean
    integer :: max_itr_poisson
    procedure(generic_sides_bc),      private, pointer :: p_sides_bc
    procedure(generic_get_projector), public,  pointer :: p_get_projector

    abstract interface

        subroutine generic_get_projector(tree)
            import
            class(af_t), intent(inout) :: tree
        end subroutine generic_get_projector


        subroutine generic_sides_bc(box, nb, iv, coords, bc_val, bc_type)
            import box_t, dp
            type(box_t), intent(in) :: box
            integer,  intent(in)    :: nb
            integer,  intent(in)    :: iv
            real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
            real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
            integer,  intent(out)   :: bc_type
        end subroutine generic_sides_bc

    end interface



contains



    subroutine buid_projector (config, tree)
        class(af_t),     intent(inout) :: tree
        class(t_config), intent(inout) :: config
        character(len=:), allocatable  :: time_method, time_scheme

        test_case = trim(adjustl(test_case))
        select case (test_case)
            case ("Cavity")
                p_sides_bc => cavity_boundary_condition

            case ("Box")
                p_sides_bc => box_boundary_condition

            case ("Channel")
                p_sides_bc => channel_boundary_condition

            case ("Taylor")
                p_sides_bc => taylor_boundary_condition
        end select

        mg%i_phi    = i_phi     ! Solution variable
        mg%i_rhs    = i_prhs    ! Right-hand side variable
        mg%i_tmp    = i_tmp     ! Variable for temporary space
        mg%sides_bc => p_sides_bc
        mg%subtract_mean = .true.

        !> Initializes the multigrid solver using the given field (and its tree)
        call mg_init(tree, mg)

        allocate(character(50) :: time_method)
        allocate(character(50) :: time_scheme)
        call config%get(section_name='solver_time', option_name='method', &
        &   val=time_method)
        call config%get(section_name='solver_time', option_name='scheme', &
        &   val=time_scheme)
        time_method = trim(adjustl(time_method))
        time_scheme = trim(adjustl(time_scheme))

        call config%get(section_name='poisson', option_name='max_iteratoin', &
        &   val=max_itr_poisson)

        select case(time_method)
            case ("Adams_Bashforth_exp")
                select case (time_scheme)
                    case("smac")
                        p_get_projector => get_smac_projector

                    case("mac")
                        p_get_projector => get_mac_projector

                    case default
                        write(*,*) "Error: The <<", " ", time_scheme,">> ", &
                        "time integration scheme is not acceptable for projection."
                        stop
                    end select

            case ("Adams_Bashforth_simp")
                select case (time_scheme)
                    case("smac")
                        p_get_projector => get_smac_projector

                    case("mac")
                        p_get_projector => get_mac_projector

                    case default
                        write(*,*) "Error: The <<", " ", time_scheme,">> ", &
                        "time integration scheme is not acceptable for projection"
                        stop
                end select

                case default
                write(*,*) "Error: The <<", " ", time_method,">> ", &
                "time integration method is not acceptable for projection."
                stop
        end select
    end subroutine buid_projector



    subroutine get_mac_projector (tree)
        class(af_t), intent(inout) :: tree
        integer :: mg_iter

        do mg_iter = 1, max_itr_poisson
            call mg_fas_fmg(tree, mg, set_residual=.true., have_guess=(mg_iter>1))
        end do

        call af_tree_copy_cc (tree, i_phi, i_p)
    end subroutine get_mac_projector



    subroutine get_smac_projector (tree)
        class(af_t), intent(inout) :: tree
        integer :: mg_iter

        do mg_iter = 1, max_itr_poisson
            call mg_fas_fmg (tree, mg, set_residual=.true., have_guess=(mg_iter>1))
        end do
!        mean = 0.0_dp
!        call af_loop_box (tree, get_mean_ker, leaves_only = .true.)
!        call af_loop_box (tree, mean_free_ker)
    end subroutine get_smac_projector



    subroutine mean_free_ker (box)
        type(box_t), intent(inout) :: box
        integer  :: nc
        real(kind=dp) :: mean

        nc = box%n_cell
!        mean = sum(box%cc(1:nc, 1:nc, i_phi) ) / real(nc**2, kind=dp)
        box%cc(1:nc, 1:nc, i_phi) = box%cc(1:nc, 1:nc, i_phi) - ( mean )
        print*, mean
    end subroutine mean_free_ker



    subroutine get_mean_ker (box)
        type(box_t), intent(inout) :: box
        integer :: nc, i, j
        !real(kind=dp) :: mean
        !$omp critical
!        nc = box%n_cell
!        do j = 1, nc
!            do i = 1, nc
!                mean = mean + box%cc(i, j, i_phi) !* box%dr(1) * box%dr(2)
!            end do
!        end do
        mean = sum(box%cc(1:nc, 1:nc, i_phi)) * real(nc**2, kind=dp)
        !$omp end critical
        !box%cc(1:nc, 1:nc, i_phi) = box%cc(1:nc, 1:nc, i_phi) - mean
    end subroutine get_mean_ker

end module m_projection
