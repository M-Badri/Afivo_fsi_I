#include "cpp_macros.h"
module m_boundary_condition
    use m_af_all
    use m_solver_parameters

    implicit none
    private
    public :: taylor_boundary_condition
    public :: channel_boundary_condition
    public :: cavity_boundary_condition
    public :: box_boundary_condition
    public :: generic_bndr_ker


    abstract interface
        subroutine generic_bndr_ker(box, nb, iv, coords, bc_val, bc_type)
            import
            type(box_t),      intent(in)  :: box  !< Box that needs b.c.
            integer(kind=I4), intent(in)  :: nb   !< Direction
            integer(kind=I4), intent(in)  :: iv   !< Index of variable
            real(kind=R8),    intent(in)  :: coords(NDIM, box%n_cell**(NDIM-1))  !< Coordinates of boundary
            real(kind=R8),    intent(out) :: bc_val(box%n_cell**(NDIM-1))  !< Boundary values
            integer(kind=I4), intent(out) :: bc_type  !< Type of b.c.
        end subroutine generic_bndr_ker
    end interface


contains


    !> Boundary condition for the Taylor vortex problem
    subroutine taylor_boundary_condition(box, nb, iv, coords, bc_val, bc_type)
        type(box_t),      intent(in)  :: box  !< Box that needs b.c.
        integer(kind=I4), intent(in)  :: nb   !< Direction
        integer(kind=I4), intent(in)  :: iv   !< Index of variable
        real(kind=R8),    intent(in)  :: coords(NDIM, box%n_cell**(NDIM-1))  !< Coordinates of boundary
        real(kind=R8),    intent(out) :: bc_val(box%n_cell**(NDIM-1))  !< Boundary values
        integer(kind=I4), intent(out) :: bc_type  !< Type of b.c.
        integer(kind=I4)              :: n

        ! Below the solution is specified in the appropriate ghost cells
        select case( nb )
            case(af_neighb_lowx)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                end if

            case(af_neighb_highx)
                if (iv == i_uold .or. iv == i_unew) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                end if

            case(af_neighb_lowy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val = 0.0_R8
                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_neumann
                    bc_val = 0.0_R8
                end if

            case(af_neighb_highy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val = 0.0_R8
                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_neumann
                    bc_val = 0.0_R8
                end if
        end select

        if ( iv == i_p ) then
            bc_type = af_bc_dirichlet
            do n = 1, box%n_cell**(NDIM-1)
                bc_val(n) =  -0.25_R8 * (exp(-2.0_R8 * (field_time-dt) / rey)) *   &
                              (cos(2.0_R8 * (coords(1, n))) + cos(2.0d0 * coords(2, n)))
            end do
        end if

        if ( iv == i_phi ) then
            bc_type = af_bc_dirichlet
            do n = 1, box%n_cell**(NDIM-1)
                bc_val(n) =  -0.25_R8 * (exp(-2.0_R8 * (field_time-dt) / rey)) *   &
                    (cos(2.0_R8 * (coords(1, n))) + cos(2.0d0 * coords(2, n)))
            end do
        end if
    end subroutine taylor_boundary_condition


    !> Boundary condition for the Channel problem
    subroutine channel_boundary_condition( box, nb, iv, coords, bc_val, bc_type )
        type(box_t),      intent(in)  :: box  !< Box that needs b.c.
        integer(kind=I4), intent(in)  :: nb   !< Direction
        integer(kind=I4), intent(in)  :: iv   !< Index of variable
        real(R8),         intent(in)  :: coords(NDIM, box%n_cell**(NDIM-1))  !< Coordinates of boundary
        real(R8),         intent(out) :: bc_val(box%n_cell**(NDIM-1))  !< Boundary values
        integer(kind=I4), intent(out) :: bc_type  !< Type of b.c.

        real(kind=R8)                 :: u_inlet = 1.0_R8, v_inlet = 0.0_R8, inv_dr(NDIM)
        integer(kind=I4)              :: nc

        inv_dr(:) = 1.0_R8 / (box%dr(:))
        nc = box%n_cell

        ! Below the solution is specified in the appropriate ghost cells
        select case( nb )
            case(af_neighb_lowx)
                if ( iv == i_unew .or. iv == i_uold ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = u_inlet

                else if (iv == i_hphix) then
                    bc_type = af_bc_dirichlet
                    bc_val  = u_inlet

                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = v_inlet

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = v_inlet

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if


            case(af_neighb_highx)
                if ( iv == i_uold ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_unew ) then
                    bc_type   = af_bc_dirichlet
                    bc_val(:) = box%cc(nc, 1:nc, i_uold) +                               &
                    &  (-u_inlet * dt * inv_dr(1)) *                                     &
                    &  (box%cc(nc, 1:nc, i_uold) - box%cc(nc-1, 1:nc, i_uold))

                else if ( iv == i_hphix ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_vold ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val(:) = box%cc(nc, 1:nc, i_vold) +                               &
                    &  (-v_inlet * dt * inv_dr(1)) *                                     &
                    &  (box%cc(nc, 1:nc, i_vold) - box%cc(nc-1, 1:nc, i_vold))

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if


            case(af_neighb_lowy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphix ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if


            case(af_neighb_highy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphix ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if

        end select

        if ( iv == i_df ) then
            bc_type = af_bc_neumann
            bc_val  = 0.0_R8
        end if

        if ( iv == i_ibfx ) then
            bc_type = af_bc_dirichlet
            bc_val  = 0.0_R8
        end if

        if ( iv == i_ibfy ) then
            bc_type = af_bc_dirichlet
            bc_val  = 0.0_R8
        end if

        if ( iv == i_vortz ) then
            bc_type = af_bc_neumann
            bc_val  = 0.0_R8
        end if

!        if (iv == i_div_1) then
!            bc_type = af_bc_neumann
!            bc_val  = 0.0_R8
!        end if
!
!        if (iv == i_div_2) then
!            bc_type = af_bc_neumann
!            bc_val  = 0.0_R8
!        end if
    end subroutine channel_boundary_condition


    !> Boundary condition for the Cavity problem
    subroutine cavity_boundary_condition( box, nb, iv, coords, bc_val, bc_type )
        type(box_t),      intent(in)  :: box  !< Box that needs b.c.
        integer(kind=I4), intent(in)  :: nb   !< Direction
        integer(kind=I4), intent(in)  :: iv   !< Index of variable
        real(kind=R8),    intent(in)  :: coords(NDIM, box%n_cell**(NDIM-1))  !< Coordinates of boundary
        real(kind=R8),    intent(out) :: bc_val(box%n_cell**(NDIM-1))  !< Boundary values
        integer(kind=I4), intent(out) :: bc_type  !< Type of b.c.

        ! Below the solution is specified in the appropriate ghost cells
        select case( nb )
            case(af_neighb_lowx)
                if( iv == i_unew .or. iv == i_uold ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                else if( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if

            case(af_neighb_highx)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if

            case(af_neighb_lowy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if ( iv == i_p ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if ( iv == i_phi ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                end if

            case(af_neighb_highy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = -1.0_R8
                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8
                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if
        end select

            if ( iv == i_df ) then
                bc_type = af_bc_neumann
                bc_val  = 0.0_R8
            end if
    end subroutine cavity_boundary_condition




     !> Boundary condition for the Box problem
    subroutine box_boundary_condition( box, nb, iv, coords, bc_val, bc_type )
        type(box_t),      intent(in)  :: box  !< Box that needs b.c.
        integer(kind=I4), intent(in)  :: nb   !< Direction
        integer(kind=I4), intent(in)  :: iv   !< Index of variable
        real(R8),         intent(in)  :: coords(NDIM, box%n_cell**(NDIM-1))  !< Coordinates of boundary
        real(R8),         intent(out) :: bc_val(box%n_cell**(NDIM-1))  !< Boundary values
        integer(kind=I4), intent(out) :: bc_type  !< Type of b.c.

        ! Below the solution is specified in the appropriate ghost cells
        select case( nb )
            case(af_neighb_lowx)
                if ( iv == i_unew .or. iv == i_uold ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if( iv == i_hphix ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if


            case(af_neighb_highx)
                if ( iv == i_uold ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphix ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_vold ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if


            case(af_neighb_lowy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphix ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if


            case(af_neighb_highy)
                if ( iv == i_uold .or. iv == i_unew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphix ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_vold .or. iv == i_vnew ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_hphiy ) then
                    bc_type = af_bc_dirichlet
                    bc_val  = 0.0_R8

                else if ( iv == i_p ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8

                else if ( iv == i_phi ) then
                    bc_type = af_bc_neumann
                    bc_val  = 0.0_R8
                end if

        end select

        if ( iv == i_df ) then
            bc_type = af_bc_neumann
            bc_val  = 0.0_R8
        end if

        if ( iv == i_ibfx ) then
            bc_type = af_bc_dirichlet
            bc_val  = 0.0_R8
        end if

        if ( iv == i_ibfy ) then
            bc_type = af_bc_dirichlet
            bc_val  = 0.0_R8
        end if

        if ( iv == i_vortz ) then
            bc_type = af_bc_neumann
            bc_val  = 0.0_R8
        end if

!        if (iv == i_div_1) then
!            bc_type = af_bc_neumann
!            bc_val  = 0.0_R8
!        end if
!
!        if (iv == i_div_2) then
!            bc_type = af_bc_neumann
!            bc_val  = 0.0_R8
!        end if
    end subroutine box_boundary_condition


end module m_boundary_condition
