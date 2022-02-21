#include "cpp_macros.h"

module  m_initial_condition
    use m_af_all
    use m_solver_parameters

    implicit none
    private
    public :: generic_init_ker
    public :: taylor_ker
    public :: channel_ker
    public :: cavity_ker
    public :: box_ker

    abstract interface
        subroutine generic_init_ker(box)
            import
            type(box_t), intent(inout) :: box
        end subroutine generic_init_ker
    end interface



contains



    subroutine taylor_ker(box)
        type(box_t), intent(inout) :: box
        integer  :: IJK, nc
        real(dp) :: rr(NDIM)
        nc = box%n_cell

        do KJI_DO(1,nc)
           rr = af_r_cc(box, [IJK])
           box%cc(IJK, i_uold) = taylor(rr(1), rr(2), field_time, i_uold)
           box%cc(IJK, i_vold) = taylor(rr(1), rr(2), field_time, i_vold)
        end do; CLOSE_DO
    end subroutine taylor_ker



    !> This routine calculates the analytic solution in point rr
    elemental function taylor(x, y, t, i_var) result(sol)
        real(dp), intent(in) :: x, y, t
        integer,  intent(in) :: i_var
        real(dp)             :: sol

        if     (i_var == i_uold) then
           sol = -exp(-2.0_dp * t / rey) * cos(x) * sin(y)
        else if(i_var == i_vold) then
           sol = +exp(-2.0_dp * t / rey) * sin(x) * cos(y)
        else if(i_var == i_p)    then
           sol = -0.25_dp * (exp(-4.0_dp * t / rey)) *   &
                    (cos(2.0_dp * x) + cos(2.0_dp * y))
        end if
    end function taylor



    subroutine channel_ker(box)
        type(box_t), intent(inout) :: box
        integer  :: IJK, nc
        real(dp) :: rr(NDIM)
        nc = box%n_cell

        do KJI_DO(1,nc)
           rr = af_r_cc(box, [IJK])
           box%cc(IJK, i_uold) = channel(rr(1), rr(2), field_time, i_uold)
           box%cc(IJK, i_vold) = channel(rr(1), rr(2), field_time, i_vold)
        end do; CLOSE_DO
    end subroutine channel_ker



    !> This routine calculates the analytic solution in point rr
    elemental function channel(x, y, t, i_var) result(sol)
        real(dp), intent(in) :: x, y, t
        integer,  intent(in) :: i_var
        real(dp)             :: sol

        if     (i_var == i_uold) then
           sol = 1.0_dp
        else if(i_var == i_vold) then
           sol = 0.0_dp
        else if(i_var == i_p)    then
           sol = 0.0_dp
        end if
    end function channel



    subroutine box_ker(box)
        type(box_t), intent(inout) :: box
        integer  :: IJK, nc
        real(dp) :: rr(NDIM)
        nc = box%n_cell

        do KJI_DO(1,nc)
           rr = af_r_cc(box, [IJK])
           box%cc(IJK, i_uold) = 0.0_R8
           box%cc(IJK, i_vold) = 0.0_R8
        end do; CLOSE_DO
    end subroutine box_ker



    subroutine cavity_ker(box)
        type(box_t), intent(inout) :: box
        integer  :: IJK, nc
        real(dp) :: rr(NDIM)
        nc = box%n_cell

        do KJI_DO(1,nc)
           rr = af_r_cc(box, [IJK])
           box%cc(IJK, i_uold) = channel(rr(1), rr(2), field_time, i_uold)
           box%cc(IJK, i_vold) = channel(rr(1), rr(2), field_time, i_vold)
        end do; CLOSE_DO
    end subroutine cavity_ker



    !> This routine calculates the analytic solution in point rr
    elemental function cavity(x, y, t, i_var) result(sol)
        real(dp), intent(in) :: x, y, t
        integer,  intent(in) :: i_var
        real(dp)             :: sol

        if     (i_var == i_uold) then
           sol = 0.0_dp
        else if(i_var == i_vold) then
           sol = 0.0_dp
        else if(i_var == i_p)    then
           sol = 0.0_dp
        end if
    end function cavity

end module m_initial_condition
