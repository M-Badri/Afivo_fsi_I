#include "cpp_macros.h"
module m_write

use m_solver_parameters
use m_af_all


contains


    subroutine print_out_info (tree)
        implicit none
        class(af_t), intent (in) :: tree

        call af_print_info(tree)
        write(*, "(A16, I3)") " Thread numbers: ", af_get_max_threads()
        write(*, "(A12, E10.4)") " Time step: ", dt
        write(*, "(A20, L)") " Is nondimensional: ", is_nondimesional
        write(*, "(A11, E10.4)") " Reynolds: ", rey
        write(*, "(A8, E10.4)") " rho_f: ", rho_f
        write(*, "(A7, E10.4)") " nu_f: ", nu_f
        write(*, "(A8, E12.6)") " CFL_V: ", dt*(1.d0/rey)*(((1.d0/af_min_dr(tree))**2)*2)
        write(*, *) "***************************************"
        write(*, *) " "
    end subroutine print_out_info



    subroutine write_ext_forces_on_hard (i_obj, time, fx, fy, mom)
        real(kind=dp), intent(in) :: time, fx, fy, mom
        integer, intent(in) :: i_obj
        integer :: myunit
        logical :: exist
        character(len=200) :: fname, ext

        write(ext,*) i_obj
        fname = trim(adjustl("./output/lift_drag_object_")) //                           &
        &  trim(adjustl(ext)) // trim(adjustl(".dat"))

        inquire(file=fname, exist=exist)
        if (exist) then
            open(newunit=myunit, file=fname, status="old",                               &
            &  position="append", action="write")
        else
            open(newunit=myunit, file=fname, status="new", action="write")
        end if
        write(myunit, *) time, fx, fy, mom
        close(myunit)
    end subroutine write_ext_forces_on_hard



    subroutine write_enery_of_object_on_hard (i_obj, time, et, er)
        real(kind=dp), intent(in) :: time, er, et
        integer, intent(in) :: i_obj
        integer :: myunit
        logical :: exist
        character(len=200) :: fname, ext

        write(ext,*) i_obj
        fname = trim(adjustl("./output/er_et_object_")) //                           &
        &  trim(adjustl(ext)) // trim(adjustl(".dat"))

        inquire(file=fname, exist=exist)
        if (exist) then
            open(newunit=myunit, file=fname, status="old",                               &
            &  position="append", action="write")
        else
            open(newunit=myunit, file=fname, status="new", action="write")
        end if
        write(myunit, *) time, et, er
        close(myunit)
    end subroutine write_enery_of_object_on_hard



    subroutine write_reynolds_of_object_on_hard (i_obj, time, rey)
        real(kind=dp), intent(in) :: time, rey
        integer :: myunit
        logical :: exist
        character(len=200) :: fname, ext

        write(ext,*) i_obj
        fname = trim(adjustl("./output/rey_object_")) //                           &
        &  trim(adjustl(ext)) // trim(adjustl(".dat"))

        inquire(file=fname, exist=exist)
        if (exist) then
            open(newunit=myunit, file=fname, status="old",                               &
            &  position="append", action="write")
        else
            open(newunit=myunit, file=fname, status="new", action="write")
        end if
        write(myunit, *) time, rey
        close(myunit)
    end subroutine write_reynolds_of_object_on_hard



    subroutine write_out_location_velocity_on_hard (i_obj, time, x, x_d, om , om_d)
        real(kind=dp), intent(in) :: time
        real(kind=dp), intent(in) :: x(:), x_d(:), om(:), om_d(:)
        integer, intent(in) :: i_obj
        integer :: myunit
        logical :: exist
        character(len=200) :: fname, ext

        write(ext,*) i_obj
        fname = trim(adjustl("./output/XY_")) //                      &
        &  trim(adjustl(ext)) // trim(adjustl(".dat"))

        inquire(file=fname, exist=exist)
        if (exist) then
            open(newunit=myunit, file=fname, status="old",                               &
            &  position="append", action="write")
        else
            open(newunit=myunit, file=fname, status="new", action="write")
        end if
        write(myunit, *) time, x
        close(myunit)

        fname = trim(adjustl("./output/XY_dot_")) //                      &
        &  trim(adjustl(ext)) // trim(adjustl(".dat"))

        inquire(file=fname, exist=exist)
        if (exist) then
            open(newunit=myunit, file=fname, status="old",                               &
            &  position="append", action="write")
        else
            open(newunit=myunit, file=fname, status="new", action="write")
        end if
        write(myunit, *) time, x_d
        close(myunit)

        fname = trim(adjustl("./output/OMEGA_")) //                      &
        &  trim(adjustl(ext)) // trim(adjustl(".dat"))

        inquire(file=fname, exist=exist)
        if (exist) then
            open(newunit=myunit, file=fname, status="old",                               &
            &  position="append", action="write")
        else
            open(newunit=myunit, file=fname, status="new", action="write")
        end if
        write(myunit, *) time, om
        close(myunit)


        fname = trim(adjustl("./output/OMEGA_dot_")) //                      &
        &  trim(adjustl(ext)) // trim(adjustl(".dat"))

        inquire(file=fname, exist=exist)
        if (exist) then
            open(newunit=myunit, file=fname, status="old",                               &
            &  position="append", action="write")
        else
            open(newunit=myunit, file=fname, status="new", action="write")
        end if
        write(myunit, *) time, om_d
        close(myunit)
    end subroutine write_out_location_velocity_on_hard



    subroutine time_stamp (time, time_step, n_obj, fx, fy, cfl_c, is_on_hard, step_time)
        use face
        implicit none

        real(kind=dp), intent(in) :: time, fx(:), fy(:), cfl_c, step_time
        integer, intent(in) :: n_obj, time_step
        character(len=:), allocatable :: time_message, hard_message
        character(len=12) :: stime
        character(len=20)  :: obj_str
        integer :: i
        logical :: is_on_hard

        write(stime, "(E12.4E3)") time
        stime = trim(adjustl(stime))
        time_message =                                                                   &
        &  colorize('Time:', color_fg='blue', color_bg='white' ,style='bold_on' ) //     &
        &  colorize( stime,  color_fg='blue_intense', style='underline_on')
        write (*, '(A)') time_message
        write (*, '(A11, I6)') "Time step: ", time_step
        write (*, "(A16, E8.2)") "Convective CFL: ", cfl_c
        write (*, "(A41, E9.3)") "Elapsed time for the last time step (s): ", step_time
        if (is_on_hard) then
            hard_message =                                                               &
            & colorize ('The data was saved on hard for time: ', color_fg='green') //    &
            & colorize (stime, color_fg='green')  //                                     &
            & colorize ('s', color_fg='green')

            write (*, '(A)')  hard_message
        end if

# if NDIM == 2
        do i = 1, n_obj
            write(obj_str, *) i
            obj_str = "Objet: " // trim(adjustl(obj_str) )
            write(*,"(A)") colorize( obj_str, color_fg='cyan_intense')
            write(*, "( A, E12.4E3, A, E12.4E3)") " Fx: ", fx, ", Fy: ", fy
        end do
# endif
    end subroutine time_stamp


end module m_write
