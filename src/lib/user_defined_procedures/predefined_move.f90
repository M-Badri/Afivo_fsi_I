#include "cpp_macros.h"
module m_predefined_move
    use m_af_all
    use m_solver_parameters

    procedure(generic_predefined_linear_move),     pointer :: p_predefined_linear_move
    procedure(generic_predefined_rotation_move),   pointer :: p_predefined_rotation_move
    procedure(generic_predefined_linear_velocity), pointer :: p_predefined_linear_velocity
    procedure(generic_predefined_rotation_velocity), pointer ::                          &
    &  p_predefined_rotation_velocity

    procedure(generic_predefined_linear_acceleration), pointer ::                        &
    &  p_predefined_linear_acceleration
     procedure(generic_predefined_rotation_acceleration), pointer ::                     &
    &  p_predefined_rotation_acceleration

    real (kind=R8) :: a, f, kc, dia, umax

    abstract interface
        function generic_predefined_linear_move (xc) result (xc_new)
            import
            real(kind=R8), intent(in) :: xc(NDIM)
            real(kind=R8)             :: xc_new(NDIM)
        end function generic_predefined_linear_move


        function generic_predefined_rotation_move (om) result (om_new)
            import
            real(kind=R8), intent(in) :: om(3)
            real(kind=R8)             :: om_new(3)
        end function generic_predefined_rotation_move


        function generic_predefined_rotation_velocity (om) result (om_d)
            import
            real(kind=R8), intent(in) :: om(3)
            real(kind=R8)             :: om_d(3)
        end function generic_predefined_rotation_velocity


        function generic_predefined_linear_velocity (xc) result (xc_d)
            import
            real(kind=R8), intent(in) :: xc(NDIM)
            real(kind=R8)             :: xc_d(NDIM)
        end function generic_predefined_linear_velocity


        function generic_predefined_rotation_acceleration (om) result (om_dd)
            import
            real(kind=R8), intent(in) :: om(3)
            real(kind=R8)             :: om_dd(3)
        end function generic_predefined_rotation_acceleration


        function generic_predefined_linear_acceleration (xc) result (xc_dd)
            import
            real(kind=R8), intent(in) :: xc(NDIM)
            real(kind=R8)             :: xc_dd(NDIM)
        end function generic_predefined_linear_acceleration
    end interface


contains

   subroutine build_predefined_move ()
        implicit none
        kc = 5.0_R8
        a  = 1.0_R8
        f  = 1.0_R8 / kc
        dia  = 1.0_R8
        umax = 1.0_R8
        p_predefined_linear_move           => linear_osci
        p_predefined_rotation_move         => no_rotation
        p_predefined_linear_velocity       => linear_osci_velocity
        p_predefined_rotation_velocity     => no_rotation_velocity
        p_predefined_rotation_acceleration => no_rotation_acceleration
        p_predefined_linear_acceleration   => linear_osci_acceleration
   end subroutine build_predefined_move



    function linear_osci (xc) result (xc_new)
        implicit none
        real(kind=R8), intent(in) :: xc(NDIM)
        real(kind=R8)             :: xc_new(NDIM)

         xc_new(:) = - (a * sin (2.0_R8 * pi * f * field_time))
    end function linear_osci



    function no_rotation (om) result (om_new)
        implicit none
        real(kind=R8), intent(in) :: om(3)
        real(kind=R8)             :: om_new(3)

         om_new(:) = 0.0_R8
    end function no_rotation



    function linear_osci_velocity (xc) result (xc_d)
        implicit none
        real(kind=R8), intent(in) :: xc(NDIM)
        real(kind=R8)             :: xc_d(NDIM)

         xc_d(1) = - (2.0_R8 * pi * f * a * cos (2.0_R8 * pi * f * field_time))
         xc_d (2:) = 0.0_R8
    end function linear_osci_velocity



    function no_rotation_velocity (om) result (om_d)
        implicit none
        real(kind=R8), intent(in) :: om(3)
        real(kind=R8)             :: om_d(3)

         om_d(:) = 0.0_R8
    end function no_rotation_velocity



    function no_rotation_acceleration (om) result (om_dd)
        implicit none
        real(kind=R8), intent(in) :: om(3)
        real(kind=R8)             :: om_dd(3)

         om_dd(:) = 0.0_R8
    end function no_rotation_acceleration



    function linear_osci_acceleration (xc) result (xc_dd)
        implicit none
        real(kind=R8), intent(in) :: xc(NDIM)
        real(kind=R8)             :: xc_dd(NDIM)

         xc_dd(1) = + (4.0_R8 * (pi**2) * (f**2) * a * sin (2.0_R8 * pi * f * field_time))
         xc_dd (2:) = 0.0_R8
    end function linear_osci_acceleration


end module m_predefined_move
