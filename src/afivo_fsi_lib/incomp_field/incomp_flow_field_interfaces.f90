#include "cpp_macros.h"
module incomp_flow_fld_intefaces

    use m_config
    use m_af_all

!    interface

!        subroutine generic_c2f_ker(box)
!            import box_t
!            type(box_t), intent(inout) :: box
!!            integer, intent(in) :: i_var(NDIM)
!        end subroutine generic_c2f_ker


!        subroutine generic_bndr_ker(box, nb, iv, coords, bc_val, bc_type)
!            import box_t
!            import dp
!            type(box_t), intent(in) :: box  !< Box that needs b.c.
!            integer,  intent(in)    :: nb   !< Direction
!            integer,  intent(in)    :: iv   !< Index of variable
!            real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))  !< Coordinates of boundary
!            real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))  !< Boundary values
!            integer,  intent(out)   :: bc_type  !< Type of b.c.
!        end subroutine generic_bndr_ker

!     end interface

end module incomp_flow_fld_intefaces
