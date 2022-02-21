
module m_mls_interpolation
   use iso_fortran_env, only: I4=>int32, R8=>real64,  character_kinds
   use m_math_lib

   implicit none

   !   "ns" is the number of the basis functions.
   integer(kind=I4), parameter :: ns = 6
   real(kind=R8), parameter :: eps = 10e-10


contains

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !  Using the points position (x,y) and their values, approximates the the
   !  function at the required position.
   !  It assumes that the first given point is a point on a boundary with
   !  DIRICHLET condition.
   !
   !  Inputs:
   !  1- np_: number of interpolant points
   !  2- nds_pos_: a dimension(2, np_) array holding the position of the
   !     source points.
   !  3- nds_val_: the value of the source points
   !  4- x_c_, y_c_: the point where the function should be interpolated
   !
   !  Output:
   !  1- s_val: an array with the length of the "ns" ("ns" is the module
   !     variable) containing the calculated shape function. These are the coefficients of
   !     the following statement (a_n):
   !
   !     f(x,y) = a_0
   !            + a_1*(x-xc)*(du/dx)             + a_2*(y-yc)*(du/dy)
   !            + a_3*(x-xc)*(x-xc)*(d^2 u/dx^2) + a_4*(x-xc)*(y-yc)*(d^2 u/(d_x*d_y))
   !            + a_5*(y-yc)*(y-yc)*(d^2 u/dy^2)
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function mls_get_scalar_dir (np_, nds_pos_, nds_val_, c_x_, c_y_) result (s_val)
      implicit none
      !   Arguments
      integer(kind=I4)                , intent(in) :: np_
      real(kind=R8), dimension(2, np_), intent(in) :: nds_pos_
      real(kind=R8), dimension(np_)   , intent(in) :: nds_val_
      real(kind=R8)                   , intent(in) :: c_x_, c_y_
      real(kind=R8), dimension(ns)                 :: s_val

      !   Local variables
      real(kind=R8), dimension(np_, ns)  :: p_mat
      real(kind=R8), dimension(ns, ns)   :: v
      real(kind=R8), dimension(ns, ns)   :: a
      real(kind=R8), dimension(ns)       :: w
      real(kind=R8), dimension(ns)       :: b
      real(kind=R8), dimension(ns)       :: cf
      real(kind=R8), dimension(np_, np_) :: w_mat
      integer :: i

      !   Gives the matrix of the orthogonal basis functions
       p_mat = assemble_p_mat_dir (nds_pos_, c_x_, c_y_, np_, ns)

      !   Gives the matrix of the weight functions
      w_mat = assemble_w_mat (nds_pos_, c_x_, c_y_, np_)

      !   Calculates the coefficient matrix of the typical MLS.
      a = matmul(matmul(transpose(p_mat), w_mat), p_mat)

      !   Calculates the right hand side of the typical MLS
      b = matmul(matmul(transpose(p_mat), w_mat), nds_val_)

      !   SVD decomposition of the coefficient matrix.
      call svdcmp (a, ns, ns, ns, ns, w, v)

      !   Ensures that the problem is neither singular nor close to singular.
      do i = 1, ns
         if(abs (w(i)) < eps) then
            w(i) = 0.0d0
         end if
      end do

      !   Using the decomposed coefficients, calculates "cf" which is the
      !   desired shape functions.
      call svbksb(a, w, v, ns, ns, ns, ns, b, cf)

      s_val = cf
   end function mls_get_scalar_dir


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !  Using the points position (x,y) and their values or derivatives,
   !  approximates the the function at the required position.
   !  It assumes that the first given point is a point on a boundary with
   !  NEUMANN condition. Therefore, the first given value must be a directional
   !  derivate.
   !  As an input the user should provide the direction along which the
   !  derivative has been given (at the first source point)

   !  Inputs:
   !  1- np_: number of interpolant points
   !  2- nds_pos_: a dimension(2, np_) array holding the position of the
   !     source points.
   !  3- nds_val_: the value of the source points
   !  4- x_c_, y_c_: the point where the function should be interpolated
   !  5- norm_: Two components of the normal vector along which the
   !     directional derivatives for the first source point has been calculated.
   !
   !  Output:
   !  1- s_val: an array with the length of the "ns" ("ns" is the module
   !     variable) containing the calculated shape function.
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function mls_get_scalar_neu (np_, nds_pos_, nds_val_, c_x_, c_y_, norm_) result (s_val)
      implicit none
      !   Arguments
      integer(kind=I4)                , intent(in) :: np_
      real(kind=R8), dimension(2, np_), intent(in) :: nds_pos_
      real(kind=R8), dimension(np_)   , intent(in) :: nds_val_
      real(kind=R8), dimension(2)     , intent(in) :: norm_
      real(kind=R8)                   , intent(in) :: c_x_, c_y_
      real(kind=R8), dimension(ns)                 :: s_val

      !   Local variables
      real(kind=R8), dimension(np_, ns)            :: p_mat
      real(kind=R8), dimension(ns, ns)             :: v
      real(kind=R8), dimension(ns, ns)             :: a
      real(kind=R8), dimension(ns)                 :: w
      real(kind=R8), dimension(ns)                 :: b
      real(kind=R8), dimension(ns)                 :: cf
      real(kind=R8), dimension(np_, np_)           :: w_mat
      integer :: i

      !   Gives the matrix of the orthogonal basis functions for Neumann boundary
      p_mat = assemble_p_mat_neu (nds_pos_, c_x_, c_y_, norm_, np_, ns)

      !   Gives the matrix of the weight functions. The same weights as
      !   Dirichlet case
      w_mat = assemble_w_mat (nds_pos_, c_x_, c_y_, np_)

      !   Calculates the coefficient matrix of the typical MLS.
      a = matmul(matmul(transpose(p_mat), w_mat), p_mat)

      !   Calculates the right hand side of the typical MLS
      b = matmul(matmul(transpose(p_mat), w_mat), nds_val_)

      !   SVD decomposition of the coefficient matrix.
      call svdcmp (a, ns, ns, ns, ns, w, v)

      !   Ensures that the problem is neither singular nor close to singular.
      do i = 1, ns
         if(abs (w(i)) < eps ) w(i) = 0.d0
      end do

      !   Using the decomposed coefficients, calculates "cf" which is the
      !   desired shape functions.
      call svbksb(a, w, v, ns, ns, ns, ns, b, cf)

      s_val = cf
   end function mls_get_scalar_neu

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !   By the given value at the source points, returns the matrix of the
   !   basis function.
   !   In this function it is assumed that on the all of the given points the
   !   the interpolant function values are given. (Dirichlet boundary condition).
   !   The matrix is calculated based on the local coordinate reference.
   !   The origin is shifted to the position of the required point.
   !
   !   Output: the following matrix
   !    _                                                              _
   !   |  1  (x1-xc)  (y1-yc)  (x1-xc)**2  (x1-xc)*(y1-yc)  (y1-yc)**2  |
   !   |  1  (x2-xc)  (y2-yc)  (x2-xc)**2  (x2-xc)*(y2-yc)  (y2-yc)**2  |
   !   |                                                                |
   !   |  ...                                                           |
   !   |                                                                |
   !   |_ 1  (xn-xc)  (yn-yc)  (xn-xc)**2  (xn-xc)*(yn-yc)  (yn-yc)**2 _|
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function assemble_p_mat_dir (nds_pos_, c_x_, c_y_, np_, ns_) result (a)
      implicit none
      integer(kind=I4)                   :: i, np_, ns_
      real(kind=R8), dimension(2, np_)   :: nds_pos_
      real(kind=R8), dimension(np_, ns_) :: a
      real(kind=R8)                      :: c_x_, c_y_

      do i = 1, np_
         a(i,1) = 1.d0
         a(i,2) = (nds_pos_(1,i) - c_x_)
         a(i,3) = (nds_pos_(2,i) - c_y_)
         a(i,4) = 0.5d0 * (nds_pos_(1,i) - c_x_) * (nds_pos_(1,i) - c_x_)
         a(i,5) = 1.0d0 * (nds_pos_(1,i) - c_x_) * (nds_pos_(2,i) - c_y_)
         a(i,6) = 0.5d0 * (nds_pos_(2,i) - c_y_) * (nds_pos_(2,i) - c_y_)
      end do
   end function assemble_p_mat_dir


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !   By the given value at the source points, returns the matrix of the
   !   basis function.
   !   In this function it is assumed that on the all of the given points,
   !   except the first point the interpolant function values are given.
   !   On the first point the directional derivatives is given (Neumann
   !   boundary condition).
   !   The matrix is calculated based on the local coordinate reference. The
   !   origin is shifted to the position of the required point.
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function assemble_p_mat_neu  (nds_pos_, c_x_, c_y_, nrm_, np_, ns_) result (a)
      implicit none
      integer(kind=I4)                   :: np_, ns_
      real(kind=R8), dimension(2, np_)   :: nds_pos_
      real(kind=R8), dimension(np_, ns_) :: a
      real(kind=R8), dimension(2)        :: nrm_
      real(kind=R8)                      :: c_x_, c_y_
      integer(kind=I4)                   :: i

      a(1,1) = 0.d0
      a(1,2) = nrm_(1)
      a(1,3) = nrm_(2)
      a(1,4) = 2.d0 * (nds_pos_(1,1) - c_x_) * nrm_(1)
      a(1,5) = (nds_pos_(1,1) - c_x_) * nrm_(2) +  &
             & (nds_pos_(2,1) - c_y_) * nrm_(1)
      a(1,6) = 2.d0 * (nds_pos_(2,1) - c_y_) * nrm_(2)

      do i = 2, np_
         a(i,1) = 1.d0
         a(i,2) = (nds_pos_(1,i) - c_x_)
         a(i,3) = (nds_pos_(2,i) - c_y_)
         a(i,4) = 0.5d0 * (nds_pos_(1,i) - c_x_) * (nds_pos_(1,i) - c_x_)
         a(i,5) = 1.d0  * (nds_pos_(1,i) - c_x_) * (nds_pos_(2,i) - c_y_)
         a(i,6) = 0.5d0 * (nds_pos_(2,i) - c_y_) * (nds_pos_(2,i) - c_y_)
      end do
   end function assemble_p_mat_neu


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !    This function provides a diagonal matrix containing weight functions
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function assemble_w_mat (nds_pos_, c_x_, c_y_, np_) result(weights)
      implicit none
      integer(kind=I4),                 intent(in) :: np_
      real(kind=R8), dimension(2, np_), intent(in) :: nds_pos_
      real(kind=R8),                    intent(in) :: c_x_, c_y_
      real(kind=R8), dimension(np_, np_)           :: weights
      integer(kind=4)                             :: i
      real(kind=R8), dimension(np_)                :: dis
      real(kind=R8)                                :: dmi, r

      do i = 1, np_
         dis(i) = distance(nds_pos_(1,i), nds_pos_(2,i), c_x_, c_y_)
      end do

      dmi = maxval(dis)

      weights = 0.d0
      do i = 1, np_
         r = dis(i) / dmi
         if (r <= 0.5d0) then
            weights(i,i) = (2./3.) - (4.*(r**2)) + (4.*(r**3))
         else if (r <= 1.d0 .and. r > 0.5d0) then
            weights(i,i) = (4./3.) - (4.*r) + (4.*(r**2)) - ((4.*(r**3.))/3.)
         else
            weights(i,i) = 1.d0
         end if
      end do
   end function assemble_w_mat


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !  Using the points position (x,y) and their values, approximates the the
   !  function
   !  at the required position.
   !  It assumes that the first given point is a point on a boundary with
   !  DIRICHLET condition.
   !
   !  Inputs:
   !  1- np_: number of interpolant points
   !  2- nds_pos_: a dimension(2, np_) array holding the position of the
   !     source points.
   !  3- nds_val_: both components of vector values of the source points
   !  4- x_c_, y_c_: the point where the function should be interpolated
   !
   !  Output:
   !  1- s_val: an array with the length of the "ns" ("ns" is the module
   !     variable)  containing the calculated shape function.
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function mls_get_vec_dir (np_, nds_pos_, nds_val_, c_x_, c_y_) result (s_val)
      !TODO the second derivatives are incorrect. Assess the effect of 0.5
      !      coefficients in "assemble_p_mat_dir_vec" which was useful for scalar
      !      version, but without clear reason.
      implicit none
      !   Arguments
      integer(kind=I4)               , intent(in) :: np_
      real(kind=R8), dimension(2,np_), intent(in) :: nds_pos_
      real(kind=R8), dimension(2*np_), intent(in) :: nds_val_
      real(kind=R8)                  , intent(in) :: c_x_, c_y_
      real(kind=R8), dimension(2*ns)              :: s_val, cf
      !   Local variables
      real(kind=R8), dimension(2*np_+1, 2*ns)     :: p_mat
      real(kind=R8), dimension(2*ns, 2*ns)        ::  a, v
      real(kind=R8), dimension(2*np_+1)           ::  b_temp
      real(kind=R8), dimension(2*ns)              :: w
      real(kind=R8), dimension(2*ns)              :: b
      real(kind=R8), dimension(2*np_+1, 2*np_+1)  :: w_mat
      integer :: i

      !   Gives the matrix of the orthogonal basis functions
      p_mat = assemble_p_mat_dir_vec (nds_pos_, c_x_, c_y_, ns, np_)

      !   Gives the matrix of the weight functions. The same weights as
      !   Dirichlet case
      w_mat = assemble_w_mat_vec (nds_pos_, c_x_, c_y_, np_)

      !   Calculates the coefficient matrix of the typical MLS
      a = matmul(matmul(transpose(p_mat), w_mat), p_mat)

      b_temp(1:2*np_) = nds_val_(1:2*np_)
      b_temp(2*np_+1) = 0.d0

      !   Calculates the right hand side of the typical MLS
      b = matmul(matmul(transpose(p_mat), w_mat), b_temp)

      !   SVD decomposition of the coefficient matrix.
      call svdcmp (a, 2*ns, 2*ns, 2*ns, 2*ns, w, v)

      !   Ensures that the problem is neither singular nor close to singular.

      do i = 1, 2*ns
!         print*,  w(i)
!         if(abs (w(1) / w(i)) > 10e+10) w(i) = w(1) / 10e+10
!         print*, w(i)
         if (abs(w(i))  < eps)  w(i) = 0.0d0
      end do
      !   Using the decomposed coefficients, calculates "cf" which is the
      !   desired shape functions.
      call svbksb (a, w, v, 2*ns, 2*ns, 2*ns, 2*ns, b, cf)
      s_val = cf
   end function mls_get_vec_dir


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !   By the given value at the source points, returns the matrix of the
   !   basis function.
   !   The matrix is calculated based on the local coordinate reference. The
   !   origin is shifted to the position of the required point.
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function assemble_p_mat_dir_vec (nds_pos_, c_x_, c_y_, ns_, np_) result(a)
      implicit none
      integer(kind=I4)                  , intent(in) :: np_, ns_
      real(kind=R8), dimension(2, 2*np_), intent(in) :: nds_pos_
      real(kind=R8)                     , intent(in) :: c_x_, c_y_

      integer(kind=I4)                               :: i, j
      real(kind=R8), dimension(2*np_+1, 2*ns_)       :: a

      a = 0.d0

      do i = 1, 2*np_, 2
         j = (i+1) / 2
         a(i,1)    = 1.d0
         a(i,3)    = (nds_pos_(1,j) - c_x_)
         a(i,5)    = (nds_pos_(2,j) - c_y_)
         a(i,7)    = 0.5d0 * (nds_pos_(1,j) - c_x_) * (nds_pos_(1,j) - c_x_)
         a(i,9)    = 1.d0  * (nds_pos_(1,j) - c_x_) * (nds_pos_(2,j) - c_y_)
         a(i,11)   = 0.5d0 * (nds_pos_(2,j) - c_y_) * (nds_pos_(2,j) - c_y_)

         a(i+1,2)  = a(i,1)
         a(i+1,4)  = a(i,3)
         a(i+1,6)  = a(i,5)
         a(i+1,8)  = a(i,7)
         a(i+1,10) = a(i,9)
         a(i+1,12) = a(i,11)
      end do

         a(2*np_+1,1)  = 0.d0
         a(2*np_+1,2)  = 0.d0
         a(2*np_+1,3)  = 1.d0
         a(2*np_+1,4)  = 0.d0 ! 2.d0*x
         a(2*np_+1,5)  = 0.d0
         a(2*np_+1,6)  = 1.d0 ! y
         a(2*np_+1,7)  = 0.d0
         a(2*np_+1,8)  = 0.d0
         a(2*np_+1,9)  = 0.d0
         a(2*np_+1,10) = 0.d0
         a(2*np_+1,11) = 0.d0 ! 2.d0*y
         a(2*np_+1,12) = 0.d0 ! x
   end function assemble_p_mat_dir_vec


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function assemble_w_mat_vec (nds_pos_, c_x_, c_y_, np_) result(weights)
      implicit none
      integer(kind=I4)               , intent(in) :: np_
      real(kind=R8), dimension(2,np_), intent(in) :: nds_pos_
      real(kind=R8)                  , intent(in) :: c_x_, c_y_
      real(kind=R8), dimension(2*np_+1, 2*np_+1)  :: weights
      integer(kind=I4)                            :: i
      real(kind=R8), dimension(np_)               :: dis
      real(kind=R8)                               :: dmi, r

      do i = 1, np_
         dis(i) = distance(nds_pos_(1,i), nds_pos_(2,i), c_x_, c_y_)
      end do

      dmi = maxval(dis)

      weights = 0.d0

      do i = 1, np_

         r = dis(i) / dmi
         if (r <= 0.5d0) then
            weights(((2*i)-1),((2*i)-1)) =        &
!            (4.d0/3.d0) - (4.d0*r) + (4.d0*(r**2)) - ((4.d0*(r**3.d0))/3.d0)
            (2.d0/3.d0) - (4.d0*(r**2)) + (4.d0*(r**3))
         else if (r <= 1.d0 .and. r > 0.5d0) then
            weights(((2*i)-1),((2*i)-1)) =        &
!            (2.d0/3.d0) - (4.d0*(r**2)) + (4.d0*(r**3))
            (4.d0/3.d0) - (4.d0*r) + (4.d0*(r**2)) - ((4.d0*(r**3.d0))/3.d0)
         else
            weights(((2*i)-1),((2*i)-1)) = 1.d0
         end if

         weights((2*i),(2*i)) = weights(((2*i)-1),((2*i)-1))
      end do

      weights((2*np_)+1, (2*np_)+1) = 1.d0

!      weights(:, :) = 1.d0
   end function assemble_w_mat_vec



end module m_mls_interpolation
