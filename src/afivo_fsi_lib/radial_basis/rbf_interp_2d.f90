module m_rbf_iterp_2d
    use m_math_lib


contains


    subroutine rbf_interp ( m, nd, xd, w, ni, xi, fi )
    !*****************************************************************************80
    !
    !! RBF_INTERP evaluates a radial basis function interpolant.
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !
    !    Input, integer ( kind = 4 ) ND, the number of data points.
    !
    !    Input, real ( kind = 8 ) XD(M,ND), the data points.
    !
    !    Input, real ( kind = 8 ) W(ND), the weights, as computed by RBF_WEIGHTS.
    !
    !    Input, integer ( kind = 4 ) NI, the number of interpolation points.
    !
    !    Input, real ( kind = 8 ) XI(M,NI), the interpolation points.
    !
    !    Output, real ( kind = 8 ) FI(NI), the interpolated values.
    !
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) nd
      integer ( kind = 4 ) ni
      real ( kind = 8 ) fi(ni)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      real ( kind = 8 ) r(nd)
      real ( kind = 8 ) v(nd)
      real ( kind = 8 ) w(nd+3)
      real ( kind = 8 ) xd(m,nd)
      real ( kind = 8 ) xi(m,ni)

      do i = 1, ni
        do j = 1, nd
          r(j) = sqrt ( sum ( ( xi(1:m,i) - xd(1:m,j) )**2 ) )
        end do
        fi(i) = dot_product ( r, w(1:nd) ) + w(nd+1) + (w(nd+2)*xi(1, i)) + (w(nd+3)*xi(2, i))
      end do
    end subroutine rbf_interp


    subroutine rbf_weight ( m, nd, xd, fd, w )
    !*****************************************************************************80
    !
    !! RBF_WEIGHT computes weights for radial basis function interpolation.
    !
    !  Discussion:
    !
    !    We assume that there are N (nonsingular) equations in N unknowns.
    !
    !    However, it should be clear that, if we are willing to do some kind
    !    of least squares calculation, we could allow for singularity,
    !    inconsistency, or underdetermine systems.  This could be associated
    !    with data points that are very close or repeated, a smaller number
    !    of data points than function values, or some other ill-conditioning
    !    of the system arising from a peculiarity in the point spacing.
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !
    !    Input, integer ( kind = 4 ) ND, the number of data points.
    !
    !    Input, real ( kind = 8 ) XD(M,ND), the data points.
    !
    !    Input, real ( kind = 8 ) FD(ND), the function values at the data points.
    !
    !    Output, real ( kind = 8 ) W(ND), the weights.
    !
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) nd

      real ( kind = 8 ) a(nd+3,nd+3)
      real ( kind = 8 ) fd(nd+3)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      real ( kind = 8 ) r(nd)
      real ( kind = 8 ) v(nd)
      real ( kind = 8 ) w(nd+3)
      real ( kind = 8 ) xd(m,nd)

      do i = 1, nd
        do j = 1, nd
          r(j) = sqrt ( sum ( ( xd(1:m,i) - xd(1:m,j) ) ** 2 ) )
        end do
        a(i,1:nd) = r(1:nd)
      end do

      a(nd+1, 1:nd)      = 1
      a(nd+1, nd+1:nd+3) = 0

      a(nd+2, 1:nd)      = xd(1,1:nd)
      a(nd+2, nd+1:nd+3) = 0.0d0

      a(nd+3, 1:nd)      = xd(2,1:nd)
      a(nd+3, nd+1:nd+3) = 0.0d0

      a(1:nd, nd+1) = 1.0d0
      a(1:nd, nd+2) = xd(1,1:nd)
      a(1:nd, nd+3) = xd(2,1:nd)

      fd (nd+1:nd+3) = 0.0d0
    !
    !  Solve for the weights.
    !
      call r8mat_solve_svd ( nd+3, nd+3, a, fd, w )

    end subroutine rbf_weight
end module m_rbf_iterp_2d
