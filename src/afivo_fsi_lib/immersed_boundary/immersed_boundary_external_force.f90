#include "cpp_macros.h"
submodule (m_immersed_boundary) ext_force


contains


    module procedure get_ext_forces
        implicit none

        integer(kind=I4) :: n, nn, nnn
        real(kind=R8)    :: ave(n_smoothed_var), domination

        call af_loop_box( tree, nullify_ext_forces_ker )

        distance_matrix_glb = max_distance
        do n = 1, this%nobjs
            this%objects(n)%ext_force(1) = 0.0_R8
            this%objects(n)%ext_force(2) = 0.0_R8
            call this%objects(n)%get_object_ext_force( tree, this%nobjs,                 &
            &  [i_unew, i_vnew, n, i_df] )

!            data_holder_2(1:3, 1:n_smoothed_var-1)=data_holder_1(1:3, 1:n_smoothed_var-1)
!            data_holder_1(1:3, 2:n_smoothed_var)  =data_holder_2(1:3, 1:n_smoothed_var-1)
            data_holder_1(1:n_smoothed_var, 2:n_time_ave) =                              &
            &  data_holder_1(1:n_smoothed_var, 1:n_time_ave-1)

            data_holder_1(1, 1)  = this%objects(n)%ext_force(1)
            data_holder_1(2, 1)  = this%objects(n)%ext_force(2)
            data_holder_1(3, 1)  = this%objects(n)%om_d(3)
            ave(:) = 0.0d0
            domination = 0.0d0
            do nn = 1, n_time_ave
                ave(:) = ave(:) +                                                        &
                &  (data_holder_1(:, nn) * real(n_time_ave + 1 - nn, kind=R8))
                domination = domination + nn
            end do

            ave(:) = ave(:) / domination
!            this%objects(n)%ext_force(1) = ave(1)
!            this%objects(n)%ext_force(2) = ave(2)
!            this%objects(n)%om_d(3)      = ave(3)
        end do
    end procedure get_ext_forces



    module procedure get_object_ext_force
        implicit none

        this%ext_force(1) = 0.0_R8
        this%ext_force(2) = 0.0_R8
        this%moment       = 0.0_R8
        call af_loop_tree (tree, get_object_ib_ext_force_ker)

    contains

        subroutine  get_object_ib_ext_force_ker (tree, id)
            implicit none
            type(af_t),       intent(inout) :: tree
            integer(kind=I4), intent(in)    :: id
            integer(kind=I4)                :: nc, i, j, ii, n
            integer(kind=I4), allocatable   :: tags(DTIMES(:))
            integer(kind=I4)                :: next_point(NDIM), ehat(NDIM)
            integer(kind=I4)                :: n_s, n_f, n_p
            logical                         :: are_all_available
            logical                         :: is_added
            real(kind=R8),    allocatable   :: cc (DTIMES(:), :)
            real(kind=R8)                   :: x_forcing(NDIM, 1), x_gama(NDIM), xyz(NDIM, 1)
            real(kind=R8)                   :: x_point(NDIM, 1), r(NDIM)
            real(kind=R8)                   :: inv_2dr(NDIM), normal(NDIM)
            real(kind=R8)                   :: poly(ns), poly_vec(2*ns), poly_p(ns) ! MLS
            real(kind=R8)                   :: inv_dt
            real(kind=R8)                   :: tmp_scord(NDIM, 8), tmp_svel(NDIM, 8)
            real(kind=R8)                   :: tmp_fcord(NDIM, 8), tmp_fvel(NDIM, 8)
            real(kind=R8)                   :: tmp_p_scord(NDIM, 8), tmp_p_sval(8)
            real(kind=R8)                   :: tmp_p_fcord(NDIM, 8), tmp_p_fval(8)
            real(kind=R8)                   :: base_normal(NDIM)
            real(kind=R8),     allocatable  :: cordinates(:, :), velocities(:, :)
            real(kind=R8),     allocatable  :: cordinatesP(:, :), valuesP(:), tmp(:)
            real(kind=R8)                   :: dis
            real(kind=R8)                   :: da
            real(kind=R8)                   :: beta, prov_vel(NDIM)
            real(kind=R8),      allocatable :: distance_matrix_loc(:, :)
            real(kind=R8)                   :: ex_force(NDIM), ex_momont, actual_vel(NDIM)

            nc = tree%boxes(id)%n_cell
            allocate( cordinates(NDIM, 1), velocities(NDIM, 1) )
            allocate( cordinatesP(NDIM,1), valuesP(1), tmp(1) )
            allocate( cc(DTIMES(-1:nc+2), size(i_var)+2 ) )
            allocate( tags(DTIMES(-1:nc+2)) )
            call af_gc2_box( tree, id, [i_df+i_var(3)-1, i_var(1), i_var(2), i_p], cc )


            tags = outside
            is_added = .false.
            inv_dt = 1.0_dp / dt
            dis = norm2( tree%boxes(id)%dr )
            da = tree%boxes(id)%dr(1) * tree%boxes(id)%dr(2)

            ex_force  = 0.0_R8
            ex_momont = 0.0_R8

!            if (allocated(distance_matrix_loc)) deallocate (distance_matrix_loc)
            allocate( distance_matrix_loc, source = distance_matrix_glb )
!            distance_matrix_loc =  max_distance

            do j = 1, nc
                do i = 1, nc
                    if (cc(i, j, 1) < 0.0_dp) then
                        tags(i, j) = inside
                        is_added = .true.
                    end if

                    if (                                 &
                         (cc(i, j, 1) > 0.0_dp) .and.    &
                         (                               &
                         cc(i-1, j,   1) < 0.0_dp .or.   &
                         cc(i+1, j,   1) < 0.0_dp .or.   &
                         cc(i,   j-1, 1) < 0.0_dp .or.   &
                         cc(i,   j+1, 1) < 0.0_dp .or.   &
                         cc(i-1, j-1, 1) < 0.0_dp .or.   &
                         cc(i-1, j+1, 1) < 0.0_dp .or.   &
                         cc(i+1, j-1, 1) < 0.0_dp .or.   &
                         cc(i+1, j+1, 1) < 0.0_dp        &
                         ) .and. (abs(cc(i, j, 1)) <= tree%boxes(id)%dr(1) )             &
                       )  then

                        tags(i, j) = forcing
                        is_added = .true.
                    end if
                end do
            end do

            if (is_added) then
                do j = 1, nc
                    do i = 1, nc
                        if (tags(i,j) == inside) then
                            xyz(:, 1) = af_r_cc (tree%boxes(id), [IJK])
                            r(:) = xyz(:, 1) - this%xc(:)
                            actual_vel(1) = this%xc_d(1) - (this%om_d(3) * r(2))
                            actual_vel(2) = this%xc_d(2) + (this%om_d(3) * r(1))
                            tree%boxes(id)%cc(i, j, i_ibfx) = &
                            ((actual_vel(1) - tree%boxes(id)%cc(i, j, i_var(1)) )*inv_dt)
                            tree%boxes(id)%cc(i, j, i_ibfy) = &
                            ((actual_vel(2) - tree%boxes(id)%cc(i, j, i_var(2)) )*inv_dt)

                            !> Adds the forces obtained from inside the object to the drag
                            !  and lift of the object.
                            if (tree%boxes(id)%lvl == tree%highest_lvl) then
                                ex_force(1) = ex_force(1) -                              &
                                &  tree%boxes(id)%cc(i, j, i_ibfx)*da
                                ex_force(2) = ex_force(2) -                              &
                                &  tree%boxes(id)%cc(i, j, i_ibfy)*da

                                ex_momont = ex_momont                                    &
                                &  - (tree%boxes(id)%cc(i, j, i_ibfy)*da) * r(1)         &
                                &  + (tree%boxes(id)%cc(i, j, i_ibfx)*da) * r(2)
                            end if
                        end if
                    end do
                end do

                inv_2dr = 1.0_dp / (2.0_dp * tree%boxes(id)%dr)
                do j = 1, nc
                    do i = 1, nc
                        if(tags(i,j) == forcing) then
                            are_all_available = .true.
                            n_s = 0
                            n_f = 0
                            n_p = 0
                            tmp_scord(:, :)   = 0.0_dp; tmp_svel(:, :)  = 0.0_dp
                            tmp_fcord(:, :)   = 0.0_dp; tmp_fvel(:, :)  = 0.0_dp
                            tmp_p_scord(:, :) = 0.0_dp; tmp_fcord(:, :) = 0.0_dp

                            !> The closest point on the physical boundary
                            normal(1) = (cc(i+1, j, 1) - cc(i-1, j, 1)) * inv_2dr(1)
                            normal(2) = (cc(i, j+1, 1) - cc(i, j-1, 1)) * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal))
                            base_normal(:) = normal(:)
                            x_forcing(:, 1) = af_r_cc (tree%boxes(id), [IJK])
                            x_gama(:) = x_forcing(:, 1) - normal*cc(i, j, 1)

                            n_s = n_s + 1
                            tmp_scord  (:, n_s) = x_gama(:)
                            tmp_svel (1, n_s) = this%xc_d(1) -                           &
                            &  (this%om_d(3) * (x_gama(2) - this%xc(2)))
                            tmp_svel (2, n_s) = this%xc_d(2) +                           &
                            &  (this%om_d(3) * (x_gama(1) - this%xc(1)))

                            tmp_p_scord(:, 1)   = x_gama(:)
                            tmp_p_sval (1) =                                             &
                            &  base_normal(1)*this%xc_dd(1) + base_normal(2)*this%xc_dd(2)

                            !> A vector indicating the direction for searching the
                            !  other interpolating points
                            ehat(:) = ceiling( normal(:) / abs(normal(:)))

                            !> The second point on the physical boundary
                            next_point(:) = [i,j] + [ehat(1), 0]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])

                            !> The fluid point corresponding to the second point on
                            !  the physical boundary
                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 n_f = n_f + 1
                                 n_p = n_p + 1
                                 tmp_fcord(:, n_f) = x_point(:, 1)
                                 tmp_fvel(1, n_f) = cc(next_point(1), next_point(2), 2)
                                 tmp_fvel(2, n_f) = cc(next_point(1), next_point(2), 3)
                                 tmp_p_fcord(:, n_p) = x_point(:, 1)
                                 tmp_p_fval(n_p)     = cc(next_point(1), next_point(2), 4)
                                 are_all_available = .false.
                            end if

                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) -         &
                                         cc(next_point(1)-1, next_point(2), 1))          &
                                         * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) -         &
                                         cc(next_point(1), next_point(2)-1, 1))          &
                                         * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal))
                            x_gama(:) =                                                  &
                            &    x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)
                            n_s = n_s + 1
                            tmp_scord (:, n_s) = x_gama(:)
                            tmp_svel (1, n_s) = this%xc_d(1) -                           &
                            &  (this%om_d(3) * (x_gama(2) - this%xc(2)))
                            tmp_svel (2, n_s) = this%xc_d(2) +                           &
                            &  (this%om_d(3) * (x_gama(1) - this%xc(1)))

                            !> The third point on the physical boundary
                            next_point(:) = [i,j] + [0, ehat(2)]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) -         &
                            &           cc(next_point(1)-1, next_point(2), 1))           &
                            &           * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) -         &
                            &           cc(next_point(1), next_point(2)-1, 1))           &
                            &           * inv_2dr(2)
                            normal(:) = normal(:) / norm2 (normal)
                            x_gama(:) =                                                  &
                            &  x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)

                            n_s = n_s + 1
                            tmp_scord (:, n_s) = x_gama(:)
                            tmp_svel (1, n_s) = this%xc_d(1) -                           &
                            &  (this%om_d(3) * (x_gama(2) - this%xc(2)))
                            tmp_svel (2, n_s) = this%xc_d(2) +                           &
                            &  (this%om_d(3) * (x_gama(1) - this%xc(1)))

                            !> The fluid point corresponding to the third point on
                            !  on the physical boundary
                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 n_f = n_f + 1
                                 n_p = n_p + 1
                                 tmp_fcord(:, n_f) = x_point(:, 1)
                                 tmp_fvel(1, n_f)  = cc(next_point(1), next_point(2), 2)
                                 tmp_fvel(2, n_f)  = cc(next_point(1), next_point(2), 3)
                                 tmp_p_fcord(:, n_p) = x_point(:, 1)
                                 tmp_p_fval(n_p)     = cc(next_point(1), next_point(2), 4)
                            end if

                            !> The first fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [ehat(1), ehat(2)]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])

                            if (tags(next_point(1), next_point(2)) .ne. forcing ) then
                                n_f = n_f + 1
                                n_p = n_p + 1
                                tmp_fcord(:, n_f) = x_point(:, 1)
                                tmp_fvel(1, n_f) = cc(next_point(1), next_point(2), 2)
                                tmp_fvel(2, n_f) = cc(next_point(1), next_point(2), 3)
                                tmp_p_fcord(:, n_p) = x_point(:, 1)
                                tmp_p_fval(n_p) = cc(next_point(1), next_point(2), 4)
                            end if

                            !> The second fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [2*ehat(1), 0]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])

                            if (tags(next_point(1), next_point(2)) .ne. forcing ) then
                                n_f = n_f + 1
                                n_p = n_p + 1
                                tmp_fcord(:, n_f) = x_point(:, 1)
                                tmp_fvel(1, n_f)  = cc(next_point(1), next_point(2), 2)
                                tmp_fvel(2, n_f)  = cc(next_point(1), next_point(2), 3)
                                tmp_p_fcord(:, n_p) = x_point(:, 1)
                                tmp_p_fval(n_p) = cc(next_point(1), next_point(2), 4)
                            end if

                            !> The third fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [0, 2*ehat(2)]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])

                            if (tags(next_point(1), next_point(2)) .ne. forcing ) then
                                n_f = n_f + 1
                                n_p = n_p + 1
                                tmp_fcord(:, n_f) = x_point(:, 1)
                                tmp_fvel(1, n_f) = cc(next_point(1), next_point(2), 2)
                                tmp_fvel(2, n_f) = cc(next_point(1), next_point(2), 3)
                                tmp_p_fcord(:, n_p) = x_point(:, 1)
                                tmp_p_fval(n_p) = cc(next_point(1), next_point(2), 4)
                            end if

                            deallocate( cordinates, velocities, cordinatesP, valuesP, tmp )
                            allocate( velocities(NDIM, n_s+n_f) )
                            allocate( cordinates(NDIM, n_s+n_f) )
                            allocate( cordinatesP(NDIM, 1+n_p) )
                            allocate( valuesP(1+n_p) )
                            allocate( tmp(2*(n_s+n_f)) )
                            velocities(:, 1:n_s)         = tmp_svel(:, 1:n_s)
                            velocities(:, n_s+1:n_f+n_s) = tmp_fvel(:, 1:n_f)
                            cordinates(:, 1:n_s)         = tmp_scord(:, 1:n_s)
                            cordinates(:, n_s+1:n_f+n_s) = tmp_fcord(:, 1:n_f)
                            valuesP(1)                   = tmp_p_sval(1)
                            valuesP(2:n_p+1)             = tmp_p_fval(1:n_p)
                            cordinatesP(:, 1)            = tmp_p_scord(:, 1)
                            cordinatesP(:, 2:n_p+1)      = tmp_p_fcord(:, 1:n_p)

!                            do ii = 1, size(tmp) / 2
!                                tmp((2*ii))   = velocities(2, ii)
!                                tmp((2*ii)-1) = velocities(1, ii)
!                            end do

!                            poly_vec(:) = mls_get_vec_dir (n_s+n_f, cordinates,          &
!                            &  tmp, x_forcing(1,1), x_forcing(2,1))
                            !>>>
!                            poly_p(:) = mls_get_scalar_neu (1+n_p, cordinatesP, valuesP, &
!                            &  x_forcing(1, 1),  x_forcing(2, 1) , base_normal)
!                            tree%boxes(id)%cc(i, j, i_p) = poly_p(1)
                            !<<<

                            beta = cc(i, j, 1) / tree%boxes(id)%dr(1)

                            poly(:) = &
                            mls_get_scalar_dir( n_s+n_f, cordinates, velocities(1, :),   &
                                                x_forcing(1,1), x_forcing(2,1) )

                            prov_vel(1) = (beta * tree%boxes(id)%cc(i, j, i_var(1))) +   &
                                          ((1.d0 - beta) * poly(1))

                            tree%boxes(id)%cc(i, j, i_ibfx) =                            &
                            ((prov_vel(1) - tree%boxes(id)%cc(i, j, i_var(1))) * inv_dt) !&
!                            + tree%boxes(id)%cc(i, j, i_resx)

                            poly(:) = &
                            mls_get_scalar_dir( n_s+n_f, cordinates, velocities(2, :),   &
                                                x_forcing(1,1), x_forcing(2,1) )

                            prov_vel(2) = (beta * tree%boxes(id)%cc(i, j, i_var(2))) +   &
                            &             ((1.d0 - beta) * poly(1))

                            tree%boxes(id)%cc(i, j, i_ibfy) =                            &
                                 ((prov_vel(2) - tree%boxes(id)%cc(i, j, i_var(2))) * inv_dt)  !&
!                            + tree%boxes(id)%cc(i, j, i_resy)

                             !> Gets the distance of this object with the others as well
                             !  as the walls.
!                            do n = 1, nobjects + nwalls
!                                if (n==this%id) cycle
!                            end do
                            if ( abs(x_forcing(1,1) - 0.0_R8 )  <=                       &
                               distance_matrix_loc(this%id, nobjects+1)) then
                               distance_matrix_loc(this%id, nobjects+1) =                &
                               abs( x_forcing(1,1) - 0.0_R8 )
                            end if

                            if (abs( x_forcing(1,1) - domain_len(1) )  <=                &
                            & distance_matrix_loc(this%id, nobjects+2)) then
                                distance_matrix_loc(this%id, nobjects+2) =               &
                                &  abs( x_forcing(1,1) - domain_len(1) )
                            end if

                            !$omp critical
                            if ( abs(x_forcing(2,1) - 0.0_R8 )  <                        &
                            &  distance_matrix_loc(this%id, nobjects+3)) then


!                            print*, i, j, x_forcing(2,1),  distance_matrix_loc(this%id, nobjects+3), id


!                                distance_matrix_loc(this%id, nobjects+3) =               &
!                                &  abs(x_forcing(2,1) - 0.0_R8)

                            end if
                            !$omp end critical


                            if ( abs( x_forcing(2,1) - domain_len(2))  <=                &
                            &  distance_matrix_loc(this%id, nobjects+4) ) then
                                distance_matrix_loc(this%id, nobjects+4) =               &
                                &  abs( x_forcing(2,1) - domain_len(2) )
                            end if

!                             print*,  distance_matrix_loc(this%id, nobjects+3), x_forcing(2,1)

                            !> Accumulates the forces calculated on immersed surface into
                            !  the Drag and lift already contains the value from inside
                            !  the object.
                            if ( tree%boxes(id)%lvl == tree%highest_lvl) then
                                xyz(:, 1) = af_r_cc (tree%boxes(id), [IJK] )
                                r(:) = xyz(:, 1) - this%xc(:)

                                ex_force(1) = ex_force(1) -                              &
                                &  tree%boxes(id)%cc(i, j, i_ibfx)*da
                                ex_force(2) = ex_force(2) -                              &
                                &  tree%boxes(id)%cc(i, j, i_ibfy)*da

                                ex_momont = ex_momont &
                                &  - (tree%boxes(id)%cc(i, j, i_ibfy)*da) * r(1)         &
                                &  + (tree%boxes(id)%cc(i, j, i_ibfx)*da) * r(2)
                            end if
                        end if
                    end do
                end do
            end if

            !> Finds the smallest distance among all boxes containing this object.
            !$omp critical
            !> Distance to left wall
            if ( distance_matrix_loc(this%id, nobjects+1) <                              &
                &   distance_matrix_glb(this%id, nobjects+1) ) then
                distance_matrix_glb(this%id, nobjects+1) =                               &
                &   distance_matrix_loc(this%id, nobjects+1)
                distance_matrix_glb(nobjects+1, this%id) =                               &
                &   distance_matrix_loc(nobjects+1, this%id)
            end if

!            print*,  distance_matrix_glb(this%id, nobjects+3)
            !> Distance to right wall
            if ( distance_matrix_loc(this%id, nobjects+2) <                              &
            &  distance_matrix_glb(this%id, nobjects+2) ) then
                distance_matrix_glb(this%id, nobjects+2) =                               &
                &   distance_matrix_loc(this%id, nobjects+2)
                distance_matrix_glb(nobjects+2, this%id) =                               &
                &   distance_matrix_loc(nobjects+2, this%id)
            end if

            !> Distance to bottom wall
            if ( distance_matrix_loc(this%id, nobjects+3) <=                             &
            &  distance_matrix_glb(this%id, nobjects+3) ) then
                distance_matrix_glb(this%id, nobjects+3) =                               &
                &  distance_matrix_loc(this%id, nobjects+3)
!                distance_matrix_glb(, nobjects+3, this%id) =                               &
!                &  distance_matrix_loc(nobjects+3, this%id)
            end if
!            print*, distance_matrix_glb(this%id, :); print*, " "

            !> Distance to top wall
            if ( distance_matrix_loc(this%id, nobjects+4) <                              &
            &  distance_matrix_glb(this%id, nobjects+4) ) then
                distance_matrix_glb(this%id, nobjects+4) =                               &
                &  distance_matrix_loc(this%id, nobjects+4)
                distance_matrix_glb(nobjects+4, this%id) =                               &
                &  distance_matrix_loc(nobjects+4, this%id)
            end if
            distance_matrix_glb(this%id, this%id) = 0.0_R8

            this%ext_force(1) = this%ext_force(1) + ex_force(1)
            this%ext_force(2) = this%ext_force(2) + ex_force(2)
            this%moment       = this%moment + ex_momont
            !$omp end critical
        end subroutine  get_object_ib_ext_force_ker
    end procedure get_object_ext_force



    subroutine nullify_ext_forces_ker( box )
        implicit none
        type(box_t), intent(inout) :: box
        integer(kind=I4)           :: nc, i, j

        nc = box%n_cell
        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_ibfx) = 0.0_R8
                box%cc(i, j, i_ibfy) = 0.0_R8
            end do
        end do
    end subroutine nullify_ext_forces_ker



    module procedure write_out_ext_forces
        implicit none
        integer(kind=I4) :: n

        do n = 1, this%nobjs
            call this%objects(n)%write_out_ext_force (n)
        end do
    end procedure write_out_ext_forces



    module procedure write_out_ext_force
        implicit none
        call write_ext_forces_on_hard( i_obj, field_time,                                &
        &  this%ext_force(1), this%ext_force(2), this%moment )
    end procedure write_out_ext_force



    module procedure write_out_objects_energy
        implicit none
        integer(kind=I4) :: n

        do n = 1, this%nobjs
            call this%objects(n)%write_out_object_energy( n )
        end do
    end procedure write_out_objects_energy



    module procedure write_out_object_energy
        implicit none
        call write_enery_of_object_on_hard( i_obj, field_time, this%et, this%er )
    end procedure write_out_object_energy



    module procedure write_out_objects_reynolds
        implicit none
        integer(kind=I4) :: n

        do n = 1, this%nobjs
            call this%objects(n)%write_out_object_reynolds( n )
        end do
    end procedure write_out_objects_reynolds



    module procedure write_out_object_reynolds
        implicit none
        call  write_reynolds_of_object_on_hard( i_obj, field_time, this%rey )
    end procedure write_out_object_reynolds

end submodule  ext_force
