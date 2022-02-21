#include "cpp_macros.h"
module m_immersed_boundary
    use omp_lib
    use m_config
    use m_af_all
    use m_solver_parameters
    use m_incomp_flow_field
    use m_rbf_iterp_2d
    use m_mls_interpolation

    implicit none
    !> These variables are used as temporary placeholder for building each object
    real(kind=dp) :: landa_g = 1.0_dp
    real(kind=dp) :: domain_len(NDIM)
    integer, parameter :: nstencil = 8
    integer, parameter :: nnmax = 8


    type marker_t
        real(kind=dp) :: loc(NDIM)
    end type marker_t

    type object_t
        integer       :: id       ! object id
        integer       :: nmarkers
        real(kind=dp) :: center(NDIM)
        real(kind=dp) :: velocity(NDIM)
        real(kind=dp) :: a_velocity(NDIM) ! angular velocity
        real(kind=dp) :: m        ! mass
        real(kind=dp) :: im       ! mass moment of inertia
        real(kind=dp),  allocatable :: ws(:) ! Weights fro RBF representation
        type(marker_t), allocatable :: marker_list(:)
    contains
        procedure, private :: read_markers
        procedure, public :: get_object_force
!        procedure, public :: get_object_interpolated_velocity
        procedure, public :: get_object_extrapolated_velocity
!        procedure, public :: move_object
    end type object_t

    type immersed_boundary_t
        type(t_config) :: config
        type(object_t), allocatable :: objects(:)
        integer :: nobjs
        real(dp):: domain_len(NDIM)
    contains
        procedure, public :: get_forces
!        procedure, public :: get_interpolated_velocity
        procedure, public :: get_extrapolated_velocity
!        procedure, public :: move_objects
    end type immersed_boundary_t


contains


    subroutine build_immersed_boundary (tree, refine_info, config, ib)
        class(t_config), intent(in) :: config
        class(af_t), intent(inout) :: tree
        class(ref_info_t), intent(inout) :: refine_info
!        class(incomp_flow_field_t), intent(inout) :: fld
        class (immersed_boundary_t), intent(out) :: ib
        integer :: n

        ib%config = config

        !> How many objects are involved in this problem?
        call config%get(section_name='immersed_boundary', option_name='nobjects', &
        &   val=ib%nobjs)

        !> Allocates an array with length equal to number of objects
        allocate (ib%objects(ib%nobjs))

        call config%get(section_name='field', option_name='domain_len',    &
        &   val=domain_len)

        !> Builds all objects involved in this problem
        do n = 1, ib%nobjs
            call build_object (ib%objects(n), tree, refine_info, is_refined = .true.)
        end do
    end subroutine build_immersed_boundary



    subroutine build_object (obj, tree, refine_info, is_initiated, is_refined)
        class(object_t), intent(inout) :: obj
        class(af_t), intent(inout) :: tree
        class(ref_info_t), intent(inout) :: refine_info
!        class(incomp_flow_field_t), intent(inout) :: fld
        real(kind=dp), allocatable, dimension(:,:) :: coords
        real(kind=dp), allocatable, dimension(:)   :: valz
        real(kind=dp) :: e, volume, r0
        real(kind=dp) :: mincoords(NDIM), maxcoords(NDIM)
        logical, intent(in), optional :: is_initiated
        logical, intent(in), optional :: is_refined
        integer :: n
#if NDIM == 2
        integer :: ncorners = 4
        real(kind=dp) :: corner_distance(4)
#endif

        if (present(is_initiated) ) then
            if(is_initiated) then
                call obj%read_markers ()
            end if
        else
            call obj%read_markers ()
            obj%velocity(:) = predefined_velocity(time = 0.0_dp)
        end if

        !>  Calculated the distance function with RBF method
        call get_initial_df ()

        !> If "is_refined" has not been passed, the solver refines it by default.
        if (present (is_refined)) then
            if (is_refined) then
                !> Refinement based on the distance function.
                do
                    call af_gc_tree (tree, [i_p, i_uold, i_unew, i_vold, i_vnew])
                    call af_adjust_refinement (tree, ref_routine, refine_info)
                    ! If no new boxes have been added, exit the loop
                    if (refine_info%n_add == 0) exit
                end do
            end if
        else
            !> Refinement based on the distance function.
            do
                call af_gc_tree (tree, [i_p, i_uold, i_unew, i_vold, i_vnew])
                call af_adjust_refinement (tree, ref_routine, refine_info)
                ! If no new boxes have been added, exit the loop
                if (refine_info%n_add == 0) exit
            end do
        end if

!        call get_final_df ()
        call af_restrict_tree (tree, [i_df])
        call af_gc_tree (tree, [i_df])


    contains


        subroutine get_final_df()
            call af_loop_box (tree, get_final_df_ker, leaves_only=.true.)
        end subroutine get_final_df


        subroutine get_final_df_ker(box)
            implicit none
            type(box_t), intent(inout) :: box
            integer :: nc, i, j, n
            real(kind=dp) ::  xi(NDIM, 1), temp

            nc = box%n_cell
            do j = 1, nc
                do i = 1, nc
                    xi(:, 1) = af_r_cc(box, [IJK])
                    box%cc(i, j, i_df) = sign (100.d0, box%cc(i, j, i_df) )
                    do n = 1, size (obj%marker_list)
                        temp = norm2 (                       &
                        xi(:,1) - obj%marker_list(n)%loc(:))

                        box%cc(i, j, i_df) =                &
                        sign (min (temp, abs (box%cc(i, j, i_df))), box%cc(i, j, i_df))
                    end do
                end do
            end do
        end subroutine get_final_df_ker



        subroutine ref_routine(box, cell_flags)
            type(box_t), intent(in) :: box ! A list of all boxes in the tree
            integer, intent(out) :: cell_flags(DTIMES(box%n_cell))
            integer                  :: IJK, nc

            nc = box%n_cell
            do KJI_DO(1,nc)
                if( box%lvl < max_refinement                  &
                       .and.                                  &
                    abs (box%cc(i, j, i_df)) < 1.d0*box%dr(1)  &
                  ) then
                    cell_flags(IJK) = af_do_ref
                else
                    cell_flags(IJK) = af_keep_ref
                end if
            end do; CLOSE_DO
        end subroutine ref_routine


        subroutine get_initial_df()
            implicit none
            if (.not. allocated (coords) ) then
                allocate(coords(NDIM, obj%nmarkers + ncorners))
            end if

            if(.not. allocated (valz)) then
                allocate(valz(obj%nmarkers + ncorners))
            end if
            valz(1:obj%nmarkers) = landa_g

            if (.not. allocated (obj%ws)) then
                allocate(obj%ws(obj%nmarkers + ncorners))
            end if

            coords(1, 1:obj%nmarkers) = obj%marker_list(:)%loc(1)
            coords(2, 1:obj%nmarkers) = obj%marker_list(:)%loc(2)

            call get_corner_distance()

            do n = 1, ncorners
                if(n==1 .or. n==4)  then
                    coords(1, obj%nmarkers + n) = domain_len(1)
                else
                    coords(1, obj%nmarkers + n) = 0.0_dp
                end if

                if (n==1 .or. n==2) then
                    coords(2, obj%nmarkers + n) = domain_len(2)
                else
                    coords(2, obj%nmarkers + n) = 0.0_dp
                end if
                valz(obj%nmarkers + n) =  corner_distance(n) + landa_g
            end do

            mincoords = minval (coords, dim=2)
            maxcoords = maxval (coords, dim=2)
            volume = ( maxcoords(1) - mincoords(1) ) * ( maxcoords(2) - mincoords(2))
            e = 1.0D+00 / real (NDIM, dp)
            r0 = (volume / obj%nmarkers) ** e

            obj%ws = 0.d0
            call rbf_weight (2, obj%nmarkers+ncorners, coords, r0, phi3, valz, obj%ws)
            call af_loop_box (tree, get_df_ker)
            !> Fills ghost cells for distance function
            call af_gc_tree(tree, [i_df])
        end subroutine get_initial_df


        subroutine get_df_ker(box)
            type(box_t), intent(inout) :: box
            integer :: nc, i, j, nd
            real(kind=dp) ::  fi(1), xi(NDIM, 1)
            integer :: ni

            nd = size (coords, 2)
            ni = 1
            nc = box%n_cell
            do j = 1, nc
                do i = 1, nc
                    xi(:, 1) = af_r_cc(box, [IJK])
                    call rbf_interp (NDIM, nd, coords, r0, phi3, obj%ws, ni, xi, fi)
                    box%cc(i, j, i_df) = fi(1) - landa_g
                end do
            end do
        end subroutine get_df_ker


        subroutine get_corner_distance()
                implicit none
                integer :: n, corner
                real(kind=dp) :: temp_dis

                do corner = 1, ncorners
                    corner_distance(corner) =  &
                        max (domain_len(1), domain_len(2))
                    do n = 1, size (obj%marker_list)
                        temp_dis = norm2 (                                        &
                        obj%marker_list(n)%loc(:) - coords(:, obj%nmarkers + corner))

                        if(temp_dis <= corner_distance(corner)) then
                            corner_distance(corner) = temp_dis
                        end if
                    end do
                end do
        end subroutine get_corner_distance
    end subroutine build_object



    subroutine get_forces (this, tree)
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t), intent(inout) :: tree
        integer :: n

        do n = 1, this%nobjs
            call this%objects(n)%get_object_force (tree, [i_unew, i_vnew])
        end do

    end subroutine get_forces



    subroutine get_object_force (this, tree, i_var)
        class(object_t), intent(inout) :: this
        class(af_t), intent(inout) :: tree
        integer, intent(in) :: i_var(:)

!        call af_loop_tree(fld%tree, get_interpolated_ker, leaves_only=.true.)
        call af_loop_tree (tree, get_object_hybrid_force_ker) !!!!!!!!!!!!!!!!!!!!!
        call af_gc_tree (tree, [i_p, i_var(1), i_var(2)])
        call af_restrict_tree (tree, [i_p, i_var(1), i_var(2)])
!        call af_restrict_tree (tree, [i_var(1)])
!        call af_restrict_tree (tree, [i_var(2)])


    contains


        subroutine  get_object_hybrid_force_ker (tree, id)
            implicit none
            type(af_t), intent(inout) :: tree
            integer, intent(in) :: id
            integer :: nc, i, j
            logical :: is_added
            real(kind=dp), allocatable :: cc (DTIMES(:), :)
            integer, allocatable :: tags(DTIMES(:))
            real(kind=dp) :: x_forcing(NDIM, 1), x_gama(NDIM)
            real(kind=dp) :: x_point(NDIM, 1)
            real(kind=dp) :: inv_2dr(NDIM), normal(NDIM)
            ! MLS
            real(kind=dp) :: cnodes(NDIM, nstencil), poly(ns)
            real(kind=dp) :: vnodes(size(i_var), nstencil)
            integer :: nnodes, ntemp
            integer :: next_point(NDIM), ehat(NDIM)
            logical :: are_all_available
            real(kind=dp) :: inv_dt
            integer :: cntr

            nc = tree%boxes(id)%n_cell
            allocate (cc(DTIMES(-1:nc+2), size(i_var)+1 ))
            allocate (tags(DTIMES(-1:nc+2)))
            call af_gc2_box (tree, id, [i_df, i_var(1), i_var(2)], cc)

            tags = outside
            is_added = .false.
            tree%boxes(id)%cc(:, :, i_ibfx) = 0.0_dp
            tree%boxes(id)%cc(:, :, i_ibfy) = 0.0_dp
            inv_dt = 1.0_dp / dt

            do j = 1, nc
                do i = 1, nc
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
                         )                               &
                       )  then

                        tags(i, j) = forcing
                        is_added = .true.
                    end if

                    if (cc(i, j, 1) < 0.0_dp) then
                        tags(i, j) = inside
                    end if
                end do
            end do

            do j = 1, nc
                do i = 1, nc
                    if (tags(i,j) == inside) then
                        tree%boxes(id)%cc(i, j, i_ibfx) = &
                        ((this%velocity(1) - tree%boxes(id)%cc(i, j, i_var(1)) )*inv_dt) &
                        + tree%boxes(id)%cc(i, j, i_resx)
!                        tree%boxes(id)%cc(i, j, i_unew) = this%velocity(1)

                        tree%boxes(id)%cc(i, j, i_ibfy) = &
                        ((this%velocity(2) - tree%boxes(id)%cc(i, j, i_var(2)) )*inv_dt) &
                        + tree%boxes(id)%cc(i, j, i_resy)
!                        tree%boxes(id)%cc(i, j, i_vnew) = this%velocity(2)
                    end if
                end do
            end do

!            is_added = .false.


            if (is_added) then
                inv_2dr = 1.0_dp / (2.0_dp * tree%boxes(id)%dr)
                do j = 1, nc
                    do i = 1, nc

                        if(tags(i,j) == forcing) then
                            nnodes = 0
                            ntemp = 0
                            vnodes = 0.d0
                            cnodes = 0.d0
                            are_all_available = .true.

                            !> The closest point on the physical boundary
                            normal(1) = (cc(i+1, j, 1) - cc(i-1, j, 1)) * inv_2dr(1)
                            normal(2) = (cc(i, j+1, 1) - cc(i, j-1, 1)) * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal))
                            x_forcing(:, 1) = af_r_cc (tree%boxes(id), [IJK])
                            x_gama(:) = x_forcing(:, 1) - normal*cc(i, j, 1)
                            nnodes = nnodes + 1  ! nnodes = 1
                            cnodes(:, nnodes) = x_gama(:)
                            vnodes(:, nnodes) = this%velocity(:)

                            !> A vector indicating the direction for searching the
                            !  other interpolating points
                            ehat(:) = ceiling( normal(:) / abs(normal(:)))

                            !> The second point on the physical boundary
                            next_point(:) = [i,j] + [ehat(1), 0]
                            x_point(:, 1) =  &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
                                         cc(next_point(1)-1, next_point(2), 1))  &
                                         * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
                                         cc(next_point(1), next_point(2)-1, 1))  &
                                         * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal))
                            x_gama(:) = &
                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)

                            nnodes = nnodes + 1 ! nnodes = 2
                            cnodes(:, nnodes) = x_gama(:)
                            vnodes(:, nnodes) = this%velocity(:)

                            !> The fluid point corresponding to the second point on
                            !  the physical boundary
                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 if (are_all_available) then
                                    cnodes(:, nnmax-1) =  x_point(:, 1)
                                    vnodes(1, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 3)
                                 else
                                    cnodes(:, nnmax)   = x_point(:, 1)
                                    vnodes(1, nnmax)   = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax)   = &
                                        cc(next_point(1), next_point(2), 3)
                                 end if
                                 are_all_available = .false.
                                 ntemp = ntemp + 1
                            end if

                            !> The third point on the physical boundary
                            next_point(:) = [i,j] + [0, ehat(2)]
                            x_point(:, 1) =  &
                            af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
                                         cc(next_point(1)-1, next_point(2), 1))  &
                                         * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
                                         cc(next_point(1), next_point(2)-1, 1))  &
                                         * inv_2dr(2)
                            normal(:) = normal(:) / norm2 (normal)
                            x_gama(:) = &
                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)

                            nnodes = nnodes + 1
                            cnodes(:, nnodes) = x_gama
                            vnodes(:, nnodes) = this%velocity(:)

                            !> The fluid point corresponding to the third point on
                            !  on the physical boundary
                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 if (are_all_available) then
                                    cnodes(:, nnmax-1) = x_point(:, 1) ! x_gama(:)
                                    vnodes(1, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 3)
                                    are_all_available = .false.
                                    ntemp = ntemp + 1
                                 else
                                    cnodes(:, nnmax) = x_point(:, 1)! x_gama(:)
                                    vnodes(1, nnmax) = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax) = &
                                        cc(next_point(1), next_point(2), 3)
                                    are_all_available = .false.
                                    ntemp = ntemp + 1
                                 end if

                            end if

                            !> The first fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [ehat(1), ehat(2)]
                            x_point(:, 1) =   &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes)  =  x_point(:, 1)
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            !> The second fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [2*ehat(1), 0]
                            x_point(:, 1) =   &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes) = x_point(:, 1)
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            !> The third fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [0, 2*ehat(2)]
                            x_point(:, 1) =   &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes) = x_point(:, 1)

                            !> The interpolation velocity are chosen from the 2 and 3
                            !  component of the "cc" array.
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            poly(:) = &
                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(1,:),   &
                                                x_forcing(1,1), x_forcing(2,1))
!                            tree%boxes(id)%cc(i, j, i_var(1)) = poly(1)
                            tree%boxes(id)%cc(i, j, i_ibfx) = &
                            ((poly(1) - tree%boxes(id)%cc(i, j, i_var(1))) * inv_dt) &
                            + tree%boxes(id)%cc(i, j, i_resx)

                            poly(:) = &
                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(2,:),   &
                                                x_forcing(1,1), x_forcing(2,1))
!                            tree%boxes(id)%cc(i, j, i_var(2)) = poly(1)
                            tree%boxes(id)%cc(i, j, i_ibfy) = &
                            ((poly(1) - tree%boxes(id)%cc(i, j, i_var(2))) * inv_dt) &
                            + tree%boxes(id)%cc(i, j, i_resy)
                        end if
                    end do
                end do
            end if




!            cntr = 0
!            tree%boxes(id)%cc(:, :, i_ibfx) = 0.0_dp
!            tree%boxes(id)%cc(:, :, i_ibfy) = 0.0_dp
!
!            inv_dt = 1.0_dp / dt
!            tags = inside
!            is_added = .false.
!
!            do j = 1, nc
!                do i = 1, nc
!                    if (                                 &
!                         (cc(i, j, 1) < 0.0_dp) .and.    &
!                         (                               &
!                         cc(i-1, j,   1) > 0.0_dp .or.   &
!                         cc(i+1, j,   1) > 0.0_dp .or.   &
!                         cc(i,   j-1, 1) > 0.0_dp .or.   &
!                         cc(i,   j+1, 1) > 0.0_dp .or.   &
!                         cc(i-1, j-1, 1) > 0.0_dp .or.   &
!                         cc(i-1, j+1, 1) > 0.0_dp .or.   &
!                         cc(i+1, j-1, 1) > 0.0_dp .or.   &
!                         cc(i+1, j+1, 1) > 0.0_dp        &
!                         )                               &
!                       )  then
!
!                        tags(i, j) = ghost
!                        is_added = .true.
!
!                    end if
!
!                    if (cc(i, j, 1) > 0.0_dp) then
!                        tags(i, j) = outside
!                    end if
!
!                    if (                                 &
!                         (cc(i, j, 1) > 0.0_dp) .and.    &
!                         (                               &
!                         cc(i-1, j,   1) < 0.0_dp .or.   &
!                         cc(i+1, j,   1) < 0.0_dp .or.   &
!                         cc(i,   j-1, 1) < 0.0_dp .or.   &
!                         cc(i,   j+1, 1) < 0.0_dp .or.   &
!                         cc(i-1, j-1, 1) < 0.0_dp .or.   &
!                         cc(i-1, j+1, 1) < 0.0_dp .or.   &
!                         cc(i+1, j-1, 1) < 0.0_dp .or.   &
!                         cc(i+1, j+1, 1) < 0.0_dp        &
!                         )                               &
!                       )  then
!
!                        tags(i, j) = forcing
!                    end if
!                end do
!            end do
!
!            do j = 1, nc
!                do i = 1, nc
!                    if (tags(i,j) == inside) then
!                        tree%boxes(id)%cc(i, j, i_ibfx) = &
!                        ((this%velocity(1) - tree%boxes(id)%cc(i, j, i_var(1)) )*inv_dt) &
!                        - tree%boxes(id)%cc(i, j, i_resx)
!!                        tree%boxes(id)%cc(i, j, i_unew) = this%velocity(1)
!
!                        tree%boxes(id)%cc(i, j, i_ibfy) = &
!                        ((this%velocity(2) - tree%boxes(id)%cc(i, j, i_var(2)) )*inv_dt) &
!                        - tree%boxes(id)%cc(i, j, i_resy)
!!                        tree%boxes(id)%cc(i, j, i_vnew) = this%velocity(2)
!
!!                        tree%boxes(id)%cc(i, j, i_u1tmp) = tree%boxes(id)%cc(i, j, i_unew)
!!                        tree%boxes(id)%cc(i, j, i_u2tmp) = tree%boxes(id)%cc(i, j, i_vnew)
!                    end if
!                end do
!            end do
!!            print*, "ETRAP"
!!            is_added = .false.
!            if (is_added) then
!                inv_2dr = 1.0_dp / (2.0_dp * tree%boxes(id)%dr)
!                do j = 1, nc
!                    do i = 1, nc
!                        if(tags(i, j) == ghost) then
!                            cntr = cntr + 1
!!
!
!
!                            nnodes = 0
!                            ntemp  = 0
!                            vnodes = 0.d0
!                            cnodes = 0.d0
!                            are_all_available = .true.
!
!                            !> The closest point on the physical boundary
!                            normal(1) = (cc(i+1, j, 1) - cc(i-1, j, 1)) * inv_2dr(1)
!                            normal(2) = (cc(i, j+1, 1) - cc(i, j-1, 1)) * inv_2dr(2)
!                            normal(:) = normal(:) / (norm2(normal)*1.0_dp)
!                            x_forcing(:, 1) = af_r_cc (tree%boxes(id), [IJK])
!                            x_gama(:) = x_forcing(:, 1) - normal*cc(i, j, 1)
!
!!                            print*, norm2(x_gama - [acos(-1.d0), acos(-1.d0)]) - 0.5d0
!
!                            nnodes = nnodes + 1  ! nnodes = 1
!                            cnodes(:, nnodes) = x_gama(:)
!                            vnodes(:, nnodes) = this%velocity(:)
!
!                            !> A vector indicating the direction for searching the
!                            !  other interpolating points
!                            ehat(:) = ceiling( normal(:) / abs(normal(:)))
!
!                            !> The second point on the physical boundary
!                            next_point(:) = [i,j] + [ehat(1), 0]
!                            x_point(:, 1) =  &
!                            af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!
!                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
!                                         cc(next_point(1)-1, next_point(2), 1))  &
!                                         * inv_2dr(1)
!                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
!                                         cc(next_point(1), next_point(2)-1, 1))  &
!                                         * inv_2dr(2)
!                            normal(:) = normal(:) / (norm2(normal))
!                            x_gama(:) = &
!                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)
!
!                            nnodes = nnodes + 1 ! nnodes = 2
!                            cnodes(:, nnodes) = x_gama(:)
!                            vnodes(:, nnodes) = this%velocity(:)
!
!                            !> The fluid point corresponding to the second point on
!                            !  the physical boundary
!!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
!                                 if (are_all_available) then
!                                    cnodes(:, nnmax-1) =  x_point(:, 1)
!                                    vnodes(1, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                    are_all_available = .false.
!                                    ntemp = ntemp + 1
!                                 else
!                                    cnodes(:, nnmax)  =  x_point(:, 1)
!                                    vnodes(1, nnmax) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                    are_all_available = .false.
!                                    ntemp = ntemp + 1
!                                 end if
!
!!                            end if
!
!                            !> The third point on the physical boundary
!                            next_point(:) = [i,j] + [0, ehat(2)]
!                            x_point(:, 1) =  &
!                            af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
!                                         cc(next_point(1)-1, next_point(2), 1))  &
!                                         * inv_2dr(1)
!                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
!                                         cc(next_point(1), next_point(2)-1, 1))  &
!                                         * inv_2dr(2)
!                            normal(:) = normal(:) / norm2(normal)
!                            x_gama(:) = &
!                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)
!
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes) = x_gama
!                            vnodes(:, nnodes) = this%velocity(:)
!
!                            !> The fluid point corresponding to the third point on
!                            !  on the physical boundary
!!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
!                                 if (are_all_available) then
!                                    cnodes(:, nnmax-1) = x_point(:, 1)!x_gama(:)
!                                    vnodes(1, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                    are_all_available = .false.
!                                    ntemp = ntemp + 1
!                                 else
!                                    cnodes(:, nnmax) = x_point(:, 1) !x_gama(:)
!                                    vnodes(1, nnmax) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                    are_all_available = .false.
!                                    ntemp = ntemp + 1
!                                 end if
!
!!                            end if
!
!                            !> The first fluid point that can not be a ghost point
!                            next_point(:) = [i,j] + [ehat(1), ehat(2)]
!                            x_point(:, 1) =   &
!                            af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes)  =  x_point(:, 1)
!                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
!                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)
!
!                            !> The second fluid point that can not be a ghost point
!                            next_point(:) = [i,j] + [2*ehat(1), 0]
!                            x_point(:, 1) =   &
!                            af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes)  =  x_point(:, 1)
!                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
!                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)
!
!                            !> The third fluid point that can not be a ghost point
!                            next_point(:) = [i,j] + [0, 2*ehat(2)]
!                            x_point(:, 1) =   &
!                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes)  =  x_point(:, 1)
!
!                            !> The interpolation velocity are chosen from the 2 and 3
!                            !  component of the "cc" array.
!                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
!                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)
!!
!
!
!                            block
!                            real(kind=dp) :: utemp, vtemp
!
!!                            poly(:) = &
!!                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(1, :),      &
!!                                                cnodes(1, 1), cnodes(2, 1))
!!                            utemp = poly(1) + &
!!                                    poly(2) * ( x_forcing(1, 1) - cnodes(1, 1) ) + &
!!                                    poly(3) * ( x_forcing(2, 1) - cnodes(2, 1) )
!
!                            poly(:) = &
!                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(1, :),      &
!                                                x_forcing(1, 1), x_forcing(2, 1))
!
!!                            print*, poly(1), utemp
!!                            tree%boxes(id)%cc(i, j, i_var(1)) = poly(1)
!
!!                            tree%boxes(id)%cc(i, j, i_ibfx) = &
!!                            ((utemp - tree%boxes(id)%cc(i, j, i_var(1))) * inv_dt) &
!!                            - tree%boxes(id)%cc(i, j, i_resx)
!
!                            tree%boxes(id)%cc(i, j, i_ibfx) = &
!                            ((poly(1) - tree%boxes(id)%cc(i, j, i_var(1))) * inv_dt) &
!                            - tree%boxes(id)%cc(i, j, i_resx)
!
!
!
!!                            poly(:) = &
!!                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(2, :),   &
!!                                                cnodes(1, 1), cnodes(2, 1))
!
!!                            vtemp = poly(1) + &
!!                                    poly(2) * ( x_forcing(1, 1) - cnodes(1, 1) ) + &
!!                                    poly(3) * ( x_forcing(2, 1) - cnodes(2, 1) )
!
!                            poly(:) = &
!                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(2,:),   &
!                                                x_forcing(1,1), x_forcing(2,1))
!!                            tree%boxes(id)%cc(i, j, i_var(2)) = poly(1)
!!                            tree%boxes(id)%cc(i, j, i_ibfy) = &
!!                            ((vtemp - tree%boxes(id)%cc(i, j, i_var(2))) * inv_dt) &
!!                            - tree%boxes(id)%cc(i, j, i_resy)
!                            tree%boxes(id)%cc(i, j, i_ibfy) = &
!                            ((poly(1) - tree%boxes(id)%cc(i, j, i_var(2))) * inv_dt) &
!                            - tree%boxes(id)%cc(i, j, i_resy)
!                            end block
!
!                        end if
!                    end do
!                end do
!            end if
!
        end subroutine  get_object_hybrid_force_ker
    end subroutine get_object_force




!    subroutine get_interpolated_velocity (this, fld, i_var)
!        class(immersed_boundary_t), intent(inout) :: this
!        class(incomp_flow_field_t), intent(inout) :: fld
!        integer, intent(in) :: i_var(:)
!        integer :: n
!
!        do n = 1, this%nobjs
!            call this%objects(n)%get_object_interpolated_velocity(fld, i_var)
!        end do
!    end subroutine get_interpolated_velocity



!    subroutine get_object_interpolated_velocity (this, fld, i_var)
!        class(object_t), intent(inout) :: this
!        class(incomp_flow_field_t), intent(inout) :: fld
!        integer, intent(in) :: i_var(:)
!
!!        call af_loop_tree(fld%tree, get_interpolated_ker, leaves_only=.true.)
!        call af_loop_tree(fld%tree, get_interpolated_ker) !!!!!!!!!!!!!!!!!!!!!
!        call af_gc_tree(fld%tree, [i_p, i_var(1), i_var(2)])
!        call af_restrict_tree(fld%tree, i_p)
!        call af_restrict_tree(fld%tree, i_var(1))
!        call af_restrict_tree(fld%tree, i_var(2))


!    contains
!
!
!        subroutine  get_interpolated_ker(tree, id)
!            implicit none
!            type(af_t), intent(inout) :: tree
!            integer, intent(in) :: id
!            integer :: nc, i, j
!            logical :: is_added
!            real(kind=dp), allocatable :: cc (DTIMES(:), :)
!            integer, allocatable :: tags(DTIMES(:))
!            real(kind=dp) :: x_forcing(NDIM, 1), x_gama(NDIM)
!            real(kind=dp) :: x_point(NDIM, 1)
!            real(kind=dp) :: inv_2dr(NDIM), normal(NDIM)
!            ! MLS
!            real(kind=dp) :: cnodes(NDIM, nstencil), poly(ns)
!            real(kind=dp) :: vnodes(size(i_var), nstencil)
!            integer :: nnodes, ntemp
!            integer :: next_point(NDIM), ehat(NDIM)
!            logical :: are_all_available
!
!            nc = tree%boxes(id)%n_cell
!            allocate(cc(DTIMES(-1:nc+2), size(i_var)+1 ))
!            allocate(tags(DTIMES(-1:nc+2)))
!            call af_gc2_box(tree, id, [i_df, i_var(1), i_var(2)], cc)
!
!            tags = outside
!            is_added = .false.
!            do j = 1, nc
!                do i = 1, nc
!                    if (                                 &
!                         (cc(i, j, 1) > 0.0_dp) .and.    &
!                         (                               &
!                         cc(i-1, j,   1) < 0.0_dp .or.   &
!                         cc(i+1, j,   1) < 0.0_dp .or.   &
!                         cc(i,   j-1, 1) < 0.0_dp .or.   &
!                         cc(i,   j+1, 1) < 0.0_dp .or.   &
!                         cc(i-1, j-1, 1) < 0.0_dp .or.   &
!                         cc(i-1, j+1, 1) < 0.0_dp .or.   &
!                         cc(i+1, j-1, 1) < 0.0_dp .or.   &
!                         cc(i+1, j+1, 1) < 0.0_dp        &
!                         )                               &
!                       )  then
!
!                        tags(i, j) = forcing
!                        is_added = .true.
!                    end if
!
!                    if (cc(i, j, 1) < 0.0_dp) then
!                        tags(i, j) = inside
!                    end if
!                end do
!            end do
!
!            do j = 1, nc
!                do i = 1, nc
!                    if (tags(i,j) == inside) then
!                        tree%boxes(id)%cc(i, j, i_unew) = 0.0_dp
!                        tree%boxes(id)%cc(i, j, i_vnew) = 0.0_dp
!                    end if
!                end do
!            end do
!
!!            is_added = .false.
!
!            if (is_added) then
!                inv_2dr = 1.0_dp / (2.0_dp * tree%boxes(id)%dr)
!                do j = 1, nc
!                    do i = 1, nc
!                        if(tags(i,j) == forcing) then
!                            nnodes = 0
!                            ntemp = 0
!                            vnodes = 0.d0
!                            cnodes = 0.d0
!                            are_all_available = .true.
!
!                            !> The closest point on the physical boundary
!                            normal(1) = (cc(i+1, j, 1) - cc(i-1, j, 1)) * inv_2dr(1)
!                            normal(2) = (cc(i, j+1, 1) - cc(i, j-1, 1)) * inv_2dr(2)
!                            normal(:) = normal(:) / (norm2(normal)*1.0_dp)
!                            x_forcing(:, 1) = af_r_cc(tree%boxes(id), [IJK])
!                            x_gama(:) = x_forcing(:, 1) - normal*cc(i, j, 1)
!                            nnodes = nnodes + 1  ! nnodes = 1
!                            cnodes(:, nnodes) = x_gama(:)
!                            vnodes(:, nnodes) = this%velocity(:)
!
!                            !> A vector indicating the direction for searching the
!                            !  other interpolating points
!                            ehat(:) = ceiling( normal(:) / abs(normal(:)))
!
!                            !> The second point on the physical boundary
!                            next_point(:) = [i,j] + [ehat(1), 0]
!                            x_point(:, 1) =  &
!                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
!                                         cc(next_point(1)-1, next_point(2), 1))  &
!                                         * inv_2dr(1)
!                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
!                                         cc(next_point(1), next_point(2)-1, 1))  &
!                                         * inv_2dr(2)
!                            normal(:) = normal(:) / (norm2(normal))
!                            x_gama(:) = &
!                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)
!
!                            nnodes = nnodes + 1 ! nnodes = 2
!                            cnodes(:, nnodes) = x_gama(:)
!                            vnodes(:, nnodes) = this%velocity(:)
!
!                            !> The fluid point corresponding to the second point on
!                            !  the physical boundary
!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
!                                 if (are_all_available) then
!                                    cnodes(:, nnmax-1) =  x_point(:, 1)
!                                    vnodes(1, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                 else
!                                    cnodes(:, nnmax)  =  x_point(:, 1)
!                                    vnodes(1, nnmax) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                 end if
!                                 are_all_available = .false.
!                                 ntemp = ntemp + 1
!                            end if
!
!                            !> The third point on the physical boundary
!                            next_point(:) = [i,j] + [0, ehat(2)]
!                            x_point(:, 1) =  &
!                            af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
!                                         cc(next_point(1)-1, next_point(2), 1))  &
!                                         * inv_2dr(1)
!                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
!                                         cc(next_point(1), next_point(2)-1, 1))  &
!                                         * inv_2dr(2)
!                            normal(:) = normal(:) / norm2(normal)
!                            x_gama(:) = &
!                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)
!
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes) = x_gama
!                            vnodes(:, nnodes) = this%velocity(:)
!
!                            !> The fluid point corresponding to the third point on
!                            !  on the physical boundary
!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
!                                 if (are_all_available) then
!                                    cnodes(:, nnmax-1) = x_gama(:)
!                                    vnodes(1, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax-1) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                 else
!                                    cnodes(:, nnmax) = x_gama(:)
!                                    vnodes(1, nnmax) = &
!                                        cc(next_point(1), next_point(2), 2)
!                                    vnodes(2, nnmax) = &
!                                        cc(next_point(1), next_point(2), 3)
!                                 end if
!                                 are_all_available = .false.
!                                 ntemp = ntemp + 1
!                            end if
!
!                            !> The first fluid point that can not be a ghost point
!                            next_point(:) = [i,j] + [ehat(1), ehat(2)]
!                            x_point(:, 1) =   &
!                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes)  =  x_point(:, 1)
!                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
!                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)
!
!                            !> The second fluid point that can not be a ghost point
!                            next_point(:) = [i,j] + [2*ehat(1), 0]
!                            x_point(:, 1) =   &
!                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes)  =  x_point(:, 1)
!                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
!                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)
!
!                            !> The third fluid point that can not be a ghost point
!                            next_point(:) = [i,j] + [0, 2*ehat(2)]
!                            x_point(:, 1) =   &
!                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
!                            nnodes = nnodes + 1
!                            cnodes(:, nnodes)  =  x_point(:, 1)
!
!                            !> The interpolation velocity are chosen from the 2 and 3
!                            !  component of the "cc" array.
!                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
!                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)
!
!                            poly(:) = &
!                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(1,:),        &
!                                                x_forcing(1,1), x_forcing(2,1))
!                            tree%boxes(id)%cc(i, j, i_var(1)) = poly(1)
!
!                            poly(:) = &
!                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(2,:),        &
!                                                x_forcing(1,1), x_forcing(2,1))
!                            tree%boxes(id)%cc(i, j, i_var(2)) = poly(1)
!                        end if
!                    end do
!                end do
!            end if
!        end subroutine  get_interpolated_ker
!    end subroutine get_object_interpolated_velocity
!
!
!
!
    subroutine get_extrapolated_velocity (this, tree, i_var)
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t), intent(inout) :: tree
        integer, intent(in) :: i_var(:)
        integer :: n

        do n = 1, this%nobjs
            call this%objects(n)%get_object_extrapolated_velocity(tree, i_var)
        end do
    end subroutine get_extrapolated_velocity



    subroutine get_object_extrapolated_velocity (this, tree, i_var)
        class(object_t), intent(inout) :: this
        class(af_t), intent(inout) :: tree
        integer, intent(in) :: i_var(:)

!        call af_loop_tree(fld%tree, get_interpolated_ker, leaves_only=.true.)
        call af_loop_tree(tree, get_extrapolated_ker) !!!!!!!!!!!!!!!!!!!!!
!        call af_gc_tree(tree, [i_p, i_var(1), i_var(2)])
!        call af_restrict_tree(tree, i_p)
!        call af_restrict_tree(tree, i_var(1))
!        call af_restrict_tree(tree, i_var(2))


    contains


        subroutine  get_extrapolated_ker(tree, id)
            implicit none
            type(af_t), intent(inout) :: tree
            integer, intent(in) :: id
            integer :: nc, i, j
            logical :: is_added
            real(kind=dp), allocatable :: cc (DTIMES(:), :)
            integer, allocatable :: tags(DTIMES(:))
            real(kind=dp) :: x_forcing(NDIM, 1), x_gama(NDIM)
            real(kind=dp) :: x_point(NDIM, 1)
            real(kind=dp) :: inv_2dr(NDIM), normal(NDIM)
            ! MLS
            real(kind=dp) :: cnodes(NDIM, nstencil), poly(ns)
            real(kind=dp) :: vnodes(size(i_var), nstencil)
            integer :: nnodes, ntemp
            integer :: next_point(NDIM), ehat(NDIM)
            logical :: are_all_available

            nc = tree%boxes(id)%n_cell
            allocate(cc(DTIMES(-1:nc+2), size(i_var)+1 ))
            allocate(tags(DTIMES(-1:nc+2)))
            call af_gc2_box(tree, id, [i_df, i_var(1), i_var(2)], cc)

            tags = inside
            is_added = .false.
            do j = 1, nc
                do i = 1, nc
                    if (                                 &
                         (cc(i, j, 1) < 0.0_dp) .and.    &
                         (                               &
                         cc(i-1, j,   1) > 0.0_dp .or.   &
                         cc(i+1, j,   1) > 0.0_dp .or.   &
                         cc(i,   j-1, 1) > 0.0_dp .or.   &
                         cc(i,   j+1, 1) > 0.0_dp .or.   &
                         cc(i-1, j-1, 1) > 0.0_dp .or.   &
                         cc(i-1, j+1, 1) > 0.0_dp .or.   &
                         cc(i+1, j-1, 1) > 0.0_dp .or.   &
                         cc(i+1, j+1, 1) > 0.0_dp        &
                         )                               &
                       )  then

                        tags(i, j) = ghost
                        is_added = .true.
                    end if

                    if (cc(i, j, 1) > 0.0_dp) then
                        tags(i, j) = outside
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
                         )                               &
                       )  then

                        tags(i, j) = forcing
                    end if
                end do
            end do

            do j = 1, nc
                do i = 1, nc
                    if (tags(i,j) == inside) then
                        tree%boxes(id)%cc(i, j, i_var(1)) = 0.0_dp
                        tree%boxes(id)%cc(i, j, i_var(2)) = 0.0_dp
                    end if
                end do
            end do
!            print*, "ETRAP"
!            is_added = .false.
            if (is_added) then
                inv_2dr = 1.0_dp / (2.0_dp * tree%boxes(id)%dr)
                do j = 1, nc
                    do i = 1, nc
                        if(tags(i, j) == ghost) then
                            nnodes = 0
                            ntemp = 0
                            vnodes = 0.d0
                            cnodes = 0.d0
                            are_all_available = .true.

                            !> The closest point on the physical boundary
                            normal(1) = (cc(i+1, j, 1) - cc(i-1, j, 1)) * inv_2dr(1)
                            normal(2) = (cc(i, j+1, 1) - cc(i, j-1, 1)) * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal)*1.0_dp)
                            x_forcing(:, 1) = af_r_cc(tree%boxes(id), [IJK])
                            x_gama(:) = x_forcing(:, 1) - normal*cc(i, j, 1)

!                            print*, norm2(x_gama - [acos(-1.d0), acos(-1.d0)]) -0.5d0

                            nnodes = nnodes + 1  ! nnodes = 1
                            cnodes(:, nnodes) = x_gama(:)
                            vnodes(:, nnodes) = this%velocity(:)

                            !> A vector indicating the direction for searching the
                            !  other interpolating points
                            ehat(:) = ceiling( normal(:) / abs(normal(:)))

                            !> The second point on the physical boundary
                            next_point(:) = [i,j] + [ehat(1), 0]
                            x_point(:, 1) =  &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
                                         cc(next_point(1)-1, next_point(2), 1))  &
                                         * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
                                         cc(next_point(1), next_point(2)-1, 1))  &
                                         * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal))
                            x_gama(:) = &
                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)

                            nnodes = nnodes + 1 ! nnodes = 2
                            cnodes(:, nnodes) = x_gama(:)
                            vnodes(:, nnodes) = this%velocity(:)

                            !> The fluid point corresponding to the second point on
                            !  the physical boundary
!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 if (are_all_available) then
                                    cnodes(:, nnmax-1) =  x_point(:, 1)
                                    vnodes(1, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 3)
                                 else
                                    cnodes(:, nnmax)  =  x_point(:, 1)
                                    vnodes(1, nnmax) = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax) = &
                                        cc(next_point(1), next_point(2), 3)
                                 end if
                                 are_all_available = .false.
                                 ntemp = ntemp + 1
!                            end if

                            !> The third point on the physical boundary
                            next_point(:) = [i,j] + [0, ehat(2)]
                            x_point(:, 1) =  &
                            af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) - &
                                         cc(next_point(1)-1, next_point(2), 1))  &
                                         * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) - &
                                         cc(next_point(1), next_point(2)-1, 1))  &
                                         * inv_2dr(2)
                            normal(:) = normal(:) / norm2(normal)
                            x_gama(:) = &
                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)

                            nnodes = nnodes + 1
                            cnodes(:, nnodes) = x_gama
                            vnodes(:, nnodes) = this%velocity(:)

                            !> The fluid point corresponding to the third point on
                            !  on the physical boundary
!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 if (are_all_available) then
                                    cnodes(:, nnmax-1) = x_gama(:)
                                    vnodes(1, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax-1) = &
                                        cc(next_point(1), next_point(2), 3)
                                 else
                                    cnodes(:, nnmax) = x_gama(:)
                                    vnodes(1, nnmax) = &
                                        cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax) = &
                                        cc(next_point(1), next_point(2), 3)
                                 end if
                                 are_all_available = .false.
                                 ntemp = ntemp + 1
!                            end if

                            !> The first fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [ehat(1), ehat(2)]
                            x_point(:, 1) =   &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes)  =  x_point(:, 1)
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            !> The second fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [2*ehat(1), 0]
                            x_point(:, 1) =   &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes)  =  x_point(:, 1)
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            !> The third fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [0, 2*ehat(2)]
                            x_point(:, 1) =   &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes)  =  x_point(:, 1)

                            !> The interpolation velocity are chosen from the 2 and 3
                            !  component of the "cc" array.
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            poly(:) = &
                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(1,:),        &
                                                x_forcing(1,1), x_forcing(2,1))
                            tree%boxes(id)%cc(i, j, i_var(1)) = poly(1)

                            poly(:) = &
                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(2,:),        &
                                                x_forcing(1,1), x_forcing(2,1))
                            tree%boxes(id)%cc(i, j, i_var(2)) = poly(1)

                        end if
                    end do
                end do
            end if
        end subroutine  get_extrapolated_ker
    end subroutine get_object_extrapolated_velocity
!
!
    subroutine read_markers(this)
        class(object_t), intent(inout) :: this
        integer :: n
        integer :: myunit

        open(unit=myunit, file="./inputs/object_1.dat")
        read(myunit, *) this%nmarkers
        allocate (this%marker_list(this%nmarkers))
        do n = 1, this%nmarkers
            read(myunit ,*) this%marker_list(n)%loc(:)
        end do
        close(myunit)
    end subroutine read_markers



!    subroutine move_objects(this, fld)
!        class(immersed_boundary_t), intent(inout) :: this
!        class(incomp_flow_field_t), intent(inout) :: fld
!        integer :: n
!
!        do n = 1, this%nobjs
!            call this%objects(n)%move_object()
!        end do
!
!        do n = 1, this%nobjs
!            call build_object (this%objects(n), fld, is_initiated=.false.,          &
!                               is_refined=.false.)
!        end do
!    end subroutine move_objects



!    subroutine move_object(this)
!        class(object_t), intent(inout) :: this
!        integer :: n
!
!        !> All Lagrangian markers are moved according to predefined components of the
!        !  rigid object velocities
!        forall (n = 1: this%nmarkers)
!            this%marker_list(n)%loc(:) =  &
!            (predefined_velocity(field_time) * dt ) + this%marker_list(n)%loc(:)
!        end forall
!
!        !> It sets the velocity of the object which will be used in the boundary in
!        !  imposing condition.
!        this%velocity(:) = predefined_velocity(field_time)
!    end subroutine move_object



    pure function predefined_velocity (time) result(vel)
        real(kind=dp), intent(in) :: time
        real(kind=dp) :: vel(NDIM)
        vel = 0.0d0 * time
    end function predefined_velocity


end module m_immersed_boundary
