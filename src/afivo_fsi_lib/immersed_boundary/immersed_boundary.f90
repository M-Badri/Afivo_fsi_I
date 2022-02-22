#include "cpp_macros.h"
module m_immersed_boundary
    use omp_lib
    use m_solver_parameters
    use m_config
    use m_af_all
    use m_solver_parameters
    use m_incomp_flow_field
    use m_rbf_iterp_2d
    use m_mls_interpolation
    use m_write
    use m_predefined_move

    implicit none
    !> These variables are used as temporary placeholder for building each object
    type(t_config)              :: config_glb
    real(kind=R8)               :: landa_g = 1.0_R8
    integer(kind=I4), parameter :: nstencil = 8
    integer(kind=I4), parameter :: nnmax = 8
    integer(kind=I4), parameter :: n_time_ave = 5, n_smoothed_var = 3
    real(kind=R8), allocatable  :: distance_matrix_glb(:, :)
    real(kind=R8), allocatable  :: force_matrix_glb(:, :, :)
    integer(kind=I4)            :: nwalls
    real(kind=R8)               :: data_holder_1(n_smoothed_var, n_time_ave)
    real(kind=R8)               :: data_holder_2(n_smoothed_var, n_time_ave-1)
!    real(kind=R8) :: inv_epsi_1, inv_epsi_2, d_0 , xi, dx

    type marker_t
        real(kind=R8) :: loc(NDIM)
    end type marker_t


    type object_t
        integer(kind=I4) :: id       ! object id
        integer(kind=I4) :: nmarkers
        real(kind=R8)    :: xc(NDIM)
        real(kind=R8)    :: xc_d(NDIM)
        real(kind=R8)    :: xc_dd(NDIM)
        real(kind=R8)    :: force(NDIM)
        real(kind=R8)    :: m        ! mass
        real(kind=R8)    :: im       ! mass moment of inertia
        real(kind=R8)    :: ext_force(NDIM)
        real(kind=R8),  allocatable :: ws(:) ! Weights fro RBF representation
        type(marker_t), allocatable :: marker_list(:)
        logical                     :: predefined_move
#if NDIM==2
        real(kind=R8) :: om(3), om_d(3), om_dd(3)  ! omega, omega_dot, omega_dotdot
        real(kind=R8) :: moment
        real(kind=R8) :: et, er, rey
#endif

    contains
        procedure, private :: read_markers
        procedure, public  :: get_object_ext_force
!        procedure, public :: get_object_interpolated_velocity
        procedure, public  :: get_object_extrapolated_velocity
        procedure, public  :: update_an_object_forces
        procedure, public  :: update_an_object_acceleration
        procedure, public  :: update_an_object_velocity
        procedure, public  :: update_an_object_location
        procedure, public  :: write_out_ext_force
        procedure, public  :: write_out_object_energy
        procedure, public  :: write_out_object_reynolds
        procedure, public  :: write_out_location_velocity
        procedure, public  :: get_predefined_forces
        procedure, public  :: get_calculated_forces
        procedure, public  :: get_predefined_acceleration
        procedure, public  :: get_calculated_acceleration
        procedure, public  :: get_predefined_velocity
        procedure, public  :: get_calculated_velocity
        procedure, public  :: get_predefined_location
        procedure, public  :: get_calculated_location
        procedure, public  :: body_force
!        procedure, public :: collision_force
        procedure, public  :: update_lagrangian_markers_location
        procedure, public  :: rebuild_object
        procedure, public  :: update_an_object_er
        procedure, public  :: update_an_object_et
        procedure, public  :: update_an_object_rey
    end type object_t


    type immersed_boundary_t
        type(object_t), allocatable :: objects(:)
        integer(kind=I4)            :: nobjs
        real(kind=R8)               :: domain_len(NDIM)
    contains
        procedure, public :: get_ext_forces
!        procedure, public :: get_interpolated_velocity
        procedure, public  :: get_extrapolated_velocity
        procedure, public  :: write_out_ext_forces
        procedure, public  :: write_out_objects_reynolds
        procedure, public  :: write_out_objects_energy
        procedure, public  :: write_out_locations_velocities
        procedure, public  :: update_objects
        procedure, public  :: get_collision_forces
        procedure, private :: get_collision_force_wall_1
        procedure, private :: get_collision_force_wall_2
        procedure, private :: get_collision_force_wall_3
        procedure, private :: get_collision_force_wall_4
    end type immersed_boundary_t



    interface

        module subroutine  build_object (obj, i_obj, tree, is_initiated)
            class(object_t),  intent(inout)          :: obj
            class(af_t),      intent(inout)          :: tree
            integer(kind=I4), intent(in)             :: i_obj
            logical,          intent(in),   optional :: is_initiated
        end subroutine build_object


        module subroutine get_ext_forces (this, tree)
            class(immersed_boundary_t), intent(inout) :: this
            class(af_t),                intent(inout) :: tree
        end subroutine get_ext_forces


        module subroutine get_object_ext_force (this, tree, nobjects, i_var)
            class(object_t),  intent(inout) :: this
            class(af_t),      intent(inout) :: tree
            integer(kind=I4), intent(in)    :: i_var(:)
            integer(kind=I4), intent(in)    :: nobjects
        end subroutine


        module subroutine write_out_ext_forces( this )
            class(immersed_boundary_t), intent(inout) :: this
        end subroutine write_out_ext_forces


        module subroutine write_out_ext_force( this, i_obj )
            class(object_t),  intent(inout) :: this
            integer(kind=I4), intent(in)    :: i_obj
        end subroutine write_out_ext_force


        module subroutine write_out_objects_energy( this )
            class(immersed_boundary_t), intent(inout) :: this
        end subroutine write_out_objects_energy


        module subroutine write_out_object_energy( this, i_obj )
            class(object_t),  intent(inout) :: this
            integer(kind=I4), intent(in)    :: i_obj
        end subroutine write_out_object_energy


        module subroutine write_out_objects_reynolds( this )
            class(immersed_boundary_t), intent(inout) :: this
        end subroutine write_out_objects_reynolds


        module subroutine write_out_object_reynolds( this, i_obj )
            class(object_t),  intent(inout) :: this
            integer(kind=I4), intent(in)    :: i_obj
        end subroutine write_out_object_reynolds

    end interface


contains


    subroutine build_immersed_boundary( tree, refine_info, config, ib )
        implicit none
        class(t_config),            intent(in)     :: config
        class(af_t),                intent(inout)  :: tree
        class(ref_info_t),          intent(inout)  :: refine_info
        class(immersed_boundary_t), intent(out)    :: ib
        integer(kind=I4)                           :: n

        config_glb = config

        !> How many objects are involved in this problem?
        call config_glb%get( section_name='immersed_boundary', option_name='nobjects',    &
        &   val=ib%nobjs )

        !> Allocates an array with length equal to number of objects
        allocate( ib%objects(ib%nobjs) )

        call config%get( section_name='field', option_name='domain_len', val=domain_len )

        !> Builds all objects involved in this problem
        do n = 1, ib%nobjs
            call build_object( ib%objects(n), n, tree )
        end do

#if NDIM == 2
        nwalls = 4
#else
        nwalls = 6
#endif
        !> wall 1 => x_low; wall 2 => x_high; wall 3 => y_low; wall 4 => y_high
        allocate( distance_matrix_glb( ib%nobjs+nwalls, ib%nobjs+nwalls) )
        distance_matrix_glb = 1000.d0

        allocate( force_matrix_glb(NDIM, ib%nobjs+nwalls, ib%nobjs+nwalls) )

        call build_predefined_move( )

        data_holder_1 = 0.0_R8
        data_holder_2 = 0.0_R8
    end subroutine build_immersed_boundary



    subroutine write_out_locations_velocities( this )
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        integer(kind=I4)                          :: n

        do n = 1, this%nobjs
            call this%objects(n)%write_out_location_velocity( n )
        end do
    end subroutine write_out_locations_velocities



    subroutine write_out_location_velocity( this, i_obj )
        implicit none
        class(object_t),  intent(inout) :: this
        integer(kind=I4), intent(in)    :: i_obj

        call write_out_location_velocity_on_hard( i_obj, field_time,                       &
        &                                         this%xc, this%xc_d, this%om, this%om_d )
    end subroutine write_out_location_velocity



    subroutine get_extrapolated_velocity( this, tree, i_var )
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t),                intent(inout) :: tree
        integer(kind=I4),           intent(in)    :: i_var(:)
        integer(kind=I4)                          :: n

        do n = 1, this%nobjs
            call this%objects(n)%get_object_extrapolated_velocity( tree, i_var )
        end do
    end subroutine get_extrapolated_velocity



    subroutine get_object_extrapolated_velocity( this, tree, i_var )
        implicit none
        class(object_t),  intent(inout) :: this
        class(af_t),      intent(inout) :: tree
        integer(kind=I4), intent(in)    :: i_var(:)

        call af_loop_tree( tree, get_extrapolated_ker )


    contains


        subroutine  get_extrapolated_ker(tree, id)
            implicit none
            type(af_t),       intent(inout) :: tree
            integer(kind=I4), intent(in)    :: id
            integer(kind=I4)                :: nc, i, j
            logical                         :: is_added
            real(kind=R8),    allocatable   :: cc(DTIMES(:), :)
            integer(kind=I4), allocatable   :: tags(DTIMES(:))
            real(kind=R8)                   :: x_forcing(NDIM, 1), x_gama(NDIM)
            real(kind=R8)                   :: x_point(NDIM, 1)
            real(kind=R8)                   :: inv_2dr(NDIM), normal(NDIM)
            ! MLS
            real(kind=R8)    :: cnodes(NDIM, nstencil), poly(ns)
            real(kind=R8)    :: vnodes(size(i_var), nstencil)
            integer(kind=I4) :: nnodes, ntemp
            integer(kind=I4) :: next_point(NDIM), ehat(NDIM)
            logical          :: are_all_available

            nc = tree%boxes(id)%n_cell
            allocate( cc(DTIMES(-1:nc+2), size(i_var)+1 ) )
            allocate( tags(DTIMES(-1:nc+2)) )
            call af_gc2_box( tree, id, [i_df+i_var(3)-1, i_var(1), i_var(2)], cc )

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
                        tree%boxes(id)%cc(i, j, i_var(1)) = 0.0_R8
                        tree%boxes(id)%cc(i, j, i_var(2)) = 0.0_R8
                    end if
                end do
            end do

            if (is_added) then
                inv_2dr = 1.0_R8 / (2.0_R8 * tree%boxes(id)%dr)
                do j = 1, nc
                    do i = 1, nc
                        if(tags(i, j) == ghost) then
                            nnodes = 0
                            ntemp = 0
                            vnodes = 0.0_R8
                            cnodes = 0.0_R8
                            are_all_available = .true.

                            !> The closest point on the physical boundary
                            normal(1) = (cc(i+1, j, 1) - cc(i-1, j, 1)) * inv_2dr(1)
                            normal(2) = (cc(i, j+1, 1) - cc(i, j-1, 1)) * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal)*1.0_R8)
                            x_forcing(:, 1) = af_r_cc(tree%boxes(id), [IJK])
                            x_gama(:) = x_forcing(:, 1) - normal*cc(i, j, 1)

!                            print*, norm2(x_gama - [acos(-1.d0), acos(-1.d0)]) -0.5d0

                            nnodes = nnodes + 1  ! nnodes = 1
                            cnodes(:, nnodes) = x_gama(:)
                            vnodes(:, nnodes) = this%xc_d(:)

                            !> A vector indicating the direction for searching the
                            !  other interpolating points
                            ehat(:) = ceiling( normal(:) / abs(normal(:)))

                            !> The second point on the physical boundary
                            next_point(:) = [i,j] + [ehat(1), 0]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) -         &
                                         cc(next_point(1)-1, next_point(2), 1))          &
                                         * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) -         &
                                         cc(next_point(1), next_point(2)-1, 1))          &
                                         * inv_2dr(2)
                            normal(:) = normal(:) / (norm2(normal))
                            x_gama(:) = &
                            x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)

                            nnodes = nnodes + 1 ! nnodes = 2
                            cnodes(:, nnodes) = x_gama(:)
                            vnodes(:, nnodes) = this%xc_d(:)

                            !> The fluid point corresponding to the second point on
                            !  the physical boundary
!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 if (are_all_available) then
                                    cnodes(:, nnmax-1) =  x_point(:, 1)
                                    vnodes(1, nnmax-1) =                                 &
                                    &  cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax-1) =                                 &
                                    &  cc(next_point(1), next_point(2), 3)
                                 else
                                    cnodes(:, nnmax)  =  x_point(:, 1)
                                    vnodes(1, nnmax) =                                   &
                                    &  cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax) =                                   &
                                    &  cc(next_point(1), next_point(2), 3)
                                 end if
                                 are_all_available = .false.
                                 ntemp = ntemp + 1
!                            end if

                            !> The third point on the physical boundary
                            next_point(:) = [i,j] + [0, ehat(2)]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            normal(1) = (cc(next_point(1)+1, next_point(2), 1) -         &
                            &  cc(next_point(1)-1, next_point(2), 1)) * inv_2dr(1)
                            normal(2) = (cc(next_point(1), next_point(2)+1, 1) -         &
                            &  cc(next_point(1), next_point(2)-1, 1)) * inv_2dr(2)
                            normal(:) = normal(:) / norm2(normal)
                            x_gama(:) =                                                  &
                            &  x_point(:, 1) - normal*cc(next_point(1), next_point(2),1)

                            nnodes = nnodes + 1
                            cnodes(:, nnodes) = x_gama
                            vnodes(:, nnodes) = this%xc_d(:)

                            !> The fluid point corresponding to the third point on
                            !  on the physical boundary
!                            if( tags(next_point(1), next_point(2)) .ne. forcing) then
                                 if (are_all_available) then
                                    cnodes(:, nnmax-1) = x_gama(:)
                                    vnodes(1, nnmax-1) =                                 &
                                    &    cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax-1) =                                 &
                                    &    cc(next_point(1), next_point(2), 3)
                                 else
                                    cnodes(:, nnmax) = x_gama(:)
                                    vnodes(1, nnmax) =                                   &
                                    &    cc(next_point(1), next_point(2), 2)
                                    vnodes(2, nnmax) =                                   &
                                    &    cc(next_point(1), next_point(2), 3)
                                 end if
                                 are_all_available = .false.
                                 ntemp = ntemp + 1
!                            end if

                            !> The first fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [ehat(1), ehat(2)]
                            x_point(:, 1) =                                              &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes)  =  x_point(:, 1)
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            !> The second fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [2*ehat(1), 0]
                            x_point(:, 1) =                                              &
                              af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes)  =  x_point(:, 1)
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            !> The third fluid point that can not be a ghost point
                            next_point(:) = [i,j] + [0, 2*ehat(2)]
                            x_point(:, 1) =                                              &
                            &  af_r_cc(tree%boxes(id), [next_point(1), next_point(2)])
                            nnodes = nnodes + 1
                            cnodes(:, nnodes)  =  x_point(:, 1)

                            !> The interpolation velocity are chosen from the 2 and 3
                            !  component of the "cc" array.
                            vnodes(1, nnodes) = cc(next_point(1), next_point(2), 2)
                            vnodes(2, nnodes) = cc(next_point(1), next_point(2), 3)

                            poly(:) =                                                    &
                            &  mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(1,:),    &
                            &  x_forcing(1,1), x_forcing(2,1))
                            tree%boxes(id)%cc(i, j, i_var(1)) = poly(1)

                            poly(:) =                                                    &
                            mls_get_scalar_dir (nnodes+ntemp, cnodes, vnodes(2,:),       &
                            &  x_forcing(1,1), x_forcing(2,1))
                            tree%boxes(id)%cc(i, j, i_var(2)) = poly(1)

                        end if
                    end do
                end do
            end if
        end subroutine  get_extrapolated_ker
    end subroutine get_object_extrapolated_velocity



    subroutine read_markers( this, fname )
        implicit none
        class(object_t),  intent(inout) :: this
        character(len=*), intent(inout) :: fname
        integer(kind=I4)                :: n
        integer(kind=I4)                :: myunit


        fname = trim( adjustl( fname ) )
        open( newunit=myunit, file=fname )
        read( myunit, * ) this%nmarkers
        read( myunit, * ) this%xc(:)
        allocate( this%marker_list(this%nmarkers) )
        do n = 1, this%nmarkers
            read( myunit ,* ) this%marker_list(n)%loc(:)
        end do
        close( myunit )
    end subroutine read_markers



    subroutine update_objects( this, tree )
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t),                intent(inout) :: tree
        integer(kind=I4)                          :: i

        call this%get_collision_forces( tree )

        do i = 1, this%nobjs
            call this%objects(i)%update_an_object_forces( )
            call this%objects(i)%update_an_object_acceleration( )
            call this%objects(i)%update_an_object_velocity( )
            call this%objects(i)%update_an_object_er( )
            call this%objects(i)%update_an_object_et( )
            call this%objects(i)%update_an_object_rey( )
            call this%objects(i)%update_an_object_location( )
            call this%objects(i)%update_lagrangian_markers_location( )
            call this%objects(i)%rebuild_object( tree=tree )
        end do
    end subroutine update_objects



    subroutine update_an_object_er( this )
        class(object_t), intent(inout) :: this
        this%er = 0.5d0 * this%im * (this%om_d(3) ** 2)
    end subroutine update_an_object_er



    subroutine update_an_object_et( this )
        class(object_t), intent(inout) :: this
        this%et = 0.5d0 * this%m * ((this%xc_d(1)**2) + (this%xc_d(2)**2))
    end subroutine update_an_object_et



    subroutine update_an_object_rey( this )
        class(object_t), intent(inout) :: this
        this%rey = ( sqrt( (this%xc_d(1)**2)+(this%xc_d(2)**2) ) * 0.25d0 * 1.25 ) / nu_f
    end subroutine update_an_object_rey



    subroutine get_collision_forces (this, tree)
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t), intent(inout) :: tree
        integer(kind=I4) :: i

!        dx = tree%dr_base(1) * (0.5_R8**(tree%highest_lvl-1))
!        inv_epsi_1 = 1.0_R8 / ((tree%dr_base(1) * (0.5_R8**(tree%highest_lvl-1))) ** 2)
!        inv_epsi_2 = 1.0_R8 / ( tree%dr_base(1) * (0.5_R8**(tree%highest_lvl-1))      )
!        inv_epsi_1 = 1.0_R8 / (dx**2)
!        inv_epsi_2 = 1.0_R8 / (dx)
!        d_0 = 2.5_R8 * tree%dr_base(1) * (0.5_R8**(tree%highest_lvl-1)) ! 2.5*dr_min
!        d_0 = 0.01_R8 * dx ! 2.5*dr_min
!        xi  = 2.0_R8 * dx

        do i = 1, this%nobjs
            force_matrix_glb(:, i, this%nobjs+1) = this%get_collision_force_wall_1 (tree, i)
            force_matrix_glb(:, i, this%nobjs+2) = this%get_collision_force_wall_2 (tree, i)
            force_matrix_glb(:, i, this%nobjs+3) = this%get_collision_force_wall_3 (tree, i)
            force_matrix_glb(:, i, this%nobjs+4) = this%get_collision_force_wall_4 (tree, i)
!            do j = 1, this%nobjs
!                if (i /= j) then
!                    call get_collision_force_ij (this, tree)
!                end if
!            end do
        force_matrix_glb(:, i, i) = 0.0_R8
        end do
    end subroutine get_collision_forces



    subroutine update_an_object_forces (this)
        implicit none
        class(object_t), intent(inout) :: this

        if (this%predefined_move) then
            call this%get_predefined_forces ()
        else
            call this%get_calculated_forces ()
        end if
    end subroutine update_an_object_forces



    subroutine get_calculated_forces( this )
        implicit none
        class(object_t), intent(inout) :: this
        real(kind=R8)                  :: f(NDIM), col_f(NDIM), b_force(NDIM)
        integer(kind=I4)               :: n, cwall

        !> Total body forces
        b_force(:) = this%body_force( )

        !> Total collision force on this object
        col_f(1) = sum( force_matrix_glb(1, this%id, : ) )
        col_f(2) = sum( force_matrix_glb(2, this%id, : ) )

        !> Forces from: collision, fluid and body
        cwall = 0
        do n = 1, nwalls
            if (                                                 &
                force_matrix_glb(1, this%id, n) /= 0.0_R8 .or.   &
                force_matrix_glb(2, this%id, n) /= 0.0_R8        &
               ) then
                cwall = n
            end if
        end do

        if (cwall == 1 .or. cwall == 2) then
            this%force(1) = col_f(1)
            this%force(2) = b_force(2) + this%ext_force(2) + col_f(2)
        else if (cwall == 3 .or. cwall == 4) then
            this%force(1) = b_force(1) + this%ext_force(1) + col_f(1)
            this%force(2) = col_f(2)
!            print*, col_F(2)
!            if (this%xc_d(2) < 0.0_R8) then
!                this%xc_d(2) = 0.0_R8
!            end if
        else
            this%force(:) = b_force(:) + this%ext_force(:) + col_f(:)
        end if

    end subroutine get_calculated_forces



    function get_collision_force_wall_1( this, tree, i ) result( f )
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t),                intent(inout) :: tree
        integer(kind=I4),           intent(in)    :: i
        real(kind=R8)                             :: f(NDIM)
        f = 0.0d0
    end function get_collision_force_wall_1



    function get_collision_force_wall_2( this, tree, i ) result (f)
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t),                intent(inout) :: tree
        integer(kind=I4),           intent(in)    :: i
        real(kind=R8)                             :: f(NDIM)
        f = 0.0d0
    end function get_collision_force_wall_2



    function get_collision_force_wall_3( this, tree, i ) result (f)
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t),                intent(inout) :: tree
        integer(kind=I4),           intent(in)    :: i
        real(kind=R8)                             :: f(NDIM), r_im, r_pa, d_ijk, im_cent(NDIM), epsi
        real(kind=R8)                             :: inv_epsi_1, inv_epsi_2 , xi, dx

        r_pa = norm2(this%objects(1)%xc - this%objects(1)%marker_list(3)%loc)
        r_im = r_pa

        im_cent(1) = this%objects(i)%xc(1)
        im_cent(2) = 0.0_R8 - r_im
        dx = tree%dr_base(1) * (0.5_R8**(tree%highest_lvl-1))

!        d_ijk = norm2(im_cent(:) - this%objects(i)%xc(:)) - (r_pa + r_im)
!!        d_ijk = (this%objects(i)%xc(:) - im_cent(:) ) - (r_pa + r_im)
!        if (d_ijk > d_0) then
!            f(:) = 0.0_R8
!
!        else if (d_ijk >= 0.0_R8 .and. d_ijk <= d_0) then
!            f(1) = inv_epsi_1 * (this%objects(i)%xc(1) - im_cent(1)) * ((d_0 - d_ijk)**2)
!!            f(2) = inv_epsi_1 * (this%objects(i)%xc(2) - im_cent(2)) * ((d_0 - d_ijk)**2) * 100
!            f(2) = inv_epsi_1 * (this%objects(i)%xc(2) - im_cent(2)) * ((d_0 - d_ijk)**2) &
!                 + this%objects(i)%m * g * 0.33333333_R8
!
!        else if(d_ijk < 0.0_R8) then
!            f(1) = inv_epsi_2 * (this%objects(i)%xc(1) - im_cent(1) ) * (- d_ijk)
!!            f(2) = inv_epsi_2 * (this%objects(i)%xc(2) - im_cent(2) ) * (- d_ijk) * 100
!            f(2) = inv_epsi_2 * (this%objects(i)%xc(2) - im_cent(2) ) * (- d_ijk)         &
!                 + this%objects(i)%m * g * 0.33333333_R8
!        else
!            print*, "Unknown!"
!        end if


!        ! Wang 2008
!        d_ijk = norm2 (im_cent(:) - this%objects(i)%xc(:))
!        xi = 2.5_R8 * dx
!        inv_epsi_1 = 1.0_R8 / ((2*dx)**2)
!        inv_epsi_2 = 1.0_R8 / (dx)
!
!        if (d_ijk > (2*r_pa)+xi) then
!            f(:) = 0.0_R8
!
!        else if (d_ijk <= 2*r_pa) then
!            f(1) = 0.0_R8
!            f(2) = inv_epsi_1 * (this%objects(i)%xc(2) - im_cent(2)) * ((2*r_pa)-d_ijk) &
!                + this%objects(i)%m * g * 0.33333333_R8
!
!
!        else if(d_ijk < (2*r_pa)+xi .and. d_ijk > (2*r_pa) ) then
!            f(1) = 0.0_R8
!            f(2) = inv_epsi_2 * (this%objects(i)%xc(2) - im_cent(2)) * ( ( (2*r_pa)+xi-d_ijk ) ** 2) &
!                 + this%objects(i)%m * g * 0.33333333_R8
!
!        else
!            print*, "Unknown!"
!        end if


!        ! Glowinski
!        epsi = 1.0 !((2*r_pa*r_pa) * 0.1 / (r_pa+r_pa))**2
!        d_ijk = norm2 (im_cent(:) - this%objects(i)%xc(:))
!        if (d_ijk > (2*r_pa)+xi) then
!            f(:) = 0.0_R8
!        else
!            f(1) = 0.d0
!!            f(2) = (1.d0/epsi) * (d_ijk) * ( ((2*r_pa) + xi - abs(d_ijk))**2 ) &
!!            + this%objects(i)%m*g !* 0.33333333_R8
!            f(2) = (inv_epsi_2) * (d_ijk) * ( ((2*r_pa) + xi - abs(d_ijk))**2 ) &
!            + this%objects(i)%m*g !* 0.33333333_R8
!        end if


!        ! Glowinski new
!        xi = 2.0_R8 * dx
!        epsi  =  2.0_R8 * (10**(-3.0_R8)) !1.0_R8 * (10**(0.0_R8)) !1.0_R8 * (10**(-2.0_R8))
!        d_ijk = norm2 (im_cent(:) - this%objects(i)%xc(:))
!        if (d_ijk > (2*r_pa)+xi) then
!            f(:) = 0.0_R8
!        else
!            f(1) = 0.0_R8
!            f(2) = (this%objects(i)%m*g / epsi) * ( ( (d_ijk-(2*r_pa)-xi) / xi )**2 )
!        end if


        ! Niu
        xi = 2.0_R8 * dx
        epsi =  ((2.0_R8*r_pa*r_pa) /  (r_pa+r_pa) ) ** 2
!        d_ijk = norm2 ( im_cent(:) - this%objects(i)%xc(:) ) !+ r_pa
        d_ijk = this%objects(i)%xc(2)

        if ( d_ijk > (2*r_pa)+xi ) then
            f(:) = 0.0_R8
        else
            f(1) = 0.0_R8
            f(2) = 10.0d0*epsi * ( 2*( (r_pa/d_ijk)**14 ) - ( (r_pa/d_ijk)**8 ) ) * &
                   ( d_ijk/(r_pa**2) )
        end if

    end function get_collision_force_wall_3



    function get_collision_force_wall_4( this, tree, i ) result (f)
        implicit none
        class(immersed_boundary_t), intent(inout) :: this
        class(af_t),                intent(inout) :: tree
        integer(kind=I4),           intent(in)    :: i
        real(kind=R8)                             :: f(NDIM)
        f = 0.0d0
    end function get_collision_force_wall_4



    subroutine get_predefined_forces( this )
        implicit none
        class(object_t), intent(inout) :: this
        print*, "FORCES NOT YET!"
    end subroutine get_predefined_forces



    function body_force( this ) result (bforce)
        implicit none
        class(object_t), intent(inout) :: this
        real(kind=R8)                  :: bforce(NDIM)

        bforce(1) = 0.0_R8
        bforce(2) = - this%m * g * 0.3_R8  !> M_p * g * (1 - (rho_fluid / rho_p)
    end function body_force



    subroutine update_an_object_acceleration( this )
        implicit none
        class(object_t), intent(inout) :: this

        if (this%predefined_move) then
            call this%get_predefined_acceleration( )
        else
            call this%get_calculated_acceleration( )
        end if
    end subroutine update_an_object_acceleration



    subroutine get_calculated_acceleration( this )
        implicit none
        class(object_t), intent(inout) :: this
        real(kind=R8)                   :: f(NDIM), col_f(NDIM)

        this%xc_dd(2) = this%force(2) / this%m
        this%xc_dd(1) = this%force(1) / this%m
#if NDIM==2
        this%om_dd(3) = this%moment / this%im
#endif
    end subroutine get_calculated_acceleration



    subroutine get_predefined_acceleration( this )
        implicit none
        class(object_t), intent(inout) :: this
        real(kind=R8)                  :: acc(NDIM)

        this%xc_dd(:) = p_predefined_linear_acceleration (this%xc)
#if NDIM == 2
        this%om_dd(:) = p_predefined_rotation_acceleration(this%om)
#endif
    end subroutine get_predefined_acceleration



    subroutine update_an_object_velocity( this )
        implicit none
        class(object_t), intent(inout) :: this

        if (this%predefined_move) then
            call this%get_predefined_velocity( )
        else
            call this%get_calculated_velocity( )
        end if
    end subroutine update_an_object_velocity



    subroutine get_calculated_velocity( this )
        implicit none
        class(object_t), intent(inout) :: this
        real(kind=R8) :: temp

         this%xc_d(:) = (dt * this%xc_dd(:)) + this%xc_d(:)
#if NDIM == 2
        this%om_d(3) = (dt * this%om_dd(3)) + this%om_d(3)
#endif
    end subroutine  get_calculated_velocity



    subroutine update_an_object_location( this )
        implicit none
        class(object_t), intent(inout) :: this
        if (this%predefined_move) then
            call this%get_predefined_location()
        else
            call this%get_calculated_location()
        end if
    end subroutine update_an_object_location



    subroutine get_predefined_location( this )
        implicit none
        class(object_t), intent(inout) :: this

        this%xc(:) = p_predefined_linear_move( this%xc )
#if NDIM == 2
        this%om(:) = p_predefined_rotation_move( this%om )
#endif
    end subroutine get_predefined_location



    subroutine get_calculated_location( this )
        implicit none
        class(object_t), intent(inout) :: this
        real(kind=R8)                  :: dxc(NDIM)
#if NDIM == 2
        real(kind=R8)                  :: dom
#endif
        dxc(:)     = (dt * this%xc_d(:))
        this%xc(:) = this%xc(:) + dxc(:)
#if NDIM == 2
        dom     = (dt * this%om_d(3))
        this%om = this%om + dom
#endif
    end subroutine  get_calculated_location



    subroutine get_predefined_velocity( this )
        implicit none
        class(object_t), intent(inout) :: this
#if NDIM == 2
        this%xc_d(:) = p_predefined_linear_velocity (this%xc)
        this%om_d(:) = p_predefined_rotation_velocity (this%om)
#endif
    end subroutine get_predefined_velocity



    subroutine update_lagrangian_markers_location( this )
        implicit none
        class(object_t), intent(inout) :: this
        integer(kind=I4)               :: n
        real(kind=R8)                  :: r(NDIM)
#if NDIM == 2
        !$OMP DO private(r)
        do n = 1, this%nmarkers
            r(:) = this%marker_list(n)%loc(:) - this%xc(:)
            this%marker_list(n)%loc(2) = this%marker_list(n)%loc(2) +                    &
            &  ( (+this%om(3)*r(1) + this%xc_d(2)) * dt)
            this%marker_list(n)%loc(1) = this%marker_list(n)%loc(1) +                    &
            &  ( (-this%om(3)*r(2) + this%xc_d(1)) * dt)
        end do
        !$OMP END DO
#endif
    end subroutine update_lagrangian_markers_location



    subroutine rebuild_object( this, tree )
        implicit none
        class(object_t),            intent(inout)  :: this
        class(af_t),                intent(inout)  :: tree
        integer(kind=I4)                           :: i_obj
        real(kind=R8), allocatable, dimension(:,:) :: coords
        real(kind=R8), allocatable, dimension(:)   :: valz
        real(kind=R8), allocatable                 :: mass_arr(:), im_arr(:)
        integer(kind=I4)                           :: n

#if NDIM == 2
        integer(kind=I4) :: ncorners = 4
        real(kind=R8)    :: corner_distance(4)
#endif
        i_obj = this%id

        !>  Calculated the distance function with RBF method
        call get_initial_df ()
        call af_gc_tree( tree, [i_df+i_obj-1] )

        !>  Improved the distance function with geometric formula
        call af_loop_box( tree, get_final_df_ker, leaves_only=.true. )
        call af_gc_tree( tree, [i_df+i_obj-1] )
        call af_restrict_tree( tree, [i_df+i_obj-1] )


    contains


        subroutine get_final_df_ker( box )
            implicit none
            type(box_t), intent(inout) :: box
            integer(kind=I4)           :: nc, i, j, n
            real(kind=R8)              ::  xi(NDIM, 1), temp

            nc = box%n_cell
            do j = 1, nc
                do i = 1, nc
                    xi(:, 1) = af_r_cc(box, [IJK])
                    box%cc(i, j, i_df+i_obj-1) = sign (100.d0, box%cc(i, j, i_df+i_obj-1))
                    do n = 1, size( this%marker_list)
                        temp = norm2( xi(:,1) - this%marker_list(n)%loc(:) )
                        box%cc(i, j, i_df+i_obj-1) =                                     &
                        &  sign( min( temp, abs (box%cc(i, j, i_df+i_obj-1)) ),           &
                        &             box%cc(i, j, i_df+i_obj-1) )
                    end do
                end do
            end do
        end subroutine get_final_df_ker


        subroutine get_initial_df( )
            implicit none

            if ( .not. allocated( coords ) ) then
                allocate( coords(NDIM, this%nmarkers + ncorners) )
            end if

            if ( .not. allocated( valz ) ) then
                allocate( valz(this%nmarkers + ncorners + 3) )
            end if
            valz(1:this%nmarkers) = landa_g

            if ( .not. allocated (this%ws) ) then
                allocate( this%ws(this%nmarkers + ncorners + 3) )
            end if

            coords(1, 1:this%nmarkers) = this%marker_list(:)%loc(1)
            coords(2, 1:this%nmarkers) = this%marker_list(:)%loc(2)

            do n = 1, ncorners
                if(n==1 .or. n==4)  then
                    coords(1, this%nmarkers + n) = domain_len(1)
                else
                    coords(1, this%nmarkers + n) = 0.0_R8
                end if

                if (n==1 .or. n==2) then
                    coords(2, this%nmarkers + n) = domain_len(2)
                else
                    coords(2, this%nmarkers + n) = 0.0_R8
                end if
            end do

            call get_corner_distance( )

            do n = 1, ncorners
                valz(this%nmarkers + n) =  corner_distance(n) + landa_g
            end do

            valz(this%nmarkers + ncorners + 1 : this%nmarkers + ncorners + 3 ) = 0.0_R8
            this%ws = 0.0_R8

            call rbf_weight( 2, this%nmarkers+ncorners, coords, valz, this%ws )
            call af_loop_box( tree, get_df_ker )
        end subroutine get_initial_df


        subroutine get_df_ker( box )
            implicit none
            type(box_t), intent(inout) :: box
            integer(kind=I4)           :: nc, i, j, nd
            real(kind=R8)              :: fi(1), xi(NDIM, 1)
            integer(kind=I4)           :: ni

            nd = size( coords, 2 )
            ni = 1
            nc = box%n_cell
            do j = 1, nc
                do i = 1, nc
                    xi(:, 1) = af_r_cc(box, [IJK])
                    call rbf_interp (NDIM, nd, coords, this%ws, ni, xi, fi)
                    box%cc(i, j, i_df+i_obj-1) = fi(1) - landa_g
                end do
            end do
        end subroutine get_df_ker


        subroutine get_corner_distance( )
            implicit none
            integer(kind=I4) :: n, corner
            real(kind=R8)    :: temp_dis

            do corner = 1, ncorners
                corner_distance(corner) = max( domain_len(1), domain_len(2) )
                do n = 1, size (this%marker_list)
                    temp_dis = norm2 (  &
                    &  this%marker_list(n)%loc(:) - coords(:, this%nmarkers + corner))

                    if(temp_dis <= corner_distance(corner)) then
                        corner_distance(corner) = temp_dis
                    end if
                end do
            end do
        end subroutine get_corner_distance
    end subroutine rebuild_object


end module m_immersed_boundary
