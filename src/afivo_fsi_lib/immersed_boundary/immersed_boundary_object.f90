#include "cpp_macros.h"
submodule (m_immersed_boundary) build_obj


contains


     module subroutine build_object( obj, i_obj, tree, is_initiated )
        implicit none
        class(object_t),  intent(inout)        :: obj
        class(af_t),      intent(inout)        :: tree
        integer(kind=I4), intent(in)           :: i_obj
        logical,          intent(in), optional :: is_initiated
        real(kind=R8), allocatable, dimension(:,:) :: coords
        real(kind=R8), allocatable, dimension(:)   :: valz
        real(kind=R8), allocatable                 :: mass_arr(:), im_arr(:)
        integer(kind=I4)                           :: n, nobj_tmp
        character (kind=CK, len=100)               :: ext
        character (kind=CK, len=200)               :: object_fname
#if NDIM == 2
        integer(kind=I4) :: ncorners = 4
        real(kind=R8)    :: corner_distance(4)
#endif
        if (present(is_initiated) ) then
            if(is_initiated) then
                write(ext, *) i_obj
                object_fname = trim( adjustl( "./inputs/object_") ) // &
                               trim( adjustl( ext) ) //                &
                               trim( adjustl( ".dat") )
                call obj%read_markers( object_fname )
            end if
        else
            write(ext, *) i_obj
            object_fname = trim( adjustl( "./inputs/object_") ) // &
                           trim( adjustl( ext) ) //                &
                           trim( adjustl( ".dat" ) )
            call obj%read_markers( object_fname )
        end if

        obj%id = i_obj
        call config_glb%get( section_name ='immersed_boundary',                          &
        &  option_name = 'predefined_move', val = obj%predefined_move)
        !> Reads the mass of this object
        call config_glb%get( section_name ='immersed_boundary', option_name = 'nobjects',&
        &  val = nobj_tmp)
        if (allocated( mass_arr ) ) deallocate( mass_arr )
        allocate( mass_arr(nobj_tmp) )
        call config_glb%get( section_name ='immersed_boundary', option_name = 'mass',    &
        &  val = mass_arr )
        obj%m = mass_arr( obj%id )

        if ( allocated( im_arr ) ) deallocate( im_arr )
        allocate( im_arr(nobj_tmp) )
        call config_glb%get( section_name ='immersed_boundary',                         &
        &  option_name = 'moment_of_inertia', val = im_arr)
        obj%im = im_arr(obj%id)

        !>  Calculated the distance function with RBF method
        call get_initial_df( )
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
                    do n = 1, size (obj%marker_list)
                        temp = norm2 ( xi(:,1) - obj%marker_list(n)%loc(:) )
                        box%cc(i, j, i_df+i_obj-1) =                                     &
                        &  sign (min (temp, abs (box%cc(i, j, i_df+i_obj-1))),           &
                        &             box%cc(i, j, i_df+i_obj-1))
                    end do
                end do
            end do
        end subroutine get_final_df_ker


        subroutine get_initial_df()
            implicit none
            if (.not. allocated( coords ) ) then
                allocate( coords(NDIM, obj%nmarkers + ncorners) )
            end if

            if(.not. allocated( valz )) then
                allocate( valz(obj%nmarkers + ncorners + 3) )
            end if
            valz(1:obj%nmarkers) = landa_g

            if (.not. allocated (obj%ws)) then
                allocate( obj%ws(obj%nmarkers + ncorners + 3) )
            end if

            coords(1, 1:obj%nmarkers) = obj%marker_list(:)%loc(1)
            coords(2, 1:obj%nmarkers) = obj%marker_list(:)%loc(2)

            do n = 1, ncorners
                if(n==1 .or. n==4)  then
                    coords(1, obj%nmarkers + n) = domain_len(1)
                else
                    coords(1, obj%nmarkers + n) = 0.0_R8
                end if

                if (n==1 .or. n==2) then
                    coords(2, obj%nmarkers + n) = domain_len(2)
                else
                    coords(2, obj%nmarkers + n) = 0.0_R8
                end if
            end do

            call get_corner_distance()

            do n = 1, ncorners
                valz(obj%nmarkers + n) =  corner_distance(n) + landa_g
            end do

            valz(obj%nmarkers + ncorners + 1 : obj%nmarkers + ncorners + 3 ) = 0.0_R8
            obj%ws = 0.0_R8

            call rbf_weight (2, obj%nmarkers+ncorners, coords, valz, obj%ws)
            call af_loop_box (tree, get_df_ker)
        end subroutine get_initial_df


        subroutine get_df_ker(box)
            implicit none
            type(box_t), intent(inout) :: box
            integer(kind=I4)           :: nc, i, j, nd
            real(kind=R8)              ::  fi(1), xi(NDIM, 1)
            integer(kind=I4)           :: ni

            nd = size (coords, 2)
            ni = 1
            nc = box%n_cell
            do j = 1, nc
                do i = 1, nc
                    xi(:, 1) = af_r_cc(box, [IJK])
                    call rbf_interp (NDIM, nd, coords, obj%ws, ni, xi, fi)
                    box%cc(i, j, i_df+i_obj-1) = fi(1) - landa_g
                end do
            end do
        end subroutine get_df_ker


        subroutine get_corner_distance()
            implicit none
            integer(kind=I4) :: n, corner
            real(kind=R8)    :: temp_dis

            do corner = 1, ncorners
                corner_distance(corner) = max (domain_len(1), domain_len(2))
                do n = 1, size (obj%marker_list)
                    temp_dis = norm2 (                                                   &
                    &  obj%marker_list(n)%loc(:) - coords(:, obj%nmarkers + corner))

                    if(temp_dis <= corner_distance(corner)) then
                        corner_distance(corner) = temp_dis
                    end if
                end do
            end do
        end subroutine get_corner_distance
    end subroutine build_object

end submodule  build_obj
