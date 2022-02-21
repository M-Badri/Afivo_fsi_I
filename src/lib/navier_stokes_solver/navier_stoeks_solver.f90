#include "cpp_macros.h"
module m_navier_stokes_solver
    use m_solver_parameters
    use m_config
    use m_af_all
    use m_incomp_flow_field
    use m_projection
    use m_helmholtz
    use m_immersed_boundary
!    use m_write_silo

    implicit none
    private
    public :: build_navier_stokes_solver

    type(t_config) :: module_config
    integer(kind=I4) :: itr_min, itr_max
    integer(kind=I4) :: hard_print_frq, monitor_print_frq
    integer(kind=I4) :: refine_frq, force_cal_frq, back_up_frq
    logical :: has_immersed_objects, are_objects_moving
    logical :: is_flow_refined, is_refined, are_objects_refined
    real(kind=R8) :: cfl_max

    procedure (generic_solve_ker),      pointer, public  :: p_solve
    procedure (generic_update_vel_ker), pointer, private :: p_update_v_ker
    procedure (generic_update_pre_ker), pointer, private :: p_update_p_ker
    procedure (generic_ts_ker),         pointer, private :: p_timestep_ker
    procedure (generic_cgrad),          pointer, private :: p_cell_xgrad
    procedure (generic_cgrad),          pointer, private :: p_cell_ygrad
    procedure (generic_fgrad),          pointer, private :: p_face_xgrad
    procedure (generic_fgrad),          pointer, private :: p_face_ygrad
    procedure (generic_diffusion),      pointer, private :: p_diffusion
    procedure (generic_convective),     pointer, private :: p_convective
    procedure (generic_c2f_ker),        pointer, private :: p_c2f_ker
    procedure (generic_set_prhs_ker),   pointer, private :: p_set_prhs_ker
    procedure (generic_set_hrhs_ker),   pointer, private :: p_set_hrhs_ker
    procedure (generic_get_vort_ker),   pointer, private :: p_get_vort_ker
    procedure (generic_get_div_ker),    pointer, private :: p_get_div_ker


    abstract interface
        subroutine generic_c2f_ker(box, i_var)
            import
            type(box_t), intent(inout) :: box
            integer(kind=I4), intent(in) :: i_var(:)
        end subroutine generic_c2f_ker


        subroutine generic_solve_ker (tree, refine_info)
            import
            class(af_t), intent(inout) :: tree
            class(ref_info_t), intent(inout) :: refine_info
        end subroutine generic_solve_ker


        subroutine generic_update_vel_ker (box)
            import
            type(box_t), intent(inout) :: box
        end subroutine generic_update_vel_ker


        subroutine generic_update_pre_ker (box)
            import
            type(box_t), intent(inout) :: box
        end subroutine generic_update_pre_ker


        subroutine generic_ts_ker (box)
            import
            type(box_t), intent(inout) :: box
        end subroutine generic_ts_ker


        function generic_cgrad (box, i_var) result(grad)
            import
            class(box_t), intent(in) :: box
            integer(kind=I4), intent(in) :: i_var
            real(kind=dp) :: grad(DTIMES(box%n_cell))
        end function generic_cgrad


        function generic_fgrad (box, i_var) result(grad)
            import
            class(box_t), intent(in) :: box
            integer(kind=I4), intent(in) :: i_var
            real(kind=dp) :: grad(DTIMES(box%n_cell+1))
        end function generic_fgrad


        function generic_diffusion (box, i_var) result(diff)
            import
            class(box_t), intent(in) :: box
            integer(kind=I4), intent(in) :: i_var
            real(kind=dp) :: diff(DTIMES(box%n_cell))
        end function generic_diffusion


        function generic_convective (box, i_var) result(advec)
            import
            class(box_t), intent(in) :: box
            integer(kind=I4), intent(in) :: i_var(:)
            real(kind=dp) :: advec(DTIMES(box%n_cell))
        end function generic_convective


        subroutine generic_set_prhs_ker(box)
           import
           type(box_t), intent(inout) :: box
        end subroutine generic_set_prhs_ker


        subroutine generic_set_hrhs_ker(box)
           import
           type(box_t), intent(inout) :: box
        end subroutine generic_set_hrhs_ker


        subroutine generic_get_vort_ker(box)
           import
           type(box_t), intent(inout) :: box
        end subroutine generic_get_vort_ker


        subroutine generic_get_div_ker(box, i_var)
            import
            type(box_t), intent(inout) :: box
            integer(kind=I4), intent(in) :: i_var(:)
        end subroutine generic_get_div_ker
    end interface



contains



    subroutine build_navier_stokes_solver (config, tree, refine_info)
        implicit none
        type(t_config),    intent(inout) :: config
        class(af_t),       intent(inout) :: tree
        class(ref_info_t), intent(inout) :: refine_info

        module_config = config

        call config%get (section_name ='time_parameters', option_name = 'iteration_min', &
        &  val = itr_min)
        call config%get (section_name ='time_parameters', option_name = 'iteration_max', &
        &  val = itr_max)
        call config%get (section_name ='time_parameters', option_name = 'hard_print_frq',&
        &  val = hard_print_frq)
        call config%get (section_name ='time_parameters', option_name = 'back_up_frq',   &
        &  val = back_up_frq)
        call config%get (section_name ='refinement', option_name = 'refine_frq',         &
        &  val = refine_frq)
        call config%get (section_name ='time_parameters', option_name =                  &
        &  'monitor_print_frq', val = monitor_print_frq)
        call config%get (section_name ='time_parameters', option_name = 'force_cal_frq', &
        &  val = force_cal_frq)
        call config%get (section_name ='time_parameters', option_name ='dt', val = dt)
        call config%get (section_name ='immersed_boundary',                              &
        &  option_name = 'has_immersed_objects', val = has_immersed_objects)
        call config%get (section_name='immersed_boundary',                               &
        &  option_name ='are_objects_moving', val = are_objects_moving)
        call config%get (section_name ='problem', option_name ='cfl_max', val = cfl_max)
        call config%get (section_name='refinement',                                      &
        &  option_name ='is_refined', val = is_refined)
        call config%get (section_name='refinement',                                      &
        &  option_name ='are_objects_refined', val = are_objects_refined)
        call config%get (section_name='refinement',                                      &
        &  option_name ='is_flow_refined', val = is_flow_refined)

        !> Set the parameters to make the problem non-dimensional or not.
        call config%get (section_name ='problem', option_name ='is_nondimesional',       &
        &  val = is_nondimesional)
        if (is_nondimesional) then
            call config%get (section_name ='problem', option_name ='reynolds', val = rey)
            rho_f = 1.0_R8
            nu_f  = 1.0_R8
        else
            call config%get (section_name ='problem', option_name ='nu_f', val = nu_f)
            call config%get (section_name ='problem', option_name ='rho_f', val = rho_f)
            rey = 1.0_R8
        end if

        if (refine_frq == 0) then
            refine_frq = itr_max + 1
        end if

        if (has_immersed_objects) then
            p_solve => predictor_corrector_simp_ibm
            call buid_helmholtz (module_config, tree)
        else
            p_solve => ab_simp_solver
            call buid_helmholtz (module_config, tree)
        end if

        !> Builds an object of a Poisson solver being used for projection
        call buid_projector (module_config, tree)

        !> Makes the needed directories for output files
        call execute_command_line ("mkdir -p ./output/")
        call execute_command_line ("mkdir -p ./output/back_up/")
    end subroutine build_navier_stokes_solver



    subroutine predictor_corrector_simp_ibm (tree, refine_info)
        implicit none
        class(af_t), intent(inout) :: tree
        class(ref_info_t), intent(inout) :: refine_info
        type(immersed_boundary_t) :: ib
        character(kind=CK, len=100) :: fname
        real(kind=R8), allocatable :: fx_arr(:), fy_arr(:)
        real(kind=R8) :: cfl_c
        real(kind=R8) :: step_time
        integer(kind=I4) :: i
        integer(kind=I4) :: itr
        integer(kind=I4) :: count_rate, t_start, t_end
        logical :: is_on_hard = .false.

        p_timestep_ker => ab_exp_smac_timestep_ker
        p_update_v_ker => smac_update_velocity_ker_2o
        p_update_p_ker => smac_update_pressure_ker_2o
        p_cell_xgrad   => central_2order_cell_xgrad
        p_cell_ygrad   => central_2order_cell_ygrad
        p_face_xgrad   => central_2order_face_xgrad
        p_face_ygrad   => central_2order_face_ygrad
        p_diffusion    => central_2order_diffusion
        p_convective   => central_2order_convective
        p_c2f_ker      => central2order_c2f_p_decoupled_ker
        p_set_prhs_ker => central2order_set_prhs_ker
        p_set_hrhs_ker => central2order_smac_set_hrhs_withforce_ker
        p_get_vort_ker => central2order_vorticity_ker
        p_get_div_ker  => central2order_divegence_ker

        !> Initializes the immersed boundary
        call build_immersed_boundary (tree, refine_info, module_config, ib)

        !> Initializes the tree from a back up file on the hard and starts the iteration
        !  from the "itr_min" and its corresponding time.
        if (itr_min /= 1) then
            call af_destroy(tree)
            write(fname, "(A,I0)") "./output/back_up/ns_bu_" // DIMNAME // "_", itr_min
            fname = trim (adjustl (fname)) // ".dat"
            call af_read_tree (tree, fname)
!            field_time = (itr_min+1) * dt
!            itr_min = itr_min + 1
            itr_min = (itr_min+1)

            !> Fills ghost cells and boundaries.
            call af_gc_tree (tree, [i_uold, i_vold])
            call af_gc_tree (tree, [i_unew, i_vnew])
!            call af_gc_tree (tree, [i_uold, i_volder])
            call af_gc_tree (tree, [i_p])

        else
            field_time = 0.0_R8
            !> Changes the refinement when it is necessary only around the objects
            if (is_refined) then
                call refine_field (tree = tree, refine_info = refine_info,                   &
                &  object_ref = are_objects_refined, flow_ref = is_flow_refined, ib = ib)

                !> Fixed objects are refined only once before time-stepping.
                if (.not. are_objects_moving) then
                    is_refined = .false.
                end if
            end if

            !> For all boxes interpolates velocities from cell to face.
            call af_loop_box_iarg (tree, p_c2f_ker, [i_uold, i_vold])

            !> Fills ghost cells and boundaries.
            call af_gc_tree (tree, [i_uold, i_vold])
            call af_gc_tree (tree, [i_p])

        end if



        !> Writes the information of the problem before starting time steps.
        call print_out_info (tree)

!        !> Fills ghost cells and boundaries.
!        call af_gc_tree (tree, [i_uold, i_vold])
!        call af_gc_tree (tree, [i_p])

!        !> For all boxes interpolates velocities from cell to face.
!        call af_loop_box_iarg (tree, p_c2f_ker, [i_uold, i_vold])

        allocate (fx_arr(ib%nobjs), fy_arr(ib%nobjs))


        do itr = itr_min, itr_max
            call system_clock(t_start, count_rate)

            field_time = dt * (itr)


!            call af_gc_tree (tree, [i_uold, i_vold])
!            call af_gc_tree (tree, [i_p])

            !> Takes one time step integration for momentum equation WITHOUT IB force.
            !  Gives the predicted velocities from which the forces must be calculated.
            call af_loop_box (tree, p_timestep_ker)

            !> Boundary condition
            call af_gc_tree (tree, [i_unew, i_vnew])
!            call af_restrict_tree (tree, [i_unew, i_vnew])

            !> Calculates Immersed boundary forces.
            call ib%get_ext_forces (tree)
            call af_gc_tree (tree, [i_ibfx, i_ibfy])
!            call af_restrict_tree (tree, [i_ibfx, i_ibfy])

            !> Takes one time step integration for momentum equation WITH IB force.
            !  Gives the corrected velocities affected by the object(s).
            !  A Helmholtz equation is solved for taking a Crank-Nicolson step.
            call af_loop_box (tree, p_set_hrhs_ker)
            call get_helmholtz (tree)

!            !> Keeps the velocity after imposing boundary condition before projection
!            call af_tree_copy_cc(tree, i_unew, i_utmp)
!            call af_tree_copy_cc(tree, i_vnew, i_vtmp)

            !> Boundary condition
            call af_gc_tree (tree, [i_unew, i_vnew])
!            call af_restrict_tree (tree, [i_unew, i_vnew])

!            !> Calculates divergence before projection
!            call af_loop_box_iarg (tree, p_get_div_ker, [i_div_1, i_unew, i_vnew])
!            call af_gc_tree (tree, [i_div_1])
!            call af_restrict_tree (tree, [i_div_1])

            !> For all boxes interpolates velocities from cell to face.
            call af_loop_box_iarg (tree, p_c2f_ker, [i_unew, i_vnew])

            !> Calculates the corresponding pressure for a divergence-free projection
            call af_loop_box (tree, p_set_prhs_ker)
            call p_get_projector (tree)

            !> After projection fills the boundary for phi and p, then restricts those
            !  into coarser levels.
            call af_gc_tree (tree, [i_p, i_phi])
!            call af_restrict_tree (tree, [i_p, i_phi])

            !> Updates the velocity into a divergence-free fields
            call af_loop_box (tree, p_update_v_ker)
            call af_gc_tree (tree, [i_unew, i_vnew])
!            call af_restrict_tree (tree, [i_unew, i_vnew])

            !> Updates the pressure based on the result of the projection
            call af_loop_box (tree, p_update_p_ker)
            call af_gc_tree (tree, [i_p, i_phi])
!            call af_restrict_tree (tree, [i_p, i_phi])

            !> Calculates the vorticity
            call af_loop_box (tree, p_get_vort_ker)
            call af_gc_tree (tree, [i_vortz])
!            call af_restrict_tree (tree, [i_vortz])

            !> Calculating the acceleration and velocity of all objects, re-locates them
            !  and updates the corresponding distance functions.
            if (are_objects_moving) then
                call ib%update_objects (tree)
            end if

!            !> Calculates and writes out the forces.
!            if (mod(itr, force_cal_frq) == 0) then
!                call ib%writeout_ext_forces()
!                call ib%writeout_locations_velocities()
!                call ib%writeout_objects_reynolds ()
!                call ib%writeout_objects_energy ()
!                do i = 1, ib%nobjs
!                    fx_arr(i) = ib%objects(i)%ext_force(1)
!                    fy_arr(i) = ib%objects(i)%ext_force(2)
!                end do
!            end if

            !> Changes the refinement when it is necessary around the objects and the high
            !  vorticity spots.
            if (is_refined .and. (mod (itr, refine_frq) == 0)) then
                call refine_field (tree = tree, refine_info = refine_info,               &
                &  object_ref = are_objects_refined, flow_ref = is_flow_refined, ib = ib)
!                call af_gc_tree (tree, [i_unew, i_uold])
!                call af_gc_tree (tree, [i_vnew, i_vold])
!                call af_gc_tree (tree, [i_vortz, i_p])
!                call af_restrict_tree (tree, [i_unew, i_uold, i_uolder, &
!                                              i_vnew, i_vold, i_volder, &
!                                              i_vortz, i_p])
!                call af_loop_box_iarg (tree, p_c2f_ker, [i_unew, i_vnew])
!                call af_loop_box_iarg (tree, p_c2f_ker, [i_uold, i_vold])
            end if

!            !> Calculates divergence after projection
!            call af_loop_box_iarg (tree, p_get_div_ker, [i_div_2, i_unew, i_vnew])
!            call af_gc_tree (tree, [i_div_2])
!            call af_restrict_tree (tree, [i_div_2])

            !> Calculates the Convective CFL and stops the code if it is too large.
            call get_cfl_c (tree, cfl_c, dt)

            !> Prints out the controlling information.
            if (mod (itr, monitor_print_frq) == 0) then
                !> Writes on monitor.

                !> Writes out the result on the hard
                if (mod(itr, hard_print_frq) == 0) then
                    is_on_hard = .true.
                    write(fname, "(A,I0)") "ns_" // DIMNAME // "_", itr
                    !> Write the cell centered data of tree to a vtk unstructured file.
                    !  Only the leaves of the tree are used
                    call af_write_vtk(tree,                                              &
                    &  trim(fname), itr, field_time, dir="output",                       &
                    &  ixs_cc =                                                          &
                    &  [i_unew, i_vnew, i_p, i_df, i_vortz])

!                   call af_write_silo(tree,                                              &
!                   trim(fname), itr, field_time, dir="output",                           &
!                   ixs_cc = [i_unew, i_vnew, i_p, i_df, i_vortz])
                end if

                call time_stamp (time=field_time, time_step=itr, n_obj=ib%nobjs,         &
                &  fx=fx_arr, fy=fy_arr, cfl_c = cfl_c, is_on_hard = is_on_hard,         &
                &  step_time = step_time)
                is_on_hard = .false.
                write (*, *)
            end if

            !> Saves the whole data of the tree for in case resuming.
            if (mod (itr, back_up_frq) == 0 ) then
                write(fname, "(A,I0)") "./output/back_up/ns_bu_" // DIMNAME // "_", itr
                call af_write_tree (tree, filename=fname)
            end if

            !> Imposes BC before swapping the velocity.
            call af_gc_tree (tree, [i_uold, i_vold, i_unew, i_vnew])

            !> Swapping the data
            call af_tree_copy_cc(tree, i_uold, i_uolder)
            call af_tree_copy_cc(tree, i_vold, i_volder)

            call af_tree_copy_cc(tree, i_unew, i_uold)
            call af_tree_copy_cc(tree, i_vnew, i_vold)

            call af_tree_copy_fc(tree, i_f, i_folder)

            call system_clock(t_end, count_rate)
            step_time = (t_end-t_start) / real(count_rate, dp)
        end do
    end subroutine predictor_corrector_simp_ibm



    subroutine refine_field (tree, refine_info, object_ref, flow_ref, ib)
        implicit none
        class (af_t), intent(inout) :: tree
        class(ref_info_t), intent(inout) :: refine_info
        class (immersed_boundary_t), intent(in) :: ib
        logical, intent (in) :: object_ref, flow_ref
        integer :: i_obj

        if ( (.not.(flow_ref)) .and. object_ref) then
            do i_obj = 1, ib%nobjs
                do
                    call af_gc_tree (tree,                                               &
                    &  [i_vortz, i_unew, i_vnew, i_uold, i_vold, i_p])
                    call af_adjust_refinement (tree, ref_objects, refine_info)
                    ! If no new boxes have been added, exit the loop
                    if (refine_info%n_add == 0) exit
                end do
            end do

        else if ( (.not.(object_ref)) .and. flow_ref) then
            do
                call af_gc_tree (tree,                                                   &
                &  [i_vortz, i_unew, i_vnew, i_uold, i_vold, i_p])
                call af_adjust_refinement (tree, ref_flow, refine_info)
                ! If no new boxes have been added, exit the loop
                if (refine_info%n_add == 0) exit
            end do

        else if (flow_ref .and. object_ref) then
            do i_obj = 1, ib%nobjs
                do
                    call af_gc_tree (tree,                                               &
                    &  [i_vortz, i_unew, i_vnew, i_uold, i_vold, i_p])
                    call af_adjust_refinement (tree, ref_objects, refine_info)
                    ! If no new boxes have been added, exit the loop
                    if (refine_info%n_add == 0) exit
                end do

                do
                    call af_gc_tree (tree,                                               &
                    &  [i_vortz, i_unew, i_vnew, i_uold, i_vold, i_p])
                    call af_adjust_refinement (tree, ref_flow, refine_info)
                    ! If no new boxes have been added, exit the loop
                    if (refine_info%n_add == 0) exit
                end do
            end do

        else
            write(*, *) "The refinement condition is not defined."
            stop

        end if


    contains


        subroutine ref_flow (box, cell_flags)
            implicit none
            type(box_t), intent(in) :: box ! A list of all boxes in the tree
            integer(kind=I4), intent(out) :: cell_flags(DTIMES(box%n_cell))
            integer(kind=I4)              :: IJK, nc
            real(kind=R8)                 :: diff, dis
            real(kind=R8) :: loc(NDIM)
            nc = box%n_cell
            loc(:) = box%r_min(:) + [(nc/2) * box%dr(1) , (nc/2) * box%dr(2)]


!            do KJI_DO(1,nc)
!                diff = abs(box%cc(i, j, i_vortz))
!
!                if( box%cc(i, j, i_df+i_obj-1) <= 1.0_R8) then
!                    cell_flags(IJK) = af_keep_ref
!                else if (diff >= 55.0_R8 .and. box%lvl < (max_refinement))then
!                    cell_flags(IJK) = af_do_ref
!                else if (diff < 5.0_R8 .and. box%lvl > (max_refinement-1)) then
!                    cell_flags(IJK) = af_rm_ref
!                else
!                    cell_flags(IJK) = af_keep_ref
!                end if
!            end do; CLOSE_DO


!            dis = abs(box%r_min(2) - 20.d0)
!            do KJI_DO(1,nc)
!
!                diff = abs(box%cc(i, j, i_vortz))
!
!                if( box%cc(i, j, i_df+i_obj-1) <= 1.0_R8) then
!                    cell_flags(IJK) = af_keep_ref
!                else if (diff >= 55.0_R8 .and. box%lvl < (max_refinement))then
!                    cell_flags(IJK) = af_do_ref
!                else if (diff < 5.0_R8 .and. box%lvl > (max_refinement-1)) then
!                    cell_flags(IJK) = af_rm_ref
!                else
!                    cell_flags(IJK) = af_keep_ref
!                end if
!            end do; CLOSE_DO

            dis = abs(box%r_min(2) - 20.d0)
            do KJI_DO(1,nc)
!                diff = abs(box%dr(1) * (box%cc(i+1, j, i_vortz) +                        &
!                &  box%cc(i-1, j, i_vortz) - 2 * box%cc(i, j, i_vortz)) +                &
!                &  box%dr(2) * (box%cc(i, j+1, i_vortz) +                                &
!                &  box%cc(i, j-1, i_vortz) - 2 * box%cc(i, j, i_vortz)))
                diff = abs(box%cc(i, j, i_vortz))

                if( abs(box%cc(i, j,  i_df+i_obj-1)) <= 1.0_R8 ) then
                    cell_flags(IJK) = af_keep_ref
                else if (diff > 0.2_R8 .and. box%lvl < 3 .and. dis < 10.0d0 )then
                    cell_flags(IJK) = af_do_ref
                else if (diff < 0.1_R8 .and. box%lvl > 2) then
                    cell_flags(IJK) = af_rm_ref
                else
                    cell_flags(IJK) = af_keep_ref
                end if
            end do; CLOSE_DO
        end subroutine ref_flow


        subroutine ref_objects(box, cell_flags)
            implicit none
            type(box_t), intent(in) :: box
            integer(kind=I4), intent(out) :: cell_flags(DTIMES(box%n_cell))
            integer(kind=I4) :: IJK, nc
            real(kind=R8) :: cntr(NDIM)

!            nc = box%n_cell
!            do KJI_DO(1, nc)
!                if( box%lvl < (max_refinement) .and.                                     &
!                &  box%cc(i, j, i_df+i_obj-1) <= 1.0_R8) then
!                    cell_flags(IJK) = af_do_ref
!
!!                else if ( (box%lvl < (max_refinement-1)) .and.                             &
!!                &  (box%cc(i, j, i_df+i_obj-1) > 0.7_R8 ).and.                   &
!!                &  (box%cc(i, j, i_df+i_obj-1) <= 1.2_R8)    &
!!                ) then
!!                    cell_flags(IJK) = af_do_ref
!
!                else
!                    cell_flags(IJK) = af_rm_ref
!                end if
!            end do; CLOSE_DO

!            nc = box%n_cell
!            do KJI_DO(1, nc)
!                if( box%lvl < max_refinement .and.                                       &
!                &  abs(box%cc(i, j, i_df+i_obj-1)) <= 30.0_R8 * box%dr(1) ) then
!                    cell_flags(IJK) = af_do_ref
!                else
!                    cell_flags(IJK) = af_rm_ref
!                end if
!            end do; CLOSE_DO

            nc = box%n_cell
            do KJI_DO(1, nc)
                if( box%lvl < max_refinement .and.                                       &
                &  abs(box%cc(i, j, i_df+i_obj-1)) <= 0.25_R8 ) then
                    cell_flags(IJK) = af_do_ref
                else
                    cell_flags(IJK) = af_keep_ref
                end if
            end do; CLOSE_DO



!            nc = box%n_cell
!            cntr(:) = box%r_min(:) + (box%n_cell/2) * box%dr(:)
!            do KJI_DO(1, nc)
!                if( box%lvl < max_refinement .and.                                       &
!                & ( abs(cntr(1) - 1.d0) < 0.3_R8  .or.                                   &
!                & (cntr(2) < 0.9_R8 .and. cntr(2) > 0.0_R8))                             &
!                ) then
!                    cell_flags(IJK) = af_do_ref
!                else
!                    cell_flags(IJK) = af_keep_ref
!                end if
!            end do; CLOSE_DO

        end subroutine ref_objects
    end subroutine refine_field



    subroutine ab_exp_smac_timestep_ker (box)
        implicit none
        type(box_t), intent(inout) :: box
        real(kind=R8) :: res_x_old(DTIMES(box%n_cell))
        real(kind=R8) :: res_y_old(DTIMES(box%n_cell))
        real(kind=R8) :: res_x_older(DTIMES(box%n_cell))
        real(kind=R8) :: res_y_older(DTIMES(box%n_cell))
        integer(kind=I4):: nc, i, j
        real(kind=R8) :: inv_re, inv_rho, diff_coef

        inv_re = (1.0_R8 / rey)
        diff_coef = inv_re * nu_f
        inv_rho = 1.0_R8 / rho_f
        nc = box%n_cell

        res_x_old(:, :)   = - (1.5_dp) * p_convective (box, [i_uold, i_f])               &
                            + (1.5_dp  * diff_coef) * p_diffusion (box, i_uold)          &
                            - inv_rho  * p_cell_xgrad (box, i_p)
        res_y_old(:, :)   = - (1.5_dp) * p_convective (box, [i_vold, i_f])               &
                            + (1.5_dp  * diff_coef) * p_diffusion  (box, i_vold)         &
                            - inv_rho  * p_cell_ygrad (box, i_p)

        res_x_older(:, :) = + (0.5_dp) * p_convective (box, [i_uolder, i_folder])        &
                            - (0.5_dp  * diff_coef) * p_diffusion (box, i_uolder)
        res_y_older(:, :) = + (0.5_dp) * p_convective (box, [i_volder, i_folder])        &
                            - (0.5_dp  * diff_coef) * p_diffusion (box, i_volder)
        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_resx) = res_x_old(i, j) + res_x_older(i, j)
                box%cc(i, j, i_resy) = res_y_old(i, j) + res_y_older(i, j)

                box%cc(i, j, i_unew) = box%cc(i, j, i_uold) + (dt * box%cc(i, j, i_resx))
                box%cc(i, j, i_vnew) = box%cc(i, j, i_vold) + (dt * box%cc(i, j, i_resy))
            end do
        end do
    end subroutine ab_exp_smac_timestep_ker



    subroutine central2order_smac_set_hrhs_withforce_ker (box)
        type(box_t), intent(inout) :: box
        integer(kind=I4) :: i, j, nc
        real(kind=R8) :: res_x_old(DTIMES(box%n_cell))
        real(kind=R8) :: res_y_old(DTIMES(box%n_cell))
        real(kind=R8) :: res_x_older(DTIMES(box%n_cell))
        real(kind=R8) :: res_y_older(DTIMES(box%n_cell))
        real(kind=R8) :: inv_helmholtz_lambda
        real(kind=R8) :: inv_re, inv_rho, diff_coef

        inv_re = (1.0_R8 / rey)
        diff_coef = inv_re * nu_f
        inv_rho = 1.0_R8 / rho_f
        nc = box%n_cell

! FIXME (mehdi#1#): Why does the - sign  for the hemholtz_lamda work?
        inv_helmholtz_lambda = -2.0_R8 * rey * (1.0_R8 / nu_f) / dt

        res_x_old(:, :)   = -(1.5_dp) * p_convective (box, [i_uold, i_f])                &
!                            +(1.5_dp) * p_diffusion  (box, i_uold)                      &
                            +(0.5_dp*diff_coef) * p_diffusion  (box, i_uold)             &
                                      - inv_rho * p_cell_xgrad (box, i_p)
        res_y_old(:, :)   = -(1.5_dp) * p_convective (box, [i_vold, i_f])                &
!                            +(1.5_dp) * p_diffusion  (box, i_vold)                      &
                            +(0.5_dp*diff_coef) * p_diffusion  (box, i_vold)             &
                                      - inv_rho * p_cell_ygrad (box, i_p)

        res_x_older(:, :) =  (0.5_dp) * p_convective (box, [i_uolder, i_folder])!&
!                            -(0.5_dp) * p_diffusion  (box, i_uolder)
        res_y_older(:, :) =  (0.5_dp) * p_convective (box, [i_volder, i_folder])!&
!                            -(0.5_dp) * p_diffusion  (box, i_volder)

        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_hrhsx) =  inv_helmholtz_lambda * (                        &
                                           dt * (                                        &
                                                + res_x_old(i, j)                        &
                                                + res_x_older(i, j)                      &
                                                + box%cc(i, j, i_ibfx)                   &
                                                )                                        &
                                            + box%cc(i, j, i_uold)                       &
                                                                )

                box%cc(i, j, i_hrhsy) =  inv_helmholtz_lambda * (                        &
                                           dt * (                                        &
                                                + res_y_old(i, j)                        &
                                                + res_y_older(i, j)                      &
                                                + box%cc(i, j, i_ibfy)                   &
                                                )                                        &
                                            + box%cc(i, j, i_vold)                       &
                                                                )
                end do
            end do
    end subroutine central2order_smac_set_hrhs_withforce_ker



    subroutine smac_update_velocity_ker_2o (box)
        type(box_t), intent(inout) :: box
        real(kind=R8) :: d_px_c(DTIMES(box%n_cell)),   d_py_c(DTIMES(box%n_cell))
        real(kind=R8) :: d_px_f(DTIMES(box%n_cell+1)), d_py_f(DTIMES(box%n_cell+1))
        integer(kind=I4) :: nc, i, j, nf

        nc = box%n_cell
        nf = nc + 1

        d_px_c(:, :) = p_cell_xgrad (box, i_phi)
        d_py_c(:, :) = p_cell_ygrad (box, i_phi)
        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_unew) = box%cc(i, j, i_unew) - (dt * d_px_c(i, j))
                box%cc(i, j, i_vnew) = box%cc(i, j, i_vnew) - (dt * d_py_c(i, j))
            end do
        end do

        d_px_f(:, :) = p_face_xgrad (box, i_phi)
        d_py_f(:, :) = p_face_ygrad (box, i_phi)
        do j = 1, nf
            do i = 1, nf
                box%fc(i, j, 1, i_f) = box%fc(i, j, 1, i_f) - (dt * d_px_f(i, j))
                box%fc(i, j, 2, i_f) = box%fc(i, j, 2, i_f) - (dt * d_py_f(i, j))
            end do
        end do
    end subroutine smac_update_velocity_ker_2o



    subroutine smac_update_pressure_ker_2o (box)
        implicit none
        type(box_t), intent(inout) :: box
        integer(kind=I4) :: nc, i, j

        nc = box%n_cell
        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_p) = box%cc(i, j, i_p) + box%cc(i, j, i_phi)
            end do
        end do
    end subroutine smac_update_pressure_ker_2o



    function central_2order_cell_xgrad (box, i_var) result (grad)
        class(box_t), intent(in) :: box
        integer(kind=I4), intent(in) :: i_var
        real(kind=R8) :: grad(DTIMES(box%n_cell))
        real(kind=R8) :: inv_2dr(NDIM)
        integer(kind=I4) :: nc, i, j

        nc = box%n_cell
        inv_2dr(:) = 1.0_dp / (2.0_dp * box%dr(:))

        do j = 1, nc
            do i = 1, nc
                grad(i, j) = inv_2dr(1) * (box%cc(i+1, j, i_var) - box%cc(i-1, j, i_var))
            end do
        end do
    end function central_2order_cell_xgrad



    function central_2order_cell_ygrad (box, i_var) result (grad)
        class(box_t), intent(in) :: box
        integer(kind=I4), intent(in) :: i_var
        real(kind=R8) :: grad(DTIMES(box%n_cell))
        real(kind=R8) :: inv_2dr(NDIM)
        integer(kind=I4) :: nc, i, j

        nc = box%n_cell
        inv_2dr(:) = 1.0_dp / (2.0_dp * box%dr(:))

        do j = 1, nc
            do i = 1, nc
                grad(i, j) = inv_2dr(2) * (box%cc(i, j+1, i_var) - box%cc(i, j-1, i_var))
            end do
        end do
    end function central_2order_cell_ygrad



    function central_2order_face_xgrad (box, i_var) result (grad)
        class(box_t), intent(in) :: box
        integer(kind=4), intent(in) :: i_var
        real(kind=R8) :: grad(DTIMES(box%n_cell+1))
        real(kind=R8) :: inv_dr(NDIM)
        integer(kind=I4) :: nf, i, j

        nf = box%n_cell + 1
        inv_dr(:) = 1.0_dp / (box%dr(:))

        do j = 1, nf
            do i = 1, nf
                grad(i, j) = inv_dr(1) * (box%cc(i, j, i_var) - box%cc(i-1, j, i_var))
            end do
        end do
    end function central_2order_face_xgrad



    function central_2order_face_ygrad (box, i_var) result (grad)
        class(box_t), intent(in) :: box
        integer(kind=I4), intent(in) :: i_var
        real(kind=R8) :: grad(DTIMES(box%n_cell+1))
        real(kind=R8) :: inv_dr(NDIM)
        integer(kind=I4) :: nf, i, j

        nf = box%n_cell + 1
        inv_dr(:) = 1.0_dp / (box%dr(:))

        do j = 1, nf
            do i = 1, nf
                grad(i, j) = inv_dr(2) * (box%cc(i, j, i_var) - box%cc(i, j-1, i_var))
            end do
        end do
    end function central_2order_face_ygrad



    function central_2order_diffusion (box, i_var) result (diff)
        class(box_t), intent(in) :: box
        integer(kind=I4), intent(in)  :: i_var
        real(kind=R8) :: diff(DTIMES(box%n_cell))
        integer(kind=I4)  :: nc, i, j
        real(kind=R8) :: inv_dr(NDIM)
!        real(kind=R8) :: inv_re

        nc = box%n_cell
        inv_dr(:)   = 1.0_dp / box%dr(:)

        do j = 1, nc
            do i = 1, nc
                diff (i, j) =          (                                                 &
                &  inv_dr(1) * inv_dr(1) * (                                             &
                &  - 2.0_dp * box%cc(i, j, i_var)                                        &
                &  + box%cc(i+1, j, i_var) + box%cc(i-1, j, i_var)                       &
                &                          )                                             &
                &                  +                                                     &
                &  inv_dr(2) * inv_dr(2) * (                                             &
                &  - 2.0_dp * box%cc(i, j, i_var)                                        &
                &  + box%cc(i, j+1, i_var) + box%cc(i, j-1, i_var)                       &
                &                          )                                             &
                &                         )
            end do
        end do
    end function central_2order_diffusion



    function central_2order_convective (box, i_var)  result(advec)
        class(box_t), intent(in) :: box
        integer(kind=I4), intent(in) :: i_var(:)
        real(kind=R8) :: advec(DTIMES(box%n_cell))
        integer(kind=I4) :: nc, i, j
        real(kind=R8) :: inv_dr(NDIM)

        nc = box%n_cell
        inv_dr(:) = 1.0_dp / box%dr(:)

        do j = 1, nc
            do i = 1, nc
                advec(i, j) =                                                            &
                !  d(U*(i_var)) / dx
                inv_dr(1) *  (                                                           &
                &  ( ( 0.5d0 * (box%cc(i, j, i_var(1)) + box%cc(i+1, j, i_var(1))) )     &
                &  * box%fc(i+1, j, 1, i_var(2)) ) -                                     &
                &  ( ( 0.5d0 * (box%cc(i, j, i_var(1)) + box%cc(i-1, j, i_var(1))) )     &
                &  * box%fc(i,   j, 1, i_var(2)) )                                       &
                &           ) +                                                          &
                !  d(V*(i_var)) / dy
                &  inv_dr(2) * (                                                         &
                &  ( (0.5d0 * (box%cc(i, j, i_var(1)) + box%cc(i, j+1, i_var(1))) )      &
                &  * box%fc(i, j+1, 2, i_var(2)) ) -                                     &
                &  ( (0.5d0 * (box%cc(i, j, i_var(1)) + box%cc(i, j-1, i_var(1))) )      &
                &  * box%fc(i, j,   2, i_var(2)) )                                       &
                    )
            end do
        end do
    end function central_2order_convective



    subroutine central2order_c2f_ker(box, i_var)
        type(box_t), intent(inout) :: box
        integer(kind=I4), intent(in) :: i_var(:)
        integer(kind=I4) :: i, j, nf

        nf = box%n_cell + 1
        do i = 1, nf
            do j = 1, nf
                box%fc(i, j, 1, i_f) = 0.5_dp *                                          &
                &  (box%cc(i-1, j, i_var(1)) + box%cc(i, j, i_var(1)))
                box%fc(i, j, 2, i_f) = 0.5_dp *                                          &
                &  (box%cc(i, j-1, i_var(2)) + box%cc(i, j, i_var(2)))
            end do
        end do
    end subroutine central2order_c2f_ker



    subroutine central2order_c2f_p_decoupled_ker (box, i_var)
        type(box_t), intent(inout) :: box
        integer(kind=I4), intent(in) :: i_var(:)
        integer(kind=I4) :: i, j, nf
        real(kind=R8) :: lo(NDIM), hi(NDIM), inv_2dx(NDIM), inv_dx(NDIM)

        nf = box%n_cell + 1
        inv_2dx = 1.0_dp / (2.0_dp * box%dr(:))
        inv_dx  = 1.0_dp / ( box%dr(:))

! FIXME (mehdi#1#): The loop should be reduced to some manipulation for layer adjacent to the bondary
        do i = 1, nf
            do j = 1, nf
                box%fc(i, j, 1, i_f) = 0.5_dp *                                          &
                &  (box%cc(i-1, j, i_var(1)) + box%cc(i, j, i_var(1)))
                box%fc(i, j, 2, i_f) = 0.5_dp *                                          &
                &  (box%cc(i, j-1, i_var(2)) + box%cc(i, j, i_var(2)))
            end do
        end do

        do i = 2, nf-1
            do j = 2, nf-1
                lo(1) = box%cc(i-1, j, i_var(1))                                         &
                &  + (dt * inv_2dx(1) * (box%cc(i,   j, i_p) - box%cc(i-2, j, i_p)) )
                hi(1) = box%cc(i  , j, i_var(1))                                         &
                &  + (dt * inv_2dx(1) * (box%cc(i+1, j, i_p) - box%cc(i-1, j, i_p)) )

                lo(2) = box%cc(i, j-1, i_var(2))                                         &
                &  + (dt * inv_2dx(2) * (box%cc(i,   j, i_p) - box%cc(i, j-2, i_p)) )
                hi(2) = box%cc(i  , j, i_var(2))                                         &
                &  + (dt * inv_2dx(2) * (box%cc(i, j+1, i_p) - box%cc(i, j-1, i_p)) )

                box%fc(i, j, 1, i_f) = 0.5_dp * ( lo(1) + hi(1) )                        &
                &  - (dt * inv_dx(1) * (box%cc(i,   j, i_p) - box%cc(i-1, j, i_p)) )
                box%fc(i, j, 2, i_f) = 0.5_dp * ( lo(2) + hi(2) )                        &
                &  - (dt * inv_dx(1) * (box%cc(i,   j, i_p) - box%cc(i, j-1, i_p)) )
            end do
        end do
    end subroutine central2order_c2f_p_decoupled_ker



    subroutine central2order_set_prhs_ker (box)
        type(box_t), intent(inout) :: box
        integer(kind=I4) :: i, j, nc
        real(kind=R8) :: inv_dr(NDIM), inv_dt

        nc = box%n_cell
        inv_dr(:) = 0.5_dp / box%dr(:)
        inv_dt    = 1.0_dp / dt

        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_prhs) = (inv_dt) * (                                      &
                &  ( inv_dr(1) * (box%fc(i+1, j, 1, i_f) - box%fc(i, j, 1, i_f))) +      &
                &  ( inv_dr(2) * (box%fc(i, j+1, 2, i_f) - box%fc(i, j, 2, i_f)))        &
                                                 )
            end do
        end do

    end subroutine central2order_set_prhs_ker



    subroutine ab_simp_solver (tree, refine_info)

        implicit none
        class(af_t), intent(inout) :: tree
        class(ref_info_t), intent(inout) :: refine_info
        integer(kind=I4) :: itr
        character(kind=CK, len=100) :: fname
        character(kind=CK, len=100) :: time_scheme

        p_update_v_ker => smac_update_velocity_ker_2o
        p_update_p_ker => smac_update_pressure_ker_2o
        p_cell_xgrad   => central_2order_cell_xgrad
        p_cell_ygrad   => central_2order_cell_ygrad
        p_face_xgrad   => central_2order_face_xgrad
        p_face_ygrad   => central_2order_face_ygrad
        p_diffusion    => central_2order_diffusion
        p_convective   => central_2order_convective
        p_c2f_ker      => central2order_c2f_ker
        p_set_prhs_ker => central2order_set_prhs_ker
        p_set_hrhs_ker => central2order_smac_set_hrhs_ker

        call module_config%get (section_name='solver_time', option_name='scheme',        &
        &  val=time_scheme)

        field_time = 0.0_dp

        !> Fills ghost cells and boundaries
        call af_gc_tree (tree, [i_uold, i_vold])

        !> Interpolates velocities from cell to face
        call af_loop_box_iarg (tree, p_c2f_ker, [i_uold, i_vold])

        do itr = itr_min, itr_max
            field_time = dt * (itr)

            !> Boundary condition
            call af_gc_tree (tree, [i_uold, i_vold])
            call af_gc_tree (tree, [i_p])

            call af_loop_box (tree, p_set_hrhs_ker)

            !> Time step
            call get_helmholtz (tree)

            !> boundary condition
            call af_gc_tree (tree, [i_unew, i_vnew])

            !> rough flux
            call af_loop_box_iarg (tree, p_c2f_ker, [i_unew, i_vnew])

            !> calculates the corresponding pressure for a divergence-free projection
            call af_loop_box (tree, p_set_prhs_ker)
            call p_get_projector (tree)

            !> After projection fills the boundary for phi and p, then restricts
            !  those  into coarser levels.
            call af_gc_tree (tree, [i_p, i_phi])
            call af_restrict_tree (tree, [i_p, i_phi])

            !> updates the velocity into a divergence-free fields
            call af_loop_box (tree, p_update_v_ker)

            !> Updates the pressure based on the result of the projection
            call af_loop_box (tree, p_update_p_ker)

            !> Prints out
            if (mod(itr, hard_print_frq) == 0) then
                write(fname, "(A,I0)") "ns_" // DIMNAME // "_", itr

                !> Write the cell centered data of tree to a vtk unstructured file.
                !> Only the leaves of the tree are used
!                call af_write_silo(tree,  &
!                    trim(fname), itr, field_time, dir="output",                          &
!                    ixs_cc = [i_unew, i_vnew, i_p, i_df, i_hrhsx, i_hrhsy])

!                  call af_write_vtk(tree, trim(fname))

                  call af_write_vtk(tree,                                              &
                    &  trim(fname), itr, field_time, dir="output",                       &
                    &  ixs_cc =                                                          &
                    &  [i_unew, i_vnew, i_p, i_df, i_vortz])


            end if

            call af_tree_copy_cc(tree, i_uold, i_uolder)
            call af_tree_copy_cc(tree, i_vold, i_volder)

            call af_tree_copy_cc(tree, i_unew, i_uold)
            call af_tree_copy_cc(tree, i_vnew, i_vold)

            call af_tree_copy_fc(tree, i_f, i_folder)
        end do
    end subroutine ab_simp_solver



    subroutine central2order_smac_set_hrhs_ker (box)
        implicit none
        type(box_t), intent(inout) :: box
        integer(kind=I4):: i, j, nc
        real(kind=R8) :: inv_dr(NDIM), inv_2dr(NDIM)
        real(kind=R8) :: res_x_old(DTIMES(box%n_cell))
        real(kind=R8) :: res_y_old(DTIMES(box%n_cell))
        real(kind=R8) :: res_x_older(DTIMES(box%n_cell))
        real(kind=R8) :: res_y_older(DTIMES(box%n_cell))
        real(kind=R8) :: inv_helmholtz_lambda
        real(kind=R8) :: inv_re, inv_rho, diff_coef

        inv_re = (1.0_R8 / rey)
        diff_coef = inv_re * nu_f
        inv_rho = 1.0_R8 / rho_f

        nc = box%n_cell
        inv_dr(:)  = 0.5_dp / box%dr(:)
        inv_2dr(:) = 1.0_dp / (2.0_dp * box%dr(:))


! FIXME (mehdi#1#): Why does the - sign  for the hemholtz_lamda work?
!        inv_helmholtz_lambda  =  ( (2.0_dp * rey) / dt )
        inv_helmholtz_lambda  =  - 2.0_dp * rey * (1.0_R8 / nu_f) / dt

        res_x_old(:, :)   = -(1.5_dp) * p_convective (box, [i_uold, i_f])                &
        &                   +(0.5_dp*diff_coef) * p_diffusion  (box, i_uold)             &
                                      - inv_rho * p_cell_xgrad (box, i_p)
        res_y_old(:, :)   = -(1.5_dp) * p_convective (box, [i_vold, i_f])                &
        &                   +(0.5_dp*diff_coef) * p_diffusion  (box, i_vold)             &
        &                             - inv_rho * p_cell_ygrad (box, i_p)

        res_x_older(:, :) =  (0.5_dp) * p_convective (box, [i_uolder, i_folder])!&
                            !-(1.0_dp / 2.0_dp) * p_diffusion  (box, i_uolder)
        res_y_older(:, :) =  (0.5_dp) * p_convective (box, [i_volder, i_folder])!&
                            !-(1.0_dp / 2.0_dp) * p_diffusion  (box, i_volder)

        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_hrhsx) =  inv_helmholtz_lambda * (                        &
                &  + (dt * res_x_old(i, j))                                              &
                &   + (dt * res_x_older(i, j))                                           &
                &   + box%cc(i, j, i_uold)                                               &
                &        )

                box%cc(i, j, i_hrhsy) =  inv_helmholtz_lambda * (                        &
                        +  (dt * res_y_old(i, j))                                        &
                        +  (dt * res_y_older(i, j))                                      &
                        +  box%cc(i, j, i_vold)                                          &
                        )
            end do
        end do
    end subroutine central2order_smac_set_hrhs_ker



    subroutine central2order_vorticity_ker (box)
        type (box_t), intent(inout) :: box
        integer(kind=I4) :: i, j, nc
        real(kind=R8) :: inv2dr(NDIM)

        inv2dr(:) = 1.0_dp / (2.0_dp * box%dr(:))
        nc = box%n_cell
        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_vortz) =                                                  &
                inv2dr(1) * (box%cc(i+1, j, i_vnew) - box%cc(i-1, j, i_vnew)) -          &
                inv2dr(2) * (box%cc(i, j+1, i_unew) - box%cc(i, j-1, i_unew))
            end do
        end do
    end subroutine central2order_vorticity_ker



    subroutine central2order_divegence_ker (box, i_var)
        implicit none
        type(box_t), intent(inout) :: box
        integer(kind=I4), intent(in) :: i_var(:)
        integer :: i, j, nc
        real(kind=dp) :: inv2dr(NDIM)

        inv2dr(:) = 1.0_dp / (2.0_dp * box%dr(:))
        nc = box%n_cell

        do j = 1, nc
            do i = 1, nc
                box%cc(i, j, i_var(1)) =                                                 &
                &  inv2dr(1) * (box%cc(i+1, j, i_var(2)) - box%cc(i-1, j, i_var(2)) ) +  &
                &  inv2dr(2) * (box%cc(i, j+1, i_var(3)) - box%cc(i, j-1, i_var(3)) )
            end do
        end do
    end subroutine central2order_divegence_ker



    subroutine get_cfl_c (tree, cfl_c, dt)
        use face
        class(af_t), intent(inout) :: tree
        real(kind=R8), intent(inout) :: cfl_c
        real(kind=R8), intent(in)    :: dt
        real(kind=R8) :: umax, vmax
        character(kind=CK, len=:), allocatable :: error_message

        call af_tree_max_cc(tree, i_unew, umax)
        call af_tree_max_cc(tree, i_vnew, vmax)
        cfl_c = dt * ( (umax / af_min_dr(tree)) + (vmax / af_min_dr(tree) ) )

        if (cfl_c >= cfl_max) then ! cfl_max is a module variable
            write (*, *)
            error_message= colorize('STOP! too large CFL', color_fg='red')
            write(*, '(A)') error_message
            write (*, *)
            STOP
        end if
    end subroutine

end module m_navier_stokes_solver
