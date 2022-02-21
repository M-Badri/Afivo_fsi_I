#include "cpp_macros.h"
module m_solver_parameters
    use iso_fortran_env, only: I4=>int32, R8=>real64, character_kinds
    use m_af_all

    implicit none
    integer(kind=I4), parameter :: CK = character_kinds(1)
    integer(kind=I4), parameter :: n_stn = 6, n_phy = 1
    integer(kind=I4), parameter :: forcing = 2, inside = -1, outside = 1, ghost = -2
    integer(kind=I4)            :: i_uolder, i_uold, i_unew
    integer(kind=I4)            :: i_volder, i_vold, i_vnew
    integer(kind=I4)            :: i_vortz
    integer(kind=I4)            :: i_f, i_folder
    integer(kind=I4)            :: i_p, i_phi, i_prhs, i_tmp, i_tmp2, i_tmp3
    integer(kind=I4)            :: i_df
    integer(kind=I4)            :: i_ibfx, i_ibfy
    integer(kind=I4)            :: i_resx, i_resy
    integer(kind=I4)            :: i_hrhsx, i_hrhsy
    integer(kind=I4)            :: i_hphix, i_hphiy
    integer(kind=I4)            :: max_refinement
    real(kind=R8),    parameter :: g = 981.0_R8
    real(kind=R8),    parameter :: pi = acos(-1.0_R8)
    real(kind=R8),    parameter :: max_distance = 1000.0_R8
    real(kind=R8)               :: field_time, dt
    real(kind=R8)               :: nu_f, rho_f, rey
    real(kind=R8)               :: domain_len(NDIM)
    character(kind=CK, len=100) :: test_case
    logical                     :: is_nondimesional
end module m_solver_parameters
