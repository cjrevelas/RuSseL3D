module xdata
!-------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------!
!***********************************PARAMETERS**********************************!
integer, parameter :: maxns = 4000

!*********************************DISCRETIZATION********************************!
integer :: ns, nx_1, np
real(8) :: ds, chainlen

!*******************************INDICES AND FLAGS*******************************!
integer :: iow, ioe, count
integer :: iseed
integer :: i1, j1, k1, k2, m1, m2, m3, id
integer :: ii1, ii2, irow
integer :: kk, iter, iterations
integer :: readfield, code_name, lshow, pr_on

!*************************************MESH**************************************!
integer :: nel, ndm, numnp, numel
integer, allocatable, dimension(:,:) :: ix
real(8), allocatable, dimension(:,:) :: xc
real(8), dimension(3) :: box_lo, box_hi, box_len

!*******************************BOUNDARY CONDITIONS*****************************!
integer :: n_dirichlet_faces
integer, allocatable, dimension(:) :: ids_dirichlet_faces

!*************************SCF MODEL AND POTENTIAL DATA**************************!
real(8) :: mon_mass, massden, kapa, kappa_T, rho_0, Rgyr, Rg_2, CN, sphere_radius
real(8) :: Temp, volume, Aps, Asio2, sigma1, sigma2, Vored, wwidth, distance
real(8) :: part_func, adh_ten, nch_per_area, free_energy_old, free_energy_new
real(8) :: exact

!*************************SCF CONVERGENCE SCHEME********************************!
real(8) :: max_error_tol, max_error, std_error, fraction, mix_coef_frac, mix_coef_kapa
integer :: scheme_type

!*********************************MUMPS OPTIONS*********************************!
integer :: mumps_matrix_type

!*****************************SIMPSON COEFFICIENTS******************************!
real(8), dimension(1:maxns) :: koeff

!*****************************PERIODIC PROFILES*********************************!
integer :: print_ev

!*****************************AUXILIARY VARIABLES*******************************!
real(8) :: start, D, t1, t2, t3, t4, t5, t6
real(8) :: test1, test2, test3, test4
real(8) :: temp1, temp2, sum_final
real(8) :: alpha, beta, coef
real(8) :: sum, term1, term2, term3, term4, term5
real(8) :: part_sum1, part_sum2, part_sum3, part_sum4

!*****************************SCF ALLOCATABLE ARRAYS****************************!
real(8), allocatable, dimension(:)   :: wa, wa_new, wa_mix, Ufield, phia_new, phi
real(8), allocatable, dimension(:,:) :: qf, qf_final
logical, allocatable, dimension(:)   :: elem_in_q0_face
!-------------------------------------------------------------------------------!
end module xdata
