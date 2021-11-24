module xdata
!--------------------------------------------------------------------!  
implicit none
!--------------------------------------------------------------------!  
integer, parameter :: maxns = 4000

integer :: ns, chainlen 
integer :: ior, iow, count, count2, init 
integer :: iseed, sl, slparams, ns_middle
integer :: i1, j1, k1, k2, m1, m2, m3
integer :: ii1, ii2, nonzero, irow, nx_1, np
integer :: kk, iterations
integer :: time_step, show, code_name, lshow, pr_on

real(8) :: lx, dx, ds, diff_number(1)
real(8) :: start, t1, t2, t3, t4, t5, t6
real(8) :: volume, D, fraction
real(8) :: test1, test2, test3, test4, sphere_radius
real(8) :: temp1, temp2, sum_final
real(8) :: alpha, beta
real(8) :: sum, term1, term2, term3, term4, term5
real(8) :: part_sum1, part_sum2, part_sum3, part_sum4
real(8) :: free_energy_old, free_energy_new
real(8) :: koeff(1:maxns), coef
real(8) :: mon_mass, massden, kapa, kappa_T, rho_0, Rgyr, Rg_2, CN
real(8) :: Temp, Aps, Asio2, sigma1, sigma2, Vored, wwidth, distance
real(8) :: part_func, adh_ten, adh_ten_alt, nch_per_area
real(8) :: exact, max_error, error

double precision, allocatable, dimension(:) :: wa, wa_new, wa_in, Ufield, phia_new, phi
double precision, allocatable, dimension(:,:) :: qf, qf_final
!--------------------------------------------------------------------!      
end module xdata 
