module xdata
     ! nx      = number of intervals used in discretizing distance in the
!           x spatial direction, normal to surface.  
!           Note: In Daoulas et al. 2005, x is symbolized as z.
     !  integer nx 
	   

! maxnx   = maximum value of nx allowed 
       integer maxnx

! ns      = number of intervals used in discretizing s, the dimensionless
!           coordinate along the chain contour
       integer ns

! maxns   = maximum value of ns allowed
       integer maxns

! lx      = width, normal to surface, of the spatial domain where
!           SCF equations are solved.  Measured in units of the
!           chain radius of gyration.       
       real*8  lx

! dx      = interval along x used in the discretized representation
!           of the Edwards diffusion equation.  Measured in units of
!           the chain radius of gyration.
       real*8  dx        

! chainlen= chain length, in monomer units (methylenes/methyls)
	   integer chainlen

       parameter (maxnx = 1000) 
       parameter (maxns = 4000)
	   

!On polymer chain discretization:

!  Example:
!                                                             ns
! 0___1___2___3___4___5___6___7___8___9___10___11___12___13___14
! 0                                                           1  (s axis)
!   Make sure, since the integrations will be performed by
!   Simpson method, that ns is an even number.
    !ds          = interval along the dimensionless chain contour

!diff_number = ds R_g^2/[2 (dx)^2], a dimensionless
!combination of the intervals
!dx along the spatial coordinate normal to the surface and
!ds along the dimensionless contour length of the chain
!appearing in the discretization of the Edwards diffusion
!equation
!      common  /general1/   diff_number,ds

      real*8               diff_number(1),ds


!wa          = self-consistent field beta N W
!      common  /fields/     wa, wa_new

      real*8               wa(1:maxnx)
      real*8               wa_new(1:maxnx)


!phia        = volume fraction profile
!      common /Phis/        phia_new

      real*8               phia_new(1:maxnx)
	 
      


       real*8      qf(1:maxnx,1:2)
      
       real*8      qf_final(1:maxnx,1:maxns)

      real*8 volume, D, fraction
      real*8 test1,test2,test3,test4,sphere_radius
      real*8 temp1, temp2, Rg_2
      real*8 adiag(1:maxnx),bdiag(1:maxnx)
      real*8 cdiag(1:maxnx),rdiag(1:maxnx)
      real*8 alpha, beta, phi(1:maxnx)
     
      real*8 wa_in(1:maxnx)
      real*8 u(1:maxnx)
      real*8 sum_final
      real*8 exact,max_error
      real*8 part_func
      real*8 sum,term1,term2,term3,term4,term5
      real*8 free_energy_old, mon_mass
      real*8 free_energy_new
	  real*8 koeff(1:maxns)
	  real*8 Ufield(1:maxnx)
      real*8 kapa, Rgyr, Vored, wwidth
      real*8 distance!, U_solid!,surf_pot
    

!***********************************************
      
	  real*8 Temp
	  real*8 CN
	  real*8 massden,Aps,Asio2,sigma1,sigma2
	  real*8 rho_0, kappa_T
	  real*8 nch_per_area
	  real*8 adh_ten,adh_ten_alt
	  real*8 error
      
!     
      integer ior,iow
      integer count,count2,init
      integer iseed,sl,slparams,ns_middle!,k
      integer i1,j1,k1,m1,m2,m3,k2!,n
      integer ii1,ii2,nonzero
	  integer irow,nx_1,np
      integer kk, iterations
      integer time_step,show,code_name,lshow,pr_on
     
      
      real*8 k_B,Kgr,K_boltzmann, N_avog,pi,coef
      parameter(k_B=8.3144598)

!in jK^-1mol
      parameter(K_boltzmann=1.38064852E-23)
!in jK^-1
      parameter(N_avog=6.02214085700E+23)
      parameter(pi= 3.14159265359)
    
    end 