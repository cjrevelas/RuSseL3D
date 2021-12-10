!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module eos
!----------------------------------------------------------------------------------------------------------!
use parser_vars
use constants
use flags
!----------------------------------------------------------------------------------------------------------!
integer :: eos_type

real(8) :: hlf_kappa_T
real(8) :: T_star, P_star, V_star, rho_star
real(8) :: T_tilde, P_tilde, rho_tilde_bulk, rsl_N
!----------------------------------------------------------------------------------------------------------!
contains
    function eos_ff(phi) result(ff)
    implicit none
    real(8), intent(in) :: phi
    real(8)             :: ff, rho_tilde

    if (eos_type.eq.eos_helfand) then
        ff = 0.5d0 / hlf_kappa_T * (phi - 1.d0)**2
    elseif (eos_type.eq.eos_sl) then
        rho_tilde = rho_tilde_bulk*phi
        if (rho_tilde.gt.1.d0) then
           rho_tilde = 0.9999d0
        endif
        ff = P_star*(T_tilde*rho_tilde - rho_tilde**2 + T_tilde*(1-rho_tilde)*LOG(1-rho_tilde))
    endif
    end function eos_ff

    function eos_df_drho(phi) result(df_drho)
    implicit none
    real(8), intent(in) :: phi
    real(8)             :: df_drho, rho_tilde

    if (eos_type.eq.eos_helfand) then
        df_drho = (phi-1.d0)/(hlf_kappa_T*rho_mol_bulk*n_avog)
    elseif (eos_type.eq.eos_sl) then
        rho_tilde = rho_tilde_bulk*phi
        if (rho_tilde.gt.1.d0) then
           rho_tilde = 0.9999d0
        endif
        df_drho = boltz_const_Joule_K*T_star*(-2.d0*rsl_N*rho_tilde - T_tilde*rsl_N*LOG(1.d0-rho_tilde))
    endif
    end function eos_df_drho

    function eos_rho_tilde_0(T_tilde, P_tilde, rsl) result(rho_tilde_0)
    implicit none
    real(8), intent(in) :: T_tilde, P_tilde, rsl
    real(8)             :: rho_tilde_0, rho_tilde_0_init=0.80d0, rho_tilde_0_new=0.d0
    real(8)             :: func=0.d0, func_deriv=0.d0
    real(8)             :: tolerance=1.d-09, err_norm=2.d01

    integer :: ii, max_iter=1000

    rho_tilde_0 = rho_tilde_0_init

    ii = 0
    do while((ii.lt.max_iter).and.(err_norm.gt.tolerance))
        ii = ii + 1

        func            = rho_tilde_0**2 + P_tilde + T_tilde*(LOG(1.d0-rho_tilde_0) + (1.d0-1.d0/rsl)*rho_tilde_0)
        func_deriv      = 2.d0*rho_tilde_0 + T_tilde*((1.d0-1.d0/rsl)-1.d0/(1.d0-rho_tilde_0))
        rho_tilde_0_new = rho_tilde_0 - func / func_deriv
        err_norm        = ABS(rho_tilde_0_new-rho_tilde_0)

        if (rho_tilde_0_new.lt.1.d0) then
            rho_tilde_0 = rho_tilde_0_new
        else
            rho_tilde_0_new = 0.95d0
            rho_tilde_0     = rho_tilde_0_new
        endif

        write(6,'(I5,2(2X,F16.9),2X,E16.9)') ii, rho_tilde_0, rho_tilde_0_new, err_norm
    enddo
    if (ii.eq.max_iter) write(6,'(A40)') "slvle solver diverged.."
    end function eos_rho_tilde_0
!----------------------------------------------------------------------------------------------------------!
end module eos
