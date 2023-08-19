!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module eos_mod
!----------------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod
use constants_mod
use flags_mod
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
integer :: eosType

real(8) :: helfandCompressibility
real(8) :: tempStar, pressStar, volStar, rhoStar
real(8) :: tempTilde, pressTilde, rhoTildeBulk, rslN
!----------------------------------------------------------------------------------------------------------------------------!
  contains
    function eos_ff(phi) result(ff)
      implicit none
      real(8), intent(in) :: phi
      real(8)             :: ff, rhoTilde

      ff = 0.0d0
      if (eosType.eq.eosHelfand) then
        ff = 5.0d-1 / helfandCompressibility * (phi - 1.0d0)**2.0d0
      elseif (eosType.eq.eosSanchezLacombe) then
        rhoTilde = rhoTildeBulk*phi
        if (rhoTilde.gt.1.0d0) rhoTilde = 0.96d0

        ff = pressStar*(tempTilde*rhoTilde - rhoTilde**2.0d0 + tempTilde*(1.0d0-rhoTilde)*LOG(1.0d0-rhoTilde))
      endif

      return
    end function eos_ff

    function eos_df_drho(phi) result(df_drho)
      implicit none
      real(8), intent(in) :: phi
      real(8)             :: df_drho, rhoTilde

      df_drho = 0.0d0
      if (eosType.eq.eosHelfand) then
        df_drho = (phi-1.0d0)/(helfandCompressibility*molarBulkDensity*n_avog)
      elseif (eosType.eq.eosSanchezLacombe) then
        rhoTilde = rhoTildeBulk*phi
        if (rhoTilde.gt.1.0d0) rhoTilde = 0.96d0

        df_drho = boltz_const_Joule_K*tempStar*(-2.0d0*rslN*rhoTilde - tempTilde*rslN*LOG(1.0d0-rhoTilde))
      endif

      return
    end function eos_df_drho

    function eosRhoTildeZero(tempTilde, pressTilde, rsl) result(rhoTildeZero)
      implicit none
      real(8), intent(in) :: tempTilde, pressTilde, rsl
      real(8)             :: rhoTildeZero, rhoTildeZero_init=8.0d-1, rhoTildeZeroNew=0.0d0
      real(8)             :: func=0.0d0, funcDeriv=0.0d0
      real(8)             :: tolerance=1.0d-9, errorNorm=2.0d1

      integer :: iter, maxNumIterations=1000

      rhoTildeZero = rhoTildeZero_init

      iter = 0
      do while((iter.lt.maxNumIterations).and.(errorNorm.gt.tolerance))
        iter = iter + 1

        func            = rhoTildeZero**2.0d0 + pressTilde + tempTilde*(LOG(1.0d0-rhoTildeZero) + (1.0d0-1.0d0/rsl)*rhoTildeZero)
        funcDeriv       = 2.0d0*rhoTildeZero + tempTilde*((1.0d0-1.0d0/rsl)-1.0d0/(1.0d0-rhoTildeZero))
        rhoTildeZeroNew = rhoTildeZero - func / funcDeriv
        errorNorm       = ABS(rhoTildeZeroNew-rhoTildeZero)

        if (rhoTildeZeroNew.lt.1.0d0) then
          rhoTildeZero = rhoTildeZeroNew
        else
          rhoTildeZeroNew = 9.5d-1
          rhoTildeZero     = rhoTildeZeroNew
        endif

        write(6,'(I5,2(2X,F16.9),2X,E16.9)') iter, rhoTildeZero, rhoTildeZeroNew, errorNorm
      enddo

      if (iter.eq.maxNumIterations) write(6,'(A40)') "slvle solver diverged.."
    end function eosRhoTildeZero
!----------------------------------------------------------------------------------------------------------------------------!
end module eos_mod
