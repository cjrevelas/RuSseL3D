!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module flags_mod
    integer, parameter :: contour_uniform = 1
    integer, parameter :: contour_symm    = 2
    integer, parameter :: contour_asymm   = 3
    integer, parameter :: contour_hybrid  = 4
    integer, parameter :: sys_m           = 1
    integer, parameter :: sys_mg          = 2
    integer, parameter :: eos_helfand     = 0
    integer, parameter :: eos_sl          = 1
    integer, parameter :: eos_sl_grad     = 2
end module flags_mod
