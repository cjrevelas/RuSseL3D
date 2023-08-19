!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module flags_mod
!----------------------------------------------------------------------------------------------------------------------------!
integer, parameter :: contourUniform    = 1
integer, parameter :: contourSymmetric  = 2
integer, parameter :: contourAsymmetric = 3
integer, parameter :: contourHybrid     = 4

integer, parameter :: systemMatrix        = 1
integer, parameter :: systemMatrixGrafted = 2

integer, parameter :: eosHelfand            = 0
integer, parameter :: eosSanchezLacombe     = 1
integer, parameter :: eosSanchezLacombeGrad = 2

integer, parameter :: mumpsAsymmetric       = 0
integer, parameter :: mumpsPositiveDefinite = 1
integer, parameter :: mumpsGeneralSymmetric = 2
!----------------------------------------------------------------------------------------------------------------------------!
end module flags_mod
