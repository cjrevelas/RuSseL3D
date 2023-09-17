!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module arrays_mod
!----------------------------------------------------------------------------------------------------------------------------!
real(8), allocatable, dimension(:)   :: dsEdwMatrix, dsEdwGrafted, dsConvMatrix, dsConvGrafted
real(8), allocatable, dimension(:)   :: coeffEdwMatrix, coeffEdwGrafted, coeffConvMatrix, coeffConvGrafted
real(8), allocatable, dimension(:)   :: xsEdwMatrix, xsConvMatrix, xsEdwGrafted, xsConvGrafted
real(8), allocatable, dimension(:)   :: wwField, wwFieldNew, wwFieldMixed, uuField
real(8), allocatable, dimension(:)   :: phiMatrix, phiGrafted, phiTotal
real(8), allocatable, dimension(:)   :: nodeVolume
real(8), allocatable, dimension(:,:) :: phiGraftedIndiv
real(8), allocatable, dimension(:)   :: d2phi_dr2
real(8), allocatable, dimension(:,:) :: qqMatrix, qqMatrixFinal, qqMatrixInterp, qqMatrixInterpGrafted
real(8), allocatable, dimension(:,:) :: qqGrafted, qqGraftedFinal, qqGraftedInterp
!----------------------------------------------------------------------------------------------------------------------------!
end module arrays_mod
