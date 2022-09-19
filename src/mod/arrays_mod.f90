!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module arrays_mod
!----------------------------------------------------------------------------------------------------------------------------!
real(8), allocatable, dimension(:)   :: ds_mx_ed, ds_gr_ed, ds_mx_conv, ds_gr_conv
real(8), allocatable, dimension(:)   :: coeff_mx_ed, coeff_gr_ed, coeff_mx_conv, coeff_gr_conv
real(8), allocatable, dimension(:)   :: xs_mx_ed, xs_mx_conv, xs_gr_ed, xs_gr_conv
real(8), allocatable, dimension(:)   :: wwField, wwFieldNew, wwFieldMixed, uuField
real(8), allocatable, dimension(:)   :: phiMatrix, phiGrafted, phiTotal
real(8), allocatable, dimension(:)   :: nodeVolume
real(8), allocatable, dimension(:,:) :: phiGraftedIndiv
real(8), allocatable, dimension(:)   :: d2phi_dr2
real(8), allocatable, dimension(:,:) :: qqMatrix, qqMatrixFinal, qqMatrixInterp, qqMatrixInterpGrafted
real(8), allocatable, dimension(:,:) :: qqGrafted, qqGraftedFinal, qqGraftedInterp
!----------------------------------------------------------------------------------------------------------------------------!
end module arrays_mod
