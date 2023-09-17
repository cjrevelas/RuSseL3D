!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine InitArrays()
!-----------------------------------------------------------------------------------------------------------!
use geometry_mod,    only: numNodes
use kcw_mod,         only: rdiag1
use parser_vars_mod, only: matrixExist, numEdwPointsMatrix, numConvolPointsMatrix, &
                           graftedExist, numEdwPointsGrafted, numConvolPointsGrafted
use arrays_mod
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
allocate(rdiag1(numNodes))
allocate(d2phi_dr2(numNodes))
allocate(wwField(numNodes),wwFieldMixed(numNodes),wwFieldNew(numNodes),uuField(numNodes))
allocate(phiMatrix(numNodes),phiGrafted(numNodes),phiTotal(numNodes))
allocate(nodeVolume(numNodes))
allocate(dsEdwMatrix(numEdwPointsMatrix+1))
allocate(xsEdwMatrix(numEdwPointsMatrix+1))
allocate(coeffEdwMatrix(numEdwPointsMatrix+1))
allocate(qqMatrix(2,numNodes))
allocate(qqMatrixFinal(numEdwPointsMatrix+1,numNodes))
allocate(qqMatrixInterpGrafted(numConvolPointsGrafted+1,numNodes))
allocate(qqGraftedInterp(numConvolPointsGrafted+1,numNodes))

if (matrixExist.eq.1) then
  allocate(qqMatrixInterp(numConvolPointsMatrix+1,numNodes))
  allocate(dsConvMatrix(numConvolPointsMatrix+1))
  allocate(xsConvMatrix(numConvolPointsMatrix+1))
  allocate(coeffConvMatrix(numConvolPointsMatrix+1))

  qqMatrixInterp  = 0.0d0
  dsConvMatrix    = 0.0d0
  xsConvMatrix    = 0.0d0
  coeffConvMatrix = 0.0d0
endif

nodeVolume     = 0.0d0
d2phi_dr2      = 0.0d0
wwField        = 0.0d0
wwFieldMixed   = 0.0d0
wwFieldNew     = 0.0d0
uuField        = 0.0d0
rdiag1         = 0.0d0
phiMatrix      = 0.0d0
phiGrafted     = 0.0d0
phiTotal       = 0.0d0
dsEdwMatrix    = 0.0d0
xsEdwMatrix    = 0.0d0
coeffEdwMatrix = 0.0d0

qqMatrix              = 0.0d0
qqMatrixFinal         = 0.0d0
qqMatrixInterpGrafted = 0.0d0

if (graftedExist.eq.1) then
  allocate(dsEdwGrafted(numEdwPointsGrafted+1),dsConvGrafted(numConvolPointsGrafted+1))
  allocate(xsEdwGrafted(numEdwPointsGrafted+1),xsConvGrafted(numConvolPointsGrafted+1))
  allocate(coeffEdwGrafted(numEdwPointsGrafted+1),coeffConvGrafted(numConvolPointsGrafted+1))
  allocate(qqGrafted(2,numNodes))
  allocate(qqGraftedFinal(numEdwPointsGrafted+1,numNodes))

  dsEdwGrafted     = 0.0d0
  dsConvGrafted    = 0.0d0
  xsEdwGrafted     = 0.0d0
  xsConvGrafted    = 0.0d0
  coeffEdwGrafted  = 0.0d0
  coeffConvGrafted = 0.0d0
  qqGrafted        = 0.0d0
  qqGraftedFinal   = 0.0d0
  qqGraftedInterp  = 0.0d0
endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine InitArrays
