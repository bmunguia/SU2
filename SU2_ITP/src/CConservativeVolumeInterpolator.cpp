/*!
 * \file CConservativeVolumeInterpolator.cpp
 * \brief Implementation of the conservative solution interpolation subroutines.
 *        This file contains the conservative interpolation logic using Alauzet 2015 method.
 * \author B. Munguía
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include <queue>

#include "../include/CConservativeVolumeInterpolator.hpp"
#include "../../Common/include/fem/fem_standard_element.hpp"

CConservativeVolumeInterpolator::CConservativeVolumeInterpolator(SU2_Comm MPICommunicator)
    : CVolumeInterpolator(MPICommunicator) { }

void CConservativeVolumeInterpolator::Interpolate(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                  CSolver** solver_container_src, CSolver** solver_container_dst) {
  if (rank == MASTER_NODE) {
    cout << endl << "----------------------------- Interpolation -----------------------------" << endl;
    cout << "Performing conservative solution interpolation from source mesh to destination mesh..." << endl;
    cout << "Source mesh: " << geometry_src->GetGlobal_nPointDomain() << " points, ";
    cout << geometry_src->GetGlobal_nElemDomain() << " elements" << endl;
    cout << "Destination mesh: " << geometry_dst->GetGlobal_nPointDomain() << " points" << endl;
  }

  /*--- Build the ADTs ---*/
  InitializeADTs(config, geometry_src, geometry_dst);

  /*--- Call the conservative interpolation method ---*/
  for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {
    auto solver_src = solver_container_src[iSol];
    auto solver_dst = solver_container_dst[iSol];
    if (solver_src && solver_dst)
      ConservativeInterpolation(config, geometry_src, geometry_dst, solver_src, solver_dst);
  }
}

void CConservativeVolumeInterpolator::ConservativeInterpolation(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                                CSolver* solver_src, CSolver* solver_dst) {
  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Apply the curvature correction to the destination nodes    ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Applying curvature correction." << endl;
  vector<su2double> coorDst;
  vector<su2double> coorDstCorrected;
  for (unsigned long iPoint = 0; iPoint < nPoint_dst; iPoint++) {
    for (unsigned short k = 0; k < nDim; ++k) {
      coorDst.push_back(geometry_dst->nodes->GetCoord(iPoint, k));
    }
  }
  ApplyCurvatureCorrection(config, geometry_src, geometry_dst, nDim, coorDst, coorDstCorrected);

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Localize destination nodes on the source mesh              ---*/
  /*---         containingElems is a map from destination mesh nodes to    ---*/
  /*---         containing elements on the source mesh, and pointsFailed   ---*/ 
  /*---         is all the nodes for which no containing element was found ---*/
  /*--------------------------------------------------------------------------*/
  vector<unsigned long> containingElems;
  vector<int> containingElemRanks;
  vector<unsigned long> pointsFailed;
  PointLocalization(geometry_src, coorDstCorrected, containingElems, containingElemRanks, pointsFailed);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute solution mass and gradient on source mesh          ---*/
  /*--------------------------------------------------------------------------*/
  vector<vector<su2double> > srcElemMass(nElem_src, vector<su2double>(nVar, 0.0));
  vector<vector<su2double> > srcElemGrad(nElem_src, vector<su2double>(nVar * nDim, 0.0));
  ComputeSolutionMass(geometry_src, solver_src, srcElemMass, srcElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Compute the intersection of elements K_dst with elements   ---*/
  /*---         K_src it overlaps                                          ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing element intersections." << endl;
  map<unsigned long, vector<unsigned long>> overlappingElements;
  ComputeOverlappingElements(geometry_src, geometry_dst, containingElems, overlappingElements);

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Mesh the intersection polygon/polyhedron of each pair      ---*/ 
  /*---         (K_dst, K_src_i)                                           ---*/
  /*--------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------*/
  /*--- Step 6: Compute destination mesh mass and gradient using Gauss     ---*/
  /*---         quadrature                                                 ---*/
  /*--------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------*/
  /*--- Step 7: Correct the gradient to enforce the maximum principle      ---*/
  /*--------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------*/
  /*--- Step 8: Perform averaging to get solution at vertices.             ---*/
  /*--------------------------------------------------------------------------*/
}

void CConservativeVolumeInterpolator::PointLocalization(CGeometry* geometry_src, 
                                                        const vector<su2double>& coor_corrected, 
                                                        vector<unsigned long>& containingElems,
                                                        vector<int>& containingElemRanks, 
                                                        vector<unsigned long>& pointsFailed) {
  /*--- Search for containing elements for the given coordinates ---*/
  CADTElemClass& volumeADT = GetSourceVolumeADT();

  const unsigned long nDOFsDst = coor_corrected.size() / nDim;

  /*--- Loop over the DOFs to be interpolated ---*/
  containingElems.clear();
  containingElemRanks.clear();
  pointsFailed.clear();

  /*--- Initialize with invalid values ---*/
  containingElems.resize(nDOFsDst, ULONG_MAX);
  containingElemRanks.resize(nDOFsDst, -1);

  for (unsigned long l = 0; l < nDOFsDst; ++l) {
    /*--- Set a pointer to the coordinates to be searched ---*/
    const su2double* coor = coor_corrected.data() + l * nDim;

    /*--- Carry out the containment search and check if it was successful ---*/
    unsigned short subElemID;
    unsigned long elemID;
    int rankID;
    su2double parCoor[3], weightsInterpol[8];

    bool foundElement = volumeADT.DetermineContainingElement(coor, subElemID, elemID, rankID, parCoor, weightsInterpol);

    if (foundElement) {
      /*--- Store element information ---*/
      containingElems[l] = elemID;
      containingElemRanks[l] = rankID;
    } else {
      /*--- Containment search failed - store the index ---*/
      pointsFailed.push_back(l);
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Volume search finished. " << pointsFailed.size() << " points failed." << endl << flush;
  }
}

void CConservativeVolumeInterpolator::ComputeSolutionMass(CGeometry* geometry, 
                                                          CSolver* solver,
                                                          vector<vector<su2double> >& elemMass,
                                                          vector<vector<su2double> >& elemGrad) {

  const unsigned short nVar = solver_src->GetnVar();
  for (unsigned long elemID = 0; elemID < geometry->GetnElem(); ++elemID) {
      auto* elem = geometry->elem[elemID];
      unsigned short nNodes = elem->GetnNodes();
      unsigned short VTK_Type = elem->GetVTK_Type();

      /*--- Get solution at nodes ---*/
      vector<vector<su2double> > solAtNodes(nNodes, vector<su2double>(nVar, 0.0));
      vector<su2double> nodeCoords(nNodes * nDim, 0.0);

      for (unsigned short iNode = 0; iNode < nNodes; ++iNode) {
          unsigned long nodeID = elem->GetNode(iNode);
          for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
              solAtNodes[iNode][iVar] = solver->GetNodes()->GetSolution(nodeID, iVar);
          }
          for (unsigned short k = 0; k < nDim; ++k) {
              nodeCoords[iNode * nDim + k] = geometry->nodes->GetCoord(nodeID, k);
          }
      }

      /*--- Create a standard element for integration ---*/
      /*--- For FVM, we typically use linear elements (nPoly = 1) ---*/
      unsigned short nPoly = 1;
      bool constJac = false;
      unsigned short orderExact = 2 * nPoly;

      CFEMStandardElement stdElement(VTK_Type, nPoly, constJac, nullptr, orderExact);

      /*--- Get integration points and weights ---*/
      const unsigned short nInt = stdElement.GetNIntegration();
      const su2double* weights = stdElement.GetWeightsIntegration();

      /*--- Get basis functions and derivatives at integration points ---*/
      const su2double* lagBasis = stdElement.GetBasisFunctionsIntegration();
      const su2double* drLagBasis = stdElement.GetDrBasisFunctionsIntegration();
      const su2double* dsLagBasis = stdElement.GetDsBasisFunctionsIntegration();
      const su2double* dtLagBasis = (nDim == 3) ? stdElement.GetDtBasisFunctionsIntegration() : nullptr;

      /*--- Initialize mass and gradient ---*/
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
          elemMass[elemID][iVar] = 0.0;
          for (unsigned short k = 0; k < nDim; ++k) {
              elemGrad[elemID][iVar * nDim + k] = 0.0;
          }
      }

      /*--- Loop over integration points ---*/
      for (unsigned short iInt = 0; iInt < nInt; ++iInt) {

          /*--- Compute Jacobian of transformation ---*/
          su2double dxdr = 0.0, dydr = 0.0, dzdr = 0.0;
          su2double dxds = 0.0, dyds = 0.0, dzds = 0.0;
          su2double dxdt = 0.0, dydt = 0.0, dzdt = 0.0;

          for (unsigned short iNode = 0; iNode < nNodes; ++iNode) {
              unsigned short ind = iInt * nNodes + iNode;

              dxdr += nodeCoords[iNode * nDim + 0] * drLagBasis[ind];
              dydr += nodeCoords[iNode * nDim + 1] * drLagBasis[ind];
              if (nDim == 3) dzdr += nodeCoords[iNode * nDim + 2] * drLagBasis[ind];

              if (dsLagBasis) {
                  dxds += nodeCoords[iNode * nDim + 0] * dsLagBasis[ind];
                  dyds += nodeCoords[iNode * nDim + 1] * dsLagBasis[ind];
                  if (nDim == 3) dzds += nodeCoords[iNode * nDim + 2] * dsLagBasis[ind];
              }

              if (dtLagBasis) {
                  dxdt += nodeCoords[iNode * nDim + 0] * dtLagBasis[ind];
                  dydt += nodeCoords[iNode * nDim + 1] * dtLagBasis[ind];
                  dzdt += nodeCoords[iNode * nDim + 2] * dtLagBasis[ind];
              }
          }

          /*--- Compute Jacobian determinant ---*/
          su2double jacobian;
          if (nDim == 2) {
              jacobian = dxdr * dyds - dydr * dxds;
          } else {
              jacobian = dxdr * (dyds * dzdt - dzds * dydt) -
                        dydr * (dxds * dzdt - dzds * dxdt) +
                        dzdr * (dxds * dydt - dyds * dxdt);
          }
          jacobian = abs(jacobian);

          /*--- Compute inverse Jacobian for gradient transformation ---*/
          su2double drdx, drdy, drdz = 0.0;
          su2double dsdx, dsdy, dsdz = 0.0;
          su2double dtdx = 0.0, dtdy = 0.0, dtdz = 0.0;

          if (nDim == 2) {
              su2double jacInv = 1.0 / jacobian;
              drdx = dyds * jacInv;
              drdy = -dxds * jacInv;
              dsdx = -dydr * jacInv;
              dsdy = dxdr * jacInv;
          } else {
              su2double jacInv = 1.0 / jacobian;
              drdx = (dyds * dzdt - dzds * dydt) * jacInv;
              drdy = (dzds * dxdt - dxds * dzdt) * jacInv;
              drdz = (dxds * dydt - dyds * dxdt) * jacInv;
              dsdx = (dzdr * dydt - dydr * dzdt) * jacInv;
              dsdy = (dxdr * dzdt - dzdr * dxdt) * jacInv;
              dsdz = (dydr * dxdt - dxdr * dydt) * jacInv;
              dtdx = (dydr * dzds - dzdr * dyds) * jacInv;
              dtdy = (dzdr * dxds - dxdr * dzds) * jacInv;
              dtdz = (dxdr * dyds - dydr * dxds) * jacInv;
          }

          /*--- Integration weight including Jacobian ---*/
          su2double intWeight = weights[iInt] * jacobian;

          /*--- Loop over variables ---*/
          for (unsigned short iVar = 0; iVar < nVar; ++iVar) {

              /*--- Interpolate solution and its gradient at integration point ---*/
              su2double solVal = 0.0;
              su2double dudr = 0.0, duds = 0.0, dudt = 0.0;

              for (unsigned short iNode = 0; iNode < nNodes; ++iNode) {
                  unsigned short ind = iInt * nNodes + iNode;
                  su2double nodeVal = solAtNodes[iNode][iVar];

                  solVal += lagBasis[ind] * nodeVal;
                  dudr += drLagBasis[ind] * nodeVal;
                  if (dsLagBasis) duds += dsLagBasis[ind] * nodeVal;
                  if (dtLagBasis) dudt += dtLagBasis[ind] * nodeVal;
              }

              /*--- Add contribution to mass integral ---*/
              elemMass[elemID][iVar] += solVal * intWeight;

              /*--- Transform gradient to physical coordinates and add to gradient integral ---*/
              su2double dudx = dudr * drdx + duds * dsdx + dudt * dtdx;
              su2double dudy = dudr * drdy + duds * dsdy + dudt * dtdy;
              su2double dudz = 0.0;
              if (nDim == 3) dudz = dudr * drdz + duds * dsdz + dudt * dtdz;

              /*--- For FVM, gradient is constant per element, so we can average ---*/
              elemGrad[elemID][iVar * nDim + 0] += dudx * intWeight;
              elemGrad[elemID][iVar * nDim + 1] += dudy * intWeight;
              if (nDim == 3) elemGrad[elemID][iVar * nDim + 2] += dudz * intWeight;
          }
      }

      /*--- For FVM with constant gradients, normalize by element volume ---*/
      su2double elemVolume = elem->GetVolume();
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
          for (unsigned short k = 0; k < nDim; ++k) {
              elemGrad[elemID][iVar * nDim + k] /= elemVolume;
          }
      }
  }
}

void CConservativeVolumeInterpolator::ComputeOverlappingElements(CGeometry* geometry_src, 
                                                                 CGeometry* geometry_dst,
                                                                 const vector<unsigned long>& containingElems,
                                                                 map<unsigned long, vector<unsigned long>>& overlappingElements) {
  overlappingElements.clear();

  /*--- Only handle triangular elements for now ---*/
  if (nDim != 2) {
    if (rank == MASTER_NODE) {
      cout << "Warning: Element intersection currently only implemented for 2D. ";
      cout << "Skipping intersection computation for mesh." << endl;
    }
    return;
  }

  unsigned long totalOverlaps = 0;

  /*--- Loop over all destination elements ---*/
  for (unsigned long dstElemID = 0; dstElemID < geometry_dst->GetnElem(); ++dstElemID) {
    auto* dstElem = geometry_dst->elem[dstElemID];

    /*--- Skip non-triangular elements ---*/
    if (dstElem->GetVTK_Type() != TRIANGLE) {
      continue;
    }

    /*--- Step 1: Build initial list from elements K_src containing vertices of K_dst ---*/
    set<unsigned long> candidateList;

    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      /*--- Find corresponding point index in corrected coordinates ---*/
      if (nodeID < containingElems.size() && containingElems[nodeID] != ULONG_MAX) {
        candidateList.insert(containingElems[nodeID]);

        /*--- TODO: Handle degenerate cases - add neighbors for vertices on edges/vertices ---*/
        /*--- For now, just add the containing element ---*/
      }
    }

    /*--- Get destination triangle vertices ---*/
    su2double dstTri[6];
    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      dstTri[iNode * 2 + 0] = geometry_dst->nodes->GetCoord(nodeID, 0);
      dstTri[iNode * 2 + 1] = geometry_dst->nodes->GetCoord(nodeID, 1);
    }

    /*--- Step 2: Process candidate list, adding new elements during intersection ---*/
    set<unsigned long> processedElems;
    queue<unsigned long> toProcess;

    /*--- Initialize queue with initial candidate list ---*/
    for (unsigned long srcElemID : candidateList) {
      toProcess.push(srcElemID);
    }

    /*--- Process elements, potentially adding new ones during intersection ---*/
    while (!toProcess.empty()) {
      unsigned long srcElemID = toProcess.front();
      toProcess.pop();

      /*--- Skip if already processed ---*/
      if (processedElems.count(srcElemID) > 0) continue;
      processedElems.insert(srcElemID);

      auto* srcElem = geometry_src->elem[srcElemID];
      if (srcElem->GetVTK_Type() != TRIANGLE) continue;

      /*--- Get source triangle vertices ---*/
      su2double srcTri[6];
      for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        unsigned long nodeID = srcElem->GetNode(iNode);
        srcTri[iNode * 2 + 0] = geometry_src->nodes->GetCoord(nodeID, 0);
        srcTri[iNode * 2 + 1] = geometry_src->nodes->GetCoord(nodeID, 1);
      }

      /*--- Check for intersection and detect new candidates ---*/
      vector<su2double> intersectionPoints;
      set<unsigned long> newCandidates;

      if (TriangleTriangleIntersection(geometry_src, srcTri, dstTri, srcElemID,
                                       intersectionPoints, newCandidates)) {
        overlappingElements[dstElemID].push_back(srcElemID);
        totalOverlaps++;

        /*--- Add new candidates to processing queue ---*/
        for (unsigned long newElemID : newCandidates) {
          if (processedElems.count(newElemID) == 0) {
            toProcess.push(newElemID);
          }
        }
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Found " << totalOverlaps << " total overlapping element pairs." << endl;
    cout << "Number of destination elements with overlaps: " << overlappingElements.size() << endl;
  }
}

bool CConservativeVolumeInterpolator::TriangleTriangleIntersection(CGeometry* geometry_src, 
                                                                   const su2double dstTri[6], 
                                                                   const su2double srcTri[6],
                                                                   unsigned long srcElemID,
                                                                   vector<su2double>& intersectionPoints,
                                                                   set<unsigned long>& newCandidates) {
  intersectionPoints.clear();
  newCandidates.clear();

  const su2double EPS = 1e-12;

  /*--- Triangle 0 (destination) vertices ---*/
  su2double t0v0[2] = {dstTri[0], dstTri[1]};
  su2double t0v1[2] = {dstTri[2], dstTri[3]};
  su2double t0v2[2] = {dstTri[4], dstTri[5]};

  /*--- Triangle 1 (source) vertices ---*/
  su2double t1v0[2] = {srcTri[0], srcTri[1]};
  su2double t1v1[2] = {srcTri[2], srcTri[3]};
  su2double t1v2[2] = {srcTri[4], srcTri[5]};

  /*--- Step 1: Compute vertex powers (signed distances) ---*/
  vector<su2double> cloudPoints;

  /*--- Powers of destination vertices w.r.t. source edges ---*/
  su2double power_t0v0_e01 = ComputeSignedDistance(t0v0, t1v0, t1v1);
  su2double power_t0v0_e12 = ComputeSignedDistance(t0v0, t1v1, t1v2);
  su2double power_t0v0_e20 = ComputeSignedDistance(t0v0, t1v2, t1v0);

  su2double power_t0v1_e01 = ComputeSignedDistance(t0v1, t1v0, t1v1);
  su2double power_t0v1_e12 = ComputeSignedDistance(t0v1, t1v1, t1v2);
  su2double power_t0v1_e20 = ComputeSignedDistance(t0v1, t1v2, t1v0);

  su2double power_t0v2_e01 = ComputeSignedDistance(t0v2, t1v0, t1v1);
  su2double power_t0v2_e12 = ComputeSignedDistance(t0v2, t1v1, t1v2);
  su2double power_t0v2_e20 = ComputeSignedDistance(t0v2, t1v2, t1v0);

  /*--- Powers of source vertices w.r.t. destination edges ---*/
  su2double power_t1v0_e01 = ComputeSignedDistance(t1v0, t0v0, t0v1);
  su2double power_t1v0_e12 = ComputeSignedDistance(t1v0, t0v1, t0v2);
  su2double power_t1v0_e20 = ComputeSignedDistance(t1v0, t0v2, t0v0);

  su2double power_t1v1_e01 = ComputeSignedDistance(t1v1, t0v0, t0v1);
  su2double power_t1v1_e12 = ComputeSignedDistance(t1v1, t0v1, t0v2);
  su2double power_t1v1_e20 = ComputeSignedDistance(t1v1, t0v2, t0v0);

  su2double power_t1v2_e01 = ComputeSignedDistance(t1v2, t0v0, t0v1);
  su2double power_t1v2_e12 = ComputeSignedDistance(t1v2, t0v1, t0v2);
  su2double power_t1v2_e20 = ComputeSignedDistance(t1v2, t0v2, t0v0);

  /*--- Check if destination vertices are inside source triangle ---*/
  if (power_t0v0_e01 >= -EPS && power_t0v0_e12 >= -EPS && power_t0v0_e20 >= -EPS) {
    cloudPoints.push_back(t0v0[0]); cloudPoints.push_back(t0v0[1]);
  }
  if (power_t0v1_e01 >= -EPS && power_t0v1_e12 >= -EPS && power_t0v1_e20 >= -EPS) {
    cloudPoints.push_back(t0v1[0]); cloudPoints.push_back(t0v1[1]);
  }
  if (power_t0v2_e01 >= -EPS && power_t0v2_e12 >= -EPS && power_t0v2_e20 >= -EPS) {
    cloudPoints.push_back(t0v2[0]); cloudPoints.push_back(t0v2[1]);
  }

  /*--- Check if source vertices are inside destination triangle ---*/
  /*--- If so, add their vertex balls to candidates ---*/
  if (power_t1v0_e01 >= -EPS && power_t1v0_e12 >= -EPS && power_t1v0_e20 >= -EPS) {
    cloudPoints.push_back(t1v0[0]); cloudPoints.push_back(t1v0[1]);
    AddVertexBallToCandidates(geometry_src, srcElemID, 0, newCandidates);
  }
  if (power_t1v1_e01 >= -EPS && power_t1v1_e12 >= -EPS && power_t1v1_e20 >= -EPS) {
    cloudPoints.push_back(t1v1[0]); cloudPoints.push_back(t1v1[1]);
    AddVertexBallToCandidates(geometry_src, srcElemID, 1, newCandidates);
  }
  if (power_t1v2_e01 >= -EPS && power_t1v2_e12 >= -EPS && power_t1v2_e20 >= -EPS) {
    cloudPoints.push_back(t1v2[0]); cloudPoints.push_back(t1v2[1]);
    AddVertexBallToCandidates(geometry_src, srcElemID, 2, newCandidates);
  }

  /*--- Add vertices that lie on edges (degenerate cases) ---*/
  if (abs(power_t0v0_e01) < EPS) { cloudPoints.push_back(t0v0[0]); cloudPoints.push_back(t0v0[1]); }
  if (abs(power_t0v0_e12) < EPS) { cloudPoints.push_back(t0v0[0]); cloudPoints.push_back(t0v0[1]); }
  if (abs(power_t0v0_e20) < EPS) { cloudPoints.push_back(t0v0[0]); cloudPoints.push_back(t0v0[1]); }

  if (abs(power_t0v1_e01) < EPS) { cloudPoints.push_back(t0v1[0]); cloudPoints.push_back(t0v1[1]); }
  if (abs(power_t0v1_e12) < EPS) { cloudPoints.push_back(t0v1[0]); cloudPoints.push_back(t0v1[1]); }
  if (abs(power_t0v1_e20) < EPS) { cloudPoints.push_back(t0v1[0]); cloudPoints.push_back(t0v1[1]); }

  if (abs(power_t0v2_e01) < EPS) { cloudPoints.push_back(t0v2[0]); cloudPoints.push_back(t0v2[1]); }
  if (abs(power_t0v2_e12) < EPS) { cloudPoints.push_back(t0v2[0]); cloudPoints.push_back(t0v2[1]); }
  if (abs(power_t0v2_e20) < EPS) { cloudPoints.push_back(t0v2[0]); cloudPoints.push_back(t0v2[1]); }

  if (abs(power_t1v0_e01) < EPS) { cloudPoints.push_back(t1v0[0]); cloudPoints.push_back(t1v0[1]); }
  if (abs(power_t1v0_e12) < EPS) { cloudPoints.push_back(t1v0[0]); cloudPoints.push_back(t1v0[1]); }
  if (abs(power_t1v0_e20) < EPS) { cloudPoints.push_back(t1v0[0]); cloudPoints.push_back(t1v0[1]); }

  if (abs(power_t1v1_e01) < EPS) { cloudPoints.push_back(t1v1[0]); cloudPoints.push_back(t1v1[1]); }
  if (abs(power_t1v1_e12) < EPS) { cloudPoints.push_back(t1v1[0]); cloudPoints.push_back(t1v1[1]); }
  if (abs(power_t1v1_e20) < EPS) { cloudPoints.push_back(t1v1[0]); cloudPoints.push_back(t1v1[1]); }

  if (abs(power_t1v2_e01) < EPS) { cloudPoints.push_back(t1v2[0]); cloudPoints.push_back(t1v2[1]); }
  if (abs(power_t1v2_e12) < EPS) { cloudPoints.push_back(t1v2[0]); cloudPoints.push_back(t1v2[1]); }
  if (abs(power_t1v2_e20) < EPS) { cloudPoints.push_back(t1v2[0]); cloudPoints.push_back(t1v2[1]); }

  /*--- Step 2: Check edge-edge intersections and add neighbors when intersected ---*/
  vector<su2double> edgeIntersections;

  /*--- Check all destination edges vs source edges ---*/
  if (LineSegmentIntersection(t0v0, t0v1, t1v0, t1v1, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 0, newCandidates); // edge 0-1
  }
  if (LineSegmentIntersection(t0v0, t0v1, t1v1, t1v2, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 1, newCandidates); // edge 1-2
  }
  if (LineSegmentIntersection(t0v0, t0v1, t1v2, t1v0, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 2, newCandidates); // edge 2-0
  }

  if (LineSegmentIntersection(t0v1, t0v2, t1v0, t1v1, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 0, newCandidates); // edge 0-1
  }
  if (LineSegmentIntersection(t0v1, t0v2, t1v1, t1v2, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 1, newCandidates); // edge 1-2
  }
  if (LineSegmentIntersection(t0v1, t0v2, t1v2, t1v0, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 2, newCandidates); // edge 2-0
  }

  if (LineSegmentIntersection(t0v2, t0v0, t1v0, t1v1, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 0, newCandidates); // edge 0-1
  }
  if (LineSegmentIntersection(t0v2, t0v0, t1v1, t1v2, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 1, newCandidates); // edge 1-2
  }
  if (LineSegmentIntersection(t0v2, t0v0, t1v2, t1v0, edgeIntersections)) {
    for (size_t i = 0; i < edgeIntersections.size(); ++i)
      cloudPoints.push_back(edgeIntersections[i]);
    AddEdgeNeighborToCandidates(geometry_src, srcElemID, 2, newCandidates); // edge 2-0
  }

  /*--- Remove duplicate points ---*/
  for (unsigned int i = 0; i < cloudPoints.size(); i += 2) {
    for (unsigned int j = i + 2; j < cloudPoints.size(); j += 2) {
      if (abs(cloudPoints[i] - cloudPoints[j]) < EPS &&
          abs(cloudPoints[i+1] - cloudPoints[j+1]) < EPS) {
        cloudPoints.erase(cloudPoints.begin() + j, cloudPoints.begin() + j + 2);
        j -= 2;
      }
    }
  }

  /*--- Check intersection result ---*/
  unsigned int numIntersection = cloudPoints.size() / 2;

  if (numIntersection < 3) {
    /*--- No intersection or degenerate case ---*/
    intersectionPoints = cloudPoints;
    return numIntersection > 0;
  }

  /*--- Step 3: Compute convex hull ---*/
  ComputeConvexHull(cloudPoints);
  intersectionPoints = cloudPoints;

  return intersectionPoints.size() >= 6; // At least 3 points for a valid intersection polygon
}

su2double CConservativeVolumeInterpolator::ComputeSignedDistance(const su2double point[2],
                                                                 const su2double lineStart[2],
                                                                 const su2double lineEnd[2]) {
  /*--- Compute signed distance using cross product ---*/
  su2double dx = lineEnd[0] - lineStart[0];
  su2double dy = lineEnd[1] - lineStart[1];
  su2double px = point[0] - lineStart[0];
  su2double py = point[1] - lineStart[1];

  /*--- Cross product gives twice the signed area ---*/
  return dx * py - dy * px;
}

bool CConservativeVolumeInterpolator::LineSegmentIntersection(const su2double P0[2], const su2double P1[2],
                                                              const su2double Q0[2], const su2double Q1[2],
                                                              vector<su2double>& intersections) {
  const su2double EPS = 1e-12;
  intersections.clear();

  /*--- Following Alauzet method: Let e_P = [P0 P1] and e_Q = [Q0 Q1] be two edges ---*/

  /*--- Compute signed distances ---*/
  su2double distP0_eQ = ComputeSignedDistance(P0, Q0, Q1);
  su2double distP1_eQ = ComputeSignedDistance(P1, Q0, Q1);
  su2double distQ0_eP = ComputeSignedDistance(Q0, P0, P1);
  su2double distQ1_eP = ComputeSignedDistance(Q1, P0, P1);

  /*--- Count zero powers ---*/
  bool zeroP0 = abs(distP0_eQ) < EPS;
  bool zeroP1 = abs(distP1_eQ) < EPS;
  bool zeroQ0 = abs(distQ0_eP) < EPS;
  bool zeroQ1 = abs(distQ1_eP) < EPS;

  int zeroCount = int(zeroP0) + int(zeroP1) + int(zeroQ0) + int(zeroQ1);

  /*--- Handle degenerate cases ---*/
  if (zeroCount > 0) {
    if (zeroCount == 1) {
      /*--- Only one power is zero ---*/
      if (zeroP0 && (distQ0_eP * distQ1_eP < 0)) {
        intersections.push_back(P0[0]);
        intersections.push_back(P0[1]);
        return true;
      } else if (zeroP1 && (distQ0_eP * distQ1_eP < 0)) {
        intersections.push_back(P1[0]);
        intersections.push_back(P1[1]);
        return true;
      } else if (zeroQ0 && (distP0_eQ * distP1_eQ < 0)) {
        intersections.push_back(Q0[0]);
        intersections.push_back(Q0[1]);
        return true;
      } else if (zeroQ1 && (distP0_eQ * distP1_eQ < 0)) {
        intersections.push_back(Q1[0]);
        intersections.push_back(Q1[1]);
        return true;
      }
      return false;
    }
    else if (zeroCount == 2) {
      /*--- Two powers are zero, one for each edge ---*/
      if (zeroP0 && zeroQ0) {
        intersections.push_back(P0[0]); // P0 = Q0
        intersections.push_back(P0[1]);
        return true;
      } else if (zeroP0 && zeroQ1) {
        intersections.push_back(P0[0]); // P0 = Q1
        intersections.push_back(P0[1]);
        return true;
      } else if (zeroP1 && zeroQ0) {
        intersections.push_back(P1[0]); // P1 = Q0
        intersections.push_back(P1[1]);
        return true;
      } else if (zeroP1 && zeroQ1) {
        intersections.push_back(P1[0]); // P1 = Q1
        intersections.push_back(P1[1]);
        return true;
      }
      /*--- Two powers zero on same edge - no intersection ---*/
      return false;
    }
    else if (zeroCount == 4) {
      /*--- All powers are zero - edges are aligned ---*/
      /*--- Use parametric coordinates like in CADTElemClass::Dist2ToLine ---*/
      /*--- X = X0 + (r+1)*(X1-X0)/2, -1 <= r <= 1 ---*/

      /*--- Convert edge P to parametric form: P(r) = P0 + (r+1)*(P1-P0)/2 ---*/
      /*--- Convert edge Q to parametric form: Q(s) = Q0 + (s+1)*(Q1-Q0)/2 ---*/

      /*--- Use the same parametrization as CADTElemClass::Dist2ToLine ---*/
      /*--- V0 = coor - (X1+X0)/2, V1 = (X1-X0)/2 ---*/
      su2double V0_P[2], V1_P[2], V0_Q[2], V1_Q[2];
      for (unsigned short k = 0; k < 2; ++k) {
        V0_P[k] = P0[k] - (P1[k] + P0[k]) / 2.0;
        V1_P[k] = (P1[k] - P0[k]) / 2.0;
        V0_Q[k] = Q0[k] - (Q1[k] + Q0[k]) / 2.0;
        V1_Q[k] = (Q1[k] - Q0[k]) / 2.0;
      }

      /*--- Find parameter values where edges overlap ---*/
      /*--- Since edges are aligned, we can project onto the direction vector ---*/
      su2double dotV1P_V1P = V1_P[0] * V1_P[0] + V1_P[1] * V1_P[1];

      if (dotV1P_V1P < EPS * EPS) {
        /*--- Edge P is degenerate (point) ---*/
        return false;
      }

      /*--- Project Q endpoints onto P parametric space ---*/
      /*--- For Q0: find r such that P(r) is closest to Q0 ---*/
      su2double diff_Q0[2] = {Q0[0] - (P1[0] + P0[0]) / 2.0, Q0[1] - (P1[1] + P0[1]) / 2.0};
      su2double r_Q0 = (diff_Q0[0] * V1_P[0] + diff_Q0[1] * V1_P[1]) / dotV1P_V1P;

      /*--- For Q1: find r such that P(r) is closest to Q1 ---*/
      su2double diff_Q1[2] = {Q1[0] - (P1[0] + P0[0]) / 2.0, Q1[1] - (P1[1] + P0[1]) / 2.0};
      su2double r_Q1 = (diff_Q1[0] * V1_P[0] + diff_Q1[1] * V1_P[1]) / dotV1P_V1P;

      /*--- Ensure r_Q0 <= r_Q1 for easier overlap computation ---*/
      if (r_Q0 > r_Q1) swap(r_Q0, r_Q1);

      /*--- Check for overlap using parametric coordinates ---*/
      su2double overlapStart = max(-1.0, r_Q0);
      su2double overlapEnd = min(1.0, r_Q1);

      if (overlapStart <= overlapEnd + EPS) {
        /*-------------------------------------------------------------------------------------------*/
        /*--- From Alauzet: "If all powers are zero, then the two edges are aligned. There is     ---*/
        /*--- intersection if and only if the edges overlap each other. There are two sub-cases:  ---*/
        /*---   * one intersection point that is the common point of the two edges                ---*/
        /*---   * two intersection points that are the end-points of the edge included in the     ---*/
        /*---     other one or one end-point of each edge in the other case"                      ---*/
        /*-------------------------------------------------------------------------------------------*/

        if (abs(overlapStart - overlapEnd) < EPS) {
          /*--- Single point intersection - common endpoint ---*/
          su2double r = overlapStart;
          intersections.push_back(P0[0] + (r + 1.0) * (P1[0] - P0[0]) / 2.0);
          intersections.push_back(P0[1] + (r + 1.0) * (P1[1] - P0[1]) / 2.0);
          return true;
        } else {
          /*--- Segment intersection - return both endpoints of overlap ---*/
          /*--- Start point of overlap ---*/
          su2double r_start = overlapStart;
          intersections.push_back(P0[0] + (r_start + 1.0) * (P1[0] - P0[0]) / 2.0);
          intersections.push_back(P0[1] + (r_start + 1.0) * (P1[1] - P0[1]) / 2.0);

          /*--- End point of overlap ---*/
          su2double r_end = overlapEnd;
          intersections.push_back(P0[0] + (r_end + 1.0) * (P1[0] - P0[0]) / 2.0);
          intersections.push_back(P0[1] + (r_end + 1.0) * (P1[1] - P0[1]) / 2.0);

          return true;
        }
      }

      /*--- No overlap ---*/
      return false;
    }
  }

  /*--- Standard case: no degenerate cases ---*/
  /*--- Check intersection condition: opposite signs ---*/
  if (distP0_eQ * distP1_eQ >= 0 || distQ0_eP * distQ1_eP >= 0) {
    return false;
  }

  /*--- Compute intersection point ---*/
  /*--- X = P0 + signedDistance(P0, e_Q) / (signedDistance(P0, e_Q) - signedDistance(P1, e_Q)) * vec{P0P1} ---*/
  su2double t = distP0_eQ / (distP0_eQ - distP1_eQ);

  /*--- vec{P0P1} = P1 - P0 ---*/
  intersections.push_back(P0[0] + t * (P1[0] - P0[0]));
  intersections.push_back(P0[1] + t * (P1[1] - P0[1]));

  return true;
}

void CConservativeVolumeInterpolator::AddEdgeNeighborToCandidates(CGeometry* geometry_src, 
                                                                  const unsigned long srcElemID,
                                                                  const unsigned short edgeIndex, 
                                                                  set<unsigned long>& newCandidates) {
  /*--- Get the element ---*/
  auto* srcElem = geometry_src->elem[srcElemID];

  /*--- Get the two nodes of the edge ---*/
  unsigned long node0 = srcElem->GetNode(edgeIndex);
  unsigned long node1 = srcElem->GetNode((edgeIndex + 1) % 3); // Circular indexing for triangle

  /*--- Find elements that share this edge (contain both nodes) ---*/
  for (unsigned long jElem = 0u; jElem < geometry_src->nodes->GetnElem(node0); ++jElem) {
    unsigned long elem0 =  geometry_src->nodes->GetElem(node0, jElem);
    if (elem0 == srcElemID) continue; // Skip self
    for (auto kElem = 0u; kElem < geometry_src->nodes->GetnElem(node1); ++kElem) {
      unsigned long elem1 =  geometry_src->nodes->GetElem(node1, kElem);
      if (elem0 == elem1) {
        newCandidates.insert(elem0);
        break;
      }
    }
  }
}

void CConservativeVolumeInterpolator::AddVertexBallToCandidates(CGeometry* geometry_src, 
                                                                const unsigned long srcElemID,
                                                                const unsigned short vertexIndex, 
                                                                set<unsigned long>& newCandidates) {
  /*--- Get the vertex node ID ---*/
  auto* srcElem = geometry_src->elem[srcElemID];
  unsigned long nodeID = srcElem->GetNode(vertexIndex);

  /*--- Add all elements other than the current source element that contain this vertex ---*/
  for (auto jElem = 0u; jElem < geometry_src->nodes->GetnElem(nodeID); ++jElem) {
    unsigned long elemID = geometry_src->nodes->GetElem(nodeID, jElem);
    if (elemID == srcElemID) continue; // Skip self

    auto* neighborElem = geometry_src->elem[elemID];
    if (neighborElem->GetVTK_Type() != TRIANGLE) continue;

    newCandidates.insert(elemID);
  }
}

void CConservativeVolumeInterpolator::ComputeConvexHull(vector<su2double>& points) {
  unsigned int n = points.size() / 2;
  if (n < 3) return;

  /*--- Convert to point pairs for easier handling ---*/
  vector<pair<su2double, su2double>> pts;
  for (unsigned int i = 0; i < n; ++i) {
    pts.push_back(make_pair(points[i*2], points[i*2+1]));
  }

  /*--- Compute centroid ---*/
  su2double cx = 0.0, cy = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    cx += pts[i].first;
    cy += pts[i].second;
  }
  cx /= n;
  cy /= n;

  /*--- Sort points by angle from centroid ---*/
  auto compareAngle = [&](const pair<su2double, su2double>& a, const pair<su2double, su2double>& b) {
    su2double angle_a = atan2(a.second - cy, a.first - cx);
    su2double angle_b = atan2(b.second - cy, b.first - cx);
    return angle_a < angle_b;
  };

  sort(pts.begin(), pts.end(), compareAngle);

  /*--- Convert back to flat array ---*/
  points.clear();
  for (const auto& pt : pts) {
    points.push_back(pt.first);
    points.push_back(pt.second);
  }
}