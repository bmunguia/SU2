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

void CConservativeVolumeInterpolator::Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
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

  /*--- Preprocess the solution to get the primitive variables ---*/
  /*--- TODO: other solver configurations                      ---*/
  solver_container_dst[FLOW_SOL]->Preprocessing(geometry_dst, solver_container_dst, config, 0, 0, RUNTIME_FLOW_SYS, false);
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE) {
    solver_container_dst[TURB_SOL]->Postprocessing(geometry_dst, solver_container_dst, config, 0);
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
  if (rank == MASTER_NODE) cout << "Performing containment search." << endl;
  vector<unsigned long> containingElems;
  vector<int> containingElemRanks;
  vector<unsigned long> pointsFailed;
  PointLocalization(geometry_src, coorDstCorrected, containingElems, containingElemRanks, pointsFailed);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute solution mass and gradient on source mesh          ---*/
  /*--------------------------------------------------------------------------*/
  const unsigned short nVar = solver_src->GetnVar();
  vector<vector<su2double> > srcElemMass(nElem_src, vector<su2double>(nVar, 0.0));
  vector<vector<su2double> > srcElemGrad(nElem_src, vector<su2double>(nVar * nDim, 0.0));
  ComputeSolutionMass(geometry_src, solver_src, srcElemMass, srcElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Compute the intersection of elements K_dst with elements   ---*/
  /*---         K_src it overlaps                                          ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing element intersections." << endl;
  map<unsigned long, vector<unsigned long>> overlappingElements;
  map<unsigned long, vector<vector<su2double>>> intersectionMeshes;
  ComputeOverlappingElements(geometry_src, geometry_dst, containingElems, overlappingElements, intersectionMeshes);

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Compute destination mesh mass and gradient using Gauss     ---*/
  /*---         quadrature over intersection regions                       ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing destination mesh mass and gradients." << endl;
  vector<vector<su2double>> dstElemMass(nElem_dst, vector<su2double>(nVar, 0.0));
  vector<vector<su2double>> dstElemGrad(nElem_dst, vector<su2double>(nVar * nDim, 0.0));
  ComputeDestinationMassAndGradient(geometry_src, geometry_dst, solver_src,
                                    overlappingElements, intersectionMeshes,
                                    srcElemMass, srcElemGrad,
                                    dstElemMass, dstElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 6: Correct the gradient to enforce the maximum principle      ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Applying local maximum principle correction." << endl;
  ApplyMaximumPrincipleCorrection(geometry_src, geometry_dst, solver_src,
                                  overlappingElements, srcElemMass, srcElemGrad,
                                  dstElemMass, dstElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 7: Perform averaging to get solution at vertices.             ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Distributing solution to destination nodes." << endl;
  DistributeSolutionToNodes(geometry_dst, solver_dst, dstElemMass, dstElemGrad);
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

  const unsigned short nVar = solver->GetnVar();
  for (unsigned long elemID = 0; elemID < geometry->GetnElem(); ++elemID) {
    auto* elem = geometry->elem[elemID];
    unsigned short nNodes = elem->GetnNodes();
    unsigned short VTK_Type = elem->GetVTK_Type();

    /*--- Get solution at nodes ---*/
    vector<vector<su2double> > solAtNodes(nNodes, vector<su2double>(nVar, 0.0));

    for (unsigned short iNode = 0; iNode < nNodes; ++iNode) {
      unsigned long nodeID = elem->GetNode(iNode);
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
        solAtNodes[iNode][iVar] = solver->GetNodes()->GetSolution(nodeID, iVar);
      }
    }

    /*--- Use CFEMStandardElement helper function for triangular elements ---*/
    if (VTK_Type != TRIANGLE) continue;

    /*--- Get triangle vertex coordinates ---*/
    su2double vertexCoords[6];
    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        unsigned long nodeID = elem->GetNode(iNode);
        vertexCoords[iNode * 2 + 0] = geometry->nodes->GetCoord(nodeID, 0);
        vertexCoords[iNode * 2 + 1] = geometry->nodes->GetCoord(nodeID, 1);
    }

    ComputeTriangleMassAndGradient(vertexCoords, solAtNodes, nVar,
                                   elemMass[elemID], elemGrad[elemID]);
  }
}

void CConservativeVolumeInterpolator::ComputeOverlappingElements(CGeometry* geometry_src,
                                                                 CGeometry* geometry_dst,
                                                                 const vector<unsigned long>& containingElems,
                                                                 map<unsigned long, vector<unsigned long>>& overlappingElements,
                                                                 map<unsigned long, vector<vector<su2double>>>& intersectionMeshes) {
  overlappingElements.clear();
  intersectionMeshes.clear();

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

    /*--------------------------------------------------------------------------*/
    /*--- Step 1: Build initial list from elements K_src containing vertices ---*/
    /*---         of K_dst                                                   ---*/
    /*--------------------------------------------------------------------------*/
    set<unsigned long> candidateList;

    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      /*--- Find corresponding point index in corrected coordinates ---*/
      if (containingElems[nodeID] != ULONG_MAX) {
        candidateList.insert(containingElems[nodeID]);

        /*--- TODO: Handle degenerate cases - add neighbors for vertices on edges/vertices ---*/
        /*--- For now, just add the containing element ---*/
      }
    }

    /*--- Skip this destination element if no initial candidates found ---*/
    if (candidateList.empty()) continue;

    /*--- Get destination triangle vertices ---*/
    su2double dstTri[6];
    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      dstTri[iNode * 2 + 0] = geometry_dst->nodes->GetCoord(nodeID, 0);
      dstTri[iNode * 2 + 1] = geometry_dst->nodes->GetCoord(nodeID, 1);
    }

    /*--------------------------------------------------------------------------*/
    /*--- Step 2: Process candidate list, adding new elements during         ---*/
    /*---         intersection, and meshing the convex polygon formed by the ---*/
    /*---         intersection                                               ---*/
    /*--------------------------------------------------------------------------*/
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

        /*--- Mesh the intersection polygon and store it ---*/
        vector<su2double> triangulatedMesh;
        MeshConvexPolygon(intersectionPoints, triangulatedMesh);
        intersectionMeshes[dstElemID].push_back(triangulatedMesh);

        /*--- Add new candidates to processing queue ---*/
        for (auto newElemID : newCandidates) {
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

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Compute vertex powers (signed distances)                   ---*/
  /*--------------------------------------------------------------------------*/
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
  if (power_t0v0_e01 >= EPS && power_t0v0_e12 >= EPS && power_t0v0_e20 >= EPS) {
    cloudPoints.push_back(t0v0[0]); cloudPoints.push_back(t0v0[1]);
  }
  if (power_t0v1_e01 >= EPS && power_t0v1_e12 >= EPS && power_t0v1_e20 >= EPS) {
    cloudPoints.push_back(t0v1[0]); cloudPoints.push_back(t0v1[1]);
  }
  if (power_t0v2_e01 >= EPS && power_t0v2_e12 >= EPS && power_t0v2_e20 >= EPS) {
    cloudPoints.push_back(t0v2[0]); cloudPoints.push_back(t0v2[1]);
  }

  /*--- Check if source vertices are inside destination triangle ---*/
  /*--- If so, add their vertex balls to candidates ---*/
  if (power_t1v0_e01 >= EPS && power_t1v0_e12 >= EPS && power_t1v0_e20 >= EPS) {
    cloudPoints.push_back(t1v0[0]); cloudPoints.push_back(t1v0[1]);
    AddVertexBallToCandidates(geometry_src, srcElemID, 0, newCandidates);
  }
  if (power_t1v1_e01 >= EPS && power_t1v1_e12 >= EPS && power_t1v1_e20 >= EPS) {
    cloudPoints.push_back(t1v1[0]); cloudPoints.push_back(t1v1[1]);
    AddVertexBallToCandidates(geometry_src, srcElemID, 1, newCandidates);
  }
  if (power_t1v2_e01 >= EPS && power_t1v2_e12 >= EPS && power_t1v2_e20 >= EPS) {
    cloudPoints.push_back(t1v2[0]); cloudPoints.push_back(t1v2[1]);
    AddVertexBallToCandidates(geometry_src, srcElemID, 2, newCandidates);
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Check edge-edge intersections with detailed analysis       ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Check all destination edges vs source edges ---*/
  IntersectionResult result;

  result = LineSegmentIntersection(t0v0, t0v1, t1v0, t1v1);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 0, newCandidates); // edge 0-1
    }
  }

  result = LineSegmentIntersection(t0v0, t0v1, t1v1, t1v2);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 1, newCandidates); // edge 1-2
    }
  }

  result = LineSegmentIntersection(t0v0, t0v1, t1v2, t1v0);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 2, newCandidates); // edge 2-0
    }
  }

  result = LineSegmentIntersection(t0v1, t0v2, t1v0, t1v1);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 0, newCandidates); // edge 0-1
    }
  }

  result = LineSegmentIntersection(t0v1, t0v2, t1v1, t1v2);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 1, newCandidates); // edge 1-2
    }
  }

  result = LineSegmentIntersection(t0v1, t0v2, t1v2, t1v0);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 2, newCandidates); // edge 2-0
    }
  }

  result = LineSegmentIntersection(t0v2, t0v0, t1v0, t1v1);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 0, newCandidates); // edge 0-1
    }
  }

  result = LineSegmentIntersection(t0v2, t0v0, t1v1, t1v2);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 1, newCandidates); // edge 1-2
    }
  }

  result = LineSegmentIntersection(t0v2, t0v0, t1v2, t1v0);
  if (result.type != NO_INTERSECTION) {
    for (size_t i = 0; i < result.points.size(); ++i)
      cloudPoints.push_back(result.points[i]);
    if (result.type == VERTEX_ON_EDGE && result.isFirstEdge) {
      AddVertexBallToCandidates(geometry_src, srcElemID, result.vertexIndex, newCandidates);
    } else if (result.type == EDGE_CROSSING) {
      AddFaceNeighborToCandidates(geometry_src, srcElemID, 2, newCandidates); // edge 2-0
    }
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

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute convex hull from intersection points               ---*/
  /*--------------------------------------------------------------------------*/
  ComputeConvexHull(cloudPoints);
  intersectionPoints = cloudPoints;

  return intersectionPoints.size() >= 6; // At least 3 points for a valid intersection polygon
}

su2double CConservativeVolumeInterpolator::ComputeLineParametricCoord(const su2double point[2],
                                                                      const su2double lineStart[2],
                                                                      const su2double lineEnd[2]) {
  /*--- Convert point to parametric coordinate on line segment ---*/
  /*--- Line: X = lineStart + (r+1)*(lineEnd-lineStart)/2, r in [-1,1] ---*/
  /*--- Solve for r: r = 2*(X-lineStart)/(lineEnd-lineStart) - 1 ---*/

  su2double dx = lineEnd[0] - lineStart[0];
  su2double dy = lineEnd[1] - lineStart[1];
  su2double px = point[0] - lineStart[0];
  su2double py = point[1] - lineStart[1];

  /*--- Project point onto line direction and normalize ---*/
  su2double dot_product = px * dx + py * dy;
  su2double line_length_sq = dx * dx + dy * dy;

  if (line_length_sq < 1e-12) {
    return 0.0; // Degenerate line
  }

  su2double t = dot_product / line_length_sq;
  return 2.0 * t - 1.0; // Convert from [0,1] to [-1,1]
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

IntersectionResult CConservativeVolumeInterpolator::LineSegmentIntersection(const su2double P0[2], const su2double P1[2],
                                                                            const su2double Q0[2], const su2double Q1[2]) {
  const su2double EPS = 1e-12;
  IntersectionResult result;
  result.type = NO_INTERSECTION;
  result.vertexIndex = -1;
  result.isFirstEdge = false;
  result.points.clear();

  /*--- Following Alauzet: Let e_P = [P0 P1] and e_Q = [Q0 Q1] be two edges ---*/

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
      /*--- Case 1: Only one power is zero - vertex lies on edge ---*/
      if (zeroP0 && (distQ0_eP * distQ1_eP < 0)) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 0;
        result.isFirstEdge = true;
        result.points.push_back(P0[0]);
        result.points.push_back(P0[1]);
        return result;
      } else if (zeroP1 && (distQ0_eP * distQ1_eP < 0)) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 1;
        result.isFirstEdge = true;
        result.points.push_back(P1[0]);
        result.points.push_back(P1[1]);
        return result;
      } else if (zeroQ0 && (distP0_eQ * distP1_eQ < 0)) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 0;
        result.isFirstEdge = false;
        result.points.push_back(Q0[0]);
        result.points.push_back(Q0[1]);
        return result;
      } else if (zeroQ1 && (distP0_eQ * distP1_eQ < 0)) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 1;
        result.isFirstEdge = false;
        result.points.push_back(Q1[0]);
        result.points.push_back(Q1[1]);
        return result;
      }
      return result; // NO_INTERSECTION
    }
    else if (zeroCount == 2) {
      /*--- Case 2: Two powers are zero - endpoint intersections ---*/
      if (zeroP0 && zeroQ0) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 0;
        result.isFirstEdge = true; // P0 = Q0
        result.points.push_back(P0[0]);
        result.points.push_back(P0[1]);
        return result;
      } else if (zeroP0 && zeroQ1) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 0;
        result.isFirstEdge = true; // P0 = Q1
        result.points.push_back(P0[0]);
        result.points.push_back(P0[1]);
        return result;
      } else if (zeroP1 && zeroQ0) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 1;
        result.isFirstEdge = true; // P1 = Q0
        result.points.push_back(P1[0]);
        result.points.push_back(P1[1]);
        return result;
      } else if (zeroP1 && zeroQ1) {
        result.type = VERTEX_ON_EDGE;
        result.vertexIndex = 1;
        result.isFirstEdge = true; // P1 = Q1
        result.points.push_back(P1[0]);
        result.points.push_back(P1[1]);
        return result;
      }
      return result; // NO_INTERSECTION
    }
    else if (zeroCount == 4) {
      /*--- Case 3: All powers are zero - edges are aligned ---*/
      result.type = EDGE_OVERLAP;

      /*--- Use parametric coordinates to find overlap ---*/
      su2double r_Q0 = ComputeLineParametricCoord(Q0, P0, P1);
      su2double r_Q1 = ComputeLineParametricCoord(Q1, P0, P1);

      if (r_Q0 > r_Q1) {
        su2double temp = r_Q0;
        r_Q0 = r_Q1;
        r_Q1 = temp;
      }

      su2double overlapStart = max(-1.0, r_Q0);
      su2double overlapEnd = min(1.0, r_Q1);

      if (overlapStart <= overlapEnd + EPS) {
        if (abs(overlapStart - overlapEnd) < EPS) {
          /*--- Single point overlap ---*/
          su2double r = overlapStart;
          result.points.push_back(P0[0] + (r + 1.0) * (P1[0] - P0[0]) / 2.0);
          result.points.push_back(P0[1] + (r + 1.0) * (P1[1] - P0[1]) / 2.0);
        } else {
          /*--- Segment overlap ---*/
          su2double r_start = overlapStart;
          result.points.push_back(P0[0] + (r_start + 1.0) * (P1[0] - P0[0]) / 2.0);
          result.points.push_back(P0[1] + (r_start + 1.0) * (P1[1] - P0[1]) / 2.0);

          su2double r_end = overlapEnd;
          result.points.push_back(P0[0] + (r_end + 1.0) * (P1[0] - P0[0]) / 2.0);
          result.points.push_back(P0[1] + (r_end + 1.0) * (P1[1] - P0[1]) / 2.0);
        }
        return result;
      }

      result.type = NO_INTERSECTION;
      return result;
    }
  }

  /*--- Case 4: Standard intersection - edges cross at interior points ---*/
  if (distP0_eQ * distP1_eQ >= 0 || distQ0_eP * distQ1_eP >= 0) {
    return result; // NO_INTERSECTION
  }

  /*--- Compute intersection point ---*/
  su2double t = distP0_eQ / (distP0_eQ - distP1_eQ);

  result.type = EDGE_CROSSING;
  result.points.push_back(P0[0] + t * (P1[0] - P0[0]));
  result.points.push_back(P0[1] + t * (P1[1] - P0[1]));

  return result;
}

void CConservativeVolumeInterpolator::AddFaceNeighborToCandidates(CGeometry* geometry,
                                                                  const unsigned long elemID,
                                                                  const unsigned short faceIndex,
                                                                  set<unsigned long>& newCandidates) {
  /*--- Get the neighbor element across this face ---*/
  long neighElemID = geometry->elem[elemID]->GetNeighbor_Elements(faceIndex);

  /*--- Only add if a valid neighbor exists (not boundary) ---*/
  if (neighElemID >= 0) {
    newCandidates.insert(static_cast<unsigned long>(neighElemID));
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

void CConservativeVolumeInterpolator::MeshConvexPolygon(const vector<su2double>& polygonPoints,
                                                        vector<su2double>& triangles) {
  triangles.clear();

  unsigned int numPoints = polygonPoints.size() / 2;

  /*--- Handle edge cases ---*/
  if (numPoints < 3) {
    return; // Cannot mesh a polygon with less than 3 points
  }

  if (numPoints == 3) {
    /*--- Already a triangle, just copy the points ---*/
    triangles = polygonPoints;
    return;
  }

  /*--- For a convex polygon with n points (n >= 4), create (n-2) triangles ---*/
  /*--- using fan triangulation from the first vertex ---*/
  /*--- This handles polygons with 4, 5, 6, or more vertices ---*/

  /*--- Get first vertex coordinates ---*/
  su2double x0 = polygonPoints[0];
  su2double y0 = polygonPoints[1];

  /*--- Create triangles by connecting first vertex to consecutive edge pairs ---*/
  for (unsigned int i = 1; i < numPoints - 1; ++i) {
    /*--- Get coordinates of the other two vertices ---*/
    su2double x1 = polygonPoints[i * 2];
    su2double y1 = polygonPoints[i * 2 + 1];
    su2double x2 = polygonPoints[(i + 1) * 2];
    su2double y2 = polygonPoints[(i + 1) * 2 + 1];

    /*--- Check orientation to ensure consistent winding ---*/
    su2double cross = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);

    if (abs(cross) > 1e-12) { // Only add non-degenerate triangles
      if (cross > 0) {
        /*--- Counter-clockwise orientation - add triangle as is ---*/
        triangles.push_back(x0); triangles.push_back(y0);
        triangles.push_back(x1); triangles.push_back(y1);
        triangles.push_back(x2); triangles.push_back(y2);
      } else {
        /*--- Clockwise orientation - reverse order to maintain consistent winding ---*/
        triangles.push_back(x0); triangles.push_back(y0);
        triangles.push_back(x2); triangles.push_back(y2);
        triangles.push_back(x1); triangles.push_back(y1);
      }
    }
  }
}

void CConservativeVolumeInterpolator::ComputeDestinationMassAndGradient(CGeometry* geometry_src,
                                                                        CGeometry* geometry_dst,
                                                                        CSolver* solver_src,
                                                                        const map<unsigned long, vector<unsigned long>>& overlappingElements,
                                                                        const map<unsigned long, vector<vector<su2double>>>& intersectionMeshes,
                                                                        const vector<vector<su2double>>& srcElemMass,
                                                                        const vector<vector<su2double>>& srcElemGrad,
                                                                        vector<vector<su2double>>& dstElemMass,
                                                                        vector<vector<su2double>>& dstElemGrad) {

  const unsigned short nVar = solver_src->GetnVar();

  /*--- Loop over all destination elements that have intersections ---*/
  for (const auto& elemPair : overlappingElements) {
    unsigned long dstElemID = elemPair.first;
    const vector<unsigned long>& srcElemIDs = elemPair.second;

    /*--- Get the corresponding intersection meshes ---*/
    const vector<vector<su2double>>& intersectionMeshList = intersectionMeshes.at(dstElemID);

    /*--- Initialize destination element mass and gradient ---*/
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      dstElemMass[dstElemID][iVar] = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        dstElemGrad[dstElemID][iVar * nDim + iDim] = 0.0;
      }
    }

    /*--- Process each intersection region T_j = intersection(K_dst, K_src_j) ---*/
    for (size_t j = 0; j < srcElemIDs.size(); ++j) {
      unsigned long srcElemID = srcElemIDs[j];
      const vector<su2double>& intersectionMesh = intersectionMeshList[j];

      /*--- Gauss quadrature over all triangles in the intersection mesh ---*/
      unsigned int numTriangles = intersectionMesh.size() / 6; // 6 coordinates per triangle

      for (unsigned int iTri = 0; iTri < numTriangles; ++iTri) {
        /*--- Get triangle vertices ---*/
        su2double x0 = intersectionMesh[iTri * 6 + 0];
        su2double y0 = intersectionMesh[iTri * 6 + 1];
        su2double x1 = intersectionMesh[iTri * 6 + 2];
        su2double y1 = intersectionMesh[iTri * 6 + 3];
        su2double x2 = intersectionMesh[iTri * 6 + 4];
        su2double y2 = intersectionMesh[iTri * 6 + 5];

        /*--- Triangle area ---*/
        su2double area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));

        if (area < 1e-12) continue; // Skip degenerate triangles

        /*--- Use 1-point Gauss quadrature (centroid rule) for exact integration ---*/
        su2double xc = (x0 + x1 + x2) / 3.0; // Triangle centroid
        su2double yc = (y0 + y1 + y2) / 3.0;

        /*--- Evaluate solution at centroid using source element data ---*/
        /*--- u(x) = u_K_src + gra(u_K_src) · (x - x_K_src) ---*/
        /*--- where x_K_src is the source element centroid ---*/

        /*--- Get source element centroid ---*/
        auto* srcElem = geometry_src->elem[srcElemID];
        su2double srcCentroid[2] = {0.0, 0.0};
        for (unsigned short iNode = 0; iNode < 3; ++iNode) {
          unsigned long nodeID = srcElem->GetNode(iNode);
          srcCentroid[0] += geometry_src->nodes->GetCoord(nodeID, 0);
          srcCentroid[1] += geometry_src->nodes->GetCoord(nodeID, 1);
        }
        srcCentroid[0] /= 3.0;
        srcCentroid[1] /= 3.0;

        /*--- Compute displacement from source centroid to integration point ---*/
        su2double dx = xc - srcCentroid[0];
        su2double dy = yc - srcCentroid[1];

        /*--- Integrate mass: int_T (u dA) ---*/
        for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
          /*--- Get solution value at source element (piecewise constant) ---*/
          su2double u_src = srcElemMass[srcElemID][iVar] / geometry_src->elem[srcElemID]->GetVolume();

          /*--- Apply gradient correction: u = u_src + gra(u_src) · (x - x_src) ---*/
          su2double gradx = srcElemGrad[srcElemID][iVar * nDim + 0];
          su2double grady = srcElemGrad[srcElemID][iVar * nDim + 1];
          su2double u_corrected = u_src + gradx * dx + grady * dy;

          /*--- Add contribution to destination mass: weight * area * u ---*/
          dstElemMass[dstElemID][iVar] += area * u_corrected;

          /*--- Add contribution to destination gradient integral ---*/
          /*--- int_T (gra(u) dA) = int_T (gra(u_src) dA) (constant gradient over source element) ---*/
          dstElemGrad[dstElemID][iVar * nDim + 0] += area * gradx;
          dstElemGrad[dstElemID][iVar * nDim + 1] += area * grady;
        }
      }
    }

    /*--- Convert gradient integral to volume average: gra(u_dst) = int_K_dst (gra(u) dA) / |K_dst| ---*/
    su2double dstVolume = geometry_dst->elem[dstElemID]->GetVolume();
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        dstElemGrad[dstElemID][iVar * nDim + iDim] /= dstVolume;
      }
    }
  }
}

void CConservativeVolumeInterpolator::ApplyMaximumPrincipleCorrection(CGeometry* geometry_src,
                                                                      CGeometry* geometry_dst,
                                                                      CSolver* solver_src,
                                                                      const map<unsigned long, vector<unsigned long>>& overlappingElements,
                                                                      const vector<vector<su2double>>& srcElemMass,
                                                                      const vector<vector<su2double>>& srcElemGrad,
                                                                      vector<vector<su2double>>& dstElemMass,
                                                                      vector<vector<su2double>>& dstElemGrad) {

  const unsigned short nVar = solver_src->GetnVar();
  const su2double EPS = 1e-12;

  /*--- Loop over all destination elements that have overlaps ---*/
  for (const auto& elemPair : overlappingElements) {
    unsigned long dstElemID = elemPair.first;
    const vector<unsigned long>& srcElemIDs = elemPair.second;

    auto* dstElem = geometry_dst->elem[dstElemID];
    if (dstElem->GetVTK_Type() != TRIANGLE) continue;

    /*--- Get destination element vertices ---*/
    su2double dstVertices[6];
    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      dstVertices[iNode * 2 + 0] = geometry_dst->nodes->GetCoord(nodeID, 0);
      dstVertices[iNode * 2 + 1] = geometry_dst->nodes->GetCoord(nodeID, 1);
    }

    /*--- Compute element centroid (barycenter G_K) ---*/
    su2double G_K[2] = {
      (dstVertices[0] + dstVertices[2] + dstVertices[4]) / 3.0,
      (dstVertices[1] + dstVertices[3] + dstVertices[5]) / 3.0
    };

    /*--- For each variable, apply Alauzet's maximum principle correction ---*/
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      /*--------------------------------------------------------------------------*/
      /*--- Step 1: Compute local bounds from overlapping source elements      ---*/
      /*--------------------------------------------------------------------------*/
      su2double u_min = 1e20;
      su2double u_max = -1e20;

      /*--- Find all vertices Q from source elements K_src that K overlaps ---*/
      for (unsigned long srcElemID : srcElemIDs) {
        auto* srcElem = geometry_src->elem[srcElemID];
        if (srcElem->GetVTK_Type() != TRIANGLE) continue;

        /*--- Get solution values at vertices of source element ---*/
        for (unsigned short iNode = 0; iNode < 3; ++iNode) {
          unsigned long nodeID = srcElem->GetNode(iNode);
          su2double u_vertex = solver_src->GetNodes()->GetSolution(nodeID, iVar);
          u_min = min(u_min, u_vertex);
          u_max = max(u_max, u_vertex);
        }
      }

      /*--- Skip if no valid bounds found ---*/
      if (u_min > 1e19 || u_max < -1e19) continue;

      /*--------------------------------------------------------------------------*/
      /*--- Step 2: Get current solution at destination element                 ---*/
      /*--------------------------------------------------------------------------*/
      su2double elemVolume = dstElem->GetVolume();
      su2double u_K_G = dstElemMass[dstElemID][iVar] / elemVolume;  // u_K(G_K)
      su2double grad_x = dstElemGrad[dstElemID][iVar * nDim + 0];   // gra(u_K) · x
      su2double grad_y = dstElemGrad[dstElemID][iVar * nDim + 1];   // gra(u_K) · y

      /*--------------------------------------------------------------------------*/
      /*--- Step 3: Compute u_K(P_i) at each vertex using Taylor expansion     ---*/
      /*--------------------------------------------------------------------------*/
      su2double u_K_P[3]; // Values at vertices P_0, P_1, P_2

      for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        /*--- Vector G_K P_i ---*/
        su2double dx = dstVertices[iNode * 2 + 0] - G_K[0];
        su2double dy = dstVertices[iNode * 2 + 1] - G_K[1];

        /*--- u_K(P_i) = u_K(G_K) + gra(u_K) · G_K P_i ---*/
        u_K_P[iNode] = u_K_G + grad_x * dx + grad_y * dy;
      }

      /*--------------------------------------------------------------------------*/
      /*--- Step 4: Check if maximum principle is violated                     ---*/
      /*--------------------------------------------------------------------------*/
      bool violatesMaxPrinciple = false;
      for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        if (u_K_P[iNode] < u_min - EPS || u_K_P[iNode] > u_max + EPS) {
          violatesMaxPrinciple = true;
          break;
        }
      }

      /*--- If no violation, skip correction ---*/
      if (!violatesMaxPrinciple) continue;

      /*--------------------------------------------------------------------------*/
      /*--- Step 5: Apply Alauzet's correction algorithm                       ---*/
      /*--------------------------------------------------------------------------*/
      /*--- Sort vertices by solution value: u_K(P_0) ≤ u_K(P_1) ≤ u_K(P_2) ---*/
      vector<pair<su2double, unsigned short>> sortedValues;
      for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        sortedValues.push_back(make_pair(u_K_P[iNode], iNode));
      }
      sort(sortedValues.begin(), sortedValues.end());

      su2double u_P0 = sortedValues[0].first;  // Smallest value
      su2double u_P1 = sortedValues[1].first;  // Middle value
      su2double u_P2 = sortedValues[2].first;  // Largest value

      /*--- Apply first correction pass ---*/
      su2double u_M_P2 = min(u_P2, u_max);
      su2double u_M_P1 = min(u_P1 + 0.5 * max(0.0, u_P2 - u_max), u_max);
      su2double u_M_P0 = 3.0 * u_K_G - u_M_P1 - u_M_P2;

      /*--- Apply second correction pass ---*/
      su2double u_tilde_P0 = max(u_M_P0, u_min);
      su2double u_tilde_P1 = max(u_M_P1 - 0.5 * max(0.0, u_min - u_M_P0), u_min);
      su2double u_tilde_P2 = 3.0 * u_K_G - u_tilde_P0 - u_tilde_P1;

      /*--------------------------------------------------------------------------*/
      /*--- Step 6: Compute corrected mass and gradient from new nodal values  ---*/
      /*--------------------------------------------------------------------------*/
      /*--- Put corrected values back in original vertex order ---*/
      su2double u_tilde[3];
      u_tilde[sortedValues[0].second] = u_tilde_P0;
      u_tilde[sortedValues[1].second] = u_tilde_P1;
      u_tilde[sortedValues[2].second] = u_tilde_P2;

      /*--- Prepare vertex solutions for single variable ---*/
      vector<vector<su2double>> solAtVertices(3, vector<su2double>(1));
      for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        solAtVertices[iNode][0] = u_tilde[iNode];
      }

      /*--- Use helper function to compute mass and gradient robustly ---*/
      vector<su2double> correctedMass, correctedGrad;
      ComputeTriangleMassAndGradient(dstVertices, solAtVertices, 1, correctedMass, correctedGrad);

      /*--- Update destination element data ---*/
      dstElemMass[dstElemID][iVar] = correctedMass[0];
      dstElemGrad[dstElemID][iVar * nDim + 0] = correctedGrad[0];
      dstElemGrad[dstElemID][iVar * nDim + 1] = correctedGrad[1];
    }
  }
}

void CConservativeVolumeInterpolator::DistributeSolutionToNodes(CGeometry* geometry_dst,
                                                                CSolver* solver_dst,
                                                                const vector<vector<su2double>>& dstElemMass,
                                                                const vector<vector<su2double>>& dstElemGrad) {

  const unsigned short nVar = solver_dst->GetnVar();
  const su2double EPS = 1e-12;

  /*--- Initialize vertex solution arrays ---*/
  vector<vector<su2double>> vertexSolution(geometry_dst->GetnPoint(), vector<su2double>(nVar, 0.0));
  vector<su2double> vertexWeight(geometry_dst->GetnPoint(), 0.0);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Accumulate contributions from each element                 ---*/
  /*--------------------------------------------------------------------------*/
  for (unsigned long dstElemID = 0; dstElemID < geometry_dst->GetnElem(); ++dstElemID) {
    auto* dstElem = geometry_dst->elem[dstElemID];

    if (dstElem->GetVTK_Type() != TRIANGLE) continue;

    /*--- Check if element has valid data ---*/
    bool hasData = false;
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      if (abs(dstElemMass[dstElemID][iVar]) > EPS) {
        hasData = true;
        break;
      }
    }
    if (!hasData) continue;

    /*--- Get element vertices and centroid ---*/
    su2double dstVertices[6];
    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      dstVertices[iNode * 2 + 0] = geometry_dst->nodes->GetCoord(nodeID, 0);
      dstVertices[iNode * 2 + 1] = geometry_dst->nodes->GetCoord(nodeID, 1);
    }

    su2double G_K[2] = {
      (dstVertices[0] + dstVertices[2] + dstVertices[4]) / 3.0,
      (dstVertices[1] + dstVertices[3] + dstVertices[5]) / 3.0
    };

    su2double elemVolume = dstElem->GetVolume();

    /*--- Loop over variables ---*/
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      /*--- Get element-centered solution ---*/
      su2double u_K_G = dstElemMass[dstElemID][iVar] / elemVolume;
      su2double grad_x = dstElemGrad[dstElemID][iVar * nDim + 0];
      su2double grad_y = dstElemGrad[dstElemID][iVar * nDim + 1];

      /*--- Interpolate to each vertex ---*/
      for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        unsigned long nodeID = dstElem->GetNode(iNode);

        /*--- Vector from centroid to vertex ---*/
        su2double dx = dstVertices[iNode * 2 + 0] - G_K[0];
        su2double dy = dstVertices[iNode * 2 + 1] - G_K[1];

        /*--- Linear reconstruction: u(P_i) = u(G_K) + gra(u) · (P_i - G_K) ---*/
        su2double vertexValue = u_K_G + grad_x * dx + grad_y * dy;

        /*--- Accumulate weighted contribution ---*/
        vertexSolution[nodeID][iVar] += vertexValue * elemVolume;
      }
    }

    /*--- Accumulate weights for all vertices of this element ---*/
    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      vertexWeight[nodeID] += elemVolume;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Compute weighted averages and set solution                 ---*/
  /*--------------------------------------------------------------------------*/
  for (unsigned long nodeID = 0; nodeID < geometry_dst->GetnPoint(); ++nodeID) {
    if (vertexWeight[nodeID] > EPS) {
      /*--- Compute average ---*/
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
        su2double avgValue = vertexSolution[nodeID][iVar] / vertexWeight[nodeID];
        solver_dst->GetNodes()->SetSolution(nodeID, iVar, avgValue);
      }
    }
  }
}

void CConservativeVolumeInterpolator::ComputeTriangleMassAndGradient(const su2double vertexCoords[6],
                                                                     const vector<vector<su2double>>& vertexSolutions,
                                                                     unsigned short nVar,
                                                                     vector<su2double>& mass,
                                                                     vector<su2double>& gradient) {
  /*--- Initialize output ---*/
  mass.assign(nVar, 0.0);
  gradient.assign(nVar * nDim, 0.0);

  /*--- Create standard triangular element ---*/
  unsigned short VTK_Type = TRIANGLE;
  unsigned short nPoly = 1;
  bool constJac = false;
  unsigned short orderExact = 2 * nPoly;

  CFEMStandardElement stdElement(VTK_Type, nPoly, constJac, nullptr, orderExact);

  /*--- Get integration points and basis functions ---*/
  const unsigned short nInt = stdElement.GetNIntegration();
  const su2double* weights = stdElement.GetWeightsIntegration();
  const su2double* lagBasis = stdElement.GetBasisFunctionsIntegration();
  const su2double* drLagBasis = stdElement.GetDrBasisFunctionsIntegration();
  const su2double* dsLagBasis = stdElement.GetDsBasisFunctionsIntegration();

  /*--- Convert vertex coordinates to node coordinates array ---*/
  su2double nodeCoords[6];
  for (unsigned short i = 0; i < 6; ++i) {
    nodeCoords[i] = vertexCoords[i];
  }

  /*--- Loop over integration points ---*/
  for (unsigned short iInt = 0; iInt < nInt; ++iInt) {

    /*--- Compute Jacobian of transformation ---*/
    su2double dxdr = 0.0, dydr = 0.0;
    su2double dxds = 0.0, dyds = 0.0;

    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned short ind = iInt * 3 + iNode;

      dxdr += nodeCoords[iNode * 2 + 0] * drLagBasis[ind];
      dydr += nodeCoords[iNode * 2 + 1] * drLagBasis[ind];
      dxds += nodeCoords[iNode * 2 + 0] * dsLagBasis[ind];
      dyds += nodeCoords[iNode * 2 + 1] * dsLagBasis[ind];
    }

    /*--- Compute Jacobian determinant ---*/
    su2double jacobian = abs(dxdr * dyds - dydr * dxds);

    /*--- Compute inverse Jacobian for gradient transformation ---*/
    su2double jacInv = 1.0 / jacobian;
    su2double drdx = dyds * jacInv;
    su2double drdy = -dxds * jacInv;
    su2double dsdx = -dydr * jacInv;
    su2double dsdy = dxdr * jacInv;

    /*--- Integration weight including Jacobian ---*/
    su2double intWeight = weights[iInt] * jacobian;

    /*--- Loop over variables ---*/
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {

      /*--- Interpolate solution and its gradient at integration point ---*/
      su2double solVal = 0.0;
      su2double dudr = 0.0, duds = 0.0;

      for (unsigned short iNode = 0; iNode < 3; ++iNode) {
        unsigned short ind = iInt * 3 + iNode;
        su2double nodeVal = vertexSolutions[iNode][iVar];

        solVal += lagBasis[ind] * nodeVal;
        dudr += drLagBasis[ind] * nodeVal;
        duds += dsLagBasis[ind] * nodeVal;
      }

      /*--- Add contribution to mass integral ---*/
      mass[iVar] += solVal * intWeight;

      /*--- Transform gradient to physical coordinates and add to gradient integral ---*/
      su2double dudx = dudr * drdx + duds * dsdx;
      su2double dudy = dudr * drdy + duds * dsdy;

      gradient[iVar * nDim + 0] += dudx * intWeight;
      gradient[iVar * nDim + 1] += dudy * intWeight;
    }
  }

  /*--- For FVM with constant gradients, normalize by element volume ---*/
  /*--- Volume is computed as the sum of integration weights * Jacobian ---*/
  su2double elemVolume = 0.0;
  for (unsigned short iInt = 0; iInt < nInt; ++iInt) {
    su2double dxdr = 0.0, dydr = 0.0, dxds = 0.0, dyds = 0.0;
    for (unsigned short iNode = 0; iNode < 3; ++iNode) {
      unsigned short ind = iInt * 3 + iNode;
      dxdr += nodeCoords[iNode * 2 + 0] * drLagBasis[ind];
      dydr += nodeCoords[iNode * 2 + 1] * drLagBasis[ind];
      dxds += nodeCoords[iNode * 2 + 0] * dsLagBasis[ind];
      dyds += nodeCoords[iNode * 2 + 1] * dsLagBasis[ind];
    }
    su2double jacobian = abs(dxdr * dyds - dydr * dxds);
    elemVolume += weights[iInt] * jacobian;
  }

  /*--- Normalize gradients by volume ---*/
  for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
    for (unsigned short k = 0; k < nDim; ++k) {
      gradient[iVar * nDim + k] /= elemVolume;
    }
  }
}