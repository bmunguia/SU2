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

#include <cmath>
#include <iomanip>

#include "../include/CConservativeVolumeInterpolator.hpp"

CConservativeVolumeInterpolator::CConservativeVolumeInterpolator(SU2_Comm MPICommunicator)
    : CVolumeInterpolator(MPICommunicator) {}

void CConservativeVolumeInterpolator::Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                  CSolver** solver_container_src, CSolver** solver_container_dst) {
  if (rank == MASTER_NODE) {
    cout << endl << "----------------------------- Interpolation -----------------------------" << endl;
    cout << "Conservative solution interpolation from source mesh to destination mesh." << endl;
    cout << "Source mesh: " << geometry_src->GetGlobal_nPointDomain() << " points, ";
    cout << geometry_src->GetGlobal_nElemDomain() << " elements." << endl;
    cout << "Destination mesh: " << geometry_dst->GetGlobal_nPointDomain() << " points, ";
    cout << geometry_dst->GetGlobal_nElemDomain() << " elements." << endl;
  }

  /*--- Build the ADTs ---*/
  InitializeADTs(config, geometry_src, geometry_dst);

  /*--- Build the R-trees ---*/
  InitializeRTrees(geometry_src, geometry_dst);

  /*--- Initialize the FEM standard element ---*/
  if (!stdElement_ptr) {
    if (rank == MASTER_NODE) cout << "Creating FEM standard element." << flush;
    stdElement_ptr = InitializeFEMStandardElement(config);
    if (rank == MASTER_NODE) cout << " Done." << endl;
  }

  /*--- Call the conservative interpolation method ---*/
  for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
    auto solver_src = solver_container_src[iSol];
    auto solver_dst = solver_container_dst[iSol];
    if (solver_src && solver_dst)
      ConservativeInterpolation(config, geometry_src, geometry_dst, solver_src, solver_dst);
  }

  /*--- Preprocess the solution to get the primitive variables ---*/
  /*--- TODO: other solver configurations                      ---*/
  solver_container_dst[FLOW_SOL]->Preprocessing(geometry_dst, solver_container_dst, config, 0, 0, RUNTIME_FLOW_SYS,
                                                false);
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE) {
    solver_container_dst[TURB_SOL]->Postprocessing(geometry_dst, solver_container_dst, config, 0);
  }
}

void CConservativeVolumeInterpolator::ConservativeInterpolation(const CConfig* config,
                                                                CGeometry* geometry_src,
                                                                CGeometry* geometry_dst,
                                                                CSolver* solver_src,
                                                                CSolver* solver_dst) {
  /*--------------------------------------------------------------------------*/
  /*--- Step 0: Initialize destination coordinate vector. If applying a    ---*/
  /*---         curvature correction, this vector will be modified.        ---*/
  /*---         Otherwise, it will just contain the original coordinates.  ---*/
  /*--------------------------------------------------------------------------*/
  nVar = solver_src->GetnVar();
  vector<su2double> coorDst;
  InitializeCoords(geometry_dst, coorDst);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Compute the intersection of elements K_dst with elements   ---*/
  /*---         K_src it overlaps, and mesh the intersection regions.      ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing element intersections." << endl;
  IntersectionMeshMap overlapMeshes;
  vector<unsigned long> incompleteOverlaps;
  set<unsigned long> uncoveredNodes;
  CreateIntersectionMeshes(geometry_src, geometry_dst, coorDst, overlapMeshes, incompleteOverlaps, uncoveredNodes);

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Compute solution mass and gradient on source mesh.         ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing source mesh mass and gradients." << endl;
  vector<vector<su2double>> srcElemMass;
  vector<vector<su2double>> srcElemGrad;
  ComputeSourceSolutionMass(geometry_src, solver_src, srcElemMass, srcElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute destination mesh mass and gradient using Gauss     ---*/
  /*---         quadrature over intersection regions.                      ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing destination mesh mass and gradients." << endl;
  vector<vector<su2double>> dstElemMass;
  vector<vector<su2double>> dstElemGrad;
  ComputeDestinationMassAndGradient(geometry_src, geometry_dst, solver_src, overlapMeshes,
                                    incompleteOverlaps, srcElemMass, srcElemGrad, dstElemMass,
                                    dstElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Correct the gradient to enforce the maximum principle.     ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Applying local maximum principle correction." << endl;
  ApplyMaximumPrincipleCorrection(geometry_src, geometry_dst, solver_src, coorDst, overlapMeshes,
                                  srcElemMass, srcElemGrad, dstElemMass, dstElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Perform averaging to get solution at vertices.             ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Distributing solution to destination nodes." << endl;
  DistributeSolutionToNodes(geometry_dst, solver_dst, coorDst, dstElemMass, dstElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 6: Handle volume nodes outside the source domain.             ---*/
  /*--------------------------------------------------------------------------*/
  if (!uncoveredNodes.empty()) {
    if (rank == MASTER_NODE) cout << "Distributing solution to nodes outside source domain." << endl;
    DistributeSolutionOutsideDomain(geometry_src, geometry_dst, solver_src, solver_dst, coorDst,
                                    overlapMeshes, uncoveredNodes, dstElemMass, dstElemGrad);
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 7: Handle surface nodes via linear interpolation.             ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Performing surface interpolation. " << endl;
  SurfaceInterpolation(geometry_src, geometry_dst, solver_src, solver_dst, coorDst);
}

void CConservativeVolumeInterpolator::PointLocalization(CGeometry* geometry_src,
                                                        const vector<su2double>& coor_dst,
                                                        vector<optional<unsigned long>>& containingElems,
                                                        vector<int>& containingElemRanks,
                                                        vector<unsigned long>& pointsFailed) {
  /*--- Search for containing elements for the given coordinates ---*/
  CADTElemClass& volumeADT = GetSourceVolumeADT();
  CADTElemClass& surfaceADT = GetSourceSurfaceADT();

  const unsigned long nDOFsDst = coor_dst.size() / nDim;

  /*--- Loop over the DOFs to be interpolated ---*/
  containingElems.clear();
  containingElemRanks.clear();
  pointsFailed.clear();

  /*--- Initialize with invalid values ---*/
  containingElems.resize(nDOFsDst, nullopt);
  containingElemRanks.resize(nDOFsDst, -1);

  for (auto l = 0u; l < nDOFsDst; ++l) {
    /*--- Set a pointer to the coordinates to be searched ---*/
    const su2double* coor = coor_dst.data() + l * nDim;

    /*--- Carry out the containment search and check if it was successful ---*/
    unsigned short subElemID;
    unsigned long elemID;
    int rankID;
    su2double parCoor[3], weightsInterpol[8];

    if (volumeADT.DetermineContainingElement(coor, subElemID, elemID, rankID, parCoor,
                                             weightsInterpol)) {
      /*--- Store element information ---*/
      containingElems[l] = elemID;
      containingElemRanks[l] = rankID;
    } else {
      /*--- No source volume element contains this node - fallback         ---*/
      /*--- The nearest volume element will be used to find intersections, ---*/
      /*--- but also add to the list of failed points, which will be       ---*/
      /*--- processed later via surface interpolation.                     ---*/

      pointsFailed.push_back(l);

      /*--- Find nearest volume element ---*/
      su2double dist;
      volumeADT.DetermineNearestElement(coor, dist, subElemID, elemID, rankID);
      containingElems[l] = elemID;
      containingElemRanks[l] = rankID;
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Volume search finished. " << pointsFailed.size() << " points failed." << endl << flush;
  }
}

void CConservativeVolumeInterpolator::ComputeSourceSolutionMass(CGeometry* geometry, CSolver* solver,
                                                                vector<vector<su2double>>& elemMass,
                                                                vector<vector<su2double>>& elemGrad) {
  elemMass.resize(nElem_src, vector<su2double>(nVar, 0.0));
  elemGrad.resize(nElem_src, vector<su2double>(nVar * nDim, 0.0));

  su2double vertexCoor[6];

  for (auto elemID = 0u; elemID < geometry->GetnElem(); ++elemID) {
    auto* elem = geometry->elem[elemID];
    const unsigned short nNodes = elem->GetnNodes();
    const unsigned short VTK_Type = elem->GetVTK_Type();
    const su2double elemVol = elem->GetVolume();

    /*--- Get solution at nodes ---*/
    vector<vector<su2double>> vertexSol(nNodes, vector<su2double>(nVar, 0.0));

    for (auto iNode = 0u; iNode < nNodes; ++iNode) {
      unsigned long nodeID = elem->GetNode(iNode);
      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        vertexSol[iNode][iVar] = solver->GetNodes()->GetSolution(nodeID, iVar);
      }
    }

    if (VTK_Type != TRIANGLE) continue;

    /*--- Get triangle vertex coordinates ---*/
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = elem->GetNode(iNode);
      vertexCoor[iNode * 2 + 0] = geometry->nodes->GetCoord(nodeID, 0);
      vertexCoor[iNode * 2 + 1] = geometry->nodes->GetCoord(nodeID, 1);
    }

    /*--- Compute mass and gradient ---*/
    ComputeTriangleMassAndGradient(vertexCoor, vertexSol, elemVol, nVar,
                                   elemMass[elemID], elemGrad[elemID]);
  }
}

void CConservativeVolumeInterpolator::CreateIntersectionMeshes(CGeometry* geometry_src,
                                                               CGeometry* geometry_dst,
                                                               const vector<su2double> &coor_dst,
                                                               IntersectionMeshMap& overlapMeshes,
                                                               vector<unsigned long>& incompleteOverlaps,
                                                               set<unsigned long>& uncoveredNodes) {
  overlapMeshes.clear();
  incompleteOverlaps.clear();
  uncoveredNodes.clear();

  /*--- Initialize all nodes as potentially uncovered ---*/
  for (auto l = 0u; l < nPoint_dst; ++l) {
    if (!geometry_dst->nodes->GetPhysicalBoundary(l))
      uncoveredNodes.insert(l);
  }

  /*--- Get references to both R-trees ---*/
  CRTreeSearchBase& srcRTree = GetSourceRTree();
  CRTreeSearchBase& dstRTree = GetDestRTree();

  /*--- Storage for triangle vertex coordinates ---*/
  su2double srcTri[6], dstTri[6];

  /*--- Storage for the R-tree searches ---*/
  set<unsigned long> containedNodes;
  set<unsigned long> candidateElems;

  /*--- Storage for result of intersection test ---*/
  vector<su2double> intersectionPoints;
  vector<su2double> intersectionElemCoords;
  vector<su2double> intersectionElemVols;

  /*--- Tolerances and counts for intersection volume check ---*/
  su2double absDiffTol[5] = {1e-10, 1e-8, 1e-6, 1e-4, 1e-2};
  su2double relDiffTol[5] = {2e-2, 5e-2, 1e-1, 2e-1, 5e-1};
  const su2double overlapTol = 0.95;
  vector<unsigned long> countAbsDiff(5, 0);
  vector<unsigned long> countRelDiff(5, 0);

  /*--- Only handle triangular elements for now ---*/
  if (nDim != 2) {
    if (rank == MASTER_NODE) {
      cout << "Warning: Element intersection currently only implemented for 2D. ";
      cout << "Skipping intersection computation for mesh." << endl;
    }
    return;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Build initial candidate list by searching for destination  ---*/
  /*---         nodes in source elements (reverse search).                 ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) {
    cout << "Building reverse candidate mapping (source elements -> destination nodes)." << endl;
  }

  /*--- Storage for reverse candidate mapping: candidateSrcElems[dstElemID] = {srcElemIDs} ---*/
  vector<set<unsigned long>> candidateSrcElems(geometry_dst->GetnElem());

  /*--- Loop over all source elements to build reverse mapping ---*/
  set<unsigned long> containedDstNodes;
  for (auto srcElemID = 0u; srcElemID < geometry_src->GetnElem(); ++srcElemID) {
    auto* srcElem = geometry_src->elem[srcElemID];

    /*--- Skip non-triangular source elements ---*/
    if (srcElem->GetVTK_Type() != TRIANGLE) continue;

    /*--- Get source triangle vertices ---*/
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = srcElem->GetNode(iNode);
      for (auto iDim = 0u; iDim < nDim; ++iDim)
        srcTri[iNode * nDim + iDim] = geometry_src->nodes->GetCoord(nodeID, iDim);
    }

    /*--- Find destination nodes within this source element's bounding box ---*/
    containedDstNodes.clear();
    dstRTree.SearchNodesInElement(srcTri, containedDstNodes, 1e-9);

    /*--- Add this source element as candidate for all destination elements containing these nodes ---*/
    for (auto dstNodeID : containedDstNodes) {
      const unsigned short nElem_node = geometry_dst->nodes->GetnElem(dstNodeID);
      for (auto jElem = 0; jElem < nElem_node; ++jElem) {
        unsigned long dstElemID = geometry_dst->nodes->GetElem(dstNodeID, jElem);
        candidateSrcElems[dstElemID].insert(srcElemID);
      }
    }
  }

  if (rank == MASTER_NODE) {
    unsigned long totalReverseCandidates = 0;
    for (const auto& candidates : candidateSrcElems) {
      totalReverseCandidates += candidates.size();
    }
    cout << "Reverse search found " << totalReverseCandidates << " total candidate element pairs." << endl;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Main loop over destination elements - combine forward and  ---*/
  /*---         reverse searches, then test for intersections.             ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) {
    cout << "Generating intersection region meshes." << endl;
  }

  unsigned long totalOverlaps = 0;
  unsigned long failedIntersections = 0;
  for (auto dstElemID = 0u; dstElemID < geometry_dst->GetnElem(); ++dstElemID) {
    auto* dstElem = geometry_dst->elem[dstElemID];

    /*--- Skip non-triangular elements ---*/
    if (dstElem->GetVTK_Type() != TRIANGLE) continue;

    /*--------------------------------------------------------------------------*/
    /*--- Step 2a: Get destination triangle vertices.                        ---*/
    /*--------------------------------------------------------------------------*/
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      for (auto iDim = 0u; iDim < nDim; ++iDim)
        dstTri[iNode * nDim + iDim] = geometry_dst->nodes->GetCoord(nodeID, iDim);
    }

    /*--------------------------------------------------------------------------*/
    /*--- Step 2b: Forward search - find source nodes in dest element bbox.  ---*/
    /*--------------------------------------------------------------------------*/
    containedNodes.clear();
    srcRTree.SearchNodesInElement(dstTri, containedNodes, 1e-8);

    /*--------------------------------------------------------------------------*/
    /*--- Step 2c: Create the list of candidate source elements. Use mesh    ---*/
    /*---          connectivity data to form the list of K_src that contain  ---*/
    /*---          points in containedNodes. Also add reverse candidates.    ---*/
    /*--------------------------------------------------------------------------*/
    candidateElems.clear();

    /*--- Forward search candidates: source elements containing nodes in dest bounding box ---*/
    for (auto srcNodeID : containedNodes) {
      const unsigned short nElem_node = geometry_src->nodes->GetnElem(srcNodeID);
      for (auto jElem = 0; jElem < nElem_node; ++jElem) {
        unsigned long srcElemID = geometry_src->nodes->GetElem(srcNodeID, jElem);
        candidateElems.insert(srcElemID);
      }
    }

    /*--- Reverse search candidates: source elements from precomputed mapping ---*/
    for (auto srcElemID : candidateSrcElems[dstElemID]) {
      candidateElems.insert(srcElemID);
    }

    /*--------------------------------------------------------------------------*/
    /*--- Step 2d: Test the candidate elements. For any which intersect,     ---*/
    /*---          mesh the intersection region.                             ---*/
    /*--------------------------------------------------------------------------*/
    su2double intersectionVol = 0.0;
    set<unsigned long> covered_vertices;
    for (auto srcElemID : candidateElems) {
      auto* srcElem = geometry_src->elem[srcElemID];

      /*--- Skip non-triangular source elements ---*/
      if (srcElem->GetVTK_Type() != TRIANGLE) continue;

      /*--- Get source triangle vertices ---*/
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        unsigned long nodeID = srcElem->GetNode(iNode);
        for (auto iDim = 0u; iDim < nDim; ++iDim)
          srcTri[iNode * nDim + iDim] = geometry_src->nodes->GetCoord(nodeID, iDim);
      }

      /*--- Check for intersection and detect new candidates ---*/
      intersectionPoints.clear();
      intersectionElemCoords.clear();
      intersectionElemVols.clear();

      bool isCoveredVertexP[3] = {};
      bool isCoveredVertexQ[3] = {};
      if (TriangleTriangleIntersection(dstTri, srcTri, intersectionPoints,
                                       intersectionElemCoords, intersectionElemVols,
                                      isCoveredVertexP, isCoveredVertexQ)) {
        overlapMeshes[dstElemID].push_back({srcElemID, intersectionElemCoords, intersectionElemVols});
        intersectionVol += accumulate(intersectionElemVols.begin(), intersectionElemVols.end(), 0.0);
        totalOverlaps++;

        for (auto iNode = 0u; iNode < 3; ++iNode) {
          if (isCoveredVertexP[iNode])
            covered_vertices.insert(dstElem->GetNode(iNode));
        }
      } else {
        /*--- Triangle intersection failed ---*/
        failedIntersections++;
      }
    }

    /*--- Check for incomplete overlap and apply standard treatment ---*/
    const su2double dstElemVol = dstElem->GetVolume();
    const su2double relDiff = abs(dstElemVol - intersectionVol) / dstElemVol;

    if (relDiff > overlapTol) {
      /*--- Mark element as having incomplete overlap ---*/
      incompleteOverlaps.push_back(dstElemID);
    }

    /*--- If this element has any intersection, mark its nodes as covered ---*/
    if (intersectionVol > 0.0) {
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        unsigned long dstNodeID = dstElem->GetNode(iNode);
        uncoveredNodes.erase(dstNodeID);
      }
      // for (auto dstNodeID : covered_vertices) {
      //   uncoveredNodes.erase(dstNodeID);
      // }
    }

    /*--- Compare total triangle area with destination element volume ---*/
    const su2double absDiff = abs(dstElemVol - intersectionVol);
    const su2double relDiff_stats = absDiff / dstElemVol;
    for (auto i = 0u; i < 5; ++i) {
      if (absDiff > absDiffTol[i]) countAbsDiff[i]++;
      if (relDiff_stats > relDiffTol[i]) countRelDiff[i]++;
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Total overlapping element pairs: " << totalOverlaps << "." << endl;
    cout << "Number of destination elements with overlaps: " << overlapMeshes.size() << endl;
    cout << "Number of failed triangle intersections: " << failedIntersections << "." << endl;
    cout << "Total destination elements processed: " << geometry_dst->GetnElem() << endl;
  }

  if (rank == MASTER_NODE) {
    cout << "Area conservation (absolute difference):" << endl;
    for (auto i = 0u; i < 5; ++i) {
      cout << "  Number exceeding " << scientific << setprecision(1);
      cout << absDiffTol[i] << ": ";
      cout << countAbsDiff[i] << endl;
    }
    cout << "Area conservation (relative difference):" << endl;
    for (auto i = 0u; i < 5; ++i) {
      cout << "  Number exceeding " << fixed << setprecision(0) << setw(3);
      cout << relDiffTol[i] * 100 << "%: ";
      cout << countRelDiff[i] << endl;
    }
    cout << "Elements with incomplete overlap: " << incompleteOverlaps.size() << endl;
    cout << "Nodes outside source domain: " << uncoveredNodes.size() << endl;
  }
}

bool CConservativeVolumeInterpolator::TriangleTriangleIntersection(const su2double dstTri[6],
                                                                   const su2double srcTri[6],
                                                                   vector<su2double>& intersectionPoints,
                                                                   vector<su2double>& intersectionElemCoords,
                                                                   vector<su2double>& intersectionElemVols,
                                                                   bool isCoveredVertexP[3],
                                                                   bool isCoveredVertexQ[3]) {
  const su2double EPS = 1e-9;

  /*--- Triangle vertices: P is destination and Q is source ---*/
  su2double P[3][2] = {{dstTri[0], dstTri[1]}, {dstTri[2], dstTri[3]}, {dstTri[4], dstTri[5]}};
  su2double Q[3][2] = {{srcTri[0], srcTri[1]}, {srcTri[2], srcTri[3]}, {srcTri[4], srcTri[5]}};

  /*--- Edge definitions: edge j connects vertex j to vertex (j+1)%3 ---*/
  /*--- Edge 0: Q0-Q1, Edge 1: Q1-Q2, Edge 2: Q2-Q0 for triangle Q   ---*/
  /*--- Edge 0: P0-P1, Edge 1: P1-P2, Edge 2: P2-P0 for triangle P   ---*/

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Compute all 18 vertex-edge powers handle degenerate        ---*/
  /*---         intersection cases (1, 2, or 4 signed distances being 0).  ---*/
  /*--------------------------------------------------------------------------*/
  /*--- power_P[i][j] = power of vertex Pi w.r.t. edge j of triangle Q ---*/
  /*--- power_Q[i][j] = power of vertex Qi w.r.t. edge j of triangle P ---*/
  su2double power_P[3][3], power_Q[3][3];
  for (auto i = 0u; i < 3; ++i) {
    for (auto j = 0u; j < 3; ++j) {
      power_P[i][j] = ComputeSignedDistance(P[i], Q[j], Q[(j+1)%3]);
      power_Q[i][j] = ComputeSignedDistance(Q[i], P[j], P[(j+1)%3]);
    }
  }

  /*--- Track what we find in each step for proper inclusion test logic ---*/
  bool hasInteriorVertices = false;

  /*--- Edge definitions ---*/
  su2double* P_edges[3][2] = {{P[0], P[1]}, {P[1], P[2]}, {P[2], P[0]}};
  su2double* Q_edges[3][2] = {{Q[0], Q[1]}, {Q[1], Q[2]}, {Q[2], Q[0]}};

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Handle degenerate edge-edge intersection cases.            ---*/
  /*--------------------------------------------------------------------------*/
  bool isDegenerateVertexP[3] = {};
  bool isDegenerateVertexQ[3] = {};
  bool isDegenerateEdgePair[3][3] = {};
  ProcessDegenerateEdgeIntersections(P_edges, Q_edges, power_P, power_Q, EPS, intersectionPoints,
                                     isDegenerateVertexP, isDegenerateVertexQ, isDegenerateEdgePair);

  /*--- Track covered vertices ---*/
  for (auto i = 0u; i < 3; ++i) {
    if (isDegenerateVertexP[i]) isCoveredVertexP[i] = true;
    if (isDegenerateVertexQ[i]) isCoveredVertexQ[i] = true;
  }

  /*--- Check if KP vertices are strictly inside KQ ---*/
  for (auto i = 0u; i < 3; ++i) {
    if ((power_P[i][0] > -EPS) &&
        (power_P[i][1] > -EPS) &&
        (power_P[i][2] > -EPS) &&
        (!isDegenerateVertexP[i])) {
      intersectionPoints.push_back(P[i][0]);
      intersectionPoints.push_back(P[i][1]);
      isCoveredVertexP[i] = true;
    }
  }

  /*--- Check if KQ vertices are strictly inside KP ---*/
  for (auto i = 0u; i < 3; ++i) {
    if ((power_Q[i][0] > -EPS) &&
        (power_Q[i][1] > -EPS) &&
        (power_Q[i][2] > -EPS) &&
        (!isDegenerateVertexQ[i])) {
      intersectionPoints.push_back(Q[i][0]);
      intersectionPoints.push_back(Q[i][1]);
      isCoveredVertexQ[i] = true;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Process edge-edge intersections, ignoring any degenerate   ---*/
  /*---         cases since they were handled in Step 1.                   ---*/
  /*--------------------------------------------------------------------------*/
  for (auto jQ = 0u; jQ < 3; ++jQ) {
    bool wasIntersected = false;
    for (auto iP = 0u; iP < 3; ++iP) {
      /*--- Skip degenerate edge pairs (already handled in Step 1) ---*/
      if (isDegenerateEdgePair[iP][jQ]) {
        wasIntersected = true;
        continue;
      }

      su2double* P0 = P_edges[iP][0];
      su2double* P1 = P_edges[iP][1];
      su2double* Q0 = Q_edges[jQ][0];
      su2double* Q1 = Q_edges[jQ][1];

      /*--- Compute signed distances for this edge pair ---*/
      su2double powerP0_eQ = ComputeSignedDistance(P0, Q0, Q1);
      su2double powerP1_eQ = ComputeSignedDistance(P1, Q0, Q1);
      su2double powerQ0_eP = ComputeSignedDistance(Q0, P0, P1);
      su2double powerQ1_eP = ComputeSignedDistance(Q1, P0, P1);

      /*--- Check for edge-edge intersection ---*/
      if ((powerP0_eQ * powerP1_eQ < 0) && (powerQ0_eP * powerQ1_eP < 0)) {
        /*--- Compute intersection point ---*/
        su2double t = powerP0_eQ / (powerP0_eQ - powerP1_eQ);

        su2double intersectionX = P0[0] + t * (P1[0] - P0[0]);
        su2double intersectionY = P0[1] + t * (P1[1] - P0[1]);

        intersectionPoints.push_back(intersectionX);
        intersectionPoints.push_back(intersectionY);

        wasIntersected = true;
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Mesh the intersection polygon.                             ---*/
  /*--------------------------------------------------------------------------*/
  MeshConvexPolygon(intersectionPoints, intersectionElemCoords, intersectionElemVols);

  return (!intersectionElemCoords.empty());
}

void CConservativeVolumeInterpolator::ProcessDegenerateEdgeIntersections(su2double* P_edges[3][2], su2double* Q_edges[3][2],
                                                                         const su2double power_P[3][3], const su2double power_Q[3][3],
                                                                         const su2double EPS,
                                                                         vector<su2double>& intersectionPoints,
                                                                         bool isDegenerateVertexP[3],
                                                                         bool isDegenerateVertexQ[3],
                                                                         bool isDegenerateEdgePair[3][3]) {
  auto addPoint = [&](su2double x, su2double y, int Pi, int Qj) {
    bool alreadyAdded = false;
    if (Pi >= 0) {
      if (isDegenerateVertexP[Pi]) alreadyAdded = true;
      isDegenerateVertexP[Pi] = true;
    }
    if (Qj >= 0) {
      if (isDegenerateVertexQ[Qj]) alreadyAdded = true;
      isDegenerateVertexQ[Qj] = true;
    }
    if (!alreadyAdded) {
      intersectionPoints.push_back(x);
      intersectionPoints.push_back(y);
    }
  };

  /*--- Process all 9 edge-edge combinations for degenerate cases ---*/
  for (auto iP = 0u; iP < 3; ++iP) {
    for (auto jQ = 0u; jQ < 3; ++jQ) {
      /*--- Indices of nodes belonging to current edges ---*/
      int iP0 = iP, iP1 = (iP + 1) % 3;
      int jQ0 = jQ, jQ1 = (jQ + 1) % 3;

      /*--- Coordinates of edge nodes ---*/
      su2double* P0 = P_edges[iP][0];  // P vertex iP
      su2double* P1 = P_edges[iP][1];  // P vertex (iP+1)%3
      su2double* Q0 = Q_edges[jQ][0];  // Q vertex jQ
      su2double* Q1 = Q_edges[jQ][1];  // Q vertex (jQ+1)%3

      /*--- Use precomputed powers instead of recalculating ---*/
      /*--- Edge iP of triangle P goes from vertex iP to vertex (iP+1)%3 ---*/
      /*--- Edge jQ of triangle Q goes from vertex jQ to vertex (jQ+1)%3 ---*/
      su2double powerP0_eQ = power_P[iP0][jQ];  // Power of P vertex iP w.r.t. Q edge jQ
      su2double powerP1_eQ = power_P[iP1][jQ];  // Power of P vertex (iP+1)%3 w.r.t. Q edge jQ
      su2double powerQ0_eP = power_Q[jQ0][iP];  // Power of Q vertex jQ w.r.t. P edge iP
      su2double powerQ1_eP = power_Q[jQ1][iP];  // Power of Q vertex (jQ+1)%3 w.r.t. P edge iP

      /*--- Count how many powers are zero ---*/
      int zeroCount = 0;
      bool P0_zero = (abs(powerP0_eQ) < EPS);
      bool P1_zero = (abs(powerP1_eQ) < EPS);
      bool Q0_zero = (abs(powerQ0_eP) < EPS);
      bool Q1_zero = (abs(powerQ1_eP) < EPS);

      if (P0_zero) zeroCount++;
      if (P1_zero) zeroCount++;
      if (Q0_zero) zeroCount++;
      if (Q1_zero) zeroCount++;

      /*--- Mark edge pair as degenerate if any zero powers ---*/
      if (zeroCount > 0) isDegenerateEdgePair[iP][jQ] = true;

      /*--------------------------------------------------------------------------*/
      /*--- Case 1: Only one power is zero                                     ---*/
      /*--------------------------------------------------------------------------*/
      if (zeroCount == 1) {
        if (P0_zero) {
          /*--- P0 lies on edge Q, intersection if Q0 and Q1 on opposite sides of edge P ---*/
          if (powerQ0_eP * powerQ1_eP < 0) {
            addPoint(P0[0], P0[1], iP0, -1);
          }
        } else if (P1_zero) {
          /*--- P1 lies on edge Q, intersection if Q0 and Q1 on opposite sides of edge P ---*/
          if (powerQ0_eP * powerQ1_eP < 0) {
            addPoint(P1[0], P1[1], iP1, -1);
          }
        } else if (Q0_zero) {
          /*--- Q0 lies on edge P, intersection if P0 and P1 on opposite sides of edge Q ---*/
          if (powerP0_eQ * powerP1_eQ < 0) {
            addPoint(Q0[0], Q0[1], -1, jQ0);
          }
        } else if (Q1_zero) {
          /*--- Q1 lies on edge P, intersection if P0 and P1 on opposite sides of edge Q ---*/
          if (powerP0_eQ * powerP1_eQ < 0) {
            addPoint(Q1[0], Q1[1], -1, jQ1);
          }
        }
      } // if zeroCount == 1

      /*--------------------------------------------------------------------------*/
      /*--- Case 2: Two powers are zero, one for each edge                     ---*/
      /*--------------------------------------------------------------------------*/
      else if (zeroCount == 2) {
        if (P0_zero && Q0_zero) {
          /*--- P0 = Q0, common vertex ---*/
          addPoint(P0[0], P0[1], iP0, jQ0);
        } else if (P0_zero && Q1_zero) {
          /*--- P0 = Q1, common vertex ---*/
          addPoint(P0[0], P0[1], iP0, jQ1);
        } else if (P1_zero && Q0_zero) {
          /*--- P1 = Q0, common vertex ---*/
          addPoint(P1[0], P1[1], iP1, jQ0);
        } else if (P1_zero && Q1_zero) {
          /*--- P1 = Q1, common vertex ---*/
          addPoint(P1[0], P1[1], iP1, jQ1);
        }
      } // if zeroCount == 2

      /*--------------------------------------------------------------------------*/
      /*--- Case 3: All four powers are zero (edges are collinear)             ---*/
      /*--------------------------------------------------------------------------*/
      else if (zeroCount == 4) {
        /*--- Edges are collinear, check for overlap using parametrization ---*/
        /*--- Parametrize both edges and find overlap in parameter space ---*/

        /*--- Edge P: P0 + t * (P1 - P0), t in [0,1] ---*/
        /*--- Edge Q: Q0 + s * (Q1 - Q0), s in [0,1] ---*/

        /*--- Direction vector of edge P ---*/
        su2double dx_P = P1[0] - P0[0];
        su2double dy_P = P1[1] - P0[1];
        su2double lengthP_sq = dx_P*dx_P + dy_P*dy_P;

        /*--- Find parameter values where Q vertices lie on P edge ---*/
        /*--- Q0 = P0 + s_Q0 * (P1 - P0) ---*/
        /*--- Q1 = P0 + s_Q1 * (P1 - P0) ---*/

        su2double s_Q0 = ((Q0[0] - P0[0]) * dx_P +
                          (Q0[1] - P0[1]) * dy_P) / lengthP_sq;
        su2double s_Q1 = ((Q1[0] - P0[0]) * dx_P +
                          (Q1[1] - P0[1]) * dy_P) / lengthP_sq;

        /*--- Edge P spans parameter interval [0, 1] ---*/
        /*--- Edge Q spans parameter interval [min(s_Q0, s_Q1), max(s_Q0, s_Q1)] ---*/
        su2double s_min = min(s_Q0, s_Q1);
        su2double s_max = max(s_Q0, s_Q1);

        /*--- Find overlap of intervals [0, 1] and [s_min, s_max] ---*/
        su2double overlap_start = max(0.0, s_min);
        su2double overlap_end = min(1.0, s_max);

        if (overlap_start < overlap_end) {
          /*--- There is overlap, add intersection points ---*/
          /*--- Add start point of overlap ---*/
          su2double X_start = P0[0] + overlap_start * dx_P;
          su2double Y_start = P0[1] + overlap_start * dy_P;

          /*--- Check if start point corresponds to a known vertex ---*/
          int Pi_start = -1, Qj_start = -1;
          if (abs(overlap_start - 0.0) < EPS) {
            Pi_start = iP0;  // Start point is P0
          } else if (abs(overlap_start - 1.0) < EPS) {
            Pi_start = iP1;  // Start point is P1
          }
          if (abs(s_Q0 - overlap_start) < EPS) {
            Qj_start = jQ0;  // Start point is Q0
          } else if (abs(s_Q1 - overlap_start) < EPS) {
            Qj_start = jQ1;  // Start point is Q1
          }
          addPoint(X_start, Y_start, Pi_start, Qj_start);

          /*--- Add end point of overlap if different from start ---*/
          if (abs(overlap_end - overlap_start) > EPS) {
            su2double X_end = P0[0] + overlap_end * dx_P;
            su2double Y_end = P0[1] + overlap_end * dy_P;

            /*--- Check if end point corresponds to a known vertex ---*/
            int Pi_end = -1, Qj_end = -1;
            if (abs(overlap_end - 0.0) < EPS) {
              Pi_end = iP0;  // End point is P0
            } else if (abs(overlap_end - 1.0) < EPS) {
              Pi_end = iP1;  // End point is P1
            }
            if (abs(s_Q0 - overlap_end) < EPS) {
              Qj_end = jQ0;  // End point is Q0
            } else if (abs(s_Q1 - overlap_end) < EPS) {
              Qj_end = jQ1;  // End point is Q1
            }
            addPoint(X_end, Y_end, Pi_end, Qj_end);
          } // if multiple overlap
        } // if overlap
      } // if zeroCount == 4
    } // for jQ
  } // for iP
}

su2double CConservativeVolumeInterpolator::ComputeSignedDistance(const su2double point[2],
                                                                 const su2double lineStart[2],
                                                                 const su2double lineEnd[2]) {
  /*--- Unit normal of line ---*/
  su2double Nx = lineStart[1] - lineEnd[1];
  su2double Ny = lineEnd[0] - lineStart[0];
  const su2double mag = sqrt(Nx * Nx + Ny * Ny);
  Nx /= mag;
  Ny /= mag;

  /*--- Compute signed distance using dot product between [P P{i+1}] and N ---*/
  const su2double Px = point[0] - lineStart[0];
  const su2double Py = point[1] - lineStart[1];

  return Px * Nx + Py * Ny;
}

void CConservativeVolumeInterpolator::MeshConvexPolygon(const vector<su2double>& polygonPoints,
                                                        vector<su2double>& polygonElemCoords,
                                                        vector<su2double>& polygonElemVols) {
  polygonElemCoords.clear();

  unsigned int numPolygonPoint = polygonPoints.size() / 2;

  /*--- Handle edge cases ---*/
  if (numPolygonPoint < 3) {
    return;  // Cannot mesh a polygon with less than 3 points
  }

  /*--- Lambda to add a triangle with positive area and return orientation ---*/
  auto addPositiveTriangle = [&](su2double P0[2], su2double P1[2], su2double P2[2]) -> bool {
    /*--- Check orientation ---*/
    const su2double cross = (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P1[1] - P0[1]) * (P2[0] - P0[0]);
    const su2double vol = 0.5 * abs(cross);
    if (cross > 0) {
      /*--- Counter-clockwise orientation - add triangle as is ---*/
      polygonElemCoords.push_back(P0[0]); polygonElemCoords.push_back(P0[1]);
      polygonElemCoords.push_back(P1[0]); polygonElemCoords.push_back(P1[1]);
      polygonElemCoords.push_back(P2[0]); polygonElemCoords.push_back(P2[1]);
      polygonElemVols.push_back(vol);
      return true;  // Original order was counter-clockwise
    } else {
      /*--- Clockwise orientation - reverse order ---*/
      polygonElemCoords.push_back(P0[0]); polygonElemCoords.push_back(P0[1]);
      polygonElemCoords.push_back(P2[0]); polygonElemCoords.push_back(P2[1]);
      polygonElemCoords.push_back(P1[0]); polygonElemCoords.push_back(P1[1]);
      polygonElemVols.push_back(vol);
      return false; // Original order was clockwise, had to reverse
    }
  };

  /*--- Mesh the first triangle ---*/
  su2double P0[2] = {polygonPoints[0], polygonPoints[1]};
  su2double P1[2] = {polygonPoints[2], polygonPoints[3]};
  su2double P2[2] = {polygonPoints[4], polygonPoints[5]};
  bool isCounterClockwise = addPositiveTriangle(P0, P1, P2);

  if (numPolygonPoint == 3) {
    /*--- Already a triangle, so return ---*/
    return;
  }

  /*--- For a convex polygon with n >= 4 points, create (n-2) triangles. ---*/
  /*--- Use Alauzet's incremental triangulation algorithm               ---*/

  /*--- Track boundary edges: each edge stores {startPoint, endPoint} ---*/
  /*--- Initially, the boundary consists of the edges of the first triangle ---*/
  /*--- Store edges in consistent orientation based on first triangle ---*/
  vector<pair<int, int>> boundaryEdges;
  if (isCounterClockwise) {
    /*--- Counter-clockwise: P0->P1->P2->P0 ---*/
    boundaryEdges.push_back({0, 1});  // Edge P0->P1
    boundaryEdges.push_back({1, 2});  // Edge P1->P2
    boundaryEdges.push_back({2, 0});  // Edge P2->P0
  } else {
    /*--- Clockwise: P0->P2->P1->P0 ---*/
    boundaryEdges.push_back({0, 2});  // Edge P0->P2
    boundaryEdges.push_back({2, 1});  // Edge P2->P1
    boundaryEdges.push_back({1, 0});  // Edge P1->P0
  }

  /*--- Add remaining points one by one ---*/
  for (auto i = 3u; i < numPolygonPoint; ++i) {
    su2double P[2] = {polygonPoints[i * 2 + 0], polygonPoints[i * 2 + 1]};

    /*--- Find the unique boundary edge that "views" this point (negative signed distance) ---*/
    int jj = -1;
    for (auto j = 0u; j < boundaryEdges.size(); ++j) {
      int idx1 = boundaryEdges[j].first;
      int idx2 = boundaryEdges[j].second;

      su2double testP1[2] = {polygonPoints[idx1 * 2 + 0], polygonPoints[idx1 * 2 + 1]};
      su2double testP2[2] = {polygonPoints[idx2 * 2 + 0], polygonPoints[idx2 * 2 + 1]};

      if (ComputeSignedDistance(P, testP1, testP2) < 0) {
        jj = j;
        break;
      }
    }

    /*--- Should always find exactly one boundary edge for a convex polygon ---*/
    if (jj >= 0) {
      /*--- Swap the indices since the signed distance was negative ---*/
      int idx2 = boundaryEdges[jj].first;
      int idx1 = boundaryEdges[jj].second;

      su2double testP1[2] = {polygonPoints[idx1 * 2 + 0], polygonPoints[idx1 * 2 + 1]};
      su2double testP2[2] = {polygonPoints[idx2 * 2 + 0], polygonPoints[idx2 * 2 + 1]};

      /*--- Create triangle with the new point and the boundary edge ---*/
      isCounterClockwise = addPositiveTriangle(P, testP1, testP2);

      /*--- Update boundary: remove the used edge and add two new boundary edges ---*/
      /*--- Maintain consistent orientation based on the triangle orientation ---*/
      boundaryEdges.erase(boundaryEdges.begin() + jj);

      /*--- Triangle P->edgeStart->edgeEnd is counter-clockwise ---*/
      /*--- Order should be correct by construction ---*/
      boundaryEdges.push_back({i, idx1});  // Edge from new point to start of used edge
      boundaryEdges.push_back({idx2, i});  // Edge from end of used edge to new point
    }
  }
}

void CConservativeVolumeInterpolator::ComputeDestinationMassAndGradient(CGeometry* geometry_src,
                                                                        CGeometry* geometry_dst,
                                                                        CSolver* solver_src,
                                                                        const IntersectionMeshMap& overlapMeshes,
                                                                        const vector<unsigned long>& incompleteOverlaps,
                                                                        const vector<vector<su2double>>& srcElemMass,
                                                                        const vector<vector<su2double>>& srcElemGrad,
                                                                        vector<vector<su2double>>& dstElemMass,
                                                                        vector<vector<su2double>>& dstElemGrad) {
  dstElemMass.resize(nElem_dst, vector<su2double>(nVar, 0.0));
  dstElemGrad.resize(nElem_dst, vector<su2double>(nVar * nDim, 0.0));

  /*--- Loop over all destination elements that have intersections ---*/
  const unsigned short nCoorPerElem = (nDim == 2)? 6 : 12;

  /*--- Counters for non-matching boundary treatment ---*/
  unsigned long totalBoundaryElems = 0;
  unsigned long nonConservativeTreatment = 0;

  for (const auto& intersection : overlapMeshes) {
    unsigned long dstElemID = intersection.first;
    const IntersectionMesh& srcElemMeshes = intersection.second;

    /*--- Initialize destination element mass and gradient ---*/
    const auto* dstElem = geometry_dst->elem[dstElemID];
    auto& dstMass = dstElemMass[dstElemID];
    auto& dstGrad = dstElemGrad[dstElemID];
    fill(dstMass.begin(), dstMass.end(), 0.0);
    fill(dstGrad.begin(), dstGrad.end(), 0.0);

    /*--- Process each intersection region T_j = intersection(K_dst, K_src_j) ---*/
    su2double intersectionVol = 0.0;

    for (const auto& srcElemMesh : srcElemMeshes) {
      unsigned long srcElemID = srcElemMesh.srcElemID;
      const vector<su2double>& triElemCoords = srcElemMesh.coords;
      const vector<su2double>& triElemVols = srcElemMesh.vols;
      unsigned int numTri = triElemCoords.size() / nCoorPerElem;

      for (auto iTri = 0u; iTri < numTri; ++iTri) {
        /*--- Get triangle vertices ---*/
        const su2double* coor_tri = triElemCoords.data() + iTri * nCoorPerElem;
        const su2double x0 = coor_tri[0], y0 = coor_tri[1];
        const su2double x1 = coor_tri[2], y1 = coor_tri[3];
        const su2double x2 = coor_tri[4], y2 = coor_tri[5];

        /*--- Use precomputed area from volume data ---*/
        const su2double area = triElemVols[iTri];
        intersectionVol += area;

        /*--- Get source element data ---*/
        const auto* srcElem = geometry_src->elem[srcElemID];
        const su2double srcVol = srcElem->GetVolume();
        const su2double* G_K_src = srcElem->GetCG();

        /*--- Triangle centroid ---*/
        const su2double xi = (x0 + x1 + x2) / 3.0;
        const su2double yi = (y0 + y1 + y2) / 3.0;

        /*--- Evaluate solution and gradient at triangle centroid ---*/
        for (auto iVar = 0u; iVar < nVar; ++iVar) {
          const su2double u_src = srcElemMass[srcElemID][iVar] / srcVol;
          const su2double* grad_u = srcElemGrad[srcElemID].data() + iVar * nDim;

          /*--- Displacement from source centroid to triangle centroid ---*/
          const su2double dx = xi - G_K_src[0];
          const su2double dy = yi - G_K_src[1];

          /*--- Solution at triangle centroid ---*/
          const su2double u_tri = u_src + grad_u[0] * dx + grad_u[1] * dy;

          /*--- Accumulate mass ---*/
          dstMass[iVar] += u_tri * area;

          /*--- Accumulate gradient ---*/
          dstGrad[iVar * nDim + 0] += grad_u[0] * area;
          dstGrad[iVar * nDim + 1] += grad_u[1] * area;
        }
      }
    }

    /*--- Apply standard treatment: normalize by intersection volume, and ---*/
    /*--- scale mass to destination element volume ---*/
    const su2double dstElemVol = dstElem->GetVolume();
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      dstMass[iVar] = dstMass[iVar] * dstElemVol / intersectionVol;
      dstGrad[iVar * nDim + 0] /= intersectionVol;
      dstGrad[iVar * nDim + 1] /= intersectionVol;
    }
  }
}

void CConservativeVolumeInterpolator::ApplyMaximumPrincipleCorrection(CGeometry* geometry_src,
                                                                      CGeometry* geometry_dst,
                                                                      CSolver* solver_src,
                                                                      const vector<su2double>& coor_dst,
                                                                      const IntersectionMeshMap& overlapMeshes,
                                                                      const vector<vector<su2double>>& srcElemMass,
                                                                      const vector<vector<su2double>>& srcElemGrad,
                                                                      vector<vector<su2double>>& dstElemMass,
                                                                      vector<vector<su2double>>& dstElemGrad) {
  const su2double EPS = 1e-16;

  /*--- Loop over all destination elements that have overlaps ---*/
  su2double u_tilde[3];
  vector<su2double> correctedMass(1, 0.0);
  vector<su2double> correctedGrad(1 * nDim, 0.0);
  vector<vector<su2double>> vertexSol(3, vector<su2double>(1));
  for (const auto& intersection : overlapMeshes) {
    unsigned long dstElemID = intersection.first;
    const IntersectionMesh& srcElemMeshes = intersection.second;

    auto* dstElem = geometry_dst->elem[dstElemID];
    if (dstElem->GetVTK_Type() != TRIANGLE) continue;

    /*--- Get destination element vertices ---*/
    su2double dstVertices[6];
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      dstVertices[iNode * 2 + 0] = coor_dst[nodeID * nDim + 0];
      dstVertices[iNode * 2 + 1] = coor_dst[nodeID * nDim + 1];
    }

    /*--- Element volume and centroid (barycenter G_K) ---*/
    const su2double elemVol = dstElem->GetVolume();
    const su2double* G_K = dstElem->GetCG();

    /*--- For each variable, apply Alauzet's maximum principle correction ---*/
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      /*--------------------------------------------------------------------------*/
      /*--- Step 1: Compute local bounds from overlapping source elements.     ---*/
      /*--------------------------------------------------------------------------*/
      su2double u_min = 1e20;
      su2double u_max = -1e20;

      /*--- Find all vertices Q from source elements K_src that K overlaps ---*/
      for (const auto& srcElemMesh : srcElemMeshes) {
        unsigned long srcElemID = srcElemMesh.srcElemID;
        auto* srcElem = geometry_src->elem[srcElemID];
        if (srcElem->GetVTK_Type() != TRIANGLE) continue;

        /*--- Get solution values at vertices of source element ---*/
        for (auto iNode = 0u; iNode < 3; ++iNode) {
          unsigned long nodeID = srcElem->GetNode(iNode);
          su2double u_vertex = solver_src->GetNodes()->GetSolution(nodeID, iVar);
          u_min = min(u_min, u_vertex);
          u_max = max(u_max, u_vertex);
        }
      }

      /*--- Skip if no valid bounds found ---*/
      if (u_min > 1e19 || u_max < -1e19) continue;

      /*--------------------------------------------------------------------------*/
      /*--- Step 2: Get current solution at destination element.               ---*/
      /*--------------------------------------------------------------------------*/
      const su2double u_G = dstElemMass[dstElemID][iVar] / elemVol;
      const su2double* gradu_G = dstElemGrad[dstElemID].data() + iVar * nDim;

      /*--------------------------------------------------------------------------*/
      /*--- Step 3: Compute u_K(P_i) at each vertex using Taylor expansion.    ---*/
      /*--------------------------------------------------------------------------*/
      su2double u_K_P[3];  // Values at vertices P_0, P_1, P_2
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        /*--- Vector G_K P_i ---*/
        const su2double dx = dstVertices[iNode * 2 + 0] - G_K[0];
        const su2double dy = dstVertices[iNode * 2 + 1] - G_K[1];

        /*--- u_K(P_i) = u_K(G_K) + gra(u_K) * G_K P_i ---*/
        u_K_P[iNode] = u_G + gradu_G[0] * dx + gradu_G[1] * dy;
      }

      /*--------------------------------------------------------------------------*/
      /*--- Step 4: Check if maximum principle is violated.                    ---*/
      /*--------------------------------------------------------------------------*/
      bool violatesMaxPrinciple = false;
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        if (u_K_P[iNode] < u_min - EPS || u_K_P[iNode] > u_max + EPS) {
          violatesMaxPrinciple = true;
          break;
        }
      }

      /*--- If no violation, skip correction ---*/
      if (!violatesMaxPrinciple) continue;

      /*--------------------------------------------------------------------------*/
      /*--- Step 5: Apply Alauzet's correction algorithm.                      ---*/
      /*--------------------------------------------------------------------------*/
      /*--- Sort vertices by solution value: u_K(P_0) ≤ u_K(P_1) ≤ u_K(P_2) ---*/
      vector<pair<su2double, unsigned short>> sortedValues;
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        sortedValues.push_back(make_pair(u_K_P[iNode], iNode));
      }
      sort(sortedValues.begin(), sortedValues.end());

      const su2double u_P0 = sortedValues[0].first;  // Smallest value
      const su2double u_P1 = sortedValues[1].first;  // Middle value
      const su2double u_P2 = sortedValues[2].first;  // Largest value

      /*--- Apply first correction pass ---*/
      const su2double u_M_P2 = min(u_P2, u_max);
      const su2double u_M_P1 = min(u_P1 + 0.5 * max(0.0, u_P2 - u_max), u_max);
      const su2double u_M_P0 = 3.0 * u_G - u_M_P1 - u_M_P2;

      /*--- Apply second correction pass ---*/
      const su2double u_tilde_P0 = max(u_M_P0, u_min);
      const su2double u_tilde_P1 = max(u_M_P1 - 0.5 * max(0.0, u_min - u_M_P0), u_min);
      const su2double u_tilde_P2 = 3.0 * u_G - u_tilde_P0 - u_tilde_P1;

      /*--------------------------------------------------------------------------*/
      /*--- Step 6: Compute corrected mass and gradient from new nodal values. ---*/
      /*--------------------------------------------------------------------------*/
      /*--- Put corrected values back in original vertex order ---*/
      u_tilde[sortedValues[0].second] = u_tilde_P0;
      u_tilde[sortedValues[1].second] = u_tilde_P1;
      u_tilde[sortedValues[2].second] = u_tilde_P2;

      /*--- Prepare vertex solutions for single variable ---*/
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        vertexSol[iNode][0] = u_tilde[iNode];
      }

      /*--- Compute corrected mass and gradient ---*/
      ComputeTriangleMassAndGradient(dstVertices, vertexSol, elemVol, 1,
                                     correctedMass, correctedGrad);

      /*--- Update destination element data ---*/
      dstElemMass[dstElemID][iVar] = correctedMass[0];
      dstElemGrad[dstElemID][iVar * nDim + 0] = correctedGrad[0];
      dstElemGrad[dstElemID][iVar * nDim + 1] = correctedGrad[1];
    }
  }
}

void CConservativeVolumeInterpolator::DistributeSolutionToNodes(CGeometry* geometry_dst,
                                                                CSolver* solver_dst,
                                                                const vector<su2double> &coor_dst,
                                                                const vector<vector<su2double>>& dstElemMass,
                                                                const vector<vector<su2double>>& dstElemGrad) {
  /*--- Loop over all nodes and accumulate contributions from each element ---*/
  vector<su2double> totalValue(nVar, 0.0);
  for (auto l = 0u; l < nPoint_dst; ++l) {
    /*--- Initialize accumulation variables for this node ---*/
    su2double totalWeight = 0.0;
    fill(totalValue.begin(), totalValue.end(), 0.0);

    /*--- Get node coordinates ---*/
    const su2double* coor = coor_dst.data() + l * nDim;

    /*--- Loop over all elements that contain this node ---*/
    const unsigned short nElem_node = geometry_dst->nodes->GetnElem(l);
    for (auto j = 0u; j < nElem_node; ++j) {
      unsigned long elemID = geometry_dst->nodes->GetElem(l, j);
      auto* elem = geometry_dst->elem[elemID];

      /*--- Skip non-triangular elements ---*/
      if (elem->GetVTK_Type() != TRIANGLE) continue;

      /*--- Get element volume and centroid ---*/
      const su2double elemVol = elem->GetVolume();
      const su2double* G_K = elem->GetCG();

      /*--- Vector from centroid to node ---*/
      const su2double vec[2] = {coor[0] - G_K[0], coor[1] - G_K[1]};

      /*--- Loop over variables ---*/
      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        /*--- Get element-centered solution ---*/
        const su2double u_G = dstElemMass[elemID][iVar] / elemVol;
        const su2double* gradu_G = dstElemGrad[elemID].data() + iVar * nDim;

        /*--- Linear reconstruction: u(P_i) = u(G_K) + gra(u) * (P_i - G_K) ---*/
        const su2double vertexValue = u_G + gradu_G[0] * vec[0] + gradu_G[1] * vec[1];

        /*--- Accumulate weighted contribution ---*/
        totalValue[iVar] += vertexValue * elemVol;
      }

      /*--- Accumulate weight ---*/
      totalWeight += elemVol;
    }

    /*--- Compute weighted average and set solution for this node ---*/
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      const su2double avgValue = totalValue[iVar] / totalWeight;
      if (isnan(avgValue)) SU2_MPI::Error(string("Invalid mass at destination node ") + to_string(l) , CURRENT_FUNCTION);

      solver_dst->GetNodes()->SetSolution(l, iVar, avgValue);
      solver_dst->GetNodes()->SetSolution_Old(l, iVar, avgValue);
    }
  }
}

void CConservativeVolumeInterpolator::ComputeTriangleMassAndGradient(const su2double vertexCoor[6],
                                                                     const vector<vector<su2double>>& vertexSol,
                                                                     const su2double elemVol,
                                                                     const unsigned short nVar,
                                                                     vector<su2double>& mass,
                                                                     vector<su2double>& grad) {
  /*--- Initialize output ---*/
  fill(mass.begin(), mass.end(), 0.0);
  fill(grad.begin(), grad.end(), 0.0);

  /*--- Extract triangle vertices ---*/
  const su2double x0 = vertexCoor[0], y0 = vertexCoor[1];
  const su2double x1 = vertexCoor[2], y1 = vertexCoor[3];
  const su2double x2 = vertexCoor[4], y2 = vertexCoor[5];

  /*--- Calculate determinant and area ---*/
  const su2double det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
  const su2double area = 0.5 * abs(det);

  /*--- Loop over variables ---*/
  for (auto iVar = 0u; iVar < nVar; ++iVar) {
    /*--- Get solution values at vertices ---*/
    const su2double u0 = vertexSol[0][iVar];
    const su2double u1 = vertexSol[1][iVar];
    const su2double u2 = vertexSol[2][iVar];

    /*--- Compute mass: area * average value ---*/
    mass[iVar] = area * (u0 + u1 + u2) / 3.0;

    /*--- Compute gradient analytically ---*/
    const su2double grad_x = ((u1 - u0) * (y2 - y0) - (u2 - u0) * (y1 - y0)) / det;
    const su2double grad_y = ((u2 - u0) * (x1 - x0) - (u1 - u0) * (x2 - x0)) / det;

    grad[iVar * nDim + 0] = grad_x;
    grad[iVar * nDim + 1] = grad_y;
  }
}

void CConservativeVolumeInterpolator::ComputeTriangleMassAndGradientFEM(const su2double vertexCoor[6],
                                                                        const vector<vector<su2double>>& vertexSol,
                                                                        const su2double elemVol,
                                                                        const unsigned short nVar,
                                                                        vector<su2double>& mass,
                                                                        vector<su2double>& grad) {
  /*--- Initialize output ---*/
  fill(mass.begin(), mass.end(), 0.0);
  fill(grad.begin(), grad.end(), 0.0);

  /*--- Get integration points and basis functions from FEM standard element ---*/
  CFEMStandardElement& stdElement = GetFEMStandardElement();
  const unsigned short nInt = stdElement.GetNIntegration();
  const su2double* weights = stdElement.GetWeightsIntegration();
  const su2double* lagBasis = stdElement.GetBasisFunctionsIntegration();
  const su2double* drLagBasis = stdElement.GetDrBasisFunctionsIntegration();
  const su2double* dsLagBasis = stdElement.GetDsBasisFunctionsIntegration();

  /*--- Loop over integration points ---*/
  for (auto iInt = 0u; iInt < nInt; ++iInt) {
    /*--- Compute Jacobian of transformation ---*/
    su2double dxdr = 0.0, dydr = 0.0;
    su2double dxds = 0.0, dyds = 0.0;

    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned short ind = iInt * 3 + iNode;

      dxdr += vertexCoor[iNode * 2 + 0] * drLagBasis[ind];
      dydr += vertexCoor[iNode * 2 + 1] * drLagBasis[ind];
      dxds += vertexCoor[iNode * 2 + 0] * dsLagBasis[ind];
      dyds += vertexCoor[iNode * 2 + 1] * dsLagBasis[ind];
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
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      /*--- Interpolate solution and its gradient at integration point ---*/
      su2double solVal = 0.0;
      su2double dudr = 0.0, duds = 0.0;

      for (auto iNode = 0u; iNode < 3; ++iNode) {
        unsigned short ind = iInt * 3 + iNode;
        su2double nodeVal = vertexSol[iNode][iVar];

        solVal += lagBasis[ind] * nodeVal;
        dudr += drLagBasis[ind] * nodeVal;
        duds += dsLagBasis[ind] * nodeVal;
      }

      /*--- Add contribution to mass integral ---*/
      mass[iVar] += solVal * intWeight;

      /*--- Transform gradient to physical coordinates and add to gradient integral ---*/
      su2double dudx = dudr * drdx + duds * dsdx;
      su2double dudy = dudr * drdy + duds * dsdy;

      grad[iVar * nDim + 0] += dudx * intWeight;
      grad[iVar * nDim + 1] += dudy * intWeight;
    }
  }

  /*--- Normalize gradients by volume ---*/
  for (auto iVar = 0u; iVar < nVar; ++iVar) {
    for (auto k = 0u; k < nDim; ++k) {
      grad[iVar * nDim + k] /= elemVol;
    }
  }
}

void CConservativeVolumeInterpolator::DistributeSolutionOutsideDomain(CGeometry* geometry_src,
                                                                      CGeometry* geometry_dst,
                                                                      CSolver* solver_src,
                                                                      CSolver* solver_dst,
                                                                      const vector<su2double>& coor_dst,
                                                                      const IntersectionMeshMap& overlapMeshes,
                                                                      const set<unsigned long>& uncoveredNodes,
                                                                      const vector<vector<su2double>>& dstElemMass,
                                                                      const vector<vector<su2double>>& dstElemGrad) {
  unsigned long nVolumeExtrapolated = 0;

  for (const auto nodeID : uncoveredNodes) {
    const su2double* coor = coor_dst.data() + nodeID * nDim;

    /*--- Find nearest covered node for extrapolation ---*/
    su2double dist2Min = numeric_limits<su2double>::max();
    long nearestNodeID = -1;
    vector<su2double> grad(nVar * nDim, 0.0);

    /*--- Step 1: Search over neighbors and check if any of them are covered ---*/
    bool foundCoveredNeighbor = false;
    const unsigned short nElem_node = geometry_dst->nodes->GetnElem(nodeID);

    for (auto j = 0u; j < nElem_node; ++j) {
      unsigned long elemID = geometry_dst->nodes->GetElem(nodeID, j);
      auto* elem = geometry_dst->elem[elemID];

      for (auto iNode = 0u; iNode < elem->GetnNodes(); ++iNode) {
        unsigned long candidateNodeID = elem->GetNode(iNode);

        /*--- Skip if this candidate is also uncovered ---*/
        if (uncoveredNodes.find(candidateNodeID) != uncoveredNodes.end()) continue;

        /*--- Compute distance ---*/
        su2double dist2 = 0.0;
        for (auto k = 0u; k < nDim; ++k) {
          su2double ds = coor[k] - coor_dst[candidateNodeID * nDim + k];
          dist2 += ds * ds;
        }

        if (dist2 < dist2Min) {
          dist2Min = dist2;
          nearestNodeID = candidateNodeID;
          foundCoveredNeighbor = true;
        }
      }
    }

    /*--- Step 2: If no neighbor is covered, use R-tree search with processed set ---*/
    if (!foundCoveredNeighbor) {
      /*--- Compute initial search radius based on maximum edge length connected to this node ---*/
      su2double maxEdgeLength = 0.0;

      for (auto j = 0u; j < nElem_node; ++j) {
        unsigned long elemID = geometry_dst->nodes->GetElem(nodeID, j);
        auto* elem = geometry_dst->elem[elemID];

        /*--- Find maximum distance to any neighbor node in this element ---*/
        for (auto iNode = 0u; iNode < elem->GetnNodes(); ++iNode) {
          unsigned long neighborNodeID = elem->GetNode(iNode);
          if (neighborNodeID == nodeID) continue; // Skip self

          /*--- Compute edge length ---*/
          su2double edgeLength2 = 0.0;
          for (auto k = 0u; k < nDim; ++k) {
            su2double ds = coor[k] - coor_dst[neighborNodeID * nDim + k];
            edgeLength2 += ds * ds;
          }
          maxEdgeLength = max(maxEdgeLength, sqrt(edgeLength2));
        }
      }

      /*--- Use maximum edge length as initial search radius (fallback to small value if no edges) ---*/
      su2double searchRadius = (maxEdgeLength > 0.0) ? maxEdgeLength : 1e-3;
      const su2double radiusGrowthFactor = 2.0;

      /*--- Set to track processed nodes to avoid redundant checks ---*/
      set<unsigned long> processed;

      /*--- R-tree search with expanding radius ---*/
      CRTreeSearchBase& dstRTree = GetDestRTree();

      while (nearestNodeID == -1) {
        /*--- Create expanding bounding box around the uncovered node ---*/
        su2double elemCoor[2 * nDim];  // For 2D: [xmin, ymin, xmax, ymax]
        for (auto k = 0u; k < nDim; ++k) {
          elemCoor[k] = coor[k] - searchRadius;         // Min coordinates
          elemCoor[nDim + k] = coor[k] + searchRadius;  // Max coordinates
        }

        /*--- Search for nodes in the expanding bounding box ---*/
        set<unsigned long> candidateNodes;
        dstRTree.SearchNodesInElement(elemCoor, candidateNodes);

        /*--- Find the nearest covered node among candidates ---*/
        for (auto candidateNodeID : candidateNodes) {
          /*--- Skip if already processed ---*/
          if (processed.find(candidateNodeID) != processed.end()) continue;
          processed.insert(candidateNodeID);

          /*--- Skip if this candidate is also uncovered ---*/
          if (uncoveredNodes.find(candidateNodeID) != uncoveredNodes.end()) continue;

          /*--- Compute actual distance ---*/
          su2double dist2 = 0.0;
          for (auto k = 0u; k < nDim; ++k) {
            su2double ds = coor[k] - coor_dst[candidateNodeID * nDim + k];
            dist2 += ds * ds;
          }

          if (dist2 < dist2Min) {
            dist2Min = dist2;
            nearestNodeID = candidateNodeID;
          }
        }

        /*--- Expand search radius for next iteration ---*/
        searchRadius *= radiusGrowthFactor;
      }
    }

    /*--- Perform extrapolation if we found a suitable reference node ---*/
    if (nearestNodeID > 0) {
      /*--- Compute area-weighted average gradient for elements containing nearest node ---*/
      fill(grad.begin(), grad.end(), 0.0);
      su2double totalVol = 0.0;

      const unsigned short nElem_nearest = geometry_dst->nodes->GetnElem(nearestNodeID);
      for (auto k = 0u; k < nElem_nearest; ++k) {
        unsigned long nearestElemID = geometry_dst->nodes->GetElem(nearestNodeID, k);

        /*--- Only include elements that have overlap ---*/
        if (overlapMeshes.find(nearestElemID) != overlapMeshes.end()) {
          const su2double elemVol = geometry_dst->elem[nearestElemID]->GetVolume();
          totalVol += elemVol;

          for (auto iVar = 0u; iVar < nVar; ++iVar) {
            for (auto iDim = 0u; iDim < nDim; ++iDim) {
              grad[iVar * nDim + iDim] += dstElemGrad[nearestElemID][iVar * nDim + iDim] * elemVol;
            }
          }
        }
      }

      /*--- Normalize by total area ---*/
      if (totalVol > 0.0) {
        for (auto iVar = 0u; iVar < nVar; ++iVar) {
          for (auto iDim = 0u; iDim < nDim; ++iDim) {
            grad[iVar * nDim + iDim] /= totalVol;
          }
        }
      }

      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        /*--- Get solution value at nearest node ---*/
        su2double nearestValue = solver_dst->GetNodes()->GetSolution(nearestNodeID, iVar);

        /*--- Extrapolate: u(p_i) = u(p_nearest) + gra(u_K) * (p_i - p_nearest) ---*/
        su2double extrapolatedValue = nearestValue;
        for (auto iDim = 0u; iDim < nDim; ++iDim) {
          su2double ds = coor[iDim] - coor_dst[nearestNodeID * nDim + iDim];
          extrapolatedValue += grad[iVar * nDim + iDim] * ds;
        }

        solver_dst->GetNodes()->SetSolution(nodeID, iVar, extrapolatedValue);
        solver_dst->GetNodes()->SetSolution_Old(nodeID, iVar, extrapolatedValue);
      }

      nVolumeExtrapolated++;
    }
  }

  if (rank == MASTER_NODE)
    cout << "Volume nodes extrapolated: " << nVolumeExtrapolated << endl;
}

void CConservativeVolumeInterpolator::SurfaceInterpolation(CGeometry* geometry_src,
                                                           CGeometry* geometry_dst,
                                                           CSolver* solver_src,
                                                           CSolver* solver_dst,
                                                           const vector<su2double>& coor_dst) {
  unsigned long nSurfaceInterpolated = 0;

  CRTreeSearchBase& surfRTree = GetSourceSurfaceRTree();

  if (!surfRTree.IsSurfaceTreeBuilt()) {
    if (rank == MASTER_NODE) {
      cout << "Warning: Source R-tree not built. Cannot perform surface interpolation." << endl;
    }
    return;
  }

  /*--- Storage for the R-tree search ---*/
  unsigned long nearestNode;
  std::vector<unsigned short> markerIDs;
  std::vector<unsigned long> elemIDs;

  for (auto l = 0u; l < nPoint_dst; ++l) {
    if (!geometry_dst->nodes->GetPhysicalBoundary(l)) continue;

    /*--- Search for the nearest surface node ---*/
    const su2double* coor = coor_dst.data() + l * nDim;
    if (surfRTree.SearchNearestSurfaceNode(coor, nearestNode, markerIDs, elemIDs)) {
      unsigned short markerID;
      unsigned long elemID;
      int rankID;
      su2double dist, distMin = numeric_limits<su2double>::max();
      su2double surfCoor[3], surfCoorMin[3];

      for (auto i = 0u; i < markerIDs.size(); ++i) {
        const unsigned short iMarker = markerIDs[i];
        const unsigned long iElem = elemIDs[i];
        NearestPointOnElement(geometry_src, iMarker, iElem, coor, surfCoor,
                              dist, nDim);

        if (dist < distMin) {
          markerID = iMarker;
          elemID = iElem;
          distMin = dist;
          for (auto k = 0u; k < nDim; ++k)
            surfCoorMin[k] = surfCoor[k];
        }
      }

      /*--- Get surface element information ---*/
      unsigned short nNode = geometry_src->bound[markerID][elemID]->GetnNodes();

      /*--- Use nearest surface element nodes for interpolation ---*/
      su2double weightsInterpol[4];
      if (geometry_src->GetnDim() == 3) {
        /*--- TODO: implement 3D interpolation weights ---*/
        SU2_MPI::Error("Not implemented for 3D yet.", CURRENT_FUNCTION);
      } else {
        /*--- For 2D case (LINE elements), use inverse distance weighting ---*/
        su2double totalWeight = 0.0;

        for (auto iNode = 0u; iNode < nNode; ++iNode) {
          unsigned long nodeID = geometry_src->bound[markerID][elemID]->GetNode(iNode);

          /*--- Compute distance from interpolation point to node ---*/
          su2double dist2 = 0.0;
          for (auto k = 0u; k < nDim; ++k) {
            su2double diff = coor[k] - geometry_src->nodes->GetCoord(nodeID, k);
            dist2 += diff * diff;
          }

          /*--- Inverse distance weighting (with small epsilon to avoid division by zero) ---*/
          weightsInterpol[iNode] = 1.0 / (sqrt(dist2) + 1e-12);
          totalWeight += weightsInterpol[iNode];
        }

        /*--- Normalize weights ---*/
        for (auto iNode = 0u; iNode < nNode; ++iNode) {
          weightsInterpol[iNode] /= totalWeight;
        }
      }

      /*--- Initialize interpolated solution to zero ---*/
      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        solver_dst->GetNodes()->SetSolution(l, iVar, 0.0);
      }

      /*--- Interpolate using shape function weights ---*/
      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        unsigned long nodeID = geometry_src->bound[markerID][elemID]->GetNode(iNode);

        for (auto iVar = 0u; iVar < nVar; ++iVar) {
          su2double val = solver_src->GetNodes()->GetSolution(nodeID, iVar);
          solver_dst->GetNodes()->Add_DeltaSolution(l, iVar, weightsInterpol[iNode] * val);
        }
      }

      nSurfaceInterpolated++;
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Boundary nodes interpolated: " << nSurfaceInterpolated << endl;
  }
}
