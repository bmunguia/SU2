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
#include <iomanip>

#include "../include/CConservativeVolumeInterpolator.hpp"
#include "../../Common/include/fem/fem_standard_element.hpp"

CConservativeVolumeInterpolator::CConservativeVolumeInterpolator(SU2_Comm MPICommunicator)
    : CVolumeInterpolator(MPICommunicator) {}

void CConservativeVolumeInterpolator::Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                  CSolver** solver_container_src, CSolver** solver_container_dst) {
  if (rank == MASTER_NODE) {
    cout << endl << "----------------------------- Interpolation -----------------------------" << endl;
    cout << "Conservative solution interpolation from source mesh to destination mesh." << endl;
    cout << "Source mesh: " << geometry_src->GetGlobal_nPointDomain() << " points, ";
    cout << geometry_src->GetGlobal_nElemDomain() << " elements." << endl;
    cout << "Destination mesh: " << geometry_dst->GetGlobal_nPointDomain() << " points,";
    cout << geometry_dst->GetGlobal_nElemDomain() << " elements." << endl;
  }

  /*--- Build the ADTs ---*/
  InitializeADTs(config, geometry_src, geometry_dst);

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
  /*--- Step 0: Initialize destination coordinate vector                   ---*/
  /*--------------------------------------------------------------------------*/
  nVar = solver_src->GetnVar();
  vector<su2double> coorDst;
  InitializeCoords(geometry_dst, coorDst);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Apply the curvature correction to the destination nodes    ---*/
  /*--------------------------------------------------------------------------*/
  // if (rank == MASTER_NODE) cout << "Applying curvature correction." << endl;
  // ApplyCurvatureCorrection(config, geometry_src, geometry_dst, nDim, coorDst);

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Localize destination nodes on the source mesh              ---*/
  /*---         containingElems is a map from destination mesh nodes to    ---*/
  /*---         containing elements on the source mesh, and pointsFailed   ---*/
  /*---         is all the nodes for which no containing element was found ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Performing containment search." << endl;
  vector<optional<unsigned long>> containingElems;
  vector<int> containingElemRanks;
  vector<unsigned long> pointsFailed;
  PointLocalization(geometry_src, coorDst, containingElems, containingElemRanks,
                    pointsFailed);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute solution mass and gradient on source mesh          ---*/
  /*--------------------------------------------------------------------------*/
  vector<vector<su2double>> srcElemMass;
  vector<vector<su2double>> srcElemGrad;
  ComputeSourceSolutionMass(geometry_src, solver_src, srcElemMass, srcElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Compute the intersection of elements K_dst with elements   ---*/
  /*---         K_src it overlaps                                          ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing element intersections." << endl;
  IntersectionMesh overlappingElements;
  ComputeOverlappingElements(geometry_src, geometry_dst, coorDst, containingElems, overlappingElements);

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Compute destination mesh mass and gradient using Gauss     ---*/
  /*---         quadrature over intersection regions                       ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing destination mesh mass and gradients." << endl;
  vector<vector<su2double>> dstElemMass;
  vector<vector<su2double>> dstElemGrad;
  ComputeDestinationMassAndGradient(geometry_src, geometry_dst, solver_src, overlappingElements,
                                    srcElemMass, srcElemGrad, dstElemMass, dstElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 6: Correct the gradient to enforce the maximum principle      ---*/
  /*--------------------------------------------------------------------------*/
  // if (rank == MASTER_NODE) cout << "Applying local maximum principle correction." << endl;
  // ApplyMaximumPrincipleCorrection(geometry_src, geometry_dst, solver_src, coorDst, overlappingElements,
  //                                 srcElemMass, srcElemGrad, dstElemMass, dstElemGrad);

  /*--------------------------------------------------------------------------*/
  /*--- Step 7: Perform averaging to get solution at vertices.             ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Distributing solution to destination nodes." << endl;
  DistributeSolutionToNodes(geometry_dst, solver_dst, coorDst, dstElemMass, dstElemGrad);
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
      /*--- Find nearest surface element ---*/
      unsigned short markerID;
      su2double dist;
      su2double surfCoor[3];

      /*--- Find the closest point on the source surface mesh ---*/
      surfaceADT.DetermineNearestElement(coor, dist, markerID, elemID, rankID);

      /*--- Get surface element information ---*/
      unsigned short nNodes = geometry_src->bound[markerID][elemID]->GetnNodes();

      /*--- Get the corresponding volume element which contains the surface element ---*/
      set<unsigned long> candidateVolElems;
      bool firstNode = true;

      for (auto iNode = 0u; iNode < nNodes; ++iNode) {
        unsigned long nodeID = geometry_src->bound[markerID][elemID]->GetNode(iNode);
        unsigned short nElem_node = geometry_src->nodes->GetnElem(nodeID);

        if (firstNode) {
          /*--- For first node, add all connected volume elements ---*/
          for (auto iElem = 0u; iElem < nElem_node; ++iElem) {
            unsigned long volElemID = geometry_src->nodes->GetElem(nodeID, iElem);
            candidateVolElems.insert(volElemID);
          }
          firstNode = false;
        } else {
          /*--- For subsequent nodes, keep only elements that contain this node ---*/
          set<unsigned long> nodeElems;
          for (auto iElem = 0u; iElem < nElem_node; ++iElem) {
            unsigned long volElemID = geometry_src->nodes->GetElem(nodeID, iElem);
            nodeElems.insert(volElemID);
          }

          /*--- Intersect with previous candidates ---*/
          set<unsigned long> intersection;
          set_intersection(candidateVolElems.begin(), candidateVolElems.end(),
                          nodeElems.begin(), nodeElems.end(),
                          inserter(intersection, intersection.begin()));
          candidateVolElems = intersection;
        }
      }

      /*--- Use the first candidate volume element (should be unique) ---*/
      if (!candidateVolElems.empty()) {
        elemID = *candidateVolElems.begin();
        containingElems[l] = elemID;
        containingElemRanks[l] = rankID;
      } else {
        /*--- No volume element contains all surface nodes - fallback failed ---*/
        pointsFailed.push_back(l);
      }
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

  for (auto elemID = 0u; elemID < geometry->GetnElem(); ++elemID) {
    auto* elem = geometry->elem[elemID];
    const unsigned short nNodes = elem->GetnNodes();
    const unsigned short VTK_Type = elem->GetVTK_Type();
    const su2double elemVolume = elem->GetVolume();

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
    su2double vertexCoords[6];
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = elem->GetNode(iNode);
      vertexCoords[iNode * 2 + 0] = geometry->nodes->GetCoord(nodeID, 0);
      vertexCoords[iNode * 2 + 1] = geometry->nodes->GetCoord(nodeID, 1);
    }

    /*--- Compute mass and gradient ---*/
    ComputeTriangleMassAndGradient(vertexCoords, vertexSol, elemVolume, nVar,
                                   elemMass[elemID], elemGrad[elemID]);
  }
}

void CConservativeVolumeInterpolator::ComputeOverlappingElements(CGeometry* geometry_src,
                                                                 CGeometry* geometry_dst,
                                                                 const vector<su2double> &coor_dst,
                                                                 const vector<optional<unsigned long>>& containingElems,
                                                                 IntersectionMesh& overlappingElements) {
  overlappingElements.clear();

  /*--- Storage for triangle vertex coordinates ---*/
  su2double srcTri[6], dstTri[6];

  /*--- Only handle triangular elements for now ---*/
  if (nDim != 2) {
    if (rank == MASTER_NODE) {
      cout << "Warning: Element intersection currently only implemented for 2D. ";
      cout << "Skipping intersection computation for mesh." << endl;
    }
    return;
  }

  /*--- Loop over all destination elements ---*/
  unsigned long totalOverlaps = 0;
  unsigned long failedIntersections = 0;
  for (auto dstElemID = 0u; dstElemID < geometry_dst->GetnElem(); ++dstElemID) {
    auto* dstElem = geometry_dst->elem[dstElemID];

    /*--- Skip non-triangular elements ---*/
    if (dstElem->GetVTK_Type() != TRIANGLE) continue;

    /*--------------------------------------------------------------------------*/
    /*--- Step 1: Build initial list from elements K_src containing vertices ---*/
    /*---         of K_dst                                                   ---*/
    /*--------------------------------------------------------------------------*/
    set<unsigned long> candidateList;

    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      /*--- Find corresponding point index in corrected coordinates ---*/
      if (containingElems[nodeID].has_value()) {
        candidateList.insert(containingElems[nodeID].value());
      }
    }

    /*--- Skip this destination element if no initial candidates found ---*/
    if (candidateList.empty()) continue;

    /*--- Get destination triangle vertices ---*/
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      dstTri[iNode * 2 + 0] = coor_dst[nodeID * nDim + 0];
      dstTri[iNode * 2 + 1] = coor_dst[nodeID * nDim + 1];
    }

    /*--------------------------------------------------------------------------*/
    /*--- Step 2: Process candidate list, adding new elements during         ---*/
    /*---         intersection, and meshing the convex polygon formed by the ---*/
    /*---         intersection                                               ---*/
    /*--------------------------------------------------------------------------*/
    set<unsigned long> processedElems;
    queue<unsigned long> toProcess;

    /*--- Initialize queue with initial candidate list ---*/
    for (auto srcElemID : candidateList) {
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
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        unsigned long nodeID = srcElem->GetNode(iNode);
        srcTri[iNode * 2 + 0] = geometry_src->nodes->GetCoord(nodeID, 0);
        srcTri[iNode * 2 + 1] = geometry_src->nodes->GetCoord(nodeID, 1);
      }

      /*--- Check for intersection and detect new candidates ---*/
      vector<su2double> intersectionPoints;
      vector<su2double> meshedIntersection;
      set<unsigned long> newCandidates;

      if (TriangleTriangleIntersection(geometry_src, srcTri, dstTri, srcElemID, intersectionPoints,
                                       meshedIntersection, newCandidates)) {
        totalOverlaps++;

        overlappingElements[dstElemID].push_back(make_pair(srcElemID, meshedIntersection));
      } else {
        /*--- Triangle intersection failed ---*/
        failedIntersections++;
      }

      /*--- Add new candidates to processing queue ---*/
      for (auto newElemID : newCandidates) {
        if (processedElems.count(newElemID) == 0) {
          toProcess.push(newElemID);
        }
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Found " << totalOverlaps << " total overlapping element pairs." << endl;
    cout << "Number of destination elements with overlaps: " << overlappingElements.size() << endl;
    cout << "Number of failed triangle intersections: " << failedIntersections << "." << endl;
  }
}

bool CConservativeVolumeInterpolator::TriangleTriangleIntersection(CGeometry* geometry_src,
                                                                   const su2double dstTri[6],
                                                                   const su2double srcTri[6],
                                                                   unsigned long srcElemID,
                                                                   vector<su2double>& intersectionPoints,
                                                                   vector<su2double>& meshedIntersection,
                                                                   set<unsigned long>& newCandidates) {
  intersectionPoints.clear();
  meshedIntersection.clear();
  newCandidates.clear();

  const su2double EPS = 1e-12;

  /*--- Triangle KP (destination) vertices ---*/
  su2double P[3][2] = {{dstTri[0], dstTri[1]}, {dstTri[2], dstTri[3]}, {dstTri[4], dstTri[5]}};

  /*--- Triangle KQ (source) vertices ---*/
  su2double Q[3][2] = {{srcTri[0], srcTri[1]}, {srcTri[2], srcTri[3]}, {srcTri[4], srcTri[5]}};

  /*--- Edge definitions: edge j connects vertex j to vertex (j+1)%3 ---*/
  /*--- Edge 0: Q1-Q2, Edge 1: Q2-Q0, Edge 2: Q0-Q1 for triangle Q   ---*/
  /*--- Edge 0: P1-P2, Edge 1: P2-P0, Edge 2: P0-P1 for triangle P   ---*/

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Compute all 18 vertex-edge powers handle degenerate        ---*/
  /*---         intersection cases (1, 2, or 4 signed distances being 0)   ---*/
  /*--------------------------------------------------------------------------*/
  vector<su2double> cloudPoints;

  /*--- Compute powers: power_P[i][j] = power of vertex Pi w.r.t. edge j of triangle Q ---*/
  /*---                 power_Q[i][j] = power of vertex Qi w.r.t. edge j of triangle P ---*/
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

  /*--- Step 1: Handle degenerate edge-edge intersection cases ---*/
  bool isDegenerateEdgePair[3][3] = {};
  ProcessDegenerateEdgeIntersections(P_edges, Q_edges, power_P, power_Q, EPS, cloudPoints,
                                     isDegenerateEdgePair);

  /*--- Check if KP vertices are strictly inside KQ ---*/
  for (auto i = 0u; i < 3; ++i) {
    if ((power_P[i][0] > EPS) &&
        (power_P[i][1] > EPS) &&
        (power_P[i][2] > EPS)) {
      cloudPoints.push_back(P[i][0]);
      cloudPoints.push_back(P[i][1]);
    }
  }

  /*--- Check if KQ vertices are strictly inside KP ---*/
  for (auto i = 0u; i < 3; ++i) {
    if ((power_Q[i][0] > EPS) &&
        (power_Q[i][1] > EPS) &&
        (power_Q[i][2] > EPS)) {
      cloudPoints.push_back(Q[i][0]);
      cloudPoints.push_back(Q[i][1]);
      AddVertexBallToCandidates(geometry_src, srcElemID, i, newCandidates);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Process edge-edge intersections, ignoring any degenerate   ---*/
  /*---         cases since they were handled in Step 1                    ---*/
  /*--------------------------------------------------------------------------*/
  for (auto iP = 0u; iP < 3; ++iP) {
    for (auto jQ = 0u; jQ < 3; ++jQ) {
      /*--- Skip degenerate edge pairs (already handled in Step 1) ---*/
      if (isDegenerateEdgePair[iP][jQ]) continue;

      su2double* edgeP0 = P_edges[iP][0];
      su2double* edgeP1 = P_edges[iP][1];
      su2double* edgeQ0 = Q_edges[jQ][0];
      su2double* edgeQ1 = Q_edges[jQ][1];

      /*--- Compute signed distances for this edge pair ---*/
      su2double distP0_eQ = ComputeSignedDistance(edgeP0, edgeQ0, edgeQ1);
      su2double distP1_eQ = ComputeSignedDistance(edgeP1, edgeQ0, edgeQ1);
      su2double distQ0_eP = ComputeSignedDistance(edgeQ0, edgeP0, edgeP1);
      su2double distQ1_eP = ComputeSignedDistance(edgeQ1, edgeP0, edgeP1);

      /*--- Check for edge-edge intersection ---*/
      if ((distP0_eQ * distP1_eQ < 0) && (distQ0_eP * distQ1_eP < 0)) {
        /*--- Compute intersection point ---*/
        su2double t = distP0_eQ / (distP0_eQ - distP1_eQ);
        su2double intersectionX = edgeP0[0] + t * (edgeP1[0] - edgeP0[0]);
        su2double intersectionY = edgeP0[1] + t * (edgeP1[1] - edgeP0[1]);

        cloudPoints.push_back(intersectionX);
        cloudPoints.push_back(intersectionY);

        /*--- Add face neighbor to candidates ---*/
        AddFaceNeighborToCandidates(geometry_src, srcElemID, jQ, newCandidates);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Mesh the intersection polygon                              ---*/
  /*--------------------------------------------------------------------------*/
  MeshConvexPolygon(cloudPoints, meshedIntersection);

  return (!meshedIntersection.empty());
}

void CConservativeVolumeInterpolator::ProcessDegenerateEdgeIntersections(su2double* P_edges[3][2], su2double* Q_edges[3][2],
                                                                         const su2double power_P[3][3], const su2double power_Q[3][3],
                                                                         const su2double EPS,
                                                                         vector<su2double>& intersectionPoints,
                                                                         bool isDegenerateEdgePair[3][3]) {
  /*--- Initialize all vertices as non-degenerate ---*/
  bool isDegeneratePi[3] = {};
  bool isDegenerateQj[3] = {};

  auto addPoint = [&](su2double x, su2double y, int Pi, int Qj) {
    bool alreadyAdded = false;
    if (Pi >= 0) {
      if (isDegeneratePi[Pi]) alreadyAdded = true;
      isDegeneratePi[Pi] = true;
    }
    if (Qj >= 0) {
      if (isDegenerateQj[Qj]) alreadyAdded = true;
      isDegenerateQj[Qj] = true;
    }
    if (alreadyAdded) return;
    intersectionPoints.push_back(x);
    intersectionPoints.push_back(y);
  };

  /*--- Process all 9 edge-edge combinations for degenerate cases ---*/
  for (auto iP = 0u; iP < 3; ++iP) {
    for (auto jQ = 0u; jQ < 3; ++jQ) {
      int iP0 = iP, iP1 = (iP + 1) % 3;
      int jQ0 = jQ, jQ1 = (jQ + 1) % 3;
      su2double* edgeP0 = P_edges[iP][0];  // P vertex iP
      su2double* edgeP1 = P_edges[iP][1];  // P vertex (iP+1)%3
      su2double* edgeQ0 = Q_edges[jQ][0];  // Q vertex jQ
      su2double* edgeQ1 = Q_edges[jQ][1];  // Q vertex (jQ+1)%3

      /*--- Use precomputed powers instead of recalculating ---*/
      /*--- Edge iP of triangle P goes from vertex iP to vertex (iP+1)%3 ---*/
      /*--- Edge jQ of triangle Q goes from vertex jQ to vertex (jQ+1)%3 ---*/
      su2double distP0_eQ = power_P[iP0][jQ];  // Power of P vertex iP w.r.t. Q edge jQ
      su2double distP1_eQ = power_P[iP1][jQ];  // Power of P vertex (iP+1)%3 w.r.t. Q edge jQ
      su2double distQ0_eP = power_Q[jQ0][iP];  // Power of Q vertex jQ w.r.t. P edge iP
      su2double distQ1_eP = power_Q[jQ1][iP];  // Power of Q vertex (jQ+1)%3 w.r.t. P edge iP

      /*--- Count how many powers are zero ---*/
      int zeroCount = 0;
      bool P0_zero = (abs(distP0_eQ) < EPS);
      bool P1_zero = (abs(distP1_eQ) < EPS);
      bool Q0_zero = (abs(distQ0_eP) < EPS);
      bool Q1_zero = (abs(distQ1_eP) < EPS);

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
          if (distQ0_eP * distQ1_eP < 0) {
            addPoint(edgeP0[0], edgeP0[1], iP0, -1);
          }
        } else if (P1_zero) {
          /*--- P1 lies on edge Q, intersection if Q0 and Q1 on opposite sides of edge P ---*/
          if (distQ0_eP * distQ1_eP < 0) {
            addPoint(edgeP1[0], edgeP1[1], iP1, -1);
          }
        } else if (Q0_zero) {
          /*--- Q0 lies on edge P, intersection if P0 and P1 on opposite sides of edge Q ---*/
          if (distP0_eQ * distP1_eQ < 0) {
            addPoint(edgeQ0[0], edgeQ0[1], -1, jQ0);
          }
        } else if (Q1_zero) {
          /*--- Q1 lies on edge P, intersection if P0 and P1 on opposite sides of edge Q ---*/
          if (distP0_eQ * distP1_eQ < 0) {
            addPoint(edgeQ1[0], edgeQ1[1], -1, jQ1);
          }
        }
      }

      /*--------------------------------------------------------------------------*/
      /*--- Case 2: Two powers are zero, one for each edge                     ---*/
      /*--------------------------------------------------------------------------*/
      else if (zeroCount == 2) {
        if (P0_zero && Q0_zero) {
          /*--- P0 = Q0, common vertex ---*/
          addPoint(edgeP0[0], edgeP0[1], iP0, jQ0);
        } else if (P0_zero && Q1_zero) {
          /*--- P0 = Q1, common vertex ---*/
          addPoint(edgeP0[0], edgeP0[1], iP0, jQ1);
        } else if (P1_zero && Q0_zero) {
          /*--- P1 = Q0, common vertex ---*/
          addPoint(edgeP1[0], edgeP1[1], iP1, jQ0);
        } else if (P1_zero && Q1_zero) {
          /*--- P1 = Q1, common vertex ---*/
          addPoint(edgeP1[0], edgeP1[1], iP1, jQ1);
        }
      }

      /*--------------------------------------------------------------------------*/
      /*--- Case 3: All four powers are zero (edges are collinear)             ---*/
      /*--------------------------------------------------------------------------*/
      else if (zeroCount == 4) {
        /*--- Edges are collinear, check for overlap using parametrization ---*/
        /*--- Parametrize both edges and find overlap in parameter space ---*/

        /*--- Edge P: P_start + t * (P_end - P_start), t ∈ [0,1] ---*/
        /*--- Edge Q: Q_start + s * (Q_end - Q_start), s ∈ [0,1] ---*/

        /*--- Direction vector of edge P ---*/
        su2double dx_P = edgeP1[0] - edgeP0[0];
        su2double dy_P = edgeP1[1] - edgeP0[1];
        su2double lengthP_sq = dx_P*dx_P + dy_P*dy_P;

        if (lengthP_sq > EPS*EPS) {
          /*--- Find parameter values where Q vertices lie on P edge ---*/
          /*--- Q_start = P_start + t_Q0 * (P_end - P_start) ---*/
          /*--- Q_end   = P_start + t_Q1 * (P_end - P_start) ---*/

          su2double t_Q0 = ((edgeQ0[0] - edgeP0[0]) * dx_P +
                            (edgeQ0[1] - edgeP0[1]) * dy_P) / lengthP_sq;
          su2double t_Q1 = ((edgeQ1[0] - edgeP0[0]) * dx_P +
                            (edgeQ1[1] - edgeP0[1]) * dy_P) / lengthP_sq;

          /*--- Edge P spans parameter interval [0, 1] ---*/
          /*--- Edge Q spans parameter interval [min(t_Q0, t_Q1), max(t_Q0, t_Q1)] ---*/
          su2double t_min = min(t_Q0, t_Q1);
          su2double t_max = max(t_Q0, t_Q1);

          /*--- Find overlap of intervals [0, 1] and [t_min, t_max] ---*/
          su2double overlap_start = max(0.0, t_min);
          su2double overlap_end = min(1.0, t_max);

          if (overlap_start < overlap_end) {
            /*--- There is overlap, add intersection points ---*/
            /*--- Add start point of overlap ---*/
            su2double startX = edgeP0[0] + overlap_start * dx_P;
            su2double startY = edgeP0[1] + overlap_start * dy_P;

            /*--- Check if start point corresponds to a known vertex ---*/
            int startPi = -1, startQj = -1;
            if (abs(overlap_start - 0.0) < EPS) {
              startPi = iP0;  // Start point is P0
            } else if (abs(overlap_start - 1.0) < EPS) {
              startPi = iP1;  // Start point is P1
            }
            if (abs(t_Q0 - overlap_start) < EPS) {
              startQj = jQ0;  // Start point is Q0
            } else if (abs(t_Q1 - overlap_start) < EPS) {
              startQj = jQ1;  // Start point is Q1
            }
            addPoint(startX, startY, startPi, startQj);

            /*--- Add end point of overlap if different from start ---*/
            if (abs(overlap_end - overlap_start) > EPS) {
              su2double endX = edgeP0[0] + overlap_end * dx_P;
              su2double endY = edgeP0[1] + overlap_end * dy_P;

              /*--- Check if end point corresponds to a known vertex ---*/
              int endPi = -1, endQj = -1;
              if (abs(overlap_end - 0.0) < EPS) {
                endPi = iP0;  // End point is P0
              } else if (abs(overlap_end - 1.0) < EPS) {
                endPi = iP1;  // End point is P1
              }
              if (abs(t_Q0 - overlap_end) < EPS) {
                endQj = jQ0;  // End point is Q0
              } else if (abs(t_Q1 - overlap_end) < EPS) {
                endQj = jQ1;  // End point is Q1
              }
              addPoint(endX, endY, endPi, endQj);
            }
          }
        }
      }
    }
  }
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
                                                                const unsigned short localNodeID,
                                                                set<unsigned long>& newCandidates) {
  /*--- Get the vertex node ID ---*/
  auto* srcElem = geometry_src->elem[srcElemID];
  unsigned long nodeID = srcElem->GetNode(localNodeID);

  /*--- Add all elements other than the current source element that contain this vertex ---*/
  for (auto jElem = 0u; jElem < geometry_src->nodes->GetnElem(nodeID); ++jElem) {
    unsigned long elemID = geometry_src->nodes->GetElem(nodeID, jElem);
    if (elemID == srcElemID) continue;  // Skip self

    auto* neighborElem = geometry_src->elem[elemID];
    if (neighborElem->GetVTK_Type() != TRIANGLE) continue;

    newCandidates.insert(elemID);
  }
}

void CConservativeVolumeInterpolator::MeshConvexPolygon(const vector<su2double>& polygonPoints,
                                                        vector<su2double>& polygonMesh) {
  polygonMesh.clear();

  unsigned int numPolygonPoint = polygonPoints.size() / 2;

  /*--- Handle edge cases ---*/
  if (numPolygonPoint < 3) {
    return;  // Cannot mesh a polygon with less than 3 points
  }

  /*--- Lambda to add a triangle with positive area and return orientation ---*/
  auto addPositiveTriangle = [&](su2double P0[2], su2double P1[2], su2double P2[2]) -> bool {
    /*--- Check orientation ---*/
    su2double cross = (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P1[1] - P0[1]) * (P2[0] - P0[0]);
    if (cross > 0) {
      /*--- Counter-clockwise orientation - add triangle as is ---*/
      polygonMesh.push_back(P0[0]); polygonMesh.push_back(P0[1]);
      polygonMesh.push_back(P1[0]); polygonMesh.push_back(P1[1]);
      polygonMesh.push_back(P2[0]); polygonMesh.push_back(P2[1]);
      return true;  // Original order was counter-clockwise
    } else {
      /*--- Clockwise orientation - reverse order ---*/
      polygonMesh.push_back(P0[0]); polygonMesh.push_back(P0[1]);
      polygonMesh.push_back(P2[0]); polygonMesh.push_back(P2[1]);
      polygonMesh.push_back(P1[0]); polygonMesh.push_back(P1[1]);
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
    int edgeToConnect = -1;
    for (auto j = 0u; j < boundaryEdges.size(); ++j) {
      int idx1 = boundaryEdges[j].first;
      int idx2 = boundaryEdges[j].second;

      su2double testP1[2] = {polygonPoints[idx1 * 2 + 0], polygonPoints[idx1 * 2 + 1]};
      su2double testP2[2] = {polygonPoints[idx2 * 2 + 0], polygonPoints[idx2 * 2 + 1]};

      if (ComputeSignedDistance(P, testP1, testP2) < 0) {
        edgeToConnect = j;
        break;
      }
    }

    /*--- Should always find exactly one boundary edge for a convex polygon ---*/
    if (edgeToConnect >= 0) {
      int idx1 = boundaryEdges[edgeToConnect].first;
      int idx2 = boundaryEdges[edgeToConnect].second;

      su2double testP1[2] = {polygonPoints[idx1 * 2 + 0], polygonPoints[idx1 * 2 + 1]};
      su2double testP2[2] = {polygonPoints[idx2 * 2 + 0], polygonPoints[idx2 * 2 + 1]};

      /*--- Create triangle with the new point and the boundary edge ---*/
      isCounterClockwise = addPositiveTriangle(P, testP1, testP2);

      /*--- Update boundary: remove the used edge and add two new boundary edges ---*/
      /*--- Maintain consistent orientation based on the triangle orientation ---*/
      boundaryEdges.erase(boundaryEdges.begin() + edgeToConnect);

      if (isCounterClockwise) {
        /*--- Triangle P->edgeStart->edgeEnd is counter-clockwise ---*/
        boundaryEdges.push_back({i, idx1});  // Edge from new point to start of used edge
        boundaryEdges.push_back({idx2, i});  // Edge from end of used edge to new point
      } else {
        /*--- Triangle P->edgeEnd->edgeStart was stored (reversed) ---*/
        boundaryEdges.push_back({i, idx2});  // Edge from new point to end of used edge
        boundaryEdges.push_back({idx1, i});  // Edge from start of used edge to new point
      }
    }
  }
}

void CConservativeVolumeInterpolator::ComputeDestinationMassAndGradient(CGeometry* geometry_src,
                                                                        CGeometry* geometry_dst,
                                                                        CSolver* solver_src,
                                                                        const IntersectionMesh& overlappingElements,
                                                                        const vector<vector<su2double>>& srcElemMass,
                                                                        const vector<vector<su2double>>& srcElemGrad,
                                                                        vector<vector<su2double>>& dstElemMass,
                                                                        vector<vector<su2double>>& dstElemGrad) {
  dstElemMass.resize(nElem_dst, vector<su2double>(nVar, 0.0));
  dstElemGrad.resize(nElem_dst, vector<su2double>(nVar * nDim, 0.0));

  /*--- Loop over all destination elements that have intersections ---*/
  su2double absDiffTol[5] = {1e-10, 1e-8, 1e-6, 1e-4, 1e-2};
  vector<unsigned long> countAbsDiff(5, 0);
  su2double relDiffTol[5] = {2e-2, 5e-2, 1e-1, 2e-1, 5e-1};
  vector<unsigned long> countRelDiff(5, 0);

  for (const auto& elemPair : overlappingElements) {
    unsigned long dstElemID = elemPair.first;
    const auto& srcElemMeshPairs = elemPair.second;

    /*--- Initialize destination element mass and gradient ---*/
    const auto* dstElem = geometry_dst->elem[dstElemID];
    auto& dstMass = dstElemMass[dstElemID];
    auto& dstGrad = dstElemGrad[dstElemID];
    fill(dstMass.begin(), dstMass.end(), 0.0);
    fill(dstGrad.begin(), dstGrad.end(), 0.0);

    /*--- Track total triangle area for this destination element ---*/
    su2double totalTriangleArea = 0.0;

    /*--- Process each intersection region T_j = intersection(K_dst, K_src_j) ---*/
    for (const auto& srcMeshPair : srcElemMeshPairs) {
      unsigned long srcElemID = srcMeshPair.first;
      const vector<su2double>& intersectionMesh = srcMeshPair.second;

      /*--- Gauss quadrature over all triangles in the intersection mesh ---*/
      unsigned int numTriangles = intersectionMesh.size() / 6;  // 6 coordinates per triangle

      for (auto iTri = 0u; iTri < numTriangles; ++iTri) {
        /*--- Get triangle vertices ---*/
        const su2double* coor_tri = intersectionMesh.data() + iTri * 6;
        const su2double x0 = coor_tri[0], y0 = coor_tri[1];
        const su2double x1 = coor_tri[2], y1 = coor_tri[3];
        const su2double x2 = coor_tri[4], y2 = coor_tri[5];

        /*--- Triangle area using cross product ---*/
        const su2double cross = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        const su2double area = 0.5 * abs(cross);

        /*--- Add to total triangle area ---*/
        totalTriangleArea += area;

        /*--- Get source element properties ---*/
        const auto* srcElem = geometry_src->elem[srcElemID];
        auto srcMass = srcElemMass[srcElemID];
        auto srcGrad = srcElemGrad[srcElemID];

        const su2double srcVolume = srcElem->GetVolume();
        const su2double* G_K_src = srcElem->GetCG();

        /*--- Use 1-point Gauss quadrature (centroid rule) ---*/
        /*--- Triangle centroid coordinates ---*/
        const su2double xi = (x0 + x1 + x2) / 3.0;
        const su2double yi = (y0 + y1 + y2) / 3.0;

        /*--- Integrate mass and gradient using 1-point quadrature ---*/
        for (auto iVar = 0u; iVar < nVar; ++iVar) {
          const su2double u_src = srcMass[iVar] / srcVolume;
          const su2double* grad_u = srcGrad.data() + iVar * nDim;

          /*--- Displacement from source centroid to triangle centroid ---*/
          const su2double dx = xi - G_K_src[0];
          const su2double dy = yi - G_K_src[1];

          /*--- Evaluate solution at triangle centroid ---*/
          const su2double u_quad = u_src + grad_u[0] * dx + grad_u[1] * dy;

          /*--- Add contribution to destination mass and gradient ---*/
          dstMass[iVar] += area * u_quad;
          dstGrad[iVar * nDim + 0] += area * grad_u[0];
          dstGrad[iVar * nDim + 1] += area * grad_u[1];
        }

        /*--- Use 3-point Gauss quadrature for exact integration of linear functions ---*/
        /*--- Quadrature points in reference triangle (r,s) coordinates ---*/
        // const su2double r[3] = {0.5, 0.0, 0.5};
        // const su2double s[3] = {0.0, 0.5, 0.5};
        // const su2double w[3] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};  // Equal weights

        // /*--- Integrate mass and gradient using Gauss quadrature ---*/
        // for (auto iVar = 0u; iVar < nVar; ++iVar) {
        //   const su2double u_src = srcMass[iVar] / srcVolume;
        //   const su2double* grad_u = srcGrad.data() + iVar * nDim;

        //   su2double mass_integral = 0.0;
        //   su2double grad_integral[2] = {0.0, 0.0};

        //   /*--- Loop over quadrature points ---*/
        //   for (auto iGauss = 0u; iGauss < 3; ++iGauss) {
        //     /*--- Map reference coordinates to physical coordinates ---*/
        //     const su2double xi = (1.0 - r[iGauss] - s[iGauss]) * x0 + r[iGauss] * x1 + s[iGauss] * x2;
        //     const su2double yi = (1.0 - r[iGauss] - s[iGauss]) * y0 + r[iGauss] * y1 + s[iGauss] * y2;

        //     /*--- Displacement from source centroid to quadrature point ---*/
        //     const su2double dx = xi - G_K_src[0];
        //     const su2double dy = yi - G_K_src[1];

        //     /*--- Evaluate solution at quadrature point ---*/
        //     const su2double u_quad = u_src + grad_u[0] * dx + grad_u[1] * dy;

        //     /*--- Add weighted contribution to integrals ---*/
        //     mass_integral += w[iGauss] * u_quad;
        //     grad_integral[0] += w[iGauss] * grad_u[0];
        //     grad_integral[1] += w[iGauss] * grad_u[1];
        //   }

        //   /*--- Add contribution to destination mass and gradient ---*/
        //   dstMass[iVar] += area * mass_integral;
        //   dstGrad[iVar * nDim + 0] += area * grad_integral[0];
        //   dstGrad[iVar * nDim + 1] += area * grad_integral[1];
        // }
      }
    }

    /*--- Volume average integral: gra(u_dst) = int_K_dst (gra(u) dA) / |K_dst| ---*/
    const su2double dstVolume = dstElem->GetVolume();
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      for (auto iDim = 0u; iDim < nDim; ++iDim) {
        dstGrad[iVar * nDim + iDim] /= dstVolume;
      }
    }

    /*--- Compare total triangle area with destination element volume ---*/
    const su2double absDiff = abs(dstVolume - totalTriangleArea);
    const su2double relDiff = absDiff / dstVolume;
    for (auto i = 0u; i < 5; ++i) {
      if (absDiff > absDiffTol[i]) countAbsDiff[i]++;
      if (relDiff > relDiffTol[i]) countRelDiff[i]++;
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Area conservation check (absolute difference):" << endl;
    for (auto i = 0u; i < 5; ++i) {
      cout << "  Number exceeding " << scientific << setprecision(1);
      cout << absDiffTol[i] << ": ";
      cout << countAbsDiff[i] << endl;
    }
    cout << "Area conservation check (relative difference):" << endl;
    for (auto i = 0u; i < 5; ++i) {
      cout << "  Number exceeding " << fixed << setprecision(0) << setw(3);
      cout << relDiffTol[i] * 100 << "%: ";
      cout << countRelDiff[i] << endl;
    }
  }
}

void CConservativeVolumeInterpolator::ApplyMaximumPrincipleCorrection(CGeometry* geometry_src,
                                                                      CGeometry* geometry_dst,
                                                                      CSolver* solver_src,
                                                                      const vector<su2double>& coor_dst,
                                                                      const IntersectionMesh& overlappingElements,
                                                                      const vector<vector<su2double>>& srcElemMass,
                                                                      const vector<vector<su2double>>& srcElemGrad,
                                                                      vector<vector<su2double>>& dstElemMass,
                                                                      vector<vector<su2double>>& dstElemGrad) {
  const su2double EPS = 1e-12;

  /*--- Loop over all destination elements that have overlaps ---*/
  su2double u_tilde[3];
  vector<su2double> correctedMass(1, 0.0);
  vector<su2double> correctedGrad(1 * nDim, 0.0);
  vector<vector<su2double>> vertexSol(3, vector<su2double>(1));
  for (const auto& elemPair : overlappingElements) {
    unsigned long dstElemID = elemPair.first;
    const vector<pair<unsigned long, vector<su2double>>>& srcElemMeshPairs = elemPair.second;

    auto* dstElem = geometry_dst->elem[dstElemID];
    if (dstElem->GetVTK_Type() != TRIANGLE) continue;

    /*--- Get destination element vertices ---*/
    su2double dstVertices[6];
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      dstVertices[iNode * 2 + 0] = coor_dst[nodeID * nDim + 0];
      dstVertices[iNode * 2 + 1] = coor_dst[nodeID * nDim + 1];
    }

    /*--- Element centroid (barycenter G_K) ---*/
    const su2double* G_K = dstElem->GetCG();

    /*--- For each variable, apply Alauzet's maximum principle correction ---*/
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      /*--------------------------------------------------------------------------*/
      /*--- Step 1: Compute local bounds from overlapping source elements      ---*/
      /*--------------------------------------------------------------------------*/
      su2double u_min = 1e20;
      su2double u_max = -1e20;

      /*--- Find all vertices Q from source elements K_src that K overlaps ---*/
      for (const auto& srcMeshPair : srcElemMeshPairs) {
        unsigned long srcElemID = srcMeshPair.first;
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
      /*--- Step 2: Get current solution at destination element                ---*/
      /*--------------------------------------------------------------------------*/
      const su2double elemVolume = dstElem->GetVolume();
      const su2double u_G = dstElemMass[dstElemID][iVar] / elemVolume;
      const su2double* gradu_G = dstElemGrad[dstElemID].data() + iVar * nDim;

      /*--------------------------------------------------------------------------*/
      /*--- Step 3: Compute u_K(P_i) at each vertex using Taylor expansion     ---*/
      /*--------------------------------------------------------------------------*/
      su2double u_K_P[3];  // Values at vertices P_0, P_1, P_2
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        /*--- Vector G_K P_i ---*/
        const su2double dx = dstVertices[iNode * 2 + 0] - G_K[0];
        const su2double dy = dstVertices[iNode * 2 + 1] - G_K[1];

        /*--- u_K(P_i) = u_K(G_K) + gra(u_K) · G_K P_i ---*/
        u_K_P[iNode] = u_G + gradu_G[0] * dx + gradu_G[1] * dy;
      }

      /*--------------------------------------------------------------------------*/
      /*--- Step 4: Check if maximum principle is violated                     ---*/
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
      /*--- Step 5: Apply Alauzet's correction algorithm                       ---*/
      /*--------------------------------------------------------------------------*/
      /*--- Sort vertices by solution value: u_K(P_0) ≤ u_K(P_1) ≤ u_K(P_2) ---*/
      vector<pair<su2double, unsigned short>> sortedValues;
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        sortedValues.push_back(make_pair(u_K_P[iNode], iNode));
      }
      cout << "Unsorted: " << sortedValues[0].first << "," << sortedValues[0].second;
      cout << ", " << sortedValues[1].first << "," << sortedValues[1].second;
      cout << ", " << sortedValues[2].first << "," << sortedValues[2].second << endl;
      sort(sortedValues.begin(), sortedValues.end());
       cout << "Sorted: " << sortedValues[0].first << "," << sortedValues[0].second;
      cout << ", " << sortedValues[1].first << "," << sortedValues[1].second;
      cout << ", " << sortedValues[2].first << "," << sortedValues[2].second << endl;

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
      /*--- Step 6: Compute corrected mass and gradient from new nodal values  ---*/
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
      ComputeTriangleMassAndGradient(dstVertices, vertexSol, elemVolume, 1,
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
      const su2double elemVolume = elem->GetVolume();
      const su2double* G_K = elem->GetCG();

      /*--- Vector from centroid to node ---*/
      const su2double vec[2] = {coor[0] - G_K[0], coor[1] - G_K[1]};

      /*--- Loop over variables ---*/
      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        /*--- Get element-centered solution ---*/
        const su2double u_G = dstElemMass[elemID][iVar] / elemVolume;
        const su2double* gradu_G = dstElemGrad[elemID].data() + iVar * nDim;

        /*--- Linear reconstruction: u(P_i) = u(G_K) + gra(u) · (P_i - G_K) ---*/
        const su2double vertexValue = u_G + gradu_G[0] * vec[0] + gradu_G[1] * vec[1];

        /*--- Accumulate weighted contribution ---*/
        totalValue[iVar] += vertexValue * elemVolume;
      }

      /*--- Accumulate weight ---*/
      totalWeight += elemVolume;
    }

    /*--- Compute weighted average and set solution for this node ---*/
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      const su2double avgValue = totalValue[iVar] / totalWeight;
      solver_dst->GetNodes()->SetSolution(l, iVar, avgValue);
    }
  }
}

void CConservativeVolumeInterpolator::ComputeTriangleMassAndGradient(const su2double vertexCoords[6],
                                                                     const vector<vector<su2double>>& vertexSol,
                                                                     const su2double elemVolume,
                                                                     const unsigned short nVar,
                                                                     vector<su2double>& mass,
                                                                     vector<su2double>& grad) {
  /*--- Initialize output ---*/
  fill(mass.begin(), mass.end(), 0.0);
  fill(grad.begin(), grad.end(), 0.0);

  /*--- Extract triangle vertices ---*/
  const su2double x0 = vertexCoords[0], y0 = vertexCoords[1];
  const su2double x1 = vertexCoords[2], y1 = vertexCoords[3];
  const su2double x2 = vertexCoords[4], y2 = vertexCoords[5];

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

void CConservativeVolumeInterpolator::ComputeTriangleMassAndGradientFEM(const su2double vertexCoords[6],
                                                                        const vector<vector<su2double>>& vertexSol,
                                                                        const su2double elemVolume,
                                                                        const unsigned short nVar,
                                                                        vector<su2double>& mass,
                                                                        vector<su2double>& grad) {
  /*--- Initialize output ---*/
  fill(mass.begin(), mass.end(), 0.0);
  fill(grad.begin(), grad.end(), 0.0);

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

  /*--- Loop over integration points ---*/
  for (auto iInt = 0u; iInt < nInt; ++iInt) {
    /*--- Compute Jacobian of transformation ---*/
    su2double dxdr = 0.0, dydr = 0.0;
    su2double dxds = 0.0, dyds = 0.0;

    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned short ind = iInt * 3 + iNode;

      dxdr += vertexCoords[iNode * 2 + 0] * drLagBasis[ind];
      dydr += vertexCoords[iNode * 2 + 1] * drLagBasis[ind];
      dxds += vertexCoords[iNode * 2 + 0] * dsLagBasis[ind];
      dyds += vertexCoords[iNode * 2 + 1] * dsLagBasis[ind];
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
      grad[iVar * nDim + k] /= elemVolume;
    }
  }
}
