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
#include "../../Common/include/adt/CADTPointsOnlyClass.hpp"

CConservativeVolumeInterpolator::CConservativeVolumeInterpolator(SU2_Comm MPICommunicator)
    : CVolumeInterpolator(MPICommunicator) {}

void CConservativeVolumeInterpolator::Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                  CSolver** solver_container_src, CSolver** solver_container_dst,
                                                  bool initial_interp) {
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
    if (solver_src && solver_dst) {
      ConservativeInterpolation(config, geometry_src, geometry_dst, solver_src, solver_dst, initial_interp);
      initial_interp = false;
    }
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
                                                                CSolver* solver_dst,
                                                                bool initial_interp) {

  /*--------------------------------------------------------------------------*/
  /*--- TODO: - MPI parallel: need some sort of parallel tree search, or   ---*/
  /*---         geometric partitioning of the source mesh to align with    ---*/
  /*---         the destination mesh for overlap detection                 ---*/
  /*---       - Periodic boundaries: I currently communicate the mass      ---*/
  /*---         volume at periodic boundaries, but errors may accumulate   ---*/
  /*---         for curved boundaries, where incomplete (or nonexistent)   ---*/
  /*---         overlaps may occur; should probably use some sort of ghost ---*/
  /*---         elements                                                   ---*/
  /*--------------------------------------------------------------------------*/
  nVar = solver_src->GetnVar();

  if (initial_interp) {
    /*--------------------------------------------------------------------------*/
    /*--- Step 0: Initialize vectors.                                        ---*/
    /*--------------------------------------------------------------------------*/
    // InitializeCoords(geometry_dst, coorDst);
    srcElemMass.resize(nElem_src, vector<su2double>(nVar));
    srcElemGrad.resize(nElem_src, vector<su2double>(nVar * nDim));
    dstElemMass.resize(nElem_dst, vector<su2double>(nVar));
    dstElemGrad.resize(nElem_dst, vector<su2double>(nVar * nDim));

    /*--------------------------------------------------------------------------*/
    /*--- Step 1: Localize destination nodes on the source mesh.             ---*/
    /*---         containingElems is a map from destination mesh nodes to    ---*/
    /*---         containing elements on the source mesh, and                ---*/
    /*---         uncontainedNodes is all the nodes for which no containing  ---*/
    /*---         element was found                                          ---*/
    /*--------------------------------------------------------------------------*/
    if (rank == MASTER_NODE) cout << "Performing containment search." << endl;
    PointLocalization(geometry_src, geometry_dst);

    /*--------------------------------------------------------------------------*/
    /*--- Step 2: Compute the intersection of elements K_dst with elements   ---*/
    /*---         K_src it overlaps.                                         ---*/
    /*--------------------------------------------------------------------------*/
    if (rank == MASTER_NODE) cout << "Computing element intersections." << endl;
    CreateIntersectionMeshes(geometry_src, geometry_dst);

    /*--------------------------------------------------------------------------*/
    /*--- Step 3: Generate nearest contained nodes for uncontained           ---*/
    /*---           nodes for efficient extrapolation later.                 ---*/
    /*--------------------------------------------------------------------------*/
    if (!uncontainedNodes.empty()) {
      if (rank == MASTER_NODE) cout << "Finding nearest contained nodes for uncontained nodes." << endl;
      FindNearestContainedNodes(geometry_dst);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Compute solution mass and gradient on source mesh.         ---*/
  /*--------------------------------------------------------------------------*/
  ComputeSourceSolutionMass(geometry_src, solver_src);

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Compute destination mesh mass and gradient using Gauss     ---*/
  /*---         quadrature over intersection regions.                      ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Computing destination mesh mass and gradients." << endl;
  ComputeDestinationMassAndGradient(geometry_src, geometry_dst, solver_src);

  /*--------------------------------------------------------------------------*/
  /*--- Step 6: Correct the gradient to enforce the maximum principle.     ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Applying local maximum principle correction." << endl;
  ApplyMaximumPrincipleCorrection(geometry_src, geometry_dst, solver_src);

  /*--------------------------------------------------------------------------*/
  /*--- Step 7: Perform averaging to get solution at vertices.             ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Distributing solution to destination nodes." << endl;
  DistributeSolutionToNodes(config, geometry_dst, solver_dst);

  /*--------------------------------------------------------------------------*/
  /*--- Step 8: Handle surface nodes via linear interpolation.             ---*/
  /*--------------------------------------------------------------------------*/
  if (rank == MASTER_NODE) cout << "Performing surface interpolation." << endl;
  SurfaceInterpolation(geometry_src, geometry_dst, solver_src, solver_dst);

  /*--------------------------------------------------------------------------*/
  /*--- Step 9: Extrapolate solution to uncontained vertices.              ---*/
  /*--------------------------------------------------------------------------*/
  if (!uncontainedNodes.empty()) {
    if (rank == MASTER_NODE) cout << "Performing extrapolation to uncontained nodes." << endl;
    ExtrapolateToUncontainedNodes(geometry_dst, solver_dst);
  }
}

void CConservativeVolumeInterpolator::PointLocalization(CGeometry* geometry_src,
                                                        CGeometry* geometry_dst) {
  /*--- Search for containing elements for the given coordinates ---*/
  CADTElemClass& volumeADT = GetSourceVolumeADT();
  CADTElemClass& surfaceADT = GetSourceSurfaceADT();

  /*--- Loop over the DOFs to be interpolated ---*/
  containingElems.clear();
  containingElemRanks.clear();
  uncontainedNodes.clear();

  /*--- Initialize with invalid values ---*/
  containingElems.resize(nPoint_dst, nullopt);
  containingElemRanks.resize(nPoint_dst, -1);

  unsigned long numFailed = 0;
  for (auto l = 0u; l < nPoint_dst; ++l) {
    /*--- Set a pointer to the coordinates to be searched ---*/
    const su2double* coor = geometry_dst->nodes->GetCoord(l);

    /*--- Carry out the containment search and check if it was successful ---*/
    unsigned short subElemID;
    unsigned long elemID;
    int rankID;
    su2double parCoor[3], weightsInterpol[8];

    const bool success = volumeADT.DetermineContainingElement(coor, subElemID, elemID, rankID,
                                                              parCoor, weightsInterpol);
    if (success) {
      /*--- Store element information ---*/
      containingElems[l] = elemID;
      containingElemRanks[l] = rankID;
    } else {
      /*--- No source volume element contains this node - fallback         ---*/
      /*--- The nearest volume element will be used to find intersections  ---*/
      /*--- but also add to the list of failed points, which will be       ---*/
      /*--- processed later via extrapolation.                             ---*/
      /*--- If the node is a surface node, it will be linearly             ---*/
      /*--- interpolated anyway.                                           ---*/
      if (!geometry_dst->nodes->GetPhysicalBoundary(l)) uncontainedNodes.push_back(l);
      if (geometry_dst->nodes->GetDomain(l)) numFailed++;

      /*--- Find nearest volume element ---*/
      su2double dist;
      volumeADT.DetermineNearestElement(coor, dist, subElemID, elemID, rankID);
      containingElems[l] = elemID;
      containingElemRanks[l] = rankID;
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Point localizaion finished. " << numFailed;
    cout << " points failed containment search." << endl << flush;
  }
}

void CConservativeVolumeInterpolator::CreateIntersectionMeshes(CGeometry* geometry_src,
                                                               CGeometry* geometry_dst) {
  overlapMeshes.clear();

  /*--- Storage for triangle vertex coordinates ---*/
  su2double srcTri[6], dstTri[6];

  /*--- Storage for result of intersection test ---*/
  vector<su2double> pointCloud;
  vector<su2double> intersectionElemCoords;
  vector<su2double> intersectionElemVols;
  set<unsigned long> detectedCandidates;

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
    set<unsigned long> candidateElems;

    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      /*--- Find corresponding point index in corrected coordinates ---*/
      if (containingElems[nodeID].has_value()) {
        candidateElems.insert(containingElems[nodeID].value());
      }
    }

    /*--- Skip this destination element if no initial candidates found ---*/
    if (candidateElems.empty()) continue;

    /*--- Get destination triangle vertices ---*/
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = dstElem->GetNode(iNode);
      for (auto iDim = 0u; iDim < nDim; ++iDim)
        dstTri[iNode * nDim + iDim] = geometry_dst->nodes->GetCoord(nodeID, iDim);
    }

    /*--------------------------------------------------------------------------*/
    /*--- Step 2: Process candidate list, adding new elements during         ---*/
    /*---         intersection, and meshing the convex polygon formed by the ---*/
    /*---         intersection                                               ---*/
    /*--------------------------------------------------------------------------*/
    set<unsigned long> processedElems;
    queue<unsigned long> toProcess;

    /*--- Initialize queue with initial candidate list ---*/
    for (auto srcElemID : candidateElems) {
      toProcess.push(srcElemID);
      processedElems.insert(srcElemID);
    }

    /*--- Process elements, potentially adding new ones during intersection ---*/
    while (!toProcess.empty()) {
      unsigned long srcElemID = toProcess.front();
      toProcess.pop();

      auto* srcElem = geometry_src->elem[srcElemID];
      if (srcElem->GetVTK_Type() != TRIANGLE) continue;

      /*--- Get source triangle vertices ---*/
      for (auto iNode = 0u; iNode < 3; ++iNode) {
        unsigned long nodeID = srcElem->GetNode(iNode);
        for (auto iDim = 0u; iDim < nDim; ++iDim)
          srcTri[iNode * nDim + iDim] = geometry_src->nodes->GetCoord(nodeID, iDim);
      }

      /*--- Check for intersection and detect new candidates ---*/
      if (TriangleTriangleIntersection(geometry_src, dstTri, srcTri, srcElemID, pointCloud,
                                       intersectionElemCoords, intersectionElemVols,
                                       detectedCandidates)) {
        overlapMeshes[dstElemID].push_back({srcElemID, intersectionElemCoords, intersectionElemVols});
        totalOverlaps++;
      } else {
        /*--- Triangle intersection failed ---*/
        failedIntersections++;
      }

      /*--- Add new candidates to processing queue ---*/
      for (auto candidate : detectedCandidates) {
        if (processedElems.count(candidate) == 0) {
          toProcess.push(candidate);
          processedElems.insert(candidate);
        }
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Found " << totalOverlaps << " total overlapping element pairs." << endl;
    cout << "Number of destination elements with overlaps: " << overlapMeshes.size() << endl;
    cout << "Number of failed triangle intersections: " << failedIntersections << "." << endl;
  }
}

void CConservativeVolumeInterpolator::FindNearestContainedNodes(CGeometry* geometry_dst) {
  /*--- Clear and resize the nearestNodes vector ---*/
  nearestNodes.clear();
  nearestNodes.resize(uncontainedNodes.size(), nullopt);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Identify contained nodes that have at least one            ---*/
  /*---         uncontained neighbor (the "front" of contained nodes)      ---*/
  /*--------------------------------------------------------------------------*/
  set<unsigned long> uncontainedSet(uncontainedNodes.begin(), uncontainedNodes.end());
  vector<unsigned long> frontNodes;
  vector<su2double> frontCoords;

  /*--- Loop through all domain nodes to find front nodes ---*/
  for (auto nodeID = 0u; nodeID < nPoint_dst; ++nodeID) {
    /*--- Skip if this node is uncontained, boundary, or not domain ---*/
    if (uncontainedSet.count(nodeID) > 0 ||
        geometry_dst->nodes->GetPhysicalBoundary(nodeID) ||
        !geometry_dst->nodes->GetDomain(nodeID)) {
      continue;
    }

    /*--- Check if this contained node has any uncontained neighbors ---*/
    bool hasUncontainedNeighbor = false;
    for (auto iPoint = 0u; iPoint < geometry_dst->nodes->GetnPoint(nodeID); ++iPoint) {
      unsigned long neighborID = geometry_dst->nodes->GetPoint(nodeID, iPoint);
      if (uncontainedSet.count(neighborID) > 0) {
        hasUncontainedNeighbor = true;
        break;
      }
    }

    /*--- If this is a front node, add it to the list ---*/
    if (hasUncontainedNeighbor) {
      frontNodes.push_back(nodeID);
      const su2double* coord = geometry_dst->nodes->GetCoord(nodeID);
      for (auto iDim = 0u; iDim < nDim; ++iDim) {
        frontCoords.push_back(coord[iDim]);
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Found " << frontNodes.size() << " front nodes (contained nodes with uncontained neighbors)." << endl;
  }

  if (frontNodes.empty()) {
    /*--- No front nodes found ---*/
    if (rank == MASTER_NODE)
      cout << "No front nodes found. All uncontained nodes will have no nearest node." << endl;
    return;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Build ADT for the front nodes                              ---*/
  /*--------------------------------------------------------------------------*/
  vector<unsigned long> frontPointIDs(frontNodes.size());
  for (auto i = 0u; i < frontNodes.size(); ++i) {
    frontPointIDs[i] = i;
  }

  CADTPointsOnlyClass frontADT(nDim, frontNodes.size(), frontCoords.data(), frontPointIDs.data(), false);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: For each uncontained node, find nearest front node         ---*/
  /*--------------------------------------------------------------------------*/
  for (auto i = 0u; i < uncontainedNodes.size(); ++i) {
    unsigned long uncontainedNodeID = uncontainedNodes[i];
    const su2double* uncontainedCoord = geometry_dst->nodes->GetCoord(uncontainedNodeID);

    /*--- Find nearest front node ---*/
    su2double dist;
    unsigned long nearestFrontNodeIndex;
    int rankID;  // Not used in local ADT
    frontADT.DetermineNearestNode(uncontainedCoord, dist, nearestFrontNodeIndex, rankID);

    /*--- Store the actual node ID of the nearest front node ---*/
    nearestNodes[i] = frontNodes[nearestFrontNodeIndex];
  }

  if (rank == MASTER_NODE) {
    cout << "Generated nearest contained nodes for " << uncontainedNodes.size();
    cout << " uncontained nodes." << endl;
  }
}

void CConservativeVolumeInterpolator::ComputeSourceSolutionMass(CGeometry* geometry_src,
                                                                CSolver* solver_src) {
  /*--- Initialize source solution mass and gradient ---*/
  for (auto& v : srcElemMass) std::fill(v.begin(), v.end(), 0.0);
  for (auto& v : srcElemGrad) std::fill(v.begin(), v.end(), 0.0);

  su2double vertexCoords[6];

  for (auto elemID = 0u; elemID < geometry_src->GetnElem(); ++elemID) {
    auto* elem = geometry_src->elem[elemID];
    const unsigned short nNodes = elem->GetnNodes();
    const unsigned short VTK_Type = elem->GetVTK_Type();
    const su2double elemVolume = elem->GetVolume();

    /*--- Get solution at nodes ---*/
    vector<vector<su2double>> vertexSol(nNodes, vector<su2double>(nVar, 0.0));

    for (auto iNode = 0u; iNode < nNodes; ++iNode) {
      unsigned long nodeID = elem->GetNode(iNode);
      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        vertexSol[iNode][iVar] = solver_src->GetNodes()->GetSolution(nodeID, iVar);
      }
    }

    if (VTK_Type != TRIANGLE) continue;

    /*--- Get triangle vertex coordinates ---*/
    for (auto iNode = 0u; iNode < 3; ++iNode) {
      unsigned long nodeID = elem->GetNode(iNode);
      vertexCoords[iNode * 2 + 0] = geometry_src->nodes->GetCoord(nodeID, 0);
      vertexCoords[iNode * 2 + 1] = geometry_src->nodes->GetCoord(nodeID, 1);
    }

    /*--- Compute mass and gradient ---*/
    ComputeTriangleMassAndGradient(vertexCoords, vertexSol, elemVolume, nVar,
                                   srcElemMass[elemID], srcElemGrad[elemID]);
  }
}

void CConservativeVolumeInterpolator::ComputeDestinationMassAndGradient(CGeometry* geometry_src,
                                                                        CGeometry* geometry_dst,
                                                                        CSolver* solver_src) {
  /*--- Initialize destination solution mass and gradient ---*/
  for (auto& v : dstElemMass) std::fill(v.begin(), v.end(), 0.0);
  for (auto& v : dstElemGrad) std::fill(v.begin(), v.end(), 0.0);


  const unsigned short nCoorPerElem = (nDim == 2)? 6 : 12;

  /*--- Loop over all destination elements that have intersections ---*/
  su2double absDiffTol[5] = {1e-10, 1e-8, 1e-6, 1e-4, 1e-2};
  su2double relDiffTol[5] = {2e-2, 5e-2, 1e-1, 2e-1, 5e-1};
  vector<unsigned long> countAbsDiff(5, 0);
  vector<unsigned long> countRelDiff(5, 0);

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

    /*--- Track total triangle area for this destination element ---*/
    su2double intersectionVol = 0.0;

    /*--- Process each intersection region T_j = intersection(K_dst, K_src_j) ---*/
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

        /*--- Triangle area using cross product ---*/
        const su2double cross = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        const su2double area = 0.5 * abs(cross);

        /*--- Add to total triangle area ---*/
        intersectionVol += area;

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
          for (auto iDim = 0u; iDim < nDim; ++iDim)
            dstGrad[iVar * nDim + iDim] += area * grad_u[iDim];
        }
      }
    }

    /*--- Volume average integral: gra(u_dst) = int_K_dst (gra(u) dA) / |K_dst| ---*/
    const su2double dstVolume = dstElem->GetVolume();

    /*--- Normalize by intersection volume instead of destination volume ---*/
    /*--- For matching volumes, this should have no impact               ---*/
    /*--- For non-matching domains, this preserves P1-exactness          ---*/
    /*--- This preserves constant/linear solutions in non-matching domains ---*/
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      /*--- For mass: divide by intersection volume to get correct barycenter value ---*/
      dstMass[iVar] = dstMass[iVar] * dstVolume / intersectionVol;

      /*--- For gradient: use intersection-weighted average ---*/
      for (auto iDim = 0u; iDim < nDim; ++iDim) {
        dstGrad[iVar * nDim + iDim] /= intersectionVol;
      }
    }

    /*--- Compare total triangle area with destination element volume ---*/
    const su2double absDiff = abs(dstVolume - intersectionVol);
    const su2double relDiff = absDiff / dstVolume;
    for (auto i = 0u; i < 5; ++i) {
      if (absDiff > absDiffTol[i]) countAbsDiff[i]++;
      if (relDiff > relDiffTol[i]) countRelDiff[i]++;
    }
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
  }
}

void CConservativeVolumeInterpolator::ApplyMaximumPrincipleCorrection(CGeometry* geometry_src,
                                                                      CGeometry* geometry_dst,
                                                                      CSolver* solver_src) {
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
      dstVertices[iNode * 2 + 0] = geometry_dst->nodes->GetCoord(nodeID, 0);
      dstVertices[iNode * 2 + 1] = geometry_dst->nodes->GetCoord(nodeID, 1);
    }

    /*--- Element volume and centroid (barycenter G_K) ---*/
    const su2double elemVolume = dstElem->GetVolume();
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
      const su2double u_G = dstElemMass[dstElemID][iVar] / elemVolume;
      const su2double* gradu_G = dstElemGrad[dstElemID].data() + iVar * nDim;

      /*--------------------------------------------------------------------------*/
      /*--- Step 3: Compute u_K(P_i) at each vertex using Taylor expansion.    ---*/
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
      ComputeTriangleMassAndGradient(dstVertices, vertexSol, elemVolume, 1,
                                     correctedMass, correctedGrad);

      /*--- Update destination element data ---*/
      dstElemMass[dstElemID][iVar] = correctedMass[0];
      dstElemGrad[dstElemID][iVar * nDim + 0] = correctedGrad[0];
      dstElemGrad[dstElemID][iVar * nDim + 1] = correctedGrad[1];
    }
  }
}

void CConservativeVolumeInterpolator::DistributeSolutionToNodes(const CConfig* config,
                                                                CGeometry* geometry_dst,
                                                                CSolver* solver_dst) {
  /*--- Loop over all nodes and accumulate contributions from each element ---*/
  for (auto l = 0u; l < nPoint_dst; ++l) {
    /*--- Get pointer to node coordinates and solution mass ---*/
    const su2double* coor = geometry_dst->nodes->GetCoord(l);

    /*--- Initialize mass and volume accumulation ---*/
    for (auto iVar = 0u; iVar < nVar+1; ++iVar)
      solver_dst->GetNodes()->SetSolution_Mass(l, iVar, 0.0);

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
        solver_dst->GetNodes()->AddSolution_Mass(l, iVar, vertexValue * elemVolume);
      }

      /*--- Accumulate weight ---*/
      solver_dst->GetNodes()->AddSolution_Mass(l, nVar, elemVolume);
    }
  }

  /*--- Periodic mass and volume communication ---*/
  for (auto i = 1u; i <= config->GetnMarker_Periodic()/2; ++i) {
    solver_dst->InitiatePeriodicComms(geometry_dst, config, i, PERIODIC_INTERP);
    solver_dst->CompletePeriodicComms(geometry_dst, config, i, PERIODIC_INTERP);
  }

  /*--- Calculate conserved quantities from mass and volume ---*/
  for (auto l = 0u; l < nPoint_dst; ++l) {
    const su2double* mass = solver_dst->GetNodes()->GetSolution_Mass(l);
    const su2double vol = mass[nVar];
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      solver_dst->GetNodes()->SetSolution(l, iVar, mass[iVar] / vol);
    }
  }
  solver_dst->Set_OldSolution();
}

void CConservativeVolumeInterpolator::ExtrapolateToUncontainedNodes(CGeometry* geometry_dst,
                                                                    CSolver* solver_dst) {
  unsigned long numExtrapolated = 0;
  unsigned long numFallback = 0;
  unsigned long numFailed = 0;

  /*--- Process each uncontained node ---*/
  for (auto i = 0u; i < uncontainedNodes.size(); ++i) {
    unsigned long nodeID = uncontainedNodes[i];

    /*--- Skip boundary nodes (they should be handled by surface interpolation) ---*/
    if (geometry_dst->nodes->GetPhysicalBoundary(nodeID)) continue;

    /*--- Get nearest contained node ID ---*/
    auto nearestNodeOpt = nearestNodes[i];

    /*--- Check if we have a valid nearest node ---*/
    if (nearestNodeOpt.has_value()) {
      unsigned long nearestNodeID = nearestNodeOpt.value();
      /*--- Extrapolate from the nearest contained node ---*/
      if (ExtrapolateFromNearestNode(geometry_dst, solver_dst, nodeID, nearestNodeID)) {
        numExtrapolated++;
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Extrapolation completed: " << numExtrapolated << " nodes from nearest nodes, ";
    cout << numFallback << " nodes via fallback, " << numFailed << " nodes failed." << endl;
  }
}

bool CConservativeVolumeInterpolator::TriangleTriangleIntersection(CGeometry* geometry_src,
                                                                   const su2double dstTri[6],
                                                                   const su2double srcTri[6],
                                                                   unsigned long srcElemID,
                                                                   vector<su2double>& pointCloud,
                                                                   vector<su2double>& intersectionElemCoords,
                                                                   vector<su2double>& intersectionElemVols,
                                                                   set<unsigned long>& candidateElems) {
  pointCloud.clear();
  intersectionElemCoords.clear();
  intersectionElemVols.clear();
  candidateElems.clear();

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
  /*---         power_P[i][j] is the power of vertex Pi w.r.t. edge j of   ---*/
  /*---         triangle Q.                                                ---*/
  /*--------------------------------------------------------------------------*/
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
  const su2double* P_edges[3][2] = {{P[0], P[1]}, {P[1], P[2]}, {P[2], P[0]}};
  const su2double* Q_edges[3][2] = {{Q[0], Q[1]}, {Q[1], Q[2]}, {Q[2], Q[0]}};

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Handle degenerate edge-edge intersection cases.            ---*/
  /*--------------------------------------------------------------------------*/
  bool isDegenerateVertexP[3] = {};
  bool isDegenerateVertexQ[3] = {};
  bool isDegenerateEdgePair[3][3] = {};
  ProcessDegenerateEdgeIntersections(P_edges, Q_edges, power_P, power_Q, EPS, pointCloud,
                                     isDegenerateVertexP, isDegenerateVertexQ, isDegenerateEdgePair);

  /*--- Add vertex ball of all degenerate vertices of Q ---*/
  for (auto i = 0u; i < 3; ++i) {
    if (isDegenerateVertexQ[i]) {
      AddVertexBallToCandidates(geometry_src, srcElemID, i, candidateElems);
    }
  }

  /*--- Check if KP vertices are strictly inside KQ ---*/
  for (auto i = 0u; i < 3; ++i) {
    if ((power_P[i][0] > -EPS) &&
        (power_P[i][1] > -EPS) &&
        (power_P[i][2] > -EPS) &&
        (!isDegenerateVertexP[i])) {
      pointCloud.push_back(P[i][0]);
      pointCloud.push_back(P[i][1]);
    }
  }

  /*--- Check if KQ vertices are strictly inside KP ---*/
  for (auto i = 0u; i < 3; ++i) {
    if ((power_Q[i][0] > -EPS) &&
        (power_Q[i][1] > -EPS) &&
        (power_Q[i][2] > -EPS) &&
        (!isDegenerateVertexQ[i])) {
      pointCloud.push_back(Q[i][0]);
      pointCloud.push_back(Q[i][1]);
      AddVertexBallToCandidates(geometry_src, srcElemID, i, candidateElems);
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

      const su2double* P0 = P_edges[iP][0];
      const su2double* P1 = P_edges[iP][1];
      const su2double* Q0 = Q_edges[jQ][0];
      const su2double* Q1 = Q_edges[jQ][1];

      /*--- Compute signed distances for this edge pair ---*/
      const su2double powerP0_eQ = ComputeSignedDistance(P0, Q0, Q1);
      const su2double powerP1_eQ = ComputeSignedDistance(P1, Q0, Q1);
      const su2double powerQ0_eP = ComputeSignedDistance(Q0, P0, P1);
      const su2double powerQ1_eP = ComputeSignedDistance(Q1, P0, P1);

      /*--- Check for edge-edge intersection ---*/
      if ((powerP0_eQ * powerP1_eQ < 0) && (powerQ0_eP * powerQ1_eP < 0)) {
        /*--- Compute intersection point ---*/
        su2double t = powerP0_eQ / (powerP0_eQ - powerP1_eQ);

        su2double intersectionX = P0[0] + t * (P1[0] - P0[0]);
        su2double intersectionY = P0[1] + t * (P1[1] - P0[1]);

        pointCloud.push_back(intersectionX);
        pointCloud.push_back(intersectionY);

        wasIntersected = true;
      }
    }

    if (wasIntersected) {
      /*--- Add face neighbor to candidates ---*/
      AddFaceNeighborToCandidates(geometry_src, srcElemID, jQ, candidateElems);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Mesh the intersection polygon.                             ---*/
  /*--------------------------------------------------------------------------*/
  MeshConvexPolygon(pointCloud, intersectionElemCoords, intersectionElemVols);

  return (!intersectionElemCoords.empty());
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

void CConservativeVolumeInterpolator::ProcessDegenerateEdgeIntersections(const su2double* P_edges[3][2],
                                                                         const su2double* Q_edges[3][2],
                                                                         const su2double power_P[3][3],
                                                                         const su2double power_Q[3][3],
                                                                         const su2double EPS,
                                                                         vector<su2double>& pointCloud,
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
      pointCloud.push_back(x);
      pointCloud.push_back(y);
    }
  };

  /*--- Process all 9 edge-edge combinations for degenerate cases ---*/
  for (auto iP = 0u; iP < 3; ++iP) {
    for (auto jQ = 0u; jQ < 3; ++jQ) {
      /*--- Indices of nodes belonging to current edges ---*/
      int iP0 = iP, iP1 = (iP + 1) % 3;
      int jQ0 = jQ, jQ1 = (jQ + 1) % 3;

      /*--- Coordinates of edge nodes ---*/
      const su2double* P0 = P_edges[iP][0];  // P vertex iP
      const su2double* P1 = P_edges[iP][1];  // P vertex (iP+1)%3
      const su2double* Q0 = Q_edges[jQ][0];  // Q vertex jQ
      const su2double* Q1 = Q_edges[jQ][1];  // Q vertex (jQ+1)%3

      /*--- Use precomputed powers instead of recalculating ---*/
      /*--- Edge iP of triangle P goes from vertex iP to vertex (iP+1)%3 ---*/
      /*--- Edge jQ of triangle Q goes from vertex jQ to vertex (jQ+1)%3 ---*/
      const su2double powerP0_eQ = power_P[iP0][jQ];  // Power of P vertex iP w.r.t. Q edge jQ
      const su2double powerP1_eQ = power_P[iP1][jQ];  // Power of P vertex (iP+1)%3 w.r.t. Q edge jQ
      const su2double powerQ0_eP = power_Q[jQ0][iP];  // Power of Q vertex jQ w.r.t. P edge iP
      const su2double powerQ1_eP = power_Q[jQ1][iP];  // Power of Q vertex (jQ+1)%3 w.r.t. P edge iP

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

void CConservativeVolumeInterpolator::AddFaceNeighborToCandidates(CGeometry* geometry,
                                                                  const unsigned long elemID,
                                                                  const unsigned short faceIndex,
                                                                  set<unsigned long>& candidateElems) {
  /*--- Get the neighbor element across this face ---*/
  long neighElemID = geometry->elem[elemID]->GetNeighbor_Elements(faceIndex);

  /*--- Only add if a valid neighbor exists (not boundary) ---*/
  if (neighElemID >= 0) {
    candidateElems.insert(static_cast<unsigned long>(neighElemID));
  }
}

void CConservativeVolumeInterpolator::AddVertexBallToCandidates(CGeometry* geometry_src,
                                                                const unsigned long srcElemID,
                                                                const unsigned short localNodeID,
                                                                set<unsigned long>& candidateElems) {
  /*--- Get the vertex node ID ---*/
  auto* srcElem = geometry_src->elem[srcElemID];
  unsigned long nodeID = srcElem->GetNode(localNodeID);

  /*--- Add all elements other than the current source element that contain this vertex ---*/
  for (auto jElem = 0u; jElem < geometry_src->nodes->GetnElem(nodeID); ++jElem) {
    unsigned long elemID = geometry_src->nodes->GetElem(nodeID, jElem);
    if (elemID == srcElemID) continue;  // Skip self

    auto* neighborElem = geometry_src->elem[elemID];
    if (neighborElem->GetVTK_Type() != TRIANGLE) continue;

    candidateElems.insert(elemID);
  }
}

void CConservativeVolumeInterpolator::MeshConvexPolygon(const vector<su2double>& polygonPoints,
                                                        vector<su2double>& polygonElemCoords,
                                                        vector<su2double>& polygonElemVols) {
  polygonElemCoords.clear();
  polygonElemVols.clear();

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

bool CConservativeVolumeInterpolator::ExtrapolateFromNearestNode(CGeometry* geometry_dst,
                                                                 CSolver* solver_dst,
                                                                 unsigned long uncontainedNodeID,
                                                                 unsigned long nearestNodeID) {
  /*--- Get coordinates ---*/
  const su2double* uncontainedCoord = geometry_dst->nodes->GetCoord(uncontainedNodeID);
  const su2double* nearestCoord = geometry_dst->nodes->GetCoord(nearestNodeID);

  /*--- Displacement vector ---*/
  su2double displacement[3];
  for (auto iDim = 0u; iDim < nDim; ++iDim)
    displacement[iDim] = uncontainedCoord[0] - nearestCoord[0];

  /*--- Get solution at nearest contained node ---*/
  vector<su2double> nearestSolution(nVar);
  for (auto iVar = 0u; iVar < nVar; ++iVar) {
    nearestSolution[iVar] = solver_dst->GetNodes()->GetSolution(nearestNodeID, iVar);
  }

  /*--- Compute volume-weighted average gradient from elements containing the nearest node ---*/
  vector<su2double> avgGradient(nVar * nDim, 0.0);
  su2double totalWeight = 0.0;

  const unsigned short nElem_node = geometry_dst->nodes->GetnElem(nearestNodeID);
  for (auto j = 0u; j < nElem_node; ++j) {
    unsigned long elemID = geometry_dst->nodes->GetElem(nearestNodeID, j);
    auto* elem = geometry_dst->elem[elemID];

    /*--- Skip non-triangular elements ---*/
    if (elem->GetVTK_Type() != TRIANGLE) continue;

    /*--- Check if element has valid solution (was processed in conservative interpolation) ---*/
    bool hasValidSolution = false;
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      if (abs(dstElemMass[elemID][iVar]) > 1e-16) {
        hasValidSolution = true;
        break;
      }
    }
    if (!hasValidSolution) continue;

    /*--- Use element volume as weight ---*/
    const su2double elemVolume = elem->GetVolume();
    totalWeight += elemVolume;

    /*--- Add weighted contribution to average gradient ---*/
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      for (auto iDim = 0u; iDim < nDim; ++iDim) {
        avgGradient[iVar * nDim + iDim] += dstElemGrad[elemID][iVar * nDim + iDim] * elemVolume;
      }
    }
  }

  /*--- Normalize gradient ---*/
  for (auto iVar = 0u; iVar < nVar; ++iVar) {
    for (auto iDim = 0u; iDim < nDim; ++iDim) {
      avgGradient[iVar * nDim + iDim] /= totalWeight;
    }
  }

  /*--- Extrapolate solution using linear reconstruction ---*/
  for (auto iVar = 0u; iVar < nVar; ++iVar) {
    su2double extrapolatedValue = nearestSolution[iVar] +
                                  avgGradient[iVar * nDim + 0] * displacement[0] +
                                  avgGradient[iVar * nDim + 1] * displacement[1];

    solver_dst->GetNodes()->SetSolution(uncontainedNodeID, iVar, extrapolatedValue);
    solver_dst->GetNodes()->SetSolution_Old(uncontainedNodeID, iVar, extrapolatedValue);
  }

  return true;
}
