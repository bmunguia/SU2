/*!
 * \file CVolumeInterpolator.cpp
 * \brief Implementation of the main solution interpolation subroutines.
 *        This file contains all utility functions used for interpolation.
 * \author B. Munguía, E. van der Weide
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

 #include "../include/CVolumeInterpolator.hpp"

CVolumeInterpolator::CVolumeInterpolator(SU2_Comm MPICommunicator) {
  /*--- Set up MPI ---*/
  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
}

void CVolumeInterpolator::InitializeConfig(CConfig* driver_config, CConfig** config_container, char* zone_file_name,
                                           char* config_file_name, SU2_COMPONENT val_software, int iZone, int nZone,
                                           SU2_MPI::Comm MPICommunicator, bool isSource) {
  if (driver_config->GetnConfigFiles() > 0) {
    strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
    config_container[iZone] = new CConfig(driver_config, zone_file_name, val_software, iZone, nZone, true);
  } else {
    config_container[iZone] = new CConfig(driver_config, config_file_name, val_software, iZone, nZone, true);
  }
  config_container[iZone]->SetMPICommunicator(MPICommunicator);

  /*--- Set the multizone part of the problem. ---*/
  const bool multizone = config_container[iZone]->GetMultizone_Problem();
  if (multizone) {
    /*--- Set the interface markers for multizone ---*/
    config_container[iZone]->SetMultizone(driver_config, config_container);
  }

  /*--- Read MESH_FILENAME if source, otherwise read destination mesh ---*/
  if (!isSource) {
    /*--- Read MESH_ITP_FILENAME ---*/
    config_container[iZone]->SetMesh_FileName(config_container[iZone]->GetMesh_Itp_FileName());
  }
}

void CVolumeInterpolator::InitializeGeometry(CConfig* config, CGeometry*& geometry, int iZone, int iInst,
                                             int nZone, bool isSource) {
  /*--- Mesh initialization ---*/
  config->SetiInst(iInst);
  CGeometry* geometry_aux = nullptr;
  geometry_aux = new CPhysicalGeometry(config, iZone, nZone);
  const bool fem_solver = config->GetFEMSolver();
  if (fem_solver)
    geometry_aux->SetColorFEMGrid_Parallel(config);
  else
    geometry_aux->SetColorGrid_Parallel(config);

  /*--- Allocate the memory of the current domain, and divide the grid between the nodes ---*/
  geometry = nullptr;

  /*--- Build the grid data structures using the ParMETIS coloring. ---*/
  if (fem_solver) {
    switch (config->GetKind_FEM_Flow()) {
      case DG: {
        geometry = new CMeshFEM_DG(geometry_aux, config);
        break;
      }
    }
  } else {
    geometry = new CPhysicalGeometry(geometry_aux, config);
  }
  delete geometry_aux;

  /*--- Add the Send/Receive boundaries ---*/
  geometry->SetSendReceive(config);

  /*--- Add the Send/Receive boundaries ---*/
  geometry->SetBoundaries(config);

  /*--- Compute elements surrounding points, points surrounding points ---*/
  if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
  geometry->SetPoint_Connectivity();

  // /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/
  // if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
  // geometry->SetRCM_Ordering(config);

  // /*--- Recompute elements surrounding points, points surrounding points ---*/
  // if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
  // geometry->SetPoint_Connectivity();

  /*--- Compute elements surrounding elements ---*/
  if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
  geometry->SetElement_Connectivity();

  /*--- Check the orientation before computing geometrical quantities ---*/
  geometry->SetBoundVolume();
  if (config->GetReorientElements()) {
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry->Check_IntElem_Orientation(config);
    geometry->Check_BoundElem_Orientation(config);
  }

  /*--- Create the edge structure ---*/
  if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
  geometry->SetEdges();
  geometry->SetVertex(config);

  /*--- Compute element volumes ---*/
  if (rank == MASTER_NODE) cout << "Setting element volumes." << endl;
  geometry->SetElemVolume();

  /*--- Create the control volume structures ---*/
  if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
  SU2_OMP_PARALLEL {
    geometry->SetControlVolume(config, ALLOCATE);
    geometry->SetBoundControlVolume(config, ALLOCATE);
  }
  END_SU2_OMP_PARALLEL

  /*--- Store the global to local mapping ---*/
  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
  geometry->SetGlobal_to_Local_Point();

  /*--- Create the data structure for MPI point-to-point communications ---*/
  if (!fem_solver) geometry->PreprocessP2PComms(geometry, config);
  if (fem_solver) {
    auto* DGMesh = dynamic_cast<CMeshFEM_DG*>(geometry);
    if (rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;
    DGMesh->CreateStandardVolumeElements(config);
    if (rank == MASTER_NODE) cout << "Creating face information." << endl;
    DGMesh->CreateFaces(config);
  }

  /*--- If we have any periodic markers in this calculation, we must
       match the periodic points found on both sides of the periodic BC.
       Note that the current implementation requires a 1-to-1 matching of
       periodic points on the pair of periodic faces after the translation
       or rotation is taken into account. ---*/

  if ((config->GetnMarker_Periodic() != 0) && !fem_solver) {
    /*--- Note that we loop over pairs of periodic markers individually
          so that repeated nodes on adjacent periodic faces are properly
          accounted for in multiple places. ---*/

    for (auto iPer = 1u; iPer <= config->GetnMarker_Periodic()/2; ++iPer) {
      geometry->MatchPeriodic(config, iPer);
    }

    /*--- For Streamwise Periodic flow, find a unique reference node on the dedicated inlet marker. ---*/
    if (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)
      geometry->FindUniqueNode_PeriodicBound(config);

    /*--- Initialize the communication framework for the periodic BCs. ---*/
    geometry->PreprocessPeriodicComms(geometry, config);
  }

  if (isSource) {
    nPoint_src = geometry->GetnPoint();
    nElem_src = geometry->GetnElem();
  } else {
    nPoint_dst = geometry->GetnPoint();
    nElem_dst = geometry->GetnElem();
  }

  if (nDim != 0) assert(geometry->GetnDim() == nDim);
  else nDim = geometry->GetnDim();
}

void CVolumeInterpolator::InitializeSolver(CConfig* config, CGeometry* geometry, CSolver**& solver_container, int iZone,
                                           int iInst, int nZone) {
  switch (config->GetKindVolumeInterpolation()) {
    case VOLUME_INTERPOLATOR::LINEAR: {
      /*--- Initialize a baseline solver ---*/
      solver_container = new CSolver*[MAX_SOLS]();
      solver_container[FLOW_SOL] = new CBaselineSolver(geometry, config);
      break;
    }
    case VOLUME_INTERPOLATOR::CONSERVATIVE: {
      /*--- Only interpolating the conserved quantities, so initialize specific solvers ---*/
      MAIN_SOLVER kindSolver = config->GetKind_Solver();
      solver_container = CSolverFactory::CreateSolverContainer(kindSolver, config, geometry, 0);
      break;
    }
    default:
      SU2_MPI::Error("Volume interpolator not supported", CURRENT_FUNCTION);
  }

}

void CVolumeInterpolator::InitializeOutput(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput*& output,
                                           int iZone, int iInst, int nZone) {
  switch (config->GetKindVolumeInterpolation()) {
    case VOLUME_INTERPOLATOR::LINEAR: {
      /*--- Initialize a baseline output ---*/
      output = new CBaselineOutput(config, geometry->GetnDim(), solver_container[FLOW_SOL]);
      break;
    }
    case VOLUME_INTERPOLATOR::CONSERVATIVE: {
      /*--- Solver-specific output ---*/
      MAIN_SOLVER kindSolver = config->GetKind_Solver();
      output = COutputFactory::CreateOutput(kindSolver, config, geometry->GetnDim());
      break;
    }
    default:
      SU2_MPI::Error("Volume interpolator not supported", CURRENT_FUNCTION);
  }

  output->PreprocessVolumeOutput(config);
  output->PreprocessHistoryOutput(config, false);

}

void CVolumeInterpolator::LoadRestarts(CConfig* config, CGeometry** geometry_container, CSolver*** solver_container, int iZone,
                                       int iInst, int TimeIter, bool UpdateGeo) {
  for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
    auto solver = solver_container[iInst][iSol];
    if (solver)
      solver->LoadRestart(geometry_container, solver_container, config, TimeIter, UpdateGeo);
  }
}

void CVolumeInterpolator::InitializeADTs(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst, bool update) {
  if (!srcVolumeADT_ptr || update) {
    if (rank == MASTER_NODE) cout << "Building volume ADT." << flush;
    srcVolumeADT_ptr = BuildVolumeADT(geometry_src);
    if (rank == MASTER_NODE) cout << " Done." << endl;
  }

  if (!srcSurfaceADT_ptr || update) {
    if (rank == MASTER_NODE) cout << "Building source surface ADT." << flush;
    srcSurfaceADT_ptr = BuildSurfaceADT(config, geometry_src);
    if (rank == MASTER_NODE) cout << " Done." << endl;
  }

  if (!dstSurfaceADT_ptr || update) {
    if (rank == MASTER_NODE) cout << "Building destination surface ADT." << flush;
    dstSurfaceADT_ptr = BuildSurfaceADT(config, geometry_dst);
    if (rank == MASTER_NODE) cout << " Done." << endl;
  }
}

std::unique_ptr<CADTElemClass> CVolumeInterpolator::BuildVolumeADT(CGeometry* geometry) {
  /*--- Prepare coordinate and connectivity data for ADT ---*/
  vector<su2double> volCoor;
  vector<unsigned long> elemConn;
  vector<unsigned short> vtkType;
  vector<unsigned short> subElemID;
  vector<unsigned long> parElemID;

  /*--- Copy coordinates ---*/
  for (auto i = 0u; i <  geometry->GetnPoint(); ++i) {
    for (auto k = 0u; k < nDim; ++k) {
      volCoor.push_back(geometry->nodes->GetCoord(i, k));
    }
  }

  /*--- Copy element connectivity and metadata ---*/
  for (auto iElem = 0u; iElem < geometry->GetnElem(); ++iElem) {
    const unsigned short VTK_Type = geometry->elem[iElem]->GetVTK_Type();
    const unsigned short nDOFsPerElem = geometry->elem[iElem]->GetnNodes();

    vtkType.push_back(VTK_Type);
    subElemID.push_back(0);  // TODO: FEM subelements
    parElemID.push_back(iElem);
    for (auto iNode = 0u; iNode < nDOFsPerElem; ++iNode) {
      elemConn.push_back(geometry->elem[iElem]->GetNode(iNode));
    }
  }

  /*--- Create and return the volume ADT ---*/
  return std::unique_ptr<CADTElemClass>(
      new CADTElemClass(nDim, volCoor, elemConn, vtkType, subElemID, parElemID, false));
}

std::unique_ptr<CADTElemClass> CVolumeInterpolator::BuildSurfaceADT(const CConfig* config, CGeometry* geometry) {
  /*--- Initialize an array for the mesh points mapping ---*/
  vector<unsigned long> meshToSurface(geometry->GetnPoint(), 0);

  /*--- Define vectors for surface connectivity and metadata ---*/
  vector<unsigned long> surfaceConn;
  vector<unsigned long> elemIDs;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned short> markerIDs;

  /*--- Loop over boundary markers ---*/
  for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); ++iMarker) {
    /*--- Loop over surface elements of this marker ---*/
    for (auto iElem = 0u; iElem < geometry->GetnElem_Bound(iMarker); ++iElem) {
      const unsigned short VTK_Type = geometry->bound[iMarker][iElem]->GetVTK_Type();
      const unsigned short nDOFsPerElem = geometry->bound[iMarker][iElem]->GetnNodes();

      /*--- Set flag for mesh points on this surface ---*/
      for (auto iNode = 0u; iNode < nDOFsPerElem; ++iNode) {
        unsigned long iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
        meshToSurface[iPoint] = 1;
      }

      /*--- Store the required data ---*/
      markerIDs.push_back(iMarker);
      VTK_TypeElem.push_back(VTK_Type);
      elemIDs.push_back(iElem);
      for (auto iNode = 0u; iNode < nDOFsPerElem; ++iNode) {
        surfaceConn.push_back(geometry->bound[iMarker][iElem]->GetNode(iNode));
      }
    }
  }

  /*--- Create coordinates of points on surfaces ---*/
  vector<su2double> surfaceCoor;
  unsigned long nVertex_Surface = 0;

  for (auto i = 0u; i < geometry->GetnPoint(); ++i) {
    if (meshToSurface[i]) {
      meshToSurface[i] = nVertex_Surface++;
      for (auto k = 0u; k < nDim; ++k) {
        surfaceCoor.push_back(geometry->nodes->GetCoord(i, k));
      }
    }
  }

  /*--- Change surface connectivity to correspond to surfaceCoor indices ---*/
  for (auto i = 0u; i < surfaceConn.size(); ++i) {
    surfaceConn[i] = meshToSurface[surfaceConn[i]];
  }

  /*--- Build and return the surface ADT ---*/
  return std::unique_ptr<CADTElemClass>(
      new CADTElemClass(nDim, surfaceCoor, surfaceConn, VTK_TypeElem, markerIDs, elemIDs, true));
}

void CVolumeInterpolator::NearestPointOnElement(CGeometry* geometry, unsigned short markerID, unsigned long elemID,
                                                const su2double* coor, su2double* surfCoor, su2double& dist2Elem,
                                                const unsigned short nDim) {
  unsigned short VTK_Type = geometry->bound[markerID][elemID]->GetVTK_Type();
  unsigned short nNodes = geometry->bound[markerID][elemID]->GetnNodes();

  auto updateNearest = [&](su2double dist2Line, su2double* surfCoorLine) {
    if (dist2Line < dist2Elem) {
      dist2Elem = dist2Line;
      for (auto k = 0u; k < nDim; ++k) surfCoor[k] = surfCoorLine[k];
    }
  };

  switch(VTK_Type) {
    case LINE: {
      const unsigned long i0 = geometry->bound[markerID][elemID]->GetNode(0);
      const unsigned long i1 = geometry->bound[markerID][elemID]->GetNode(1);
      NearestPointOnLine(geometry, i0, i1, coor, surfCoor, dist2Elem, nDim);
      break;
    }
    case TRIANGLE: {
      const unsigned long i0 = geometry->bound[markerID][elemID]->GetNode(0);
      const unsigned long i1 = geometry->bound[markerID][elemID]->GetNode(1);
      const unsigned long i2 = geometry->bound[markerID][elemID]->GetNode(2);
      if (!NearestPointOnTriangle(geometry, i0, i1, i2, coor, surfCoor, dist2Elem, nDim)) {
        NearestPointOnLine(geometry, i0, i1, coor, surfCoor, dist2Elem, nDim);

        su2double dist2Line;
        su2double surfCoorLine[3];
        NearestPointOnLine(geometry, i1, i2, coor, surfCoorLine, dist2Line, nDim);
        updateNearest(dist2Line, surfCoorLine);
        NearestPointOnLine(geometry, i2, i0, coor, surfCoorLine, dist2Line, nDim);
        updateNearest(dist2Line, surfCoorLine);
      }
      break;
    }

    case QUADRILATERAL: {
      const unsigned long i0 = geometry->bound[markerID][elemID]->GetNode(0);
      const unsigned long i1 = geometry->bound[markerID][elemID]->GetNode(1);
      const unsigned long i2 = geometry->bound[markerID][elemID]->GetNode(2);
      const unsigned long i3 = geometry->bound[markerID][elemID]->GetNode(3);
      if (!NearestPointOnQuadrilateral(geometry, i0, i1, i2, i3, coor, surfCoor, dist2Elem, nDim)) {
        NearestPointOnLine(geometry, i0, i1, coor, surfCoor, dist2Elem, nDim);

        su2double dist2Line;
        su2double surfCoorLine[3];
        NearestPointOnLine(geometry, i1, i2, coor, surfCoorLine, dist2Line, nDim);
        updateNearest(dist2Line, surfCoorLine);
        NearestPointOnLine(geometry, i2, i3, coor, surfCoorLine, dist2Line, nDim);
        updateNearest(dist2Line, surfCoorLine);
        NearestPointOnLine(geometry, i3, i0, coor, surfCoorLine, dist2Line, nDim);
        updateNearest(dist2Line, surfCoorLine);
      }
      break;
    }

    default:
      /* This should not happen. */
      SU2_MPI::Error("Element type not recognized.", CURRENT_FUNCTION);

  }
}

void CVolumeInterpolator::NearestPointOnLine(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                                             const su2double* coor, su2double* surfCoor, su2double& dist2Line,
                                             const unsigned short nDim) {
  su2double x0[3], x1[3];
  for (auto k = 0u; k < nDim; ++k) {
    x0[k] = geometry->nodes->GetCoord(i0, k);
    x1[k] = geometry->nodes->GetCoord(i1, k);
  }

  /*--- Use the same parametrization as CADTElemClass::Dist2ToLine ---*/
  /*--- X = X0 + (r+1)*(X1-X0)/2, -1 <= r <= 1                     ---*/
  /*--- V0 = coor - (X1+X0)/2, V1 = (X1-X0)/2                      ---*/
  su2double V0[3], V1[3];
  for (auto k = 0u; k < nDim; ++k) {
    V0[k] = coor[k] - 0.5 * (x1[k] + x0[k]);
    V1[k] = 0.5 * (x1[k] - x0[k]);
  }

  /*--- Determine the value of r for minimum distance ---*/
  su2double dotV0V1 = 0.0, dotV1V1 = 0.0;
  for (auto k = 0u; k < nDim; ++k) {
    dotV0V1 += V0[k] * V1[k];
    dotV1V1 += V1[k] * V1[k];
  }
  su2double r = dotV0V1 / dotV1V1;
  r = max(-1.0, min(1.0, r));

  /*--- Compute the nearest point using the parametric equation ---*/
  dist2Line = 0.0;
  for (auto k = 0u; k < nDim; ++k) {
    surfCoor[k] = x0[k] + 0.5 * (r + 1.0) * (x1[k] - x0[k]);

    /*--- Also compute minimum distance squared ---*/
    const su2double ds = V0[k] - r * V1[k];
    dist2Line += ds * ds;
  }
}

bool CVolumeInterpolator::NearestPointOnTriangle(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                                                 const unsigned long i2, const su2double* coor, su2double* surfCoor,
                                                 su2double& dist2Tria, const unsigned short nDim) {
  su2double x0[3], x1[3], x2[3];
  for (auto k = 0u; k < nDim; ++k) {
    x0[k] = geometry->nodes->GetCoord(i0, k);
    x1[k] = geometry->nodes->GetCoord(i1, k);
    x2[k] = geometry->nodes->GetCoord(i2, k);
  }

  /*--- Use the same parametrization as CADTElemClass::Dist2ToTriangle   ---*/
  /*--- X = X0 + (r+1)*(X1-X0)/2 + (s+1)*(X2-X0)/2, r, s >= -1, r+s <= 0 ---*/
  /*--- V0 = coor - (X1+X2)/2, V1 = (X1-X0)/2, V2 = (X2-X0)/2            ---*/
  su2double V0[3], V1[3], V2[3];
  for (auto k = 0u; k < nDim; ++k) {
    V0[k] = coor[k] - 0.5 * (x1[k] + x2[k]);
    V1[k] = 0.5 * (x1[k] - x0[k]);
    V2[k] = 0.5 * (x2[k] - x0[k]);
  }

  /*--- Compute dot products ---*/
  su2double dotV0V1 = 0.0, dotV0V2 = 0.0, dotV1V1 = 0.0, dotV1V2 = 0.0, dotV2V2 = 0.0;
  for (auto k = 0u; k < nDim; ++k) {
    dotV0V1 += V0[k] * V1[k];
    dotV0V2 += V0[k] * V2[k];
    dotV1V1 += V1[k] * V1[k];
    dotV1V2 += V1[k] * V2[k];
    dotV2V2 += V2[k] * V2[k];
  }

  /*--- Solve for parametric coordinates ---*/
  const su2double detInv = 1.0 / (dotV1V1 * dotV2V2 - dotV1V2 * dotV1V2);
  su2double r = detInv * (dotV0V1 * dotV2V2 - dotV0V2 * dotV1V2);
  su2double s = detInv * (dotV0V2 * dotV1V1 - dotV0V1 * dotV1V2);

  /*--- Check if projection is inside triangle ---*/
  const su2double tolInsideElem = 1.e-10;
  const su2double paramLowerBound = -1.0 - tolInsideElem;

  if ((r >= paramLowerBound) && (s >= paramLowerBound) && ((r + s) <= tolInsideElem)) {
    /*--- Projection is inside triangle, compute nearest point ---*/
    dist2Tria = 0.0;
    for (auto k = 0u; k < nDim; ++k) {
      surfCoor[k] = x0[k] + 0.5 * (r + 1.0) * (x1[k] - x0[k]) + 0.5 * (s + 1.0) * (x2[k] - x0[k]);

      /*--- Also compute minimum distance squared ---*/
      const su2double ds = V0[k] - r * V1[k] - s * V2[k];
      dist2Tria += ds * ds;
    }

    return true;
  }

  /*--- The projection of the coordinate is outside the triangle. ---*/
  /*--- Return false. ---*/
  return false;
}

bool CVolumeInterpolator::NearestPointOnQuadrilateral(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                                                      const unsigned long i2, const unsigned long i3, const su2double* coor,
                                                      su2double* surfCoor, su2double& dist2Quad, const unsigned short nDim) {
  su2double x0[3], x1[3], x2[3], x3[3];
  for (auto k = 0u; k < nDim; ++k) {
    x0[k] = geometry->nodes->GetCoord(i0, k);
    x1[k] = geometry->nodes->GetCoord(i1, k);
    x2[k] = geometry->nodes->GetCoord(i2, k);
    x3[k] = geometry->nodes->GetCoord(i3, k);
  }

  /*--- TODO: Use the same parametrization as CADTElem::Dist2ToQuadrilateral ---*/
  return false;
}

void CVolumeInterpolator::ApplyCurvatureCorrection(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                   const unsigned short nDim, vector<su2double>& coor_dst,
                                                   vector<su2double>& coor_corrected) {
  /*--- Initialize corrected coordinates to original coordinates ---*/
  coor_corrected = coor_dst;

  /*--- Get references to the ADTs for both source and destination grids ---*/
  CADTElemClass& srcSurfaceADT = GetSourceSurfaceADT();
  CADTElemClass& dstSurfaceADT = GetDestinationSurfaceADT();

  /*--- Apply curvature correction if both surfaces exist ---*/
  if (!srcSurfaceADT.IsEmpty() && !dstSurfaceADT.IsEmpty()) {
    const unsigned long nDOFs = coor_dst.size() / nDim;
    for (auto l = 0u; l < nDOFs; ++l) {
      /*--- Find nearest wall points ---*/
      unsigned short srcMarkerID, dstMarkerID;
      unsigned long srcElemID, dstElemID;
      int srcRankID, dstRankID;
      su2double srcDist, dstDist;
      su2double surfCoorSrc[3], surfCoorDst[3];

      su2double* coor = coor_corrected.data() + l * nDim;

      /*--- Find the closest point on the source surface mesh ---*/
      srcSurfaceADT.DetermineNearestElement(coor, srcDist, srcMarkerID, srcElemID, srcRankID);
      NearestPointOnElement(geometry_src, srcMarkerID, srcElemID, coor, surfCoorSrc,
                            srcDist, nDim);

      /*--- Find the closest point on the destination surface mesh to the source wall point ---*/
      dstSurfaceADT.DetermineNearestElement(surfCoorSrc, dstDist, dstMarkerID, dstElemID, dstRankID);
      NearestPointOnElement(geometry_dst, dstMarkerID, dstElemID, coor, surfCoorDst,
                            dstDist, nDim);

      /*--- Determine the curvature correction, which is the vector from the  ---*/
      /*--- wall coordinate of the output grid to the wall coordinates on the ---*/
      /*--- input grid                                                        ---*/
      su2double mag = 0.0;
      for (auto k = 0u; k < nDim; ++k) mag += pow(surfCoorSrc[k] - surfCoorDst[k], 2.0);
      if (sqrt(mag) > 1e-5) {
        for (auto k = 0u; k < nDim; ++k) {
          cout << "Correction[" << k << "]: " << surfCoorSrc[k] - surfCoorDst[k];
          if (k < nDim - 1) cout << "; ";
        }
        cout << endl;
      }
      for (auto k = 0u; k < nDim; ++k) coor[k] += surfCoorSrc[k] - surfCoorDst[k];
    }
  }
}

void CVolumeInterpolator::WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput* output,
                                     unsigned long TimeIter) {
  /*--- Load history data (volume output might require some values) --- */

  output->SetHistoryOutput(geometry, solver_container, config, TimeIter, 0, 0);

  /*--- Load the data --- */

  output->LoadData(geometry, config, solver_container);

  /*--- Set the filenames ---*/

  output->SetRestartFilename(config->GetRestart_FileName());

  output->SetVolumeFilename(config->GetVolume_FileName());

  output->SetSurfaceFilename(config->GetSurfCoeff_FileName());

  auto FileFormat = config->GetVolumeOutputFiles();
  for (auto iFile = 0u; iFile < config->GetnVolumeOutputFiles(); ++iFile) {
    // if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII && FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY &&
    //     FileFormat[iFile] != OUTPUT_TYPE::CSV)
    output->WriteToFile(config, geometry, FileFormat[iFile]);
  }
}