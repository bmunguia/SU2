/*!
 * \file SU2_ITP.cpp
 * \brief Main file for the solution interpolation code (SU2_ITP).
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

#include "../include/SU2_ITP.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  unsigned short iZone, iInst;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;

  char config_file_name[MAX_STRING_SIZE];

  /*--- MPI initialization ---*/

  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Comm MPICommunicator = SU2_MPI::GetComm();

  const int rank = SU2_MPI::GetRank();
  const int size = SU2_MPI::GetSize();

  /*--- Pointer to different structures that will be used throughout the entire code ---*/

  COutput** output = nullptr;
  CGeometry*** geometry_src = nullptr;
  CGeometry*** geometry_dst = nullptr;
  CSolver*** solver_src = nullptr;
  CSolver*** solver_dst = nullptr;
  CConfig** config_src = nullptr;
  CConfig** config_dst = nullptr;
  CConfig* driver_config = nullptr;
  unsigned short* nInst = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2 || argc == 3) {
    strcpy(config_file_name, argv[1]);
  } else {
    strcpy(config_file_name, "default.cfg");
  }

  CConfig* config = nullptr;
  config = new CConfig(config_file_name, SU2_COMPONENT::SU2_ITP);

  const auto nZone = config->GetnZone();

  /*--- Definition of the containers per zones ---*/

  solver_src = new CSolver**[nZone]();
  solver_dst = new CSolver**[nZone]();
  config_src = new CConfig*[nZone]();
  config_dst = new CConfig*[nZone]();
  geometry_src = new CGeometry**[nZone]();
  geometry_dst = new CGeometry**[nZone]();
  nInst = new unsigned short[nZone];
  driver_config = nullptr;
  output = new COutput*[nZone]();

  for (iZone = 0; iZone < nZone; iZone++) {
    nInst[iZone] = 1;
  }

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_ITP, false);

  /*--- Initialize a char to store the zone filename ---*/
  char zone_file_name[MAX_STRING_SIZE];

  /*--- Store a boolean for multizone problems ---*/
  const bool multizone = config->GetMultizone_Problem();

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    InitializeConfig(driver_config, config_src, zone_file_name, config_file_name, iZone, nZone, MPICommunicator);
    InitializeConfig(driver_config, config_dst, zone_file_name, config_file_name, iZone, nZone, MPICommunicator, false);
  }

  /*--- Set the multizone part of the problem. ---*/
  if (multizone) {
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_src[iZone]->SetMultizone(driver_config, config_src);
    }
  }

  /*--- Allocate geometries and solvers ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    geometry_src[iZone] = new CGeometry*[nInst[iZone]];
    geometry_dst[iZone] = new CGeometry*[nInst[iZone]];
    solver_src[iZone] = new CSolver*[nInst[iZone]];
    solver_dst[iZone] = new CSolver*[nInst[iZone]];
    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
      geometry_src[iZone][iInst] = nullptr;
      geometry_dst[iZone][iInst] = nullptr;
      solver_src[iZone][iInst] = nullptr;
      solver_dst[iZone][iInst] = nullptr;
    }
  }

  /*--- Read the geometry for each zone ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
      InitializeGeometry(config_src[iZone], geometry_src[iZone][iInst], iZone, iInst, nZone);
      InitializeGeometry(config_dst[iZone], geometry_dst[iZone][iInst], iZone, iInst, nZone);
    }
  }

  const bool fsi = config_src[ZONE_0]->GetFSI_Simulation();
  const bool fem_solver = config_src[ZONE_0]->GetFEMSolver();

  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/

  StartTime = SU2_MPI::Wtime();

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Solution Postprocessing -----------------------" << endl;

  /*--- TODO: multizone, FSI, harmonic balance, or FEM.---*/

  if (config_src[ZONE_0]->GetTime_Domain()) {
    /*--- Unsteady simulation: merge all unsteady time steps. First,
      find the frequency and total number of files to write. ---*/

    su2double Physical_dt, Physical_t;
    unsigned long TimeIter = 0;
    bool StopCalc = false;
    bool* SolutionInstantiated = new bool[nZone];

    for (iZone = 0; iZone < nZone; iZone++) SolutionInstantiated[iZone] = false;

    /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
    if (config_src[ZONE_0]->GetTime_Domain() && config_src[ZONE_0]->GetRestart())
      TimeIter = config_src[ZONE_0]->GetRestart_Iter();

    while (TimeIter < config_src[ZONE_0]->GetnTime_Iter()) {
      /*--- Check several conditions in order to merge the correct time step files. ---*/
      Physical_dt = config_src[ZONE_0]->GetTime_Step();
      Physical_t = (TimeIter + 1) * Physical_dt;
      if (Physical_t >= config_src[ZONE_0]->GetMax_Time()) StopCalc = true;

      if ((TimeIter + 1 == config_src[ZONE_0]->GetnTime_Iter()) ||
          (TimeIter % config_src[ZONE_0]->GetVolumeOutputFrequency(0) == 0) ||
          (StopCalc)) {
        /*--- Read in the restart file for this time step ---*/
        for (iZone = 0; iZone < nZone; iZone++) {
          /*--- Set the current iteration number in the config class. ---*/
          config_src[iZone]->SetTimeIter(TimeIter);
          config_dst[iZone]->SetTimeIter(TimeIter);

          /*--- Either instantiate the solution class or load a restart file. ---*/
          if (!SolutionInstantiated[iZone] &&
              (TimeIter == 0 || (config_src[ZONE_0]->GetRestart() &&
                                 ((long)TimeIter == SU2_TYPE::Int(config_src[ZONE_0]->GetRestart_Iter()) ||
                                  TimeIter % config_src[ZONE_0]->GetVolumeOutputFrequency(0) == 0 ||
                                  TimeIter + 1 == config_src[ZONE_0]->GetnTime_Iter())))) {
            /*--- Initialize the solution classes ---*/
            solver_src[iZone][INST_0] = new CBaselineSolver(geometry_src[iZone][INST_0], config_src[iZone]);
            solver_dst[iZone][INST_0] = new CBaselineSolver(geometry_dst[iZone][INST_0], config_dst[iZone]);

            /*--- Initialize and preprocess the output ---*/
            output[iZone] = new CBaselineOutput(config_dst[iZone], geometry_dst[iZone][INST_0]->GetnDim(),
                                                solver_dst[iZone][INST_0]);
            output[iZone]->PreprocessVolumeOutput(config_dst[iZone]);
            output[iZone]->PreprocessHistoryOutput(config_dst[iZone], false);

            SolutionInstantiated[iZone] = true;
          }

          /*--- Load the solution on the source mesh ---*/
          config_src[iZone]->SetiInst(INST_0);
          config_dst[iZone]->SetiInst(INST_0);
          solver_src[iZone][INST_0]->LoadRestart(geometry_src[iZone], &solver_src[iZone], config_src[iZone], TimeIter,
                                                 true);

          /*--- Interpolate the solution ---*/
          InterpolateSolution(config_src[iZone], geometry_src[iZone][INST_0], geometry_dst[iZone][INST_0],
                              solver_src[iZone][INST_0], solver_dst[iZone][INST_0]);
        }

        if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << "." << endl;

        for (iZone = 0; iZone < nZone; iZone++) {
          WriteFiles(config_dst[iZone], geometry_dst[iZone][INST_0], &solver_dst[iZone][INST_0], output[iZone],
                     TimeIter);
        }
      }

      TimeIter++;
      if (StopCalc) break;
    }

  }

  else {
    /*--- Steady simulation: merge the single solution file. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      config_src[iZone]->SetiInst(INST_0);
      config_dst[iZone]->SetiInst(INST_0);

      /*--- Initialize the solution classes ---*/
      solver_src[iZone][INST_0] = new CBaselineSolver(geometry_src[iZone][INST_0], config_src[iZone]);
      solver_dst[iZone][INST_0] = new CBaselineSolver(geometry_dst[iZone][INST_0], config_dst[iZone]);

      /*--- Initialize and preprocess the output ---*/
      output[iZone] =
          new CBaselineOutput(config_dst[iZone], geometry_dst[iZone][INST_0]->GetnDim(), solver_dst[iZone][INST_0]);
      output[iZone]->PreprocessVolumeOutput(config_dst[iZone]);
      output[iZone]->PreprocessHistoryOutput(config_dst[iZone], false);

      /*--- Load the solution on the source mesh ---*/
      solver_src[iZone][INST_0]->LoadRestart(geometry_src[iZone], &solver_src[iZone], config_src[iZone], 0, true);

      /*--- Interpolate the solution ---*/
      InterpolateSolution(config_src[iZone], geometry_src[iZone][INST_0], geometry_dst[iZone][INST_0],
                          solver_src[iZone][INST_0], solver_dst[iZone][INST_0]);
    }
    for (iZone = 0; iZone < nZone; iZone++) {
      WriteFiles(config_dst[iZone], geometry_dst[iZone][INST_0], &solver_dst[iZone][INST_0], output[iZone], 0);
    }
  }

  delete config;
  config = nullptr;

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Finalize Solver -------------------------" << endl;

  if (geometry_src != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (geometry_src[iZone][iInst] != nullptr) {
          delete geometry_src[iZone][iInst];
        }
      }
      if (geometry_src[iZone] != nullptr) delete[] geometry_src[iZone];
    }
    delete[] geometry_src;
  }
  if (geometry_dst != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (geometry_dst[iZone][iInst] != nullptr) {
          delete geometry_dst[iZone][iInst];
        }
      }
      if (geometry_dst[iZone] != nullptr) delete[] geometry_dst[iZone];
    }
    delete[] geometry_dst;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry containers." << endl;

  if (solver_src != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (solver_src[iZone][iInst] != nullptr) {
          delete solver_src[iZone][iInst];
        }
      }
      if (solver_src[iZone] != nullptr) delete[] solver_src[iZone];
    }
    delete[] solver_src;
  }
  if (solver_dst != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (solver_dst[iZone][iInst] != nullptr) {
          delete solver_dst[iZone][iInst];
        }
      }
      if (solver_dst[iZone] != nullptr) delete[] solver_dst[iZone];
    }
    delete[] solver_dst;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSolver containers." << endl;

  if (config_src != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_src[iZone] != nullptr) {
        delete config_src[iZone];
      }
    }
    delete[] config_src;
  }
  if (config_dst != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_dst[iZone] != nullptr) {
        delete config_dst[iZone];
      }
    }
    delete[] config_dst;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig containers." << endl;

  if (output != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (output[iZone] != nullptr) {
        delete output[iZone];
      }
    }
    delete[] output;
  }
  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime - StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on " << size;
    if (size == 1)
      cout << " core." << endl;
    else
      cout << " cores." << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Exit Success (SU2_ITP) ------------------------" << endl << endl;

  /*--- Finalize MPI parallelization ---*/
  SU2_MPI::Finalize();

  return EXIT_SUCCESS;
}

void InitializeConfig(CConfig* driver_config, CConfig** config_container, char* zone_file_name, char* config_file_name,
                      int iZone, int nZone, SU2_MPI::Comm MPICommunicator, bool isSource) {
  if (driver_config->GetnConfigFiles() > 0) {
    strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
    config_container[iZone] = new CConfig(driver_config, zone_file_name, SU2_COMPONENT::SU2_ITP, iZone, nZone, true);
  } else {
    config_container[iZone] = new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_ITP, iZone, nZone, true);
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

void InitializeGeometry(CConfig* config, CGeometry*& geometry, int iZone, int iInst, int nZone) {
  int rank = SU2_MPI::GetRank();

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
}

std::unique_ptr<CADTElemClass> BuildSurfaceADT(const CConfig* config, CGeometry* geometry) {
  /*--- Build surface ADT for wall distance computation or general surface interpolation ---*/
  const unsigned short nDim = geometry->GetnDim();

  /*--- Initialize an array for the mesh points mapping ---*/
  vector<unsigned long> meshToSurface(geometry->GetnPoint(), 0);

  /*--- Define vectors for surface connectivity and metadata ---*/
  vector<unsigned long> surfaceConn;
  vector<unsigned long> elemIDs;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned short> markerIDs;

  /*--- Loop over boundary markers ---*/
  for (unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    /*--- Loop over surface elements of this marker ---*/
    for (unsigned long iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
      unsigned short VTK_Type = geometry->bound[iMarker][iElem]->GetVTK_Type();
      unsigned short nNodes = geometry->bound[iMarker][iElem]->GetnNodes();

      /*--- Set flag for mesh points on this surface ---*/
      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
        meshToSurface[geometry->bound[iMarker][iElem]->GetNode(iNode)] = 1;
      }

      /*--- Store connectivity and metadata ---*/
      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
        surfaceConn.push_back(geometry->bound[iMarker][iElem]->GetNode(iNode));
      }

      markerIDs.push_back(iMarker);
      VTK_TypeElem.push_back(VTK_Type);
      elemIDs.push_back(iElem);
    }
  }

  /*--- Create coordinates of points on surfaces ---*/
  vector<su2double> surfaceCoor;
  unsigned long nVertex_Surface = 0;

  for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (meshToSurface[iPoint]) {
      meshToSurface[iPoint] = nVertex_Surface++;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        surfaceCoor.push_back(geometry->nodes->GetCoord(iPoint, iDim));
      }
    }
  }

  /*--- Change surface connectivity to correspond to surfaceCoor indices ---*/
  for (unsigned long i = 0; i < surfaceConn.size(); i++) {
    surfaceConn[i] = meshToSurface[surfaceConn[i]];
  }

  /*--- Build and return the surface ADT ---*/
  return std::unique_ptr<CADTElemClass>(
      new CADTElemClass(nDim, surfaceCoor, surfaceConn, VTK_TypeElem, markerIDs, elemIDs, true));
}

std::unique_ptr<CADTElemClass> BuildVolumeADT(CGeometry* geometry) {
  /*--- Build volume ADT for element searching ---*/
  const unsigned short nDim = geometry->GetnDim();
  const unsigned long nPoint = geometry->GetnPoint();
  const unsigned long nElem = geometry->GetnElem();

  const unsigned short maxNodePerElem = (nDim == 2) ? 4 : 8;

  /*--- Prepare coordinate and connectivity data for ADT ---*/
  vector<su2double> volCoor;
  vector<unsigned long> elemConn;
  vector<unsigned short> vtkType;
  vector<unsigned short> subElemID;
  vector<unsigned long> parElemID;

  /*--- Copy coordinates ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      volCoor.push_back(geometry->nodes->GetCoord(iPoint, iDim));
    }
  }

  /*--- Copy element connectivity and metadata ---*/
  for (unsigned long iElem = 0; iElem < nElem; iElem++) {
    unsigned short VTK_Type = geometry->elem[iElem]->GetVTK_Type();
    unsigned short nNodes = geometry->elem[iElem]->GetnNodes();

    for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
      elemConn.push_back(geometry->elem[iElem]->GetNode(iNode));
    }

    vtkType.push_back(VTK_Type);
    subElemID.push_back(0);  // TODO: FEM subelements
    parElemID.push_back(iElem);
  }

  /*--- Create and return the volume ADT ---*/
  return std::unique_ptr<CADTElemClass>(
      new CADTElemClass(nDim, volCoor, elemConn, vtkType, subElemID, parElemID, false));
}

void InterpolateSolution(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst, CSolver* solver_src,
                         CSolver* solver_dst) {
  const int rank = SU2_MPI::GetRank();

  if (rank == MASTER_NODE) {
    cout << endl << "----------------------------- Interpolation -----------------------------" << endl;
    cout << "Interpolating solution from source mesh to destination mesh..." << endl;
  }

  const unsigned short nDim = geometry_src->GetnDim();
  const unsigned long nPoint_src = geometry_src->GetnPoint();
  const unsigned long nPoint_dst = geometry_dst->GetnPoint();
  const unsigned long nElem_src = geometry_src->GetnElem();

  if (rank == MASTER_NODE) {
    cout << "Source mesh: " << geometry_src->GetGlobal_nPointDomain() << " points, ";
    cout << geometry_src->GetGlobal_nElemDomain() << " elements" << endl;
    cout << "Destination mesh: " << geometry_dst->GetGlobal_nPointDomain() << " points" << endl;
  }

  /*--- Apply the curvature correction ---*/
  if (rank == MASTER_NODE) cout << "Applying curvature correction." << endl;
  vector<su2double> coorDst;
  vector<su2double> coorDstCorrected;
  for (unsigned long iPoint = 0; iPoint < nPoint_dst; iPoint++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      coorDst.push_back(geometry_dst->nodes->GetCoord(iPoint, iDim));
    }
  }
  ApplyCurvatureCorrection(config, geometry_src, geometry_dst, nDim, coorDst, coorDstCorrected);

  /*--- Volume interpolation ---*/
  if (rank == MASTER_NODE) cout << "Performing volume interpolation." << endl;
  vector<unsigned long> pointsFailed;
  VolumeInterpolationSolution(geometry_src, solver_src, solver_dst, coorDstCorrected, pointsFailed);

  /*--- Carry out a surface interpolation, via a minimum distance search,       */
  /*    for the points that could not be interpolated via the regular volume    */
  /*    interpolation. Print a warning about this.                           ---*/
  if (pointsFailed.size()) {
    if (rank == MASTER_NODE) {
      cout << pointsFailed.size() << " DOFs for which the containment search failed." << endl;
      cout << "A minimum distance search to the boundary of the domain is used for these points. " << endl;
    }
    unsigned long nPointsBeforeSurface = pointsFailed.size();
    SurfaceInterpolationSolution(geometry_src, solver_src, solver_dst, coorDst, pointsFailed);
  }
}

void VolumeInterpolationSolution(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                                 const vector<su2double>& coor_corrected, vector<unsigned long>& pointsFailed) {
  const int rank = SU2_MPI::GetRank();

  /*--- Step 1: Build the volume ADT for element searching ---*/
  if (rank == MASTER_NODE) cout << "Building volume ADT." << flush;

  unique_ptr<CADTElemClass> volumeADT_ptr = BuildVolumeADT(geometry_src);
  CADTElemClass& volumeADT = *volumeADT_ptr;

  if (rank == MASTER_NODE) cout << " Done." << endl;

  /*--- Step 2: Search for donor elements for the given coordinates ---*/
  const unsigned short nDim = geometry_src->GetnDim();
  const unsigned long nDOFsDst = coor_corrected.size() / nDim;
  const unsigned short nVar = solver_src->GetnVar();

  /*--- Loop over the DOFs to be interpolated ---*/
  pointsFailed.clear();
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
      /*--- Get element information ---*/
      unsigned short nNodes = geometry_src->elem[elemID]->GetnNodes();

      /*--- Initialize interpolated solution to zero ---*/
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        solver_dst->GetNodes()->SetSolution(l, iVar, 0.0);
      }

      /*--- Interpolate using shape function weights ---*/
      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
        unsigned long nodeID = geometry_src->elem[elemID]->GetNode(iNode);

        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
          su2double solValue = solver_src->GetNodes()->GetSolution(nodeID, iVar);
          su2double currentValue = solver_dst->GetNodes()->GetSolution(l, iVar);
          solver_dst->GetNodes()->SetSolution(l, iVar, currentValue + weightsInterpol[iNode] * solValue);
        }
      }
    } else {
      /*--- Containment search failed - store the index ---*/
      pointsFailed.push_back(l);
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Volume search finished. " << pointsFailed.size() << " points failed." << endl << flush;
  }
}

void SurfaceInterpolationSolution(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                                  const vector<su2double> &coor_dst, vector<unsigned long> &pointsFailed) {
  const int rank = SU2_MPI::GetRank();

  if (pointsFailed.empty()) return;

  /*--- Step 1: Build the surface ADT for minimum distance search ---*/
  if (rank == MASTER_NODE) {
    cout << "Building surface ADT for minimum distance search." << flush;
  }

  const unsigned short nDim = geometry_src->GetnDim();
  const unsigned short nVar = solver_src->GetnVar();

  /*--- Build surface ADT using existing function ---*/
  std::unique_ptr<CADTElemClass> surfaceADT_ptr = BuildSurfaceADT(nullptr, geometry_src);
  CADTElemClass& surfaceADT = *surfaceADT_ptr;

  /*--- Check if surface ADT was built successfully ---*/
  if (!surfaceADT.IsEmpty()) {

    if (rank == MASTER_NODE) cout << " Done." << endl;

    /*--- Step 2: Search for donor elements for failed points ---*/
    if (rank == MASTER_NODE) {
      cout << "Performing minimum distance search for " << pointsFailed.size()
           << " failed points." << endl << flush;
    }

    unsigned long nExtrapolated = 0;

    /*--- Loop over failed points for minimum distance search ---*/
    for (unsigned long l = 0; l < pointsFailed.size(); ++l) {
      /*--- Get coordinates of failed point ---*/
      const unsigned long pointIndex = pointsFailed[l];
      const su2double* coor = coor_dst.data() + pointIndex * nDim;

      /*--- Find nearest surface element ---*/
      unsigned short nearestMarker;
      unsigned long nearestElemID;
      int rankID;
      su2double distance;

      surfaceADT.DetermineNearestElement(coor, distance, nearestMarker, nearestElemID, rankID);

      /*--- Compute nearest point on surface element ---*/
      su2double wallCoor[3];
      ComputeNearestPointOnSurfaceElement(geometry_src, nearestMarker, nearestElemID,
                                         coor, wallCoor, nDim);

      /*--- Use nearest surface element nodes for interpolation ---*/
      /*--- This is a simplified approach - could be improved with better surface projection ---*/
      unsigned short nNodes = geometry_src->bound[nearestMarker][nearestElemID]->GetnNodes();

      /*--- Simple inverse distance weighting from surface element nodes ---*/
      vector<su2double> weights(nNodes, 0.0);
      su2double totalWeight = 0.0;

      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
        unsigned long nodeID = geometry_src->bound[nearestMarker][nearestElemID]->GetNode(iNode);

        /*--- Compute distance from interpolation point to node ---*/
        su2double dist2 = 0.0;
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          su2double diff = coor[iDim] - geometry_src->nodes->GetCoord(nodeID, iDim);
          dist2 += diff * diff;
        }

        /*--- Inverse distance weighting (with small epsilon to avoid division by zero) ---*/
        weights[iNode] = 1.0 / (sqrt(dist2) + 1e-12);
        totalWeight += weights[iNode];
      }

      /*--- Normalize weights ---*/
      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
        weights[iNode] /= totalWeight;
      }

      /*--- Interpolate solution using weighted average ---*/
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        su2double interpolatedValue = 0.0;

        for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
          unsigned long nodeID = geometry_src->bound[nearestMarker][nearestElemID]->GetNode(iNode);
          su2double solValue = solver_src->GetNodes()->GetSolution(nodeID, iVar);
          interpolatedValue += weights[iNode] * solValue;
        }

        solver_dst->GetNodes()->SetSolution(pointIndex, iVar, interpolatedValue);
      }

      nExtrapolated++;
    }

    if (rank == MASTER_NODE) {
      cout << "Surface search finished. " << nExtrapolated << " points extrapolated."
           << endl << flush;
    }

  } else {
    /*--- No surface elements found ---*/
    if (rank == MASTER_NODE) {
      cout << " No surface elements found for minimum distance search." << endl;
    }
  }
}

void ComputeNearestPointOnSurfaceElement(CGeometry* geometry, unsigned short markerID,
                                         unsigned long elemID, const su2double* coor,
                                         su2double* nearestPoint, const unsigned short nDim) {
  unsigned short VTK_Type = geometry->bound[markerID][elemID]->GetVTK_Type();
  unsigned short nNodes = geometry->bound[markerID][elemID]->GetnNodes();

  switch(VTK_Type) {
    case TRIANGLE: {
      su2double x0[3], x1[3], x2[3];
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        x0[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(0), iDim);
        x1[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(1), iDim);
        x2[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(2), iDim);
      }

      /*--- Compute nearest point on triangle using parametric approach ---*/
      /*--- Based on CADTElemClass::Dist2ToTriangle logic ---*/
      su2double V0[3], V1[3], V2[3];
      for (unsigned short k = 0; k < nDim; ++k) {
        V0[k] = coor[k] - 0.5 * (x1[k] + x2[k]);
        V1[k] = 0.5 * (x1[k] - x0[k]);
        V2[k] = 0.5 * (x2[k] - x0[k]);
      }

      /*--- Compute dot products ---*/
      su2double dotV0V1 = 0.0, dotV0V2 = 0.0, dotV1V1 = 0.0, dotV1V2 = 0.0, dotV2V2 = 0.0;
      for (unsigned short k = 0; k < nDim; ++k) {
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
        for (unsigned short k = 0; k < nDim; ++k) {
          nearestPoint[k] = 0.5 * (x1[k] + x2[k]) + r * V1[k] + s * V2[k];
        }
      } else {
        /*--- Projection is outside triangle, find nearest point on edges ---*/
        su2double minDist2 = 1e30;
        su2double bestPoint[3];

        /*--- Check each edge ---*/
        for (int iEdge = 0; iEdge < 3; ++iEdge) {
          su2double* p0 = (iEdge == 0) ? x0 : (iEdge == 1) ? x1 : x2;
          su2double* p1 = (iEdge == 0) ? x1 : (iEdge == 1) ? x2 : x0;

          /*--- Find nearest point on this edge ---*/
          su2double edge[3], toCoor[3];
          su2double edgeLength2 = 0.0, dotProduct = 0.0;

          for (unsigned short k = 0; k < nDim; ++k) {
            edge[k] = p1[k] - p0[k];
            toCoor[k] = coor[k] - p0[k];
            edgeLength2 += edge[k] * edge[k];
            dotProduct += edge[k] * toCoor[k];
          }

          su2double t = (edgeLength2 > 1e-12) ? dotProduct / edgeLength2 : 0.0;
          t = max(0.0, min(1.0, t)); // Clamp to [0,1]

          su2double edgePoint[3], dist2 = 0.0;
          for (unsigned short k = 0; k < nDim; ++k) {
            edgePoint[k] = p0[k] + t * edge[k];
            su2double diff = coor[k] - edgePoint[k];
            dist2 += diff * diff;
          }

          if (dist2 < minDist2) {
            minDist2 = dist2;
            for (unsigned short k = 0; k < nDim; ++k) {
              bestPoint[k] = edgePoint[k];
            }
          }
        }

        for (unsigned short k = 0; k < nDim; ++k) {
          nearestPoint[k] = bestPoint[k];
        }
      }
      break;
    }

    case QUADRILATERAL: {
      su2double x0[3], x1[3], x2[3], x3[3];
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        x0[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(0), iDim);
        x1[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(1), iDim);
        x2[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(2), iDim);
        x3[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(3), iDim);
      }

      /*--- Project to quad using two triangles (0-1-2 and 0-2-3), then check edges if outside ---*/
      auto nearest_on_triangle = [&](const su2double* a, const su2double* b, const su2double* c, const su2double* p, su2double* out) {
        /*--- Barycentric projection ---*/
        su2double ab[3], ac[3], ap[3];
        for (unsigned short k = 0; k < nDim; ++k) {
          ab[k] = b[k] - a[k];
          ac[k] = c[k] - a[k];
          ap[k] = p[k] - a[k];
        }
        su2double d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0;
        for (unsigned short k = 0; k < nDim; ++k) {
          d1 += ab[k] * ab[k];
          d2 += ab[k] * ac[k];
          d3 += ac[k] * ac[k];
          d4 += ab[k] * ap[k];
          d5 += ac[k] * ap[k];
        }
        su2double denom = d1 * d3 - d2 * d2;
        su2double v = (d3 * d4 - d2 * d5) / denom;
        su2double w = (d1 * d5 - d2 * d4) / denom;
        su2double u = 1.0 - v - w;
        /*--- Clamp to triangle ---*/
        if (u >= 0 && v >= 0 && w >= 0) {
          for (unsigned short k = 0; k < nDim; ++k)
            out[k] = u * a[k] + v * b[k] + w * c[k];
          return true;
        }
        return false;
      };

      su2double tri0[3], tri1[3];
      bool inTri0 = nearest_on_triangle(x0, x1, x2, coor, tri0);
      bool inTri1 = nearest_on_triangle(x0, x2, x3, coor, tri1);
      su2double dist0 = 1e30, dist1 = 1e30;
      if (inTri0) {
        dist0 = 0.0;
        for (unsigned short k = 0; k < nDim; ++k) {
          su2double d = coor[k] - tri0[k];
          dist0 += d * d;
        }
      }
      if (inTri1) {
        dist1 = 0.0;
        for (unsigned short k = 0; k < nDim; ++k) {
          su2double d = coor[k] - tri1[k];
          dist1 += d * d;
        }
      }
      if (inTri0 || inTri1) {
        if (dist0 < dist1) {
          for (unsigned short k = 0; k < nDim; ++k) nearestPoint[k] = tri0[k];
        } else {
          for (unsigned short k = 0; k < nDim; ++k) nearestPoint[k] = tri1[k];
        }
      } else {
        // If not inside, check all 4 edges
        su2double minDist2 = 1e30, best[3];
        const su2double* pts[4] = {x0, x1, x2, x3};
        for (int i = 0; i < 4; ++i) {
          const su2double* p0 = pts[i];
          const su2double* p1 = pts[(i+1)%4];
          su2double edge[3], toCoor[3];
          su2double edgeLength2 = 0.0, dotProduct = 0.0;
          for (unsigned short k = 0; k < nDim; ++k) {
            edge[k] = p1[k] - p0[k];
            toCoor[k] = coor[k] - p0[k];
            edgeLength2 += edge[k] * edge[k];
            dotProduct += edge[k] * toCoor[k];
          }
          su2double t = (edgeLength2 > 1e-12) ? dotProduct / edgeLength2 : 0.0;
          t = max(0.0, min(1.0, t));
          su2double edgePoint[3], dist2 = 0.0;
          for (unsigned short k = 0; k < nDim; ++k) {
            edgePoint[k] = p0[k] + t * edge[k];
            su2double diff = coor[k] - edgePoint[k];
            dist2 += diff * diff;
          }
          if (dist2 < minDist2) {
            minDist2 = dist2;
            for (unsigned short k = 0; k < nDim; ++k) best[k] = edgePoint[k];
          }
        }
        for (unsigned short k = 0; k < nDim; ++k) nearestPoint[k] = best[k];
      }
      break;
    }

    case LINE: {
      su2double p0[3], p1[3];
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        p0[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(0), iDim);
        p1[iDim] = geometry->nodes->GetCoord(geometry->bound[markerID][elemID]->GetNode(1), iDim);
      }
      su2double edge[3], toCoor[3];
      su2double edgeLength2 = 0.0, dotProduct = 0.0;
      for (unsigned short k = 0; k < nDim; ++k) {
        edge[k] = p1[k] - p0[k];
        toCoor[k] = coor[k] - p0[k];
        edgeLength2 += edge[k] * edge[k];
        dotProduct += edge[k] * toCoor[k];
      }
      su2double t = (edgeLength2 > 1e-12) ? dotProduct / edgeLength2 : 0.0;
      t = max(0.0, min(1.0, t));
      for (unsigned short k = 0; k < nDim; ++k) {
        nearestPoint[k] = p0[k] + t * edge[k];
      }
      break;
    }

    default:
      /* This should not happen. */
      SU2_MPI::Error("Element type not recognized.", CURRENT_FUNCTION);

  }
}

void ApplyCurvatureCorrection(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                              const unsigned short nDim, const vector<su2double>& coor_dst,
                              vector<su2double>& coor_corrected) {
  /*--- Initialize corrected coordinates to original coordinates ---*/
  coor_corrected = coor_dst;

  /*--- Build surface ADTs for both source and destination grids ---*/
  std::unique_ptr<CADTElemClass> srcSurfaceADT_ptr = BuildSurfaceADT(config, geometry_src);
  std::unique_ptr<CADTElemClass> dstSurfaceADT_ptr = BuildSurfaceADT(config, geometry_dst);

  CADTElemClass& srcSurfaceADT = *srcSurfaceADT_ptr;
  CADTElemClass& dstSurfaceADT = *dstSurfaceADT_ptr;

  /*--- Apply curvature correction if both surfaces exist ---*/
  if (!srcSurfaceADT.IsEmpty() && !dstSurfaceADT.IsEmpty()) {
    const unsigned long nDOFs = coor_dst.size() / nDim;
    for (unsigned long l = 0; l < nDOFs; ++l) {
      /*--- Find nearest wall points ---*/
      unsigned short srcMarkerID, dstMarkerID;
      unsigned long srcElemID, dstElemID;
      int srcRankID, dstRankID;
      su2double srcDist, dstDist;
      su2double wallCoorSrc[3], wallCoorDst[3];

      su2double* coor = coor_corrected.data() + l * nDim;

      /*--- Find the closest point on the source surface element ---*/
      srcSurfaceADT.DetermineNearestElement(coor, srcDist, srcMarkerID, srcElemID, srcRankID);
      ComputeNearestPointOnSurfaceElement(geometry_src, srcMarkerID, srcElemID, coor, wallCoorSrc, nDim);

      /*--- Find the closest point on the destination surface element ---*/
      dstSurfaceADT.DetermineNearestElement(coor, dstDist, dstMarkerID, dstElemID, dstRankID);
      ComputeNearestPointOnSurfaceElement(geometry_dst, dstMarkerID, dstElemID, coor, wallCoorDst, nDim);

      /*--- Determine the curvature correction, which is the vector from the     */
      /*    wall coordinate of the output grid to the wall coordinates on the    */
      /*    input grid                                                        ---*/
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) coor[iDim] += wallCoorSrc[iDim] - wallCoorDst[iDim];
    }
  }
}

void WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput* output,
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
  for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++) {
    // if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII && FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY &&
    //     FileFormat[iFile] != OUTPUT_TYPE::CSV)
    output->WriteToFile(config, geometry, FileFormat[iFile]);
  }
}
