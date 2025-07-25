/*!
 * \file SU2_ITP.cpp
 * \brief Main file for the solution interpolation code (SU2_ITP).
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
      InitializeGeometry(config_src[iZone], geometry_src[iZone][iInst], output, iZone, iInst);
      InitializeGeometry(config_dst[iZone], geometry_dst[iZone][iInst], output, iZone, iInst);
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
          ((TimeIter % config_src[ZONE_0]->GetVolumeOutputFrequency(0) == 0) && (TimeIter != 0) &&
            (config_src[ZONE_0]->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
            (config_src[ZONE_0]->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND)) ||
          (StopCalc) ||
          (((config_src[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
            (config_src[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) &&
            ((TimeIter == 0) || (TimeIter % config_src[ZONE_0]->GetVolumeOutputFrequency(0) == 0)))) {
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
          solver_src[iZone][INST_0]->LoadRestart(geometry_src[iZone], &solver_src[iZone],
                                                 config_src[iZone], TimeIter, true);

          /*--- Interpolate the solution ---*/
        }

        if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << "." << endl;

        for (iZone = 0; iZone < nZone; iZone++) {
          WriteFiles(config_dst[iZone], geometry_dst[iZone][INST_0], &solver_dst[iZone][INST_0],
                     output[iZone], TimeIter);
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
      output[iZone] = new CBaselineOutput(config_dst[iZone], geometry_dst[iZone][INST_0]->GetnDim(),
                                          solver_dst[iZone][INST_0]);
      output[iZone]->PreprocessVolumeOutput(config_dst[iZone]);
      output[iZone]->PreprocessHistoryOutput(config_dst[iZone], false);

      /*--- Load the solution on the source mesh ---*/
      solver_src[iZone][INST_0]->LoadRestart(geometry_src[iZone], &solver_src[iZone],
                                             config_src[iZone], 0, true);

      /*--- Interpolate the solution ---*/
    }
    for (iZone = 0; iZone < nZone; iZone++) {
      WriteFiles(config_dst[iZone], geometry_dst[iZone][INST_0], &solver_dst[iZone][INST_0],
                  output[iZone], 0);
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

void InitializeConfig(CConfig* driver_config, CConfig** config_container, char* zone_file_name, const char* config_file_name,
                      int iZone, int nZone, SU2_MPI::Comm MPICommunicator, bool isSource) {
  if (driver_config->GetnConfigFiles() > 0) {
    strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
    config_container[iZone] = new CConfig(driver_config, zone_file_name, SU2_COMPONENT::SU2_ITP, iZone, nZone, true);
  } else {
    config_container[iZone] =
        new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_ITP, iZone, nZone, true);
  }
  config_container[iZone]->SetMPICommunicator(MPICommunicator);

  /*--- Set the multizone part of the problem. ---*/
  const bool multizone = config_container[iZone]->GetMultizone_Problem();
  if (multizone) {
    /*--- Set the interface markers for multizone ---*/
    config_container[iZone]->SetMultizone(driver_config, config_src);
  }

  /*--- Read MESH_FILENAME if source, otherwise read destination mesh ---*/
  if (!isSource) {
    /*--- For now, just use MESH_OUT_FILENAME ---*/
    config_container[iZone]->SetMesh_FileName(config_container[iZone]->GetMesh_Out_FileName());
  }
}

void InitializeGeometry(CConfig* config, CGeometry*& geometry, COutput* output,
                        int iZone, int iInst) {
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

  /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/
  if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
  geometry->SetRCM_Ordering(config);

  /*--- Recompute elements surrounding points, points surrounding points ---*/
  if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
  geometry->SetPoint_Connectivity();

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
  const bool fem_solver = config[ZONE_0]->GetFEMSolver();
  if (!fem_solver)
    geometry->PreprocessP2PComms(geometry, config);
  if (fem_solver) {
    auto* DGMesh = dynamic_cast<CMeshFEM_DG*>(geometry);
    if (rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;
    DGMesh->CreateStandardVolumeElements(config);
    if (rank == MASTER_NODE) cout << "Creating face information." << endl;
    DGMesh->CreateFaces(config);
  }
}

void WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput* output,
                unsigned long TimeIter) {
  /*--- Load history data (volume output might require some values) --- */

  output->SetHistoryOutput(geometry, solver_container, config, TimeIter, 0, 0);

  /*--- Load the data --- */

  output->LoadData(geometry, config, solver_container);

  /*--- Set the filenames ---*/

  output->SetVolumeFilename(config->GetVolume_FileName());

  output->SetSurfaceFilename(config->GetSurfCoeff_FileName());

  for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++) {
    auto FileFormat = config->GetVolumeOutputFiles();
    if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII && FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY &&
        FileFormat[iFile] != OUTPUT_TYPE::CSV)
      output->WriteToFile(config, geometry, FileFormat[iFile]);
  }
}
