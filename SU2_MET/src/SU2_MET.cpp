/*!
 * \file SU2_MET.cpp
 * \brief Main file for the solution metric calculation code (SU2_MET).
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

#include "../include/SU2_MET.hpp"

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

  CConfig* driver_config = nullptr;
  CConfig** config = nullptr;
  CGeometry*** geometry = nullptr;
  CSolver**** solver = nullptr;
  COutput** output = nullptr;
  unsigned short* nInst = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2 || argc == 3) {
    strcpy(config_file_name, argv[1]);
  } else {
    strcpy(config_file_name, "default.cfg");
  }

  /*--- Initialize the main configuration ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_MET);

  const auto nZone = driver_config->GetnZone();

  /*--- Definition of the containers per zones ---*/
  config = new CConfig*[nZone]();
  geometry = new CGeometry**[nZone]();
  solver = new CSolver***[nZone]();
  output = new COutput*[nZone]();
  nInst = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    nInst[iZone] = 1;
  }

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_MET, false);

  /*--- Initialize a char to store the zone filename ---*/
  char zone_file_name[MAX_STRING_SIZE];

  /*--- Store a boolean for multizone problems ---*/
  const bool multizone = driver_config->GetMultizone_Problem();

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    InitializeConfig(driver_config, config, zone_file_name,
                     config_file_name, SU2_COMPONENT::SU2_MET, iZone, nZone,
                     MPICommunicator);
  }

  /*--- Set the multizone part of the problem. ---*/
  if (multizone) {
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config[iZone]->SetMultizone(driver_config, config);
    }
  }

  /*--- Allocate geometries and solvers ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    geometry[iZone] = new CGeometry*[nInst[iZone]]();
    solver[iZone] = new CSolver**[nInst[iZone]]();
  }

  /*--- Read the geometry for each zone ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
      InitializeGeometry(config[iZone], geometry[iZone][iInst], iZone, iInst, nZone);
    }
  }

  const bool fsi = config[ZONE_0]->GetFSI_Simulation();
  const bool fem_solver = config[ZONE_0]->GetFEMSolver();

  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/

  StartTime = SU2_MPI::Wtime();

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Solution Postprocessing -----------------------" << endl;

  /*--- TODO: multizone, FSI, harmonic balance, or FEM.---*/

  if (config[ZONE_0]->GetTime_Domain()) {
    /*--- Unsteady simulation: merge all unsteady time steps. First,
      find the frequency and total number of files to write. ---*/

    su2double Physical_dt, Physical_t;
    unsigned long TimeIter = 0;
    const bool dual_time_2nd = (config[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
    bool StopCalc = false;
    bool* SolutionInstantiated = new bool[nZone];

    for (iZone = 0; iZone < nZone; iZone++) SolutionInstantiated[iZone] = false;

    /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
    if (config[ZONE_0]->GetTime_Domain() && config[ZONE_0]->GetRestart())
      TimeIter = config[ZONE_0]->GetRestart_Iter();

    while (TimeIter < config[ZONE_0]->GetnTime_Iter()) {
      /*--- Check several conditions in order to merge the correct time step files. ---*/
      Physical_dt = config[ZONE_0]->GetTime_Step();
      Physical_t = (TimeIter + 1) * Physical_dt;
      if (Physical_t >= config[ZONE_0]->GetMax_Time()) StopCalc = true;

      const bool IsTime0 = (TimeIter == 0);
      const bool IsTimeWrt = (TimeIter % config[ZONE_0]->GetVolumeOutputFrequency(0) == 0);
      const bool IsTimeEnd = (TimeIter + 1 == config[ZONE_0]->GetnTime_Iter()) ||
                             (TimeIter + 2 == config[ZONE_0]->GetnTime_Iter() && dual_time_2nd);
      const bool IsTimeRestart = ((long)TimeIter == SU2_TYPE::Int(config[ZONE_0]->GetRestart_Iter()));

      if (StopCalc || IsTime0 || IsTimeWrt || IsTimeEnd || IsTimeRestart) {
        /*--- Read in the restart file for this time step ---*/
        for (iZone = 0; iZone < nZone; iZone++) {
          /*--- Set the current iteration number in the config class. ---*/
          config[iZone]->SetTimeIter(TimeIter);

          /*--- Only implemented for single-instance problems ---*/
          config[iZone]->SetiInst(INST_0);

          /*--- Either instantiate the solution class or load a restart file. ---*/
          bool initialInterp = false;
          if (!SolutionInstantiated[iZone]) {
            /*--- Initialize the solution classes ---*/
            InitializeSolver(config[iZone], geometry[iZone][INST_0], solver[iZone][INST_0],
                             iZone, INST_0, nZone);

            /*--- Initialize and preprocess the output ---*/
            InitializeOutput(config[iZone], geometry[iZone][INST_0], solver[iZone][INST_0],
                             output[iZone], iZone, INST_0, nZone);

            SolutionInstantiated[iZone] = true;
            initialInterp = true;

            /*--- Calculate the surface metric ---*/
            SurfaceMetricField(config[iZone], geometry[iZone][INST_0]);
          }

          /*--- Load the solution on the source mesh ---*/
          LoadRestarts(config[iZone], geometry[iZone], solver[iZone], iZone, INST_0, TimeIter, true);

          /*--- Normalize metric field if requested ---*/
          NormalizeMetricField(config[iZone], solver[iZone][INST_0][FLOW_SOL],
                               geometry[iZone][INST_0]);
        }

        if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << "." << endl;

        for (iZone = 0; iZone < nZone; iZone++) {
          WriteFiles(config[iZone], geometry[iZone][INST_0], solver[iZone][INST_0],
                     output[iZone], TimeIter);
        }
      }

      TimeIter++;
      if (StopCalc) break;
    }

    delete[] SolutionInstantiated;

  }

  else {
    /*--- Steady simulation: merge the single solution file. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      config[iZone]->SetiInst(INST_0);

      /*--- Initialize the solution classes ---*/
      InitializeSolver(config[iZone], geometry[iZone][INST_0], solver[iZone][INST_0],
                       iZone, INST_0, nZone);

      /*--- Initialize and preprocess the output ---*/
      InitializeOutput(config[iZone], geometry[iZone][INST_0], solver[iZone][INST_0],
                       output[iZone], iZone, INST_0, nZone);

      /*--- Calculate the surface metric ---*/
      SurfaceMetricField(config[iZone], geometry[iZone][INST_0]);

      /*--- Load the solution on the source mesh ---*/
      LoadRestarts(config[iZone], geometry[iZone], solver[iZone], iZone, INST_0, 0, true);

      /*--- Normalize metric field if requested ---*/
      NormalizeMetricField(config[iZone], solver[iZone][INST_0][FLOW_SOL],
                           geometry[iZone][INST_0]);
    }
    for (iZone = 0; iZone < nZone; iZone++) {
      WriteFiles(config[iZone], geometry[iZone][INST_0], solver[iZone][INST_0],
                 output[iZone], 0);
    }
  }

  delete config;
  config = nullptr;

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Finalize Solver -------------------------" << endl;

  if (geometry != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (geometry[iZone][iInst] != nullptr) {
          delete geometry[iZone][iInst];
        }
      }
      if (geometry[iZone] != nullptr) delete[] geometry[iZone];
    }
    delete[] geometry;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry containers." << endl;

  if (solver != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (solver[iZone][iInst] != nullptr) {
          for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {
            if (solver[iZone][iInst][iSol] != nullptr) {
              delete solver[iZone][iInst][iSol];
            }
          }
          delete solver[iZone][iInst];
        }
      }
      if (solver[iZone] != nullptr) delete[] solver[iZone];
    }
    delete[] solver;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSolver containers." << endl;

  if (config != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config[iZone] != nullptr) {
        delete config[iZone];
      }
    }
    delete[] config;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig containers." << endl;

  if (driver_config != nullptr) {
    delete driver_config;
  }
  if (rank == MASTER_NODE) cout << "Deleted driver CConfig class." << endl;

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
    cout << endl << "------------------------- Exit Success (SU2_MET) ------------------------" << endl << endl;

  /*--- Finalize MPI parallelization ---*/
  SU2_MPI::Finalize();

  return EXIT_SUCCESS;
}

void InitializeConfig(CConfig* driver_config, CConfig** config_container, char* zone_file_name,
                      char* config_file_name, SU2_COMPONENT val_software, int iZone, int nZone,
                      SU2_MPI::Comm MPICommunicator) {
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
}

void InitializeGeometry(CConfig* config, CGeometry*& geometry, int iZone, int iInst, int nZone) {
  const int rank = SU2_MPI::GetRank();
  if (rank == MASTER_NODE)
    cout << endl << "-------------------- Geometry Information ( Zone "  << iZone << " ) --------------------" << endl;

  /*--- Mesh initialization ---*/
  config->SetiInst(iInst);
  const bool fem_solver = config->GetFEMSolver();
  const bool fea = config->GetStructuralProblem();

  CGeometry* geometry_aux = nullptr;
  geometry_aux = new CPhysicalGeometry(config, iZone, nZone);

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

  /*--- Compute the surface curvature ---*/
  if (!fea) {
    if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
    geometry->ComputeSurf_Curvature(config);
  }

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
}

void InitializeSolver(CConfig* config, CGeometry* geometry, CSolver**& solver_container, int iZone,
                      int iInst, int nZone) {
  solver_container = new CSolver*[MAX_SOLS]();
  solver_container[FLOW_SOL] = new CBaselineSolver(geometry, config);
}

void InitializeOutput(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput*& output,
                                           int iZone, int iInst, int nZone) {
  output = new CMetricOutput(config, geometry->GetnDim(), solver_container[FLOW_SOL]);
  output->PreprocessVolumeOutput(config);
  output->PreprocessHistoryOutput(config, false);
}

void LoadRestarts(CConfig* config, CGeometry** geometry_container, CSolver*** solver_container, int iZone,
                                       int iInst, int TimeIter, bool UpdateGeo) {
  for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
    auto solver = solver_container[iInst][iSol];
    if (solver)
      solver->LoadRestart(geometry_container, solver_container, config, TimeIter, UpdateGeo);
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
  for (auto iFile = 0u; iFile < config->GetnVolumeOutputFiles(); ++iFile) {
    // if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII && FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY &&
    //     FileFormat[iFile] != OUTPUT_TYPE::CSV)
    output->WriteToFile(config, geometry, FileFormat[iFile]);
  }
}

vector<int> GetMetricFieldIndices(const CConfig* config, const CSolver* solver) {
  vector<int> indices;vector<string> fields = solver->GetSolutionFields();
  fields.erase(fields.begin()); // remove Point_ID
  for (size_t i = 0; i < fields.size(); ++i) {
    string field_name = fields[i].substr(1, fields[i].size() - 2);
    if (field_name.rfind("Metric_", 0) == 0) {
      indices.push_back(static_cast<int>(i));
    }
  }

  if (indices.size() == 0) {
    SU2_MPI::Error("Metric tensor not found in solution.", CURRENT_FUNCTION);
    return indices;
  }

  return indices;
}

void NormalizeMetricField(const CConfig* config, CSolver* solver, CGeometry* geometry) {
  const int rank = SU2_MPI::GetRank();

  const unsigned long nPointDomain = geometry->GetnPointDomain();
  const unsigned short nDim = geometry->GetnDim();
  const unsigned short nSymMat = solver->GetnSymMat();

  /*--- Check if normalization is enabled ---*/
  if (!config->GetNormalize_Metric()) {
    return;
  }

  /*--- Get metric field indices ---*/
  vector<int> iFields = GetMetricFieldIndices(config, solver);
  if (iFields.empty()) {
    SU2_MPI::Error("No metric fields found for normalization.", CURRENT_FUNCTION);
    return;
  } else if (iFields.size() != nSymMat) {
    if (rank == MASTER_NODE)
    SU2_MPI::Error(string("Incorrect number of tensor components. Should be ") +
                   to_string(nSymMat) + string("for ")  + to_string(nDim) + string("D."),
                   CURRENT_FUNCTION);
    return;
  }

  if (rank == MASTER_NODE) {
    cout << "Metric tensor found at indices ( ";
    for (auto i = 0u; i < nSymMat; ++i) {
      cout << iFields[i];
      if (i < nSymMat - 1) cout << ", ";
    }
    cout << " ) in solution." << endl;
  }

  /*--- Try to read integral file ---*/
  const bool tabTecplot = config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_TECPLOT;
  string filename = config->GetMetric_Integral_FileName();
  unsigned short lastindex = filename.find_last_of('.');
  filename = filename.substr(0, lastindex);
  if (tabTecplot)
    filename += ".dat";
  else
    filename += ".csv";

  ifstream integral_file(filename);
  if (!integral_file.good()) {
    if (rank == MASTER_NODE) {
      cout << "Warning: Cannot open metric integral file '" << filename
           << "'. Skipping normalization." << endl;
    }
    return;
  }

  /*--- Read and sum all integral values (space-time integral) ---*/
  su2double integral_value = 0.0;
  bool found_value = false;
  string line;
  int lines_to_skip = tabTecplot ? 2 : 1;  // Skip 2 lines for .dat, 1 line for .csv
  int lines_skipped = 0;

  while (getline(integral_file, line)) {
    /*--- Skip header lines ---*/
    if (lines_skipped < lines_to_skip) {
      lines_skipped++;
      continue;
    }

    /*--- Skip empty lines ---*/
    if (line.empty()) continue;

    /*--- Parse CSV line ---*/
    istringstream iss(line);
    string token;
    vector<string> tokens;

    /*--- Split by comma ---*/
    while (getline(iss, token, ',')) {
      /*--- Remove leading/trailing whitespace ---*/
      size_t start = token.find_first_not_of(" \t");
      if (start == string::npos) continue;
      size_t end = token.find_last_not_of(" \t");
      token = token.substr(start, end - start + 1);
      tokens.push_back(token);
    }

    /*--- We expect at least 2 values: TimeIter, Integral ---*/
    if (tokens.size() >= 2) {
      try {
        su2double integral_val = stod(tokens[1]);
        integral_value += integral_val;
        found_value = true;
      } catch (const exception& e) {
        if (rank == MASTER_NODE) {
          cout << "Warning: Could not parse line: " << line << endl;
        }
      }
    }
  }

  integral_file.close();

  if (!found_value || integral_value <= 0.0) {
    if (rank == MASTER_NODE) {
      cout << "Warning: Invalid integral value (" << integral_value
           << ") found. Skipping normalization." << endl;
    }
    return;
  }

  /*--- Create a metric container for the normalization ---*/
  su2matrix<su2double> metric_field(nPointDomain, nSymMat);

  /*--- Extract metric tensor from solution fields ---*/
  for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
    for (auto iSymMat = 0u; iSymMat < nSymMat; iSymMat++) {
      metric_field(iPoint, iSymMat) = solver->GetNodes()->GetSolution(iPoint, iFields[iSymMat]);
    }
  }

  /*--- Apply normalization using the tensor::metric interface ---*/
  const unsigned short iSensor = 0;
  normalizeMetrics<su2double, tensor::metric>(
    *geometry, *config, iSensor, integral_value, metric_field);

  /*--- Write the normalized metric back to the solution fields ---*/
  for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
    for (auto iSymMat = 0u; iSymMat < nSymMat; iSymMat++) {
      solver->GetNodes()->SetSolution(iPoint, iFields[iSymMat],
                                      metric_field(iPoint, iSymMat));
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Metric field normalization completed. Integrated determinant value: ";
    cout << setprecision(3) << scientific << integral_value << endl;
  }
}

void SurfaceMetricField(const CConfig* config, CGeometry* geometry) {
  if (!config->GetCompute_Metric_Geo()) return;

  const int rank = SU2_MPI::GetRank();

  const unsigned long nPointDomain = geometry->GetnPointDomain();
  const unsigned short nDim = geometry->GetnDim();
  const unsigned short nSymMat = 3 * (nDim - 1);

  /*--- Create a metric container  ---*/
  auto& metric_field = geometry->nodes->GetMetric();
  geometricSurfaceMetrics<su2double, tensor::metric>(
    *geometry, *config, metric_field
  );

  if (rank == MASTER_NODE)
    cout << "Surface metric calculation completed." << endl;
}