/*!
 * \file SU2_ERR.cpp
 * \brief Main file for the solution error estimation code (SU2_ERR).
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

#include "../include/SU2_ERR.hpp"

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
  CSolver*** solver_ref = nullptr;
  CConfig** config_src = nullptr;
  CConfig** config_dst = nullptr;
  CConfig** config_ref = nullptr;
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
  config = new CConfig(config_file_name, SU2_COMPONENT::SU2_ERR);

  const auto nZone = config->GetnZone();

  /*--- Definition of the containers per zones ---*/

  solver_src = new CSolver**[nZone]();
  solver_dst = new CSolver**[nZone]();
  solver_ref = new CSolver**[nZone]();
  config_src = new CConfig*[nZone]();
  config_dst = new CConfig*[nZone]();
  config_ref = new CConfig*[nZone]();
  geometry_src = new CGeometry**[nZone]();
  geometry_dst = new CGeometry**[nZone]();
  nInst = new unsigned short[nZone];
  driver_config = nullptr;
  output = new COutput*[nZone]();

  for (iZone = 0; iZone < nZone; iZone++) {
    nInst[iZone] = 1;
  }

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_ERR, false);

  /*--- Initialize a char to store the zone filename ---*/
  char zone_file_name[MAX_STRING_SIZE];

  /*--- Store a boolean for multizone problems ---*/
  const bool multizone = config->GetMultizone_Problem();

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    InitializeConfig(driver_config, config_src, zone_file_name, config_file_name, SU2_COMPONENT::SU2_ERR, iZone, nZone, MPICommunicator);
    InitializeConfig(driver_config, config_dst, zone_file_name, config_file_name, SU2_COMPONENT::SU2_ERR, iZone, nZone, MPICommunicator, false);
    InitializeConfig(driver_config, config_ref, zone_file_name, config_file_name, SU2_COMPONENT::SU2_ERR, iZone, nZone, MPICommunicator, false);

    /*--- Load the reference solution ---*/
    config_ref[iZone]->SetSolution_FileName(config_ref[iZone]->GetSolution_RefFileName());
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
    const bool dual_time_2nd = (config_src[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
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

      const bool IsTime0 = (TimeIter == 0);
      const bool IsTimeWrt = (TimeIter % config_src[ZONE_0]->GetVolumeOutputFrequency(0) == 0);
      const bool IsTimeEnd = (TimeIter + 1 == config_src[ZONE_0]->GetnTime_Iter()) ||
                             (TimeIter + 2 == config_src[ZONE_0]->GetnTime_Iter() && dual_time_2nd);
      const bool IsTimeRestart = ((long)TimeIter == SU2_TYPE::Int(config_src[ZONE_0]->GetRestart_Iter()));

      if (StopCalc || IsTime0 || IsTimeWrt || IsTimeEnd || IsTimeRestart) {
        /*--- Read in the restart file for this time step ---*/
        for (iZone = 0; iZone < nZone; iZone++) {
          /*--- Set the current iteration number in the config class. ---*/
          config_src[iZone]->SetTimeIter(TimeIter);
          config_dst[iZone]->SetTimeIter(TimeIter);

          /*--- Either instantiate the solution class or load a restart file. ---*/
          if (!SolutionInstantiated[iZone]) {
            /*--- Initialize the solution classes ---*/
            solver_src[iZone][INST_0] = new CBaselineSolver(geometry_src[iZone][INST_0], config_src[iZone]);
            solver_dst[iZone][INST_0] = new CBaselineSolver(geometry_dst[iZone][INST_0], config_dst[iZone]);
            solver_ref[iZone][INST_0] = new CBaselineSolver(geometry_dst[iZone][INST_0], config_ref[iZone]);

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
          config_ref[iZone]->SetiInst(INST_0);
          solver_src[iZone][INST_0]->LoadRestart(geometry_src[iZone], &solver_src[iZone], config_src[iZone],
                                                 TimeIter, true);

          /*--- Interpolate the solution ---*/
          InterpolateSolution(config_src[iZone], geometry_src[iZone][INST_0], geometry_dst[iZone][INST_0],
                              solver_src[iZone][INST_0], solver_dst[iZone][INST_0]);

          /*--- Load the reference solution on the destination mesh ---*/
          solver_ref[iZone][INST_0]->LoadRestart(geometry_dst[iZone], &solver_ref[iZone], config_ref[iZone],
                                                 TimeIter, true);

          /*--- Estimate the error ---*/

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

    delete[] SolutionInstantiated;

  }

  else {
    /*--- Steady simulation: merge the single solution file. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      config_src[iZone]->SetiInst(INST_0);
      config_dst[iZone]->SetiInst(INST_0);

      /*--- Initialize the solution classes ---*/
      solver_src[iZone][INST_0] = new CBaselineSolver(geometry_src[iZone][INST_0], config_src[iZone]);
      solver_dst[iZone][INST_0] = new CBaselineSolver(geometry_dst[iZone][INST_0], config_dst[iZone]);
      solver_ref[iZone][INST_0] = new CBaselineSolver(geometry_dst[iZone][INST_0], config_ref[iZone]);

      /*--- Initialize and preprocess the output ---*/
      output[iZone] =
          new CBaselineOutput(config_dst[iZone], geometry_dst[iZone][INST_0]->GetnDim(), solver_dst[iZone][INST_0]);
      output[iZone]->PreprocessVolumeOutput(config_dst[iZone]);
      output[iZone]->PreprocessHistoryOutput(config_dst[iZone], false);

      /*--- Load the solution on the source mesh ---*/
      solver_src[iZone][INST_0]->LoadRestart(geometry_src[iZone], &solver_src[iZone], config_src[iZone],
                                             0, true);

      /*--- Interpolate the solution ---*/
      InterpolateSolution(config_src[iZone], geometry_src[iZone][INST_0], geometry_dst[iZone][INST_0],
                          solver_src[iZone][INST_0], solver_dst[iZone][INST_0]);

      /*--- Load the reference solution on the destination mesh ---*/
      solver_ref[iZone][INST_0]->LoadRestart(geometry_dst[iZone], &solver_ref[iZone], config_ref[iZone],
                                             0, true);
      /*--- Estimate the error ---*/

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
  if (solver_ref != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (solver_ref[iZone][iInst] != nullptr) {
          delete solver_ref[iZone][iInst];
        }
      }
      if (solver_ref[iZone] != nullptr) delete[] solver_ref[iZone];
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
  if (config_ref != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_ref[iZone] != nullptr) {
        delete config_ref[iZone];
      }
    }
    delete[] config_ref;
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
    cout << endl << "------------------------- Exit Success (SU2_ERR) ------------------------" << endl << endl;

  /*--- Finalize MPI parallelization ---*/
  SU2_MPI::Finalize();

  return EXIT_SUCCESS;
}
