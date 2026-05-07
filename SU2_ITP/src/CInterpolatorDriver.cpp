/*!
 * \file CInterpolatorDriver.cpp
 * \brief Implementation of the interpolator driver.
 * \author B. Munguía
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

#include "../include/CInterpolatorDriver.hpp"
#include "../../SU2_CFD/include/metrics/metricUtils.hpp"

using namespace std;

CInterpolatorDriver::CInterpolatorDriver(char* confFile, SU2_Comm MPICommunicator)
    : CDriverBase(confFile, 1, MPICommunicator) {

  /*--- Read the driver configuration to determine the number of zones and
   *    interpolation method. nZone is overwritten from its initial value of 1. ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_ITP);
  nZone = driver_config->GetnZone();
  is_unsteady_ = driver_config->GetTime_Domain();

  /*--- Create the interpolator object(s) ---*/
  InitializeInterpolators();

  /*--- Read zone configs and size all CDriverBase container arrays ---*/
  InitializeConfig();

  /*--- Build source and destination geometries ---*/
  InitializeGeometry();

  /*--- Allocate init-state flags ---*/
  solution_instantiated_ = new bool[nZone]();  // zero-init -> false
  interpolation_initialized_ = new bool[nZone]();  // zero-init -> false
  metrics_initialized_   = new bool[nZone]();  // zero-init -> false

  /*--- Eagerly initialize destination/source solvers, sensor arrays, and output
   *    so the Python API can query metric sensors immediately after construction.
   *    Sensor resolution happens inside InitializeSolversAndOutput, between
   *    solver creation and output creation, so the output sees the full sensor list. ---*/
  for (auto iZone = 0u; iZone < nZone; iZone++) {
    InitializeSolversAndOutput(iZone);
  }
}

CInterpolatorDriver::~CInterpolatorDriver() {
  /*--- Finalize if the user has not already done so ---*/
  if (interpolator_ != nullptr) Finalize();
}

// ---------------------------------------------------------------------------
// Setup helpers
// ---------------------------------------------------------------------------

void CInterpolatorDriver::InitializeInterpolators() {
  interpolator_ = new CVolumeInterpolator*[nZone]();
  for (auto iZone = 0u; iZone < nZone; iZone++) {
    if (driver_config->GetKindVolumeInterpolation() == VOLUME_INTERPOLATOR::CONSERVATIVE)
      interpolator_[iZone] = new CConservativeVolumeInterpolator(SU2_MPI::GetComm());
    else
      interpolator_[iZone] = new CLinearVolumeInterpolator(SU2_MPI::GetComm());
  }
}

void CInterpolatorDriver::InitializeConfig() {
  /*--- Allocate CDriverBase container arrays (config, geometry, solver, output, etc.) ---*/
  InitializeContainers();

  /*--- Allocate source container arrays ---*/
  config_src_   = new CConfig*[nZone]();
  geometry_src_ = new CGeometry**[nZone]();
  solver_src_   = new CSolver***[nZone]();

  char zone_file_name[MAX_STRING_SIZE];
  const bool multizone = driver_config->GetMultizone_Problem();

  for (auto iZone = 0u; iZone < nZone; iZone++) {
    /*--- Source config ---*/
    interpolator_[iZone]->InitializeConfig(driver_config, config_src_, zone_file_name,
                                           config_file_name, SU2_COMPONENT::SU2_ITP,
                                           iZone, nZone, SU2_MPI::GetComm(), /*isSource=*/true);
    /*--- Destination config -> stored in CDriverBase::config_container ---*/
    interpolator_[iZone]->InitializeConfig(driver_config, config_container, zone_file_name,
                                           config_file_name, SU2_COMPONENT::SU2_ITP,
                                           iZone, nZone, SU2_MPI::GetComm(), /*isSource=*/false);
  }

  if (multizone) {
    for (auto iZone = 0u; iZone < nZone; iZone++)
      config_src_[iZone]->SetMultizone(driver_config, config_src_);
  }

  main_config = config_container[ZONE_0];

  /*--- Allocate per-zone instance-level slices for source ---*/
  for (auto iZone = 0u; iZone < nZone; iZone++) {
    geometry_src_[iZone] = new CGeometry*[nInst[iZone]]();
    solver_src_[iZone]   = new CSolver**[nInst[iZone]]();
  }
}

void CInterpolatorDriver::InitializeGeometry() {
  for (auto iZone = 0u; iZone < nZone; iZone++) {
    for (auto iInst = 0u; iInst < nInst[iZone]; iInst++) {
      /*--- Source geometry ---*/
      interpolator_[iZone]->InitializeGeometry(config_src_[iZone],
                                               geometry_src_[iZone][iInst],
                                               iZone, iInst, nZone, /*isSource=*/true);

      /*--- Destination geometry: allocate the inner [iInst][MESH_0] layers of
       *    geometry_container before passing geometry_container[iZone][iInst][MESH_0]
       *    by reference to InitializeGeometry. ---*/
      geometry_container[iZone] = new CGeometry**[nInst[iZone]]();
      geometry_container[iZone][iInst] = new CGeometry*[1]();
      interpolator_[iZone]->InitializeGeometry(config_container[iZone],
                                               geometry_container[iZone][iInst][MESH_0],
                                               iZone, iInst, nZone, /*isSource=*/false);
    }
  }
  main_geometry = geometry_container[ZONE_0][INST_0][MESH_0];
  nDim = main_geometry->GetnDim();
}

void CInterpolatorDriver::InitializeSolversAndOutput(unsigned short iZone) {
  const auto iInst = INST_0;

  config_src_[iZone]->SetiInst(iInst);
  config_container[iZone]->SetiInst(iInst);

  /*--- Allocate solver container layers for this zone in CDriverBase ---*/
  solver_container[iZone] = new CSolver***[nInst[iZone]]();
  solver_container[iZone][iInst] = new CSolver**[1]();
  solver_container[iZone][iInst][MESH_0] = nullptr;

  /*--- Initialize solvers for source and destination ---*/
  interpolator_[iZone]->InitializeSolver(config_src_[iZone],
                                         geometry_src_[iZone][iInst],
                                         solver_src_[iZone][iInst],
                                         iZone, iInst, nZone);
  interpolator_[iZone]->InitializeSolver(config_container[iZone],
                                         geometry_container[iZone][iInst][MESH_0],
                                         solver_container[iZone][iInst][MESH_0],
                                         iZone, iInst, nZone);

  /*--- Resolve metric sensor indices and allocate sensor arrays before
   *    initializing the output, so the output object sees the full sensor list
   *    and registers Hessian fields with the correct names and indices. ---*/
  if (config_container[iZone]->GetCompute_Metric()) {
    if (rank == MASTER_NODE)
      cout << "Resolving metric sensor indices." << endl;

    bool resolved = MetricUtils::ResolveSensorIndices(
      config_container[iZone],
      geometry_container[iZone][iInst][MESH_0],
      solver_container[iZone][iInst][MESH_0]
    );

    if (resolved) {
      MetricUtils::InitializeMetrics(solver_container[iZone][iInst][MESH_0]);
      const auto total = MetricUtils::TotalNumSensors(solver_container[iZone][iInst][MESH_0]);
      if (rank == MASTER_NODE && total > 0)
        cout << "Successfully resolved " << total << " metric sensors." << endl;
    } else if (rank == MASTER_NODE) {
      cout << "Warning: COMPUTE_METRIC is enabled but no valid sensors found." << endl;
    }
  }

  /*--- Initialize and preprocess the output for the destination mesh ---*/
  interpolator_[iZone]->InitializeOutput(config_container[iZone],
                                         geometry_container[iZone][iInst][MESH_0],
                                         solver_container[iZone][iInst][MESH_0],
                                         output_container[iZone],
                                         iZone, iInst, nZone);

  solution_instantiated_[iZone] = true;
}

// ---------------------------------------------------------------------------
// Driver API
// ---------------------------------------------------------------------------

void CInterpolatorDriver::Preprocess(unsigned long TimeIter) {
  this->TimeIter = TimeIter;

  for (auto iZone = 0u; iZone < nZone; iZone++) {
    config_src_[iZone]->SetTimeIter(TimeIter);
    config_container[iZone]->SetTimeIter(TimeIter);

    const bool initial_interp = !interpolation_initialized_[iZone];

    /*--- Load source restart ---*/
    interpolator_[iZone]->LoadRestarts(config_src_[iZone],
                                       geometry_src_[iZone],
                                       solver_src_[iZone],
                                       iZone, INST_0, TimeIter, /*UpdateGeo=*/true);

    /*--- Interpolate source → destination ---*/
    interpolator_[iZone]->Interpolate(config_src_[iZone],
                                      geometry_src_[iZone][INST_0],
                                      geometry_container[iZone][INST_0][MESH_0],
                                      solver_src_[iZone][INST_0],
                                      solver_container[iZone][INST_0][MESH_0],
                                      initial_interp);

    interpolation_initialized_[iZone] = true;

    /*--- Compute primitive variables so the caller can read them (e.g. to populate custom sensors) ---*/
    interpolator_[iZone]->PostprocessPrimitives(config_container[iZone],
                                                geometry_container[iZone][INST_0][MESH_0],
                                                solver_container[iZone][INST_0][MESH_0],
                                                initial_interp);
  }
}

void CInterpolatorDriver::Postprocess() {
  for (auto iZone = 0u; iZone < nZone; iZone++) {
    const bool first_metric = !metrics_initialized_[iZone];
    interpolator_[iZone]->Postprocess(config_container[iZone],
                                      geometry_container[iZone][INST_0][MESH_0],
                                      solver_container[iZone][INST_0][MESH_0],
                                      first_metric);
    metrics_initialized_[iZone] = true;
  }
}

void CInterpolatorDriver::Output(unsigned long TimeIter) {
  if (rank == MASTER_NODE)
    cout << "Writing the volume solution for time step " << TimeIter << "." << endl;

  for (auto iZone = 0u; iZone < nZone; iZone++) {
    interpolator_[iZone]->WriteFiles(config_container[iZone],
                                     geometry_container[iZone][INST_0][MESH_0],
                                     solver_container[iZone][INST_0][MESH_0],
                                     output_container[iZone],
                                     TimeIter);
  }
}

void CInterpolatorDriver::Finalize() {
  /*--- Source solvers ---*/
  if (solver_src_ != nullptr) {
    for (auto iZone = 0u; iZone < nZone; iZone++) {
      if (solver_src_[iZone] != nullptr) {
        for (auto iInst = 0u; iInst < nInst[iZone]; iInst++) {
          if (solver_src_[iZone][iInst] != nullptr) {
            for (auto iSol = 0u; iSol < MAX_SOLS; iSol++)
              delete solver_src_[iZone][iInst][iSol];
            delete[] solver_src_[iZone][iInst];
          }
        }
        delete[] solver_src_[iZone];
      }
    }
    delete[] solver_src_;
    solver_src_ = nullptr;
  }

  /*--- Source geometries ---*/
  if (geometry_src_ != nullptr) {
    for (auto iZone = 0u; iZone < nZone; iZone++) {
      if (geometry_src_[iZone] != nullptr) {
        for (auto iInst = 0u; iInst < nInst[iZone]; iInst++)
          delete geometry_src_[iZone][iInst];
        delete[] geometry_src_[iZone];
      }
    }
    delete[] geometry_src_;
    geometry_src_ = nullptr;
  }

  /*--- Source configs ---*/
  if (config_src_ != nullptr) {
    for (auto iZone = 0u; iZone < nZone; iZone++)
      delete config_src_[iZone];
    delete[] config_src_;
    config_src_ = nullptr;
  }

  /*--- Interpolators ---*/
  if (interpolator_ != nullptr) {
    for (auto iZone = 0u; iZone < nZone; iZone++)
      delete interpolator_[iZone];
    delete[] interpolator_;
    interpolator_ = nullptr;
  }

  delete[] solution_instantiated_;
  solution_instantiated_ = nullptr;
  delete[] interpolation_initialized_;
  interpolation_initialized_ = nullptr;
  delete[] metrics_initialized_;
  metrics_initialized_ = nullptr;

  /*--- Destination containers are cleaned up by CDriverBase::CommonFinalize() ---*/
  CommonFinalize();
}

// ---------------------------------------------------------------------------
// Unsteady helpers
// ---------------------------------------------------------------------------

unsigned long CInterpolatorDriver::GetStartTimeIter() const {
  const auto* cfg = config_src_[ZONE_0];
  if (cfg->GetTime_Domain() && cfg->GetRestart())
    return cfg->GetRestart_Iter();
  return 0;
}

bool CInterpolatorDriver::StopCalc(unsigned long TimeIter) const {
  const auto* cfg = config_src_[ZONE_0];
  if (TimeIter >= cfg->GetnTime_Iter()) return true;
  const su2double Physical_t = (TimeIter + 1) * cfg->GetTime_Step();
  return Physical_t >= cfg->GetMax_Time();
}

bool CInterpolatorDriver::ShouldOutput(unsigned long TimeIter) const {
  if (StopCalc(TimeIter)) return true;

  const auto* cfg = config_src_[ZONE_0];
  const bool IsTime0       = (TimeIter == 0);
  const bool IsTimeWrt     = (TimeIter % cfg->GetVolumeOutputFrequency(0) == 0);
  const bool IsTimeEnd     = (TimeIter + cfg->GetnRestartFinalIters() >= cfg->GetnTime_Iter());
  const bool IsTimeRestart = ((long)TimeIter == SU2_TYPE::Int(cfg->GetRestart_Iter()));

  return IsTime0 || IsTimeWrt || IsTimeEnd || IsTimeRestart;
}
