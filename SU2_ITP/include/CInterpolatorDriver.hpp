/*!
 * \file CInterpolatorDriver.hpp
 * \brief Driver class for solution interpolation, exposing the Python wrapper API.
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

#pragma once

#include "../../SU2_CFD/include/drivers/CDriverBase.hpp"
#include "CConservativeVolumeInterpolator.hpp"
#include "CLinearVolumeInterpolator.hpp"

/*!
 * \class CInterpolatorDriver
 * \brief Driver that wraps CVolumeInterpolator and inherits from CDriverBase so that
 *        it can be exposed via the pysu2 Python wrapper.
 *
 * Container ownership:
 *   CDriverBase::config_container   [iZone]              -> destination configs
 *   CDriverBase::geometry_container [iZone][INST_0][MESH_0] -> destination geometry
 *   CDriverBase::solver_container   [iZone][INST_0][MESH_0] -> destination solver array (CSolver**)
 *   CDriverBase::output_container   [iZone]              -> destination output
 *   CDriverBase::main_config                             -> config_container[ZONE_0]
 *   CDriverBase::main_geometry                           -> geometry_container[ZONE_0][INST_0][MESH_0]
 *
 *   config_src_    [iZone]           -> source configs (owned here, not via CDriverBase)
 *   geometry_src_  [iZone][iInst]    -> source geometry (owned here)
 *   solver_src_    [iZone][iInst]    -> source solver arrays CSolver** (owned here)
 *
 * Python workflow (steady, custom sensors):
 * \code
 *   driver = pysu2.CInterpolatorDriver("interp.cfg", comm)
 *   custom_sensors = CustomSensorRegistry({"MACH": mach_fn})
 *   custom_sensors.initialize(driver)   # GetMetricSensorIndex — from CDriverBase
 *   driver.Preprocess(0)               # LoadRestarts + Interpolate + PostprocessPrimitives
 *   custom_sensors.populate(driver)    # SetSensorAdapt per node
 *   driver.Postprocess()               # metric field computation
 *   driver.Output(0)
 *   driver.Finalize()
 * \endcode
 *
 * \author B. Munguía
 */
class CInterpolatorDriver : public CDriverBase {
 protected:
  /*--- Source-mesh containers (not part of CDriverBase) ---*/
  CConfig**    config_src_   = nullptr;  /*!< \brief Source configs [nZone]. */
  CGeometry*** geometry_src_ = nullptr;  /*!< \brief Source geometry [nZone][nInst]. */
  CSolver****  solver_src_   = nullptr;  /*!< \brief Source solvers [nZone][nInst] -> CSolver**. */

  /*--- Interpolator instances, one per zone ---*/
  CVolumeInterpolator** interpolator_ = nullptr;  /*!< \brief Volume interpolator [nZone]. */

  /*--- Lazy-init and first-call tracking ---*/
  bool* solution_instantiated_ = nullptr;  /*!< \brief Whether solvers/output are initialized per zone. */
  bool* metrics_initialized_   = nullptr;  /*!< \brief Whether first Postprocess has run per zone. */

  bool is_unsteady_ = false;  /*!< \brief True for time-domain simulations. */

 public:
  /*!
   * \brief Constructor. Reads config, builds all source and destination containers,
   *        and initialises geometries. Solvers and output are initialized lazily on
   *        the first Preprocess() call.
   * \param[in] confFile        - Configuration file name.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CInterpolatorDriver(char* confFile, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor. Calls Finalize() if not already called.
   */
  ~CInterpolatorDriver() override;

  /*!
   * \brief Load source-mesh restarts, interpolate to destination, and compute
   *        primitive variables on the destination mesh. After this returns the
   *        caller may inject custom metric sensors before calling Postprocess().
   * \param[in] TimeIter - Current time iteration index.
   */
  void Preprocess(unsigned long TimeIter);

  /*!
   * \brief No-op. Retained for API symmetry with CSinglezoneDriver.
   */
  void Run() override {}

  /*!
   * \brief Compute the metric field (sensor resolution, gradients, Hessians) on the
   *        destination mesh. Custom sensors must have been populated between
   *        Preprocess() and this call.
   */
  void Postprocess();

  /*!
   * \brief Write output files for all zones at the given time iteration.
   * \param[in] TimeIter - Current time iteration index.
   */
  void Output(unsigned long TimeIter);

  /*!
   * \brief Deallocate all owned memory (source + destination containers, interpolators).
   */
  void Finalize() override;

  /// \addtogroup PySU2
  /// @{

  /*!
   * \brief Return true when the problem is time-domain (unsteady).
   */
  bool IsUnsteady() const { return is_unsteady_; }

  /*!
   * \brief Return the first time iteration index (0, or restart iteration when restarting).
   */
  unsigned long GetStartTimeIter() const;

  /*!
   * \brief Return true when the time loop should stop for the given iteration.
   * \param[in] TimeIter - Current time iteration index.
   */
  bool StopCalc(unsigned long TimeIter) const;

  /*!
   * \brief Return true when output should be written for the given time iteration.
   * \param[in] TimeIter - Current time iteration index.
   */
  bool ShouldOutput(unsigned long TimeIter) const;

  /// @}

 protected:
  /*--- Internal setup helpers called from the constructor ---*/

  /*!
   * \brief Create the CVolumeInterpolator instances (one per zone).
   */
  void InitializeInterpolators();

  /*!
   * \brief Read zone configs for source and destination meshes and call
   *        CDriverBase::InitializeContainers() to size the base arrays.
   */
  void InitializeConfig();

  /*!
   * \brief Build source and destination geometries for all zones.
   */
  void InitializeGeometry();

  /*!
   * \brief Initialize solvers and output for a single zone (called lazily on first Preprocess).
   * \param[in] iZone - Zone index.
   */
  void InitializeSolversAndOutput(unsigned short iZone);
};
