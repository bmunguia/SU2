/*!
 * \file SU2_MET.hpp
 * \brief Headers of the main subroutines of the code SU2_MET.
 *        The subroutines and functions are in the <i>SU2_MET.cpp</i> file.
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

#include <cmath>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>

#include "../../SU2_CFD/include/solvers/CBaselineSolver.hpp"
#include "../../SU2_CFD/include/solvers/CBaselineSolver_FEM.hpp"
#include "../../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../../SU2_CFD/include/output/CBaselineOutput.hpp"
#include "../../SU2_CFD/include/output/COutputFactory.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../SU2_CFD/include/metrics/computeMetrics.hpp"

void InitializeConfig(CConfig* driver_config, CConfig** config_container, char* zone_file_name,
                      char* config_file_name, SU2_COMPONENT val_software, int iZone, int nZone,
                      SU2_MPI::Comm MPICommunicator);

void InitializeGeometry(CConfig* config, CGeometry*& geometry, int iZone, int iInst, int nZone);

void InitializeSolver(CConfig* config, CGeometry* geometry, CSolver**& solver_container, int iZone,
                      int iInst, int nZone);

void InitializeOutput(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput*& output,
                      int iZone, int iInst, int nZone);

void LoadRestarts(CConfig* config, CGeometry** geometry_container, CSolver*** solver_container, int iZone,
                  int iInst, int TimeIter, bool UpdateGeo);

void WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container,
                COutput* output, unsigned long TimeIter);

/*!
 * \brief Get the index of the metric tensor field for the current solver time step.
 *
 * \param[in] config - Definition of the particular problem.
 * \param[in] solver - Container with the solution.
 * \return Indices of the fields in the solution fields vector corresponding to the metric tensor.
 * \throws SU2_MPI::Error if the sensor field is not found in the solution fields.
 */
vector<int> GetMetricFieldIndices(const CConfig* config, const CSolver* solver);

/*!
 * \brief Normalize the metric tensor field for the current solver time step.
 *
 * \param[in] config - Definition of the particular problem.
 * \param[in] solver - Container with the solution.
 * \param[in] geometry - Geometrical definition of the problem.
 */
void NormalizeMetricField(const CConfig* config, CSolver* solver, CGeometry* geometry);