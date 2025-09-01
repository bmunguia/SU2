/*!
 * \file SU2_ERR.hpp
 * \brief Headers of the main subroutines of the code SU2_ERR.
 *        The subroutines and functions are in the <i>SU2_ERR.cpp</i> file.
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

#include "../../SU2_ITP/include/CConservativeVolumeInterpolator.hpp"
#include "../../SU2_ITP/include/CLinearVolumeInterpolator.hpp"

/*!
 * \brief Get the index of the solution field corresponding to the metric sensor.
 *
 * \param[in] config - Definition of the particular problem.
 * \param[in] solver - Container with the solution.
 * \return Index of the field in the solution fields vector corresponding to the metric sensor.
 * \throws SU2_MPI::Error if the sensor field is not found in the solution fields.
 */
int GetSensorFieldIndex(const CConfig* config, const CSolver* solver);

/*!
 * \brief Estimate the error between two solutions based on the specified metric sensor and norm.
 * \param[in] config - Definition of the particular problem.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] solver_dst - Container with the destination (interpolated) solution.
 * \param[in] solver_ref - Container with the reference (fine-mesh) solution.
 * \param[in] iFieldDst - Index of sensor in destination solution.
 * \param[in] iFieldRef - Index of sensor in reference solution.
 * \return Lp-norm error for the specified field.
 */
su2double EstimateFieldError(const CConfig* config, CGeometry* geometry,
                             CSolver* solver_dst, CSolver* solver_ref,
                             int iFieldDst, int iFieldRef);

