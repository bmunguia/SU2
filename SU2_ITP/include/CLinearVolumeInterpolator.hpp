/*!
 * \file CLinearVolumeInterpolator.hpp
 * \brief Headers of the main linear solution interpolation subroutines.
 *        The subroutines and functions are in the <i>CLinearVolumeInterpolator.cpp</i> file.
 * \author B. Munguía, E. van der Weide
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

#include "CVolumeInterpolator.hpp"

/*!
 * \class CLinearVolumeInterpolator
 * \brief Performs linear solution interpolation between meshes.
 * \author B. Munguía, E. van der Weide
 */
class CLinearVolumeInterpolator : public CVolumeInterpolator {
  public:
    /*!
     * \brief Constructor of the class.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CLinearVolumeInterpolator(SU2_Comm MPICommunicator);

    /*!
     * \brief Default destructor.
     */
    ~CLinearVolumeInterpolator() override = default;

    /*!
     * \brief Main interpolation routine using linear interpolation.
     * \param[in] config - Configuration object
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_container_src - Source mesh solver
     * \param[in] solver_container_dst - Destination mesh solver
     */
    void Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                     CSolver** solver_container_src, CSolver** solver_container_dst) override;

  private:
    void LinearInterpolation(const CConfig *config, CGeometry* geometry_src, CGeometry* geometry_dst,
                             CSolver* solver_src, CSolver* solver_dst);

    void VolumeInterpolation(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                             const vector<su2double> &coor_corrected, vector<unsigned long> &pointsFailed);

    void SurfaceInterpolation(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                              const vector<su2double> &coor_dst, vector<unsigned long> &pointsFailed);
};
