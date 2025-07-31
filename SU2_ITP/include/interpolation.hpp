/*!
 * \file interpolation.hpp
 * \brief Headers of the main solution interpolation subroutines.
 *        The subroutines and functions are in the <i>interpolation.cpp</i> file.
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

#include <cmath>
#include <memory>

#include "../../Common/include/parallelization/mpi_structure.hpp"

#include "../../SU2_CFD/include/solvers/CBaselineSolver.hpp"
#include "../../SU2_CFD/include/solvers/CBaselineSolver_FEM.hpp"
#include "../../SU2_CFD/include/output/CBaselineOutput.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/adt/CADTElemClass.hpp"
#include "../../Common/include/CConfig.hpp"

void InitializeConfig(CConfig* driver_config, CConfig** config_container, char* zone_file_name,
                      char* config_file_name, SU2_COMPONENT val_software, int iZone, int nZone,
                      SU2_MPI::Comm MPICommunicator, bool isSource = true);

void InitializeGeometry(CConfig* config, CGeometry*& geometry, int iZone, int iInst, int nZone);

std::unique_ptr<CADTElemClass> BuildSurfaceADT(const CConfig* config, CGeometry* geometry);

std::unique_ptr<CADTElemClass> BuildVolumeADT(CGeometry* geometry);

void InterpolateSolution(const CConfig *config, CGeometry* geometry_src, CGeometry* geometry_dst,
                         CSolver* solver_src, CSolver* solver_dst);

void VolumeInterpolationSolution(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                                 const vector<su2double> &coor_corrected, vector<unsigned long> &pointsFailed);

void SurfaceInterpolationSolution(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                                  const vector<su2double> &coor_dst, vector<unsigned long> &pointsFailed);

void NearestPointOnElement(CGeometry* geometry, unsigned short markerID, unsigned long elemID,
                           const su2double* coor, su2double* surfCoor, su2double& dist2Elem,
                           const unsigned short nDim);

void NearestPointOnLine(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                        const su2double* coor, su2double* surfCoor, su2double& dist2Line,
                        const unsigned short nDim);

bool NearestPointOnTriangle(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                            const unsigned long i2, const su2double* coor, su2double* surfCoor,
                            su2double& dist2Tria, const unsigned short nDim);

bool NearestPointOnQuadrilateral(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                                 const unsigned long i2, const unsigned long i3, const su2double* coor,
                                 su2double* surfCoor, su2double& dist2Quad, const unsigned short nDim);

void ApplyCurvatureCorrection(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                              const unsigned short nDim, const vector<su2double> &coor_dst,
                              vector<su2double> &coor_corrected);

void WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container,
                COutput* output, unsigned long TimeIter);
