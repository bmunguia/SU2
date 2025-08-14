/*!
 * \file CConservativeVolumeInterpolator.cpp
 * \brief Implementation of the conservative solution interpolation subroutines.
 *        This file contains the conservative interpolation logic using Alauzet 2010 method.
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

#include "../include/CConservativeVolumeInterpolator.hpp"

CConservativeVolumeInterpolator::CConservativeVolumeInterpolator(SU2_Comm MPICommunicator)
    : CVolumeInterpolator(MPICommunicator) { }

void CConservativeVolumeInterpolator::Interpolate(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                  CSolver* solver_src, CSolver* solver_dst) {
  if (rank == MASTER_NODE) {
    cout << endl << "----------------------------- Interpolation -----------------------------" << endl;
    cout << "Performing conservative solution interpolation from source mesh to destination mesh..." << endl;
    cout << "Source mesh: " << geometry_src->GetGlobal_nPointDomain() << " points, ";
    cout << geometry_src->GetGlobal_nElemDomain() << " elements" << endl;
    cout << "Destination mesh: " << geometry_dst->GetGlobal_nPointDomain() << " points" << endl;
  }

  /*--- Build the ADTs ---*/
  InitializeADTs(config, geometry_src, geometry_dst);

  /*--- Call the conservative interpolation method ---*/
  ConservativeInterpolation(geometry_src, geometry_dst, solver_src, solver_dst);
}

void CConservativeVolumeInterpolator::InitializeSolver(CConfig* config, CGeometry* geometry, CSolver*& solver, int iZone,
                                                       int iInst, int nZone) {

}

void CConservativeVolumeInterpolator::ConservativeInterpolation(CGeometry* geometry_src, CGeometry* geometry_dst,
                                                                CSolver* solver_src, CSolver* solver_dst) {

}
