/*!
 * \file CPerturbedCylinderSolution.cpp
 * \brief Implementations of the member functions of CPerturbedCylinderSolution.
 * \author B. Munguía
 * \version 8.4.0 "Harrier"
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

#include "../../../include/toolboxes/MMS/CPerturbedCylinderSolution.hpp"

CPerturbedCylinderSolution::CPerturbedCylinderSolution() : CVerificationSolution() {}

CPerturbedCylinderSolution::CPerturbedCylinderSolution(unsigned short val_nDim, unsigned short val_nVar,
                                                       unsigned short val_iMesh, CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the perturbed cylinder case!!!" << endl;
    cout << endl << flush;
  }

  Gamma = config->GetGamma();
  ovGm1 = 1.0 / (Gamma - 1.0);

  rhoInf = config->GetDensity_FreeStreamND();
  pInf   = config->GetPressure_FreeStreamND();

  const su2double* vel = config->GetVelocity_FreeStreamND();
  uInf = vel[0];
  vInf = vel[1];
  wInf = (nDim == 3) ? vel[2] : 0.0;

  perturbationFactor = config->GetCylinderVelPerturbFactor();
  cylX0 = config->GetCylinderVelPerturbX0();
  cylY0 = config->GetCylinderVelPerturbY0();
  cylR  = config->GetCylinderVelPerturbR();

  if ((config->GetTime_Marching() != TIME_MARCHING::TIME_STEPPING) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the perturbed cylinder case", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible NS or RANS equations must be selected for the perturbed cylinder case",
                   CURRENT_FUNCTION);
}

CPerturbedCylinderSolution::~CPerturbedCylinderSolution() = default;

void CPerturbedCylinderSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                             su2double* val_solution) const {
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  const su2double dx = x - cylX0;
  const su2double dy = y - cylY0;
  const su2double dist = sqrt(dx * dx + dy * dy);
  const su2double factor = (y > 0.0 && dist < cylR) ? perturbationFactor : 1.0;

  const su2double u = uInf * factor;
  const su2double v = vInf * factor;
  const su2double w = (nDim == 3) ? wInf * factor : 0.0;

  val_solution[0]        = rhoInf;
  val_solution[1]        = rhoInf * u;
  val_solution[2]        = rhoInf * v;
  val_solution[3]        = (nDim == 3) ? rhoInf * w : 0.0;
  val_solution[nVar - 1] = pInf * ovGm1 + 0.5 * rhoInf * (u * u + v * v + w * w);
}

bool CPerturbedCylinderSolution::ExactSolutionKnown() const { return false; }
