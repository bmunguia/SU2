/*!
 * \file CKelvinHelmholtzSolution.cpp
 * \brief Implementations of the member functions of CKelvinHelmholtzSolution.
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

#include "../../../include/toolboxes/MMS/CKelvinHelmholtzSolution.hpp"

CKelvinHelmholtzSolution::CKelvinHelmholtzSolution() : CVerificationSolution() {}

CKelvinHelmholtzSolution::CKelvinHelmholtzSolution(unsigned short val_nDim, unsigned short val_nVar,
                                                   unsigned short val_iMesh, CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for the
   KH instability test case. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the KH instability case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Store the KH instability parameters here. ---*/
  dyInstability = 0.05;
  w0Instability = 0.1;
  nInstability = 1.0;

  rho1Instability = 2.0;
  rho2Instability = 1.0;
  u1Instability = 0.5;
  u2Instability = -0.5;

  p0Instability = 2.5;

  /*--- Useful coefficients in which Gamma is present. ---*/
  Gamma = config->GetGamma();
  ovGm1 = 1.0 / (Gamma - 1.0);

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if ((config->GetTime_Marching() != TIME_MARCHING::TIME_STEPPING) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the KH instability", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the KH instability", CURRENT_FUNCTION);

  if ((Kind_Solver != MAIN_SOLVER::EULER) && (Kind_Solver != MAIN_SOLVER::FEM_EULER))
    SU2_MPI::Error("Euler equations must be selected for the KH instability", CURRENT_FUNCTION);

  if ((config->GetKind_FluidModel() != STANDARD_AIR) && (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the KH instability", CURRENT_FUNCTION);
}

CKelvinHelmholtzSolution::~CKelvinHelmholtzSolution() = default;

void CKelvinHelmholtzSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                          su2double* val_solution) const {
  /*--- For the case that the KH instability is run with boundary
        conditions (other possibility is with periodic conditions),
        the exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CKelvinHelmholtzSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                           su2double* val_solution) const {
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /* Compute the ramp function. */
  const su2double Rden1 = 1.0 + exp(2.0 * (y - 0.25) / dyInstability);
  const su2double Rden2 = 1.0 + exp(2.0 * (y - 0.75) / dyInstability);
  const su2double R = 1.0 + 1.0 / Rden1 - 1.0 / Rden2;

  /* Compute the density and x-velocity from the ramp function. */
  const su2double rho = rho1Instability + R * (rho2Instability - rho1Instability);
  const su2double u = u1Instability + R * (u2Instability - u1Instability);

  /* Compute the sinusoidal y-velocity perturbation. */
  const su2double vden1 = 1.0 + exp(1.25 * (y - 1.75) / dyInstability);
  const su2double vden2 = 1.0 + exp(1.25 * (y + 0.75) / dyInstability);
  const su2double v = w0Instability * sin(2.0 * nInstability * PI_NUMBER * x) * (1.0 / vden1 - 1.0 / vden2);

  /* Compute the conservative variables. Note that both 2D and 3D
     cases are treated correctly. */
  val_solution[0] = rho;
  val_solution[1] = rho * u;
  val_solution[2] = rho * v;
  val_solution[3] = 0.0;
  val_solution[nVar - 1] = p0Instability * ovGm1 + 0.5 * rho * (u * u + v * v);
}

bool CKelvinHelmholtzSolution::ExactSolutionKnown() const { return false; }