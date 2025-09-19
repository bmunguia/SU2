/*!
 * \file CSodShockTubeSolution.cpp
 * \brief Implementations of the member functions of CSodShockTubeSolution.
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

#include "../../../include/toolboxes/MMS/CSodShockTubeSolution.hpp"

CSodShockTubeSolution::CSodShockTubeSolution() : CVerificationSolution() {}

CSodShockTubeSolution::CSodShockTubeSolution(unsigned short val_nDim, unsigned short val_nVar,
                                       unsigned short val_iMesh, CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*  Write a message that the solution is initialized for the
   Sod shock tube test case. */
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Sod shock tube case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Store the Sod shock tube initial conditions here. ---*/
  rhoL = 1.0;    uL = 0.0;    pL = 1.0;     // Left state
  rhoR = 0.125;  uR = 0.0;    pR = 0.1;     // Right state
  x0 = 0.0;                                 // Interface position (middle of domain from -1 to 1)

  /* Useful coefficients in which Gamma is present. */
  Gamma = config->GetGamma();
  ovGm1 = 1.0 / (Gamma - 1.0);

  /* Perform some sanity and error checks for this solution here. */
  if ((config->GetTime_Marching() != TIME_MARCHING::TIME_STEPPING) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the Sod shock tube", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the Sod shock tube", CURRENT_FUNCTION);

  if ((config->GetKind_FluidModel() != STANDARD_AIR) && (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the Sod shock tube", CURRENT_FUNCTION);
}

CSodShockTubeSolution::~CSodShockTubeSolution() = default;

void CSodShockTubeSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                    su2double* val_solution) const {
  /*--- For the case that the Sod shock tube is run with boundary
        conditions, the exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CSodShockTubeSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                     su2double* val_solution) const {
  /* Get the x-coordinate (assuming 1D problem along x-axis) */
  const su2double x = val_coords[0];
  
  /* If time is essentially zero, return initial conditions */
  su2double rho, u, p;
  if (val_t < 1e-12) {
    if (x < x0) {
      rho = rhoL; u = uL; p = pL;
    } else {
      rho = rhoR; u = uR; p = pR;
    }
  } else {
    /* Compute sound speeds */
    const su2double cL = sqrt(Gamma * pL / rhoL);
    const su2double cR = sqrt(Gamma * pR / rhoR);
    
    /* Solve for pressure ratio P */
    const su2double P = SolvePressureRatio(pL, pR, cL, cR);
    
    /* Compute wave speeds and positions */
    const su2double c_shock = uR + cR * sqrt((Gamma - 1.0 + P * (Gamma + 1.0)) / (2.0 * Gamma));
    const su2double x_shock = x0 + c_shock * val_t;
    
    const su2double alpha = (Gamma + 1.0) / (Gamma - 1.0);
    const su2double c_contact = uL + 2.0 * cL / (Gamma - 1.0) * (1.0 - pow(P * pR / pL, (Gamma - 1.0) / 2.0 / Gamma));
    const su2double x_contact = x0 + c_contact * val_t;
    
    const su2double rho3 = rhoL * pow(P * pR / pL, 1.0 / Gamma);
    const su2double p3 = P * pR;
    const su2double c_fanright = c_contact - sqrt(Gamma * p3 / rho3);
    const su2double x_fanright = x0 + c_fanright * val_t;
    
    const su2double c_fanleft = -cL;
    const su2double x_fanleft = x0 + c_fanleft * val_t;
    
    /* Determine which region we're in and set solution accordingly */
    if (x >= x_shock) {
      /* Region R (unshocked right state) */
      rho = rhoR;
      u = uR;
      p = pR;
    } else if (x >= x_contact) {
      /* Region 2 (shocked right state) */
      rho = (1.0 + alpha * P) / (alpha + P) * rhoR;
      u = c_contact;
      p = P * pR;
    } else if (x >= x_fanright) {
      /* Region 3 (contact region) */
      rho = rho3;
      u = c_contact;
      p = p3;
    } else if (x >= x_fanleft) {
      /* Region 4 (rarefaction fan) */
      const su2double u4 = 2.0 / (Gamma + 1.0) * (cL + (x - x0) / val_t);
      rho = rhoL * pow(1.0 - (Gamma - 1.0) / 2.0 * u4 / cL, 2.0 / (Gamma - 1.0));
      u = u4;
      p = pL * pow(1.0 - (Gamma - 1.0) / 2.0 * u4 / cL, 2.0 * Gamma / (Gamma - 1.0));
    } else {
      /* Region L (unshocked left state) */
      rho = rhoL;
      u = uL;
      p = pL;
    }
  }

  /* Total energy */
  const su2double rhoE = p * ovGm1 + 0.5 * rho * (u * u);

  /* Store the solution */
  val_solution[0] = rho;
  val_solution[1] = rho * u;
  val_solution[2] = 0.0;
  val_solution[3] = 0.0;
  val_solution[nVar - 1] = rhoE;
}

su2double CSodShockTubeSolution::SolvePressureRatio(su2double pL, su2double pR, su2double cL, su2double cR) const {
  /* Use Newton's method to solve for pressure ratio P */
  su2double P = 0.5; // Initial guess
  const su2double tol = 1e-12, eps = 1e-8;
  const int max_iter = 100;
  
  for (int iter = 0; iter < max_iter; iter++) {
    const su2double f = PressureFunction(P, pL, pR, cL, cR);
    
    if (fabs(f) < tol) break;
    
    /* Compute derivative numerically for simplicity */
    const su2double df_dP = (PressureFunction(P + eps, pL, pR, cL, cR) - f) / eps;
    
    /* Newton update */
    P = P - f / df_dP;
    
    /* Ensure P stays positive */
    P = std::max(P, 1e-6);
  }
  
  return P;
}

su2double CSodShockTubeSolution::PressureFunction(su2double P, su2double pL, su2double pR, 
                                               su2double cL, su2double cR) const {
  const su2double a = (Gamma - 1.0) * (cR / cL) * (P - 1.0);
  const su2double b = sqrt(2.0 * Gamma * (2.0 * Gamma + (Gamma + 1.0) * (P - 1.0)));
  return P - pL / pR * pow(1.0 - a / b, 2.0 * Gamma / (Gamma - 1.0));
}

bool CSodShockTubeSolution::ExactSolutionKnown() const { return true; }
