/*!
 * \file CSodShockTubeSolution.hpp
 * \brief Header file for the class CSodShockTubeSolution.
 *        The implementations are in the <i>CSodShockTubeSolution.cpp</i> file.
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

#pragma once

#include <cmath>
#include "CVerificationSolution.hpp"

/*!
 * \class CSodShockTubeSolution
 * \brief Class to define the required data for the Sod shock tube.
 * \author B. Munguía
 */
class CSodShockTubeSolution final : public CVerificationSolution {
 protected:
  /*--- Sod shock tube initial conditions. ---*/
  su2double rhoL, uL, pL;   /*!< \brief Left state: density, velocity, pressure. */
  su2double rhoR, uR, pR;   /*!< \brief Right state: density, velocity, pressure. */
  su2double x0;             /*!< \brief Initial interface position. */

  /*--- Variables involving gamma. ---*/
  su2double Gamma;    /*!< \brief Gamma */
  su2double ovGm1;    /*!< \brief 1 over Gamma minus 1 */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CSodShockTubeSolution(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CSodShockTubeSolution(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_iMesh, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSodShockTubeSolution(void) override;

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double* val_coords, const su2double val_t, su2double* val_solution) const override;

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double* val_coords, const su2double val_t, su2double* val_solution) const override;

  /*!
   * \brief Whether or not the exact solution is known for this verification solution.
   * \return  - True, because the exact solution is known for the Sod shock tube case.
   */
  bool ExactSolutionKnown(void) const override;

 private:
  /*!
   * \brief Solves for the pressure ratio P in the Riemann problem using Newton's method.
   * \param[in] pL - Left pressure
   * \param[in] pR - Right pressure  
   * \param[in] cL - Left sound speed
   * \param[in] cR - Right sound speed
   * \return The pressure ratio P
   */
  su2double SolvePressureRatio(su2double pL, su2double pR, su2double cL, su2double cR) const;

  /*!
   * \brief Function whose root gives the pressure ratio P.
   * \param[in] P - Pressure ratio
   * \param[in] pL - Left pressure
   * \param[in] pR - Right pressure
   * \param[in] cL - Left sound speed  
   * \param[in] cR - Right sound speed
   * \return Function value
   */
  su2double PressureFunction(su2double P, su2double pL, su2double pR, su2double cL, su2double cR) const;
};
