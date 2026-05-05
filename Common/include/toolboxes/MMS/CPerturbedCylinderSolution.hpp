/*!
 * \file CPerturbedCylinderSolution.hpp
 * \brief Header file for the class CPerturbedCylinderSolution.
 *        The implementations are in the <i>CPerturbedCylinderSolution.cpp</i> file.
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

#pragma once

#include <cmath>
#include "CVerificationSolution.hpp"

/*!
 * \class CPerturbedCylinderSolution
 * \brief Class to define initial conditions for cylinder vortex-shedding studies.
 *        Cells in the upper half of the domain (y > 0) that lie within the radial
 *        distance CYLINDER_VEL_PERTURB_R from the cylinder center are initialized
 *        with the freestream x-velocity scaled by CYLINDER_VEL_PERTURB_FACTOR to
 *        break symmetry and produce a consistent shedding pattern. The exact
 *        solution is not known.
 * \author B. Munguía
 */
class CPerturbedCylinderSolution final : public CVerificationSolution {
 protected:
  su2double rhoInf;            /*!< \brief Freestream density (non-dimensional). */
  su2double pInf;              /*!< \brief Freestream pressure (non-dimensional). */
  su2double uInf;              /*!< \brief Freestream x-velocity (non-dimensional). */
  su2double vInf;              /*!< \brief Freestream y-velocity (non-dimensional). */
  su2double wInf;              /*!< \brief Freestream z-velocity (non-dimensional). */
  su2double perturbationFactor; /*!< \brief Velocity scaling factor applied in the perturbed zone. */
  su2double cylX0;              /*!< \brief x-coordinate of the cylinder center. */
  su2double cylY0;              /*!< \brief y-coordinate of the cylinder center. */
  su2double cylR;               /*!< \brief Outer radius of the perturbation zone. */

  su2double Gamma;  /*!< \brief Ratio of specific heats. */
  su2double ovGm1;  /*!< \brief 1 / (Gamma - 1). */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CPerturbedCylinderSolution(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nVar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CPerturbedCylinderSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                             CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CPerturbedCylinderSolution(void) override;

  /*!
   * \brief Get the solution (initial condition) at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the solution is stored.
   */
  void GetSolution(const su2double* val_coords, const su2double val_t, su2double* val_solution) const override;

  /*!
   * \brief Whether or not the exact solution is known for this verification solution.
   * \return  - False, because no analytical solution exists for the perturbed cylinder case.
   */
  bool ExactSolutionKnown(void) const override;
};
