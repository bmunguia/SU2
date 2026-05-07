/*!
 * \file computeGradientsL2Projection.hpp
 * \brief Generic implementation of lumped-mass L2-projection gradient/Hessian computation.
 * \note Only implemented for triangles (2D) and tetrahedra (3D). Boundary nodes accumulate
 *       contributions solely from interior elements, avoiding contamination from BC states.
 * \author B. Munguía
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "correctGradientsSymmetry.hpp"

namespace detail {

/*!
 * \brief Compute the inward normal to the face opposite node iNode in a simplex element.
 * \note For triangles: inward normal to edge opposite vertex iNode (magnitude = edge length).
 *       For tetrahedra: inward normal to face opposite vertex iNode (magnitude = face area).
 */
template<size_t nDim, class ElemType, class NodeType>
void GetInwardNormal(const ElemType& elem, const NodeType& nodes,
                     size_t iNode, su2double* normal) {
  if constexpr (nDim == 2) {
    const size_t p0 = elem->GetNode((iNode + 1) % 3);
    const size_t p1 = elem->GetNode((iNode + 2) % 3);
    normal[0] = nodes->GetCoord(p0, 1) - nodes->GetCoord(p1, 1);
    normal[1] = nodes->GetCoord(p1, 0) - nodes->GetCoord(p0, 0);
  } else {
    const size_t p0 = elem->GetNode((iNode + 1) % 4);
    const size_t p1 = elem->GetNode((iNode + 2) % 4);
    const size_t p2 = elem->GetNode((iNode + 3) % 4);

    const su2double v1[3] = {
      nodes->GetCoord(p1, 0) - nodes->GetCoord(p0, 0),
      nodes->GetCoord(p1, 1) - nodes->GetCoord(p0, 1),
      nodes->GetCoord(p1, 2) - nodes->GetCoord(p0, 2)
    };
    const su2double v2[3] = {
      nodes->GetCoord(p2, 0) - nodes->GetCoord(p0, 0),
      nodes->GetCoord(p2, 1) - nodes->GetCoord(p0, 1),
      nodes->GetCoord(p2, 2) - nodes->GetCoord(p0, 2)
    };

    normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
    normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
    normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

    /*--- Flip to inward if necessary. ---*/
    const size_t pOpp = elem->GetNode(iNode);
    const su2double toOpp[3] = {
      nodes->GetCoord(pOpp, 0) - (nodes->GetCoord(p0, 0) + nodes->GetCoord(p1, 0) + nodes->GetCoord(p2, 0)) / 3.0,
      nodes->GetCoord(pOpp, 1) - (nodes->GetCoord(p0, 1) + nodes->GetCoord(p1, 1) + nodes->GetCoord(p2, 1)) / 3.0,
      nodes->GetCoord(pOpp, 2) - (nodes->GetCoord(p0, 2) + nodes->GetCoord(p1, 2) + nodes->GetCoord(p2, 2)) / 3.0
    };
    if (normal[0]*toOpp[0] + normal[1]*toOpp[1] + normal[2]*toOpp[2] < 0.0) {
      normal[0] = -normal[0]; normal[1] = -normal[1]; normal[2] = -normal[2];
    }
  }
}

/*!
 * \brief Compute gradients using the lumped-mass L2-projection on simplicial elements.
 * \note Gradient at jNode accumulates: (1/(nDim+1)!) * u(iNode) * n_iNode / V_jNode
 *       for all (iNode, jNode) pairs in each element, where n_iNode is the inward normal
 *       to the face opposite iNode. Boundary nodes only receive interior-element contributions.
 */
template<size_t nDim, class FieldType, class GradientType>
void computeGradientsL2Projection(CSolver* solver,
                                  MPI_QUANTITIES kindMpiComm,
                                  PERIODIC_QUANTITIES kindPeriodicComm,
                                  CGeometry& geometry,
                                  const CConfig& config,
                                  const FieldType& field,
                                  const size_t varBegin,
                                  const size_t varEnd,
                                  const int idxVel,
                                  GradientType& gradient)
{
  const size_t nPointDomain = geometry.GetnPointDomain();
  const size_t nElem = geometry.GetnElem();
  const size_t nNode = (nDim == 2) ? 3 : 4;

  /*--- factor = 1/(nDim+1)! : 1/6 for triangles, 1/24 for tetrahedra. ---*/
  const su2double factor = (nDim == 2) ? 1.0/6.0 : 1.0/24.0;

  auto nodes = geometry.nodes;

  /*--- Zero out gradient for all domain points. ---*/
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        gradient(iPoint, iVar, iDim) = 0.0;

  su2double normal[nDim];

  /*--- Accumulate element contributions. ---*/
  for (size_t iElem = 0; iElem < nElem; ++iElem) {
    auto elem = geometry.elem[iElem];

    for (size_t iNode = 0; iNode < nNode; ++iNode) {
      const size_t iPoint = elem->GetNode(iNode);
      GetInwardNormal<nDim>(elem, nodes, iNode, normal);

      for (size_t jNode = 0; jNode < nNode; ++jNode) {
        const size_t jPoint = elem->GetNode(jNode);
        if (!nodes->GetDomain(jPoint)) continue;

        const su2double invVol = 1.0 / (nodes->GetVolume(jPoint) + nodes->GetPeriodicVolume(jPoint));
        for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
          const su2double val = factor * field(iPoint, iVar) * invVol;
          for (size_t iDim = 0; iDim < nDim; ++iDim)
            gradient(jPoint, iVar, iDim) += val * normal[iDim];
        }
      }
    }
  }

  /*--- Symmetry-plane and Euler-wall corrections. ---*/
  correctGradientsSymmetry<nDim>(geometry, config, varBegin, varEnd, idxVel, gradient);

  if (solver == nullptr) return;

  for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic) {
    solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
  }

  solver->InitiateComms(&geometry, &config, kindMpiComm);
  solver->CompleteComms(&geometry, &config, kindMpiComm);
}

/*!
 * \brief Compute Hessians using the lumped-mass L2-projection on simplicial elements.
 * \note Applies the same projection to the gradient field. Uses lower-triangular symmetric
 *       storage with index ind = row*(row+1)/2 + col (row >= col), matching Green-Gauss.
 *       Off-diagonal entries are symmetrized: H_kl += 0.5*(g_l*n[k] + g_k*n[l]).
 */
template<size_t nDim, class GradientType>
void computeHessiansL2Projection(CSolver* solver,
                                 MPI_QUANTITIES kindMpiComm,
                                 PERIODIC_QUANTITIES kindPeriodicComm,
                                 CGeometry& geometry,
                                 const CConfig& config,
                                 const GradientType& gradient,
                                 const size_t varBegin,
                                 const size_t varEnd,
                                 const int idxVel,
                                 GradientType& hessian)
{
  const size_t nPointDomain = geometry.GetnPointDomain();
  const size_t nElem = geometry.GetnElem();
  const size_t nNode = (nDim == 2) ? 3 : 4;
  const size_t nSymMat = 3 * (nDim - 1);

  const su2double factor = (nDim == 2) ? 1.0/6.0 : 1.0/24.0;

  auto nodes = geometry.nodes;

  /*--- Zero out Hessian for all domain points. ---*/
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iMat = 0; iMat < nSymMat; ++iMat)
        hessian(iPoint, iVar, iMat) = 0.0;

  su2double normal[nDim];

  /*--- Accumulate element contributions. ---*/
  for (size_t iElem = 0; iElem < nElem; ++iElem) {
    auto elem = geometry.elem[iElem];

    for (size_t iNode = 0; iNode < nNode; ++iNode) {
      const size_t iPoint = elem->GetNode(iNode);
      GetInwardNormal<nDim>(elem, nodes, iNode, normal);

      for (size_t jNode = 0; jNode < nNode; ++jNode) {
        const size_t jPoint = elem->GetNode(jNode);
        if (!nodes->GetDomain(jPoint)) continue;

        const su2double invVol = 1.0 / (nodes->GetVolume(jPoint) + nodes->GetPeriodicVolume(jPoint));
        for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {

          /*--- Diagonal entries: H_kk += factor * g_k(iPoint) * n[k] / V ---*/
          for (size_t k = 0; k < nDim; ++k) {
            const size_t ind = k * (k + 1) / 2 + k;
            hessian(jPoint, iVar, ind) += factor * gradient(iPoint, iVar, k) * normal[k] * invVol;
          }

          /*--- Off-diagonal entries (lower triangle, k > l):
           *    H_kl += 0.5 * factor * (g_l(iPoint)*n[k] + g_k(iPoint)*n[l]) / V ---*/
          for (size_t k = 1; k < nDim; ++k) {
            for (size_t l = 0; l < k; ++l) {
              const size_t ind = k * (k + 1) / 2 + l;
              hessian(jPoint, iVar, ind) += factor * 0.5 *
                (gradient(iPoint, iVar, l) * normal[k] +
                 gradient(iPoint, iVar, k) * normal[l]) * invVol;
            }
          }
        }
      }
    }
  }

  if (solver == nullptr) return;

  for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic) {
    solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
  }

  solver->InitiateComms(&geometry, &config, kindMpiComm);
  solver->CompleteComms(&geometry, &config, kindMpiComm);
}

} // namespace detail

/*--- Public instantiations for 2D and 3D. ---*/

template<class FieldType, class GradientType>
void computeGradientsL2Projection(CSolver* solver,
                                  MPI_QUANTITIES kindMpiComm,
                                  PERIODIC_QUANTITIES kindPeriodicComm,
                                  CGeometry& geometry,
                                  const CConfig& config,
                                  const FieldType& field,
                                  const size_t varBegin,
                                  const size_t varEnd,
                                  const int idxVel,
                                  GradientType& gradient) {
  switch (geometry.GetnDim()) {
    case 2:
      detail::computeGradientsL2Projection<2>(solver, kindMpiComm, kindPeriodicComm,
                                              geometry, config, field, varBegin, varEnd, idxVel, gradient);
      break;
    case 3:
      detail::computeGradientsL2Projection<3>(solver, kindMpiComm, kindPeriodicComm,
                                              geometry, config, field, varBegin, varEnd, idxVel, gradient);
      break;
    default:
      SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
  }
}

template<class GradientType>
void computeHessiansL2Projection(CSolver* solver,
                                 MPI_QUANTITIES kindMpiComm,
                                 PERIODIC_QUANTITIES kindPeriodicComm,
                                 CGeometry& geometry,
                                 const CConfig& config,
                                 const GradientType& gradient,
                                 const size_t varBegin,
                                 const size_t varEnd,
                                 const int idxVel,
                                 GradientType& hessian) {
  switch (geometry.GetnDim()) {
    case 2:
      detail::computeHessiansL2Projection<2>(solver, kindMpiComm, kindPeriodicComm,
                                             geometry, config, gradient, varBegin, varEnd, idxVel, hessian);
      break;
    case 3:
      detail::computeHessiansL2Projection<3>(solver, kindMpiComm, kindPeriodicComm,
                                             geometry, config, gradient, varBegin, varEnd, idxVel, hessian);
      break;
    default:
      SU2_MPI::Error("Too many dimensions to compute Hessians.", CURRENT_FUNCTION);
  }
}
