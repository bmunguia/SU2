/*!
 * \file common.hpp
 * \brief Common convection-related methods.
 * \author P. Gomes, F. Palacios, T. Economon
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

#include "../../CNumericsSIMD.hpp"
#include "../../util.hpp"
#include "../variables.hpp"
#include "../../../variables/CNSVariable.hpp"

/*!
 * \brief Blended difference for U-MUSCL reconstruction.
 * \param[in] gradProj - Gradient projection at point i: dot(grad_i, vector_ij).
 * \param[in] delta - Centered difference: V_j - V_i.
 * \param[in] kappa - Blending parameter.
 * \return Blended difference for reconstruction from point i.
 */
FORCEINLINE Double umusclProjection(Double gradProj,
                                    Double delta,
                                    Double kappa) {
  /*-------------------------------------------------------------------*/
  /*--- The MUSCL kappa-scheme reconstruction is typically written: ---*/
  /*---     V_L = V_i + 0.25 * dV_ij^kap, where                     ---*/
  /*---     dV_ij^kap = (1-kappa) dV_ij^upw + (1+kappa) dV_ij^cen,  ---*/
  /*---     dV_ij^cen = V_j - V_i,                                  ---*/
  /*---     dV_ij^upw = 2 grad(Vi) dot vector_ij - dV_ij^cen.       ---*/
  /*--- To maintain proper scaling for edge limiters, the result of ---*/
  /*--- this function is 0.5 * dV_ij^kap.                           ---*/
  /*-------------------------------------------------------------------*/
  return (1.0 - kappa) * gradProj + kappa * delta;
}

/*!
 * \brief MUSCL reconstruction of the specified variable.
 * \note The result should be halved when added to i (or subtracted from j).
 * \param[in] grad_i - Gradient vector at point i.
 * \param[in] vector_ij - Distance vector from i to j.
 * \param[in] delta - Centered difference: V_j - V_i.
 * \param[in] iVar - Variable index.
 * \param[in] kappa - Blending coefficient.
 * \param[in] umusclRamp - MUSCL 1st-2nd order ramp times Newton-Krylov relaxation.
 * \return Variable reconstructed from point i.
 */
template<class GradType, size_t nDim>
FORCEINLINE Double musclReconstruction(const GradType& grad,
                                       const VectorDbl<nDim>& vector_ij,
                                       const Double delta,
                                       size_t iVar,
                                       Double kappa,
                                       Double umusclRamp) {
  const Double proj = dot(grad[iVar], vector_ij);
  return umusclRamp * umusclProjection(proj, delta, kappa);
}

/*!
 * \brief Unlimited reconstruction.
 */
template<size_t nVarGrad_ = 0, size_t nDim, class VarType, class Gradient_t>
FORCEINLINE void musclUnlimited(Int iPoint,
                                Int jPoint,
                                const VectorDbl<nDim>& vector_ij,
                                const Gradient_t& gradient,
                                CPair<VarType>& V,
                                Double kappa,
                                Double umusclRamp) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  auto grad_i = gatherVariables<nVarGrad,nDim>(iPoint, gradient);
  auto grad_j = gatherVariables<nVarGrad,nDim>(jPoint, gradient);

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Centered difference, needed for U-MUSCL projection ---*/
    const Double delta_ij = V.j.all(iVar) - V.i.all(iVar);

    /*--- U-MUSCL reconstructed variables ---*/
    const Double proj_i = musclReconstruction(grad_i, vector_ij, delta_ij, iVar, kappa, umusclRamp);
    const Double proj_j = musclReconstruction(grad_j, vector_ij, delta_ij, iVar, kappa, umusclRamp);

    /*--- Apply reconstruction: V_L = V_i + 0.5 * dV_ij^kap ---*/
    V.i.all(iVar) += 0.5 * proj_i;
    V.j.all(iVar) -= 0.5 * proj_j;
  }
}

/*!
 * \brief Beta scheme reconstruction for (beta = 1/3 for third-order accuracy).
 */
template<size_t nVarGrad_ = 0, size_t nDim, class VarType, class Gradient_t>
FORCEINLINE void betaUnlimited(Int iPoint,
                               Int jPoint,
                               const VectorDbl<nDim>& vector_ij,
                               const Gradient_t& gradient,
                               CPair<VarType>& V,
                               Double beta) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  auto grad_i = gatherVariables<nVarGrad,nDim>(iPoint, gradient);
  auto grad_j = gatherVariables<nVarGrad,nDim>(jPoint, gradient);

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Centered difference: dV_ij^cent = Vj - Vi ---*/
    const Double delta_cent = V.j.all(iVar) - V.i.all(iVar);

    /*--- Upwind difference: dV_ij^upwind = 2∇V(Ci)·(SiSj) - dV_ij^cent ---*/
    const Double grad_proj_i = 2.0 * dot(grad_i[iVar], vector_ij);
    const Double grad_proj_j = 2.0 * dot(grad_j[iVar], vector_ij);
    const Double delta_upwind_i = grad_proj_i - delta_cent;
    const Double delta_upwind_j = grad_proj_j - delta_cent;

    /*--- Beta-weighted difference: dV_ij^beta = (1-beta)ΔV_ij^cent + beta ΔV_ij^upwind ---*/
    const Double delta_beta_i = (1.0 - beta) * delta_cent + beta * delta_upwind_i;
    const Double delta_beta_j = (1.0 - beta) * delta_cent + beta * delta_upwind_j;

    /*--- Apply reconstruction: V_ij = Vi + 0.5 * ΔV_ij^β ---*/
    V.i.all(iVar) += 0.5 * delta_beta_i;
    V.j.all(iVar) -= 0.5 * delta_beta_j;
  }
}

/*!
 * \brief Limited reconstruction with point-based limiter.
 */
template<size_t nVarGrad_ = 0, size_t nDim, class VarType, class Limiter_t, class Gradient_t>
FORCEINLINE void musclPointLimited(Int iPoint,
                                   Int jPoint,
                                   const VectorDbl<nDim>& vector_ij,
                                   const Limiter_t& limiter,
                                   const Gradient_t& gradient,
                                   CPair<VarType>& V,
                                   Double kappa,
                                   Double umusclRamp) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  auto lim_i = gatherVariables<nVarGrad>(iPoint, limiter);
  auto lim_j = gatherVariables<nVarGrad>(jPoint, limiter);

  auto grad_i = gatherVariables<nVarGrad,nDim>(iPoint, gradient);
  auto grad_j = gatherVariables<nVarGrad,nDim>(jPoint, gradient);

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Centered difference, needed for U-MUSCL projection ---*/
    const Double delta_ij = V.j.all(iVar) - V.i.all(iVar);

    /*--- U-MUSCL reconstructed variables ---*/
    const Double proj_i = musclReconstruction(grad_i, vector_ij, delta_ij, iVar, kappa, umusclRamp);
    const Double proj_j = musclReconstruction(grad_j, vector_ij, delta_ij, iVar, kappa, umusclRamp);

    /*--- Apply reconstruction: V_L = V_i + 0.5 * lim * dV_ij^kap ---*/
    V.i.all(iVar) += 0.5 * lim_i(iVar) * proj_i;
    V.j.all(iVar) -= 0.5 * lim_j(iVar) * proj_j;
  }
}

/*!
 * \brief Limited reconstruction with edge-based limiter.
 */
template<size_t nVarGrad_ = 0, size_t nDim, class VarType, class Gradient_t>
FORCEINLINE void musclEdgeLimited(Int iPoint,
                                  Int jPoint,
                                  const VectorDbl<nDim>& vector_ij,
                                  const Gradient_t& gradient,
                                  CPair<VarType>& V,
                                  Double kappa,
                                  Double umusclRamp) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  auto grad_i = gatherVariables<nVarGrad,nDim>(iPoint, gradient);
  auto grad_j = gatherVariables<nVarGrad,nDim>(jPoint, gradient);

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Centered difference, needed for U-MUSCL projection and limiter ---*/
    const Double delta_ij = V.j.all(iVar) - V.i.all(iVar);
    const Double delta_ij_2 = pow(delta_ij, 2) + 1e-6;

    /*--- U-MUSCL reconstructed variables ---*/
    const Double proj_i = musclReconstruction(grad_i, vector_ij, delta_ij, iVar, kappa, umusclRamp);
    const Double proj_j = musclReconstruction(grad_j, vector_ij, delta_ij, iVar, kappa, umusclRamp);

    /// TODO: Customize the limiter function.
    const Double lim_i = (delta_ij_2 + proj_i*delta_ij) / (pow(proj_i,2) + delta_ij_2);
    const Double lim_j = (delta_ij_2 + proj_j*delta_ij) / (pow(proj_j,2) + delta_ij_2);

    /*--- Apply reconstruction: V_L = V_i + 0.5 * lim * dV_ij^kap ---*/
    V.i.all(iVar) += 0.5 * lim_i * proj_i;
    V.j.all(iVar) -= 0.5 * lim_j * proj_j;
  }
}

/*!
 * \brief Configurable Piperno limiter function φ(r, k) where r = 1/R (inverse slope ratio).
 * \param[in] r - Inverse slope ratio r = Δu_{i-1/2} / Δu_{i+1/2}.
 * \param[in] k - Piperno coefficient (k >= 1). k=1 recovers the standard Piperno limiter.
 * \note Uses branchless SIMD-friendly evaluation across all three regions.
 */
FORCEINLINE Double pipernoLimiterFunction(Double r, su2double k) {
  const Double pos_r = fmax(r, 0.0);

  /*--- Region R > k ---*/
  /*--- s = r / (1 + r*(1-k)); φ = 1 + (3/2 s + 1)(s-1)³ ---*/
  const Double s = pos_r / fmax(1.0 + pos_r * (1.0 - k), 1e-6);
  const Double s_minus_1 = s - 1.0;
  const Double phi_mild = 1.0 + (1.5 * s + 1.0) * pow(s_minus_1, 3);

  /*--- Region 1/k <= R <= k ---*/
  const Double phi_flat = 1.0;

  /*--- Region R < 1/k ---*/
  /*--- φ = (3r²-6r+19) / ((r-k)³+3r²-6r+19) ---*/
  const Double r_minus_k = pos_r - k;
  const Double r_sq = pow(pos_r, 2);
  const Double numerator = 3.0 * r_sq - 6.0 * pos_r + 19.0;
  const Double denominator = pow(r_minus_k, 3) + numerator;
  const Double phi_strong = numerator / fmax(denominator, 1e-6);

  /*--- Select region using branchless evaluation ---*/
  const Double inv_k = 1.0 / k;
  const Double phi = (pos_r < inv_k) * phi_mild + (pos_r >= inv_k) * ((pos_r <= k) * phi_flat + (pos_r > k) * phi_strong);

  return (r > 0.0) * phi;
}

/*!
 * \brief Configurable Piperno slope limiter reconstruction (edge formulation).
 * \param[in] k - Piperno coefficient (k >= 1). k=1 recovers the standard Piperno limiter.
 */
template<size_t nVarGrad_ = 0, size_t nDim, class VarType, class Gradient_t>
FORCEINLINE void musclPiperno(Int iPoint,
                              Int jPoint,
                              const VectorDbl<nDim>& vector_ij,
                              const Gradient_t& gradient,
                              CPair<VarType>& V,
                              su2double k) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  auto grad_i = gatherVariables<nVarGrad,nDim>(iPoint, gradient);
  auto grad_j = gatherVariables<nVarGrad,nDim>(jPoint, gradient);

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Compute differences for Piperno scheme ---*/
    const Double delta_ij = V.j.all(iVar) - V.i.all(iVar);  // Δu_{i+1/2}
    const Double proj_i = dot(grad_i[iVar], vector_ij);     // ∇u_i · Δx_{ij}
    const Double proj_j = dot(grad_j[iVar], vector_ij);     // ∇u_j · Δx_{ij}

    /*--- Compute upwind differences consistent with beta scheme ---*/
    /*--- Δu_{i-1/2} = 2∇u_i·Δx - Δu_{i+1/2} (upwind difference) ---*/
    /*--- Δu_{i+3/2} = 2∇u_j·Δx - Δu_{i+1/2} (upwind difference) ---*/
    const Double delta_imhalf = 2.0 * proj_i - delta_ij;
    const Double delta_jphalf = 2.0 * proj_j - delta_ij;

    /*--- Compute slope ratios r_i = 1/R_i and r_j = 1/R_j ---*/
    /*--- r_i = Δu_{i-1/2} / Δu_{i+1/2}, r_j = Δu_{i+3/2} / Δu_{i+1/2} ---*/
    const Double sign_delta_ij = (delta_ij >= 0.0) - (delta_ij < 0.0);
    const Double inv_delta_ij = sign_delta_ij / fmax(abs(delta_ij), 1e-6);

    const Double inv_R_i = delta_imhalf * inv_delta_ij;
    const Double inv_R_j = delta_jphalf * inv_delta_ij;

    /*--- Compute Piperno limiter functions φ(r, k) ---*/
    /*--- ψ(R) = (1/3 Δu_{i-1/2} + 2/3 Δu_{i+1/2}) φ(r) ---*/
    const Double phi_inv_R_i = pipernoLimiterFunction(inv_R_i, k);
    const Double phi_inv_R_j = pipernoLimiterFunction(inv_R_j, k);

    const Double proj_lim_i = (ONE3 * delta_imhalf + TWO3 * delta_ij) * phi_inv_R_i;
    const Double proj_lim_j = (ONE3 * delta_jphalf + TWO3 * delta_ij) * phi_inv_R_j;

    /*--- Apply Piperno reconstruction ---*/
    /*--- u_{i+1/2,L}^{lim} = u_i + 0.5 ψ(R_i) Δu_{i-1/2} ---*/
    /*--- u_{i+1/2,R}^{lim} = u_{i+1} - 0.5 ψ(R_j) Δu_{i+3/2} ---*/
    V.i.all(iVar) += 0.5 * proj_lim_i;
    V.j.all(iVar) -= 0.5 * proj_lim_j;
  }
}

/*!
 * \brief Retrieve primitive variables for points i/j, reconstructing them if needed.
 * \note Density and enthalpy are recomputed from ideal gas EOS.
 * \param[in] iEdge, iPoint, jPoint - Edge and its nodes.
 * \param[in] gamma - Heat capacity ratio.
 * \param[in] gasConst - Specific gas constant.
 * \param[in] muscl - If true, reconstruct, else simply fetch.
 * \param[in] kappa - Blending coefficient for MUSCL reconstruction.
 * \param[in] umusclRamp - MUSCL 1st-2nd order ramp times Newton-Krylov relaxation.
 * \param[in] limiterType - Type of flux limiter.
 * \param[in] pipernoK - Piperno limiter k coefficient (only used when limiterType == PIPERNO).
 * \param[in] V1st - Pair of compressible flow primitives for nodes i,j.
 * \param[in] vector_ij - Distance vector from i to j.
 * \param[in] solution - Entire solution container (a derived CVariable).
 * \return Pair of primitive variables.
 */
template<class ReconVarType, class PrimVarType, size_t nDim, class VariableType>
FORCEINLINE CPair<ReconVarType> reconstructPrimitives(Int iEdge, Int iPoint, Int jPoint,
                                                      const su2double& gamma,
                                                      const su2double& gasConst,
                                                      bool muscl,
                                                      const su2double& kappa,
                                                      const su2double& umusclRamp,
                                                      LIMITER limiterType,
                                                      const su2double& pipernoK,
                                                      const CPair<PrimVarType>& V1st,
                                                      const VectorDbl<nDim>& vector_ij,
                                                      const VariableType& solution) {
  static_assert(ReconVarType::nVar <= PrimVarType::nVar);

  const auto& gradients = solution.GetGradient_Reconstruction();
  const auto& limiters = solution.GetLimiter_Primitive();

  CPair<ReconVarType> V;

  for (size_t iVar = 0; iVar < ReconVarType::nVar; ++iVar) {
    V.i.all(iVar) = V1st.i.all(iVar);
    V.j.all(iVar) = V1st.j.all(iVar);
  }

  if (muscl) {
    /*--- Reconstruct density and enthalpy without using their gradients. ---*/
    constexpr auto nVarGrad = ReconVarType::nVar - 2;
    switch (limiterType) {
    case LIMITER::NONE:
      musclUnlimited<nVarGrad>(iPoint, jPoint, vector_ij, gradients, V, kappa, umusclRamp);
      break;
    case LIMITER::VAN_ALBADA_EDGE:
      musclEdgeLimited<nVarGrad>(iPoint, jPoint, vector_ij, gradients, V, kappa, umusclRamp);
      break;
    case LIMITER::PIPERNO:
      musclPiperno<nVarGrad>(iPoint, jPoint, vector_ij, gradients, V, pipernoK);
      break;
    default:
      musclPointLimited<nVarGrad>(iPoint, jPoint, vector_ij, limiters, gradients, V, kappa, umusclRamp);
      break;
    }
    /*--- Recompute density using the reconstructed pressure and temperature. ---*/
    V.i.density() = V.i.pressure() / (gasConst * V.i.temperature());
    V.j.density() = V.j.pressure() / (gasConst * V.j.temperature());

    /*--- Reconstruct enthalpy using dH/dT = Cp and dH/dv = v. Recomputing enthalpy would cause
     * stability issues because we use rho E = rho H - P, which loses its relation to temperature
     * if both rho and H are recomputed. This only seems to be an issue for wall-function meshes.
     * NOTE: This "one-sided" reconstruction does not lose much of the U-MUSCL benefit, because
     * the static enthalpy is linear, and for the KE term we do the equivalent of using the
     * average dH/dv, which is exact for quadratic functions. ---*/
    const su2double cp = gasConst * gamma / (gamma - 1);
    V.i.enthalpy() += cp * (V.i.temperature() - V1st.i.temperature());
    V.j.enthalpy() += cp * (V.j.temperature() - V1st.j.temperature());
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      V.i.enthalpy() += 0.5 * (pow(V.i.velocity(iDim), 2) - pow(V1st.i.velocity(iDim), 2));
      V.j.enthalpy() += 0.5 * (pow(V.j.velocity(iDim), 2) - pow(V1st.j.velocity(iDim), 2));
    }

    /*--- Detect a non-physical reconstruction based on negative pressure or density. ---*/
    const Double neg_p_or_rho = fmax(fmin(V.i.pressure(), V.j.pressure()) < 0.0,
                                     fmin(V.i.density(), V.j.density()) < 0.0);
    /*--- Test the sign of the Roe-averaged speed of sound. ---*/
    const Double R = sqrt(abs(V.j.density() / V.i.density()));
    /*--- Delay dividing by R+1 until comparing enthalpy and velocity magnitude. ---*/
    const Double enthalpy = R*V.j.enthalpy() + V.i.enthalpy();
    Double v_squared = 0.0;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      v_squared += pow(R*V.j.velocity(iDim) + V.i.velocity(iDim), 2);
    }
    /*--- Multiply enthalpy by R+1 since v^2 was not divided by (R+1)^2.
     * Note: a = sqrt((gamma-1) * (H - 0.5 * v^2)) ---*/
    const Double neg_sound_speed = enthalpy * (R+1) < 0.5 * v_squared;

    /*--- Revert to first order if the state is non-physical. ---*/
    Double bad_recon = fmax(neg_p_or_rho, neg_sound_speed);
    /*--- Handle SIMD dimensions 1 by 1. ---*/
    for (size_t k = 0; k < Double::Size; ++k) {
      bad_recon[k] = solution.UpdateNonPhysicalEdgeCounter(iEdge[k], bad_recon[k]);
    }
    for (size_t iVar = 0; iVar < ReconVarType::nVar; ++iVar) {
      V.i.all(iVar) = bad_recon * V1st.i.all(iVar) + (1-bad_recon) * V.i.all(iVar);
      V.j.all(iVar) = bad_recon * V1st.j.all(iVar) + (1-bad_recon) * V.j.all(iVar);
    }
  }
  return V;
}

/*!
 * \brief Compute and return the P tensor (compressible flow, ideal gas).
 */
template<size_t nDim, class RandomAccessIterator>
FORCEINLINE MatrixDbl<nDim+2> pMatrix(Double gamma, Double density, const RandomAccessIterator& velocity,
                                      Double projVel, Double speedSound, const VectorDbl<nDim>& normal) {
  MatrixDbl<nDim+2> pMat;
  const Double vel2 = 0.5*squaredNorm<nDim>(velocity);

  if (nDim == 2) {
    pMat(0,0) = 1.0;
    pMat(0,1) = 0.0;

    pMat(1,0) = velocity[0];
    pMat(1,1) = density*normal(1);

    pMat(2,0) = velocity[1];
    pMat(2,1) = -density*normal(0);

    pMat(3,0) = vel2;
    pMat(3,1) = density*(velocity[0]*normal(1) - velocity[1]*normal(0));
  }
  else {
    pMat(0,0) = normal(0);
    pMat(0,1) = normal(1);
    pMat(0,2) = normal(2);

    pMat(1,0) = velocity[0]*normal(0);
    pMat(1,1) = velocity[0]*normal(1) - density*normal(2);
    pMat(1,2) = velocity[0]*normal(2) + density*normal(1);

    pMat(2,0) = velocity[1]*normal(0) + density*normal(2);
    pMat(2,1) = velocity[1]*normal(1);
    pMat(2,2) = velocity[1]*normal(2) - density*normal(0);

    pMat(3,0) = velocity[2]*normal(0) - density*normal(1);
    pMat(3,1) = velocity[2]*normal(1) + density*normal(0);
    pMat(3,2) = velocity[2]*normal(2);

    pMat(4,0) = vel2*normal(0) + density*(velocity[1]*normal(2) - velocity[2]*normal(1));
    pMat(4,1) = vel2*normal(1) - density*(velocity[0]*normal(2) - velocity[2]*normal(0));
    pMat(4,2) = vel2*normal(2) + density*(velocity[0]*normal(1) - velocity[1]*normal(0));
  }

  /*--- Last two columns. ---*/

  const Double rhoOn2 = 0.5*density;
  const Double rhoOnTwoC = rhoOn2 / speedSound;
  pMat(0,nDim) = rhoOnTwoC;
  pMat(0,nDim+1) = rhoOnTwoC;

  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    pMat(iDim+1,nDim) = rhoOnTwoC * velocity[iDim] + rhoOn2 * normal(iDim);
    pMat(iDim+1,nDim+1) = rhoOnTwoC * velocity[iDim] - rhoOn2 * normal(iDim);
  }

  pMat(nDim+1,nDim) = rhoOnTwoC * vel2 + rhoOn2 * (projVel + speedSound/(gamma-1));
  pMat(nDim+1,nDim+1) = rhoOnTwoC * vel2 - rhoOn2 * (projVel - speedSound/(gamma-1));

  return pMat;
}

/*!
 * \brief Compute and return the inverse P tensor (compressible flow, ideal gas).
 */
template<size_t nDim, class RandomAccessIterator>
FORCEINLINE MatrixDbl<nDim+2> pMatrixInv(Double gamma, Double density, const RandomAccessIterator& velocity,
                                         Double projVel, Double speedSound, const VectorDbl<nDim>& normal) {
  MatrixDbl<nDim+2> pMatInv;

  const Double c2 = pow(speedSound,2);
  const Double vel2 = 0.5*squaredNorm<nDim>(velocity);
  const Double oneOnRho = 1 / density;

  if (nDim == 2) {
    Double tmp = (gamma-1)/c2;
    pMatInv(0,0) = 1.0 - tmp*vel2;
    pMatInv(0,1) = tmp*velocity[0];
    pMatInv(0,2) = tmp*velocity[1];
    pMatInv(0,3) = -tmp;

    pMatInv(1,0) = (normal(0)*velocity[1]-normal(1)*velocity[0])*oneOnRho;
    pMatInv(1,1) = normal(1)*oneOnRho;
    pMatInv(1,2) = -normal(0)*oneOnRho;
    pMatInv(1,3) = 0.0;
  }
  else {
    Double tmp = (gamma-1)/c2 * normal(0);
    pMatInv(0,0) = normal(0) - tmp*vel2 - (normal(2)*velocity[1]-normal(1)*velocity[2])*oneOnRho;
    pMatInv(0,1) = tmp*velocity[0];
    pMatInv(0,2) = tmp*velocity[1] + normal(2)*oneOnRho;
    pMatInv(0,3) = tmp*velocity[2] - normal(1)*oneOnRho;
    pMatInv(0,4) = -tmp;

    tmp = (gamma-1)/c2 * normal(1);
    pMatInv(1,0) = normal(1) - tmp*vel2 + (normal(2)*velocity[0]-normal(0)*velocity[2])*oneOnRho;
    pMatInv(1,1) = tmp*velocity[0] - normal(2)*oneOnRho;
    pMatInv(1,2) = tmp*velocity[1];
    pMatInv(1,3) = tmp*velocity[2] + normal(0)*oneOnRho;
    pMatInv(1,4) = -tmp;

    tmp = (gamma-1)/c2 * normal(2);
    pMatInv(2,0) = normal(2) - tmp*vel2 - (normal(1)*velocity[0]-normal(0)*velocity[1])*oneOnRho;
    pMatInv(2,1) = tmp*velocity[0] + normal(1)*oneOnRho;
    pMatInv(2,2) = tmp*velocity[1] - normal(0)*oneOnRho;
    pMatInv(2,3) = tmp*velocity[2];
    pMatInv(2,4) = -tmp;
  }

  /*--- Last two rows. ---*/

  const Double gamma_minus_1_on_rho_times_c = (gamma-1) / (density*speedSound);

  for (size_t iVar = nDim; iVar < nDim+2; ++iVar) {
    Double sign = (iVar==nDim)? 1 : -1;
    pMatInv(iVar,0) = -sign*projVel*oneOnRho + gamma_minus_1_on_rho_times_c * vel2;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      pMatInv(iVar,iDim+1) = sign*normal(iDim)*oneOnRho - gamma_minus_1_on_rho_times_c * velocity[iDim];
    }
    pMatInv(iVar,nDim+1) = gamma_minus_1_on_rho_times_c;
  }

  return pMatInv;
}

/*!
 * \brief Convective projected (onto normal) flux (compressible flow).
 */
template<class PrimVarType, class ConsVarType, size_t nDim>
FORCEINLINE VectorDbl<nDim+2> inviscidProjFlux(const PrimVarType& V,
                                               const ConsVarType& U,
                                               const VectorDbl<nDim>& normal) {
  static_assert(ConsVarType::nVar == nDim+2);
  Double mdot = dot(U.momentum(), normal);
  VectorDbl<nDim+2> flux;
  flux(0) = mdot;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    flux(iDim+1) = mdot*V.velocity(iDim) + normal(iDim)*V.pressure();
  }
  flux(nDim+1) = mdot*V.enthalpy();
  return flux;
}

/*!
 * \brief Jacobian of the convective flux (compressible flow, ideal gas).
 */
template<size_t nDim, class RandomAccessIterator>
FORCEINLINE MatrixDbl<nDim+2> inviscidProjJac(Double gamma, RandomAccessIterator velocity,
                                              Double energy, const VectorDbl<nDim>& normal,
                                              Double scale) {
  MatrixDbl<nDim+2> jac;

  Double projVel = dot(velocity, normal);
  Double gamma_m_1 = gamma-1;
  Double phi = 0.5*gamma_m_1*squaredNorm<nDim>(velocity);
  Double a1 = gamma*energy - phi;

  jac(0,0) = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    jac(0,iDim+1) = scale * normal(iDim);
  }
  jac(0,nDim+1) = 0.0;

  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    jac(iDim+1,0) = scale * (normal(iDim)*phi - velocity[iDim]*projVel);
    for (size_t jDim = 0; jDim < nDim; ++jDim) {
      jac(iDim+1,jDim+1) = scale * (normal(jDim)*velocity[iDim] - gamma_m_1*normal(iDim)*velocity[jDim]);
    }
    jac(iDim+1,iDim+1) += scale * projVel;
    jac(iDim+1,nDim+1) = scale * gamma_m_1 * normal(iDim);
  }

  jac(nDim+1,0) = scale * projVel * (phi-a1);
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    jac(nDim+1,iDim+1) = scale * (normal(iDim)*a1 - gamma_m_1*velocity[iDim]*projVel);
  }
  jac(nDim+1,nDim+1) = scale * gamma * projVel;

  return jac;
}

/*!
 * \brief (Low) Dissipation coefficient for Roe schemes.
 */
template<class VariableType>
FORCEINLINE Double roeDissipation(Int iPoint,
                                  Int jPoint,
                                  ENUM_ROELOWDISS type,
                                  const VariableType& solution) {
  if (type == NO_ROELOWDISS) {
    return 1.0;
  }

  const auto& sol = static_cast<const CNSVariable&>(solution);
  const auto& sensor = sol.GetSensor();
  const auto& dissip = sol.GetRoe_Dissipation();

  const Double si = gatherVariables(iPoint, sensor);
  const Double sj = gatherVariables(jPoint, sensor);
  const Double avgSensor = 0.5 * (si + sj);

  const Double di = gatherVariables(iPoint, dissip);
  const Double dj = gatherVariables(jPoint, dissip);
  const Double avgDissip = 0.5 * (di + dj);

  /*--- A minimum level of upwinding is used to enhance stability. ---*/
  constexpr passivedouble minDissip = 0.05;

  switch (type) {
    case FD:
    case FD_DUCROS: {
      Double d = fmax(minDissip, 1.0 - avgDissip);

      if (type == FD_DUCROS) {
        /*--- See Jonhsen et al. JCP 229 (2010) pag. 1234 ---*/
        d = fmax(d, 0.05 + 0.95*(avgSensor > 0.65));
      }
      return d;
    }
    case NTS:
      return fmax(minDissip, avgDissip);

    case NTS_DUCROS:
      /*--- See Xiao et al. INT J HEAT FLUID FL 51 (2015) pag. 141
       * https://doi.org/10.1016/j.ijheatfluidflow.2014.10.007 ---*/
      return fmax(minDissip, avgSensor+avgDissip - avgSensor*avgDissip);

    default:
      assert(false);
      return 1.0;
  }
}

/*!
 * \brief Correct spectral radius (avgLambda) for stretching.
 */
template<class VariableType, class T>
FORCEINLINE Double correctedSpectralRadius(Int iPoint,
                                           Int jPoint,
                                           Double avgLambda,
                                           T stretchParam,
                                           const VariableType& solution) {

  const auto lambda_i = gatherVariables(iPoint, solution.GetLambda());
  const Double phi_i = pow(0.25*lambda_i/avgLambda, stretchParam);

  const auto lambda_j = gatherVariables(jPoint, solution.GetLambda());
  const Double phi_j = pow(0.25*lambda_j/avgLambda, stretchParam);

  return 4*phi_i*phi_j / (phi_i + phi_j) * avgLambda;
}

/*!
 * \brief Update of a flux Jacobian due to a scalar dissipation term.
 */
template<class VariableType, size_t nVar>
FORCEINLINE void scalarDissipationJacobian(const VariableType& V,
                                           Double gamma,
                                           Double dissipConst,
                                           MatrixDbl<nVar>& jac) {
  /*--- Diagonal entries. ---*/
  for (size_t iVar = 0; iVar < nVar-1; ++iVar) {
    jac(iVar,iVar) += dissipConst;
  }
  jac(nVar-1,nVar-1) += dissipConst * gamma;

  /*--- N-1 columns of last row. ---*/
  dissipConst *= gamma-1.0;
  for (size_t iDim = 0; iDim < VariableType::nDim; ++iDim) {
    jac(nVar-1,iDim+1) -= dissipConst * V.velocity(iDim);
    jac(nVar-1,0) += dissipConst * pow(V.velocity(iDim), 2);
  }
}
