/*!
 * \file CLimiterDetails.hpp
 * \brief A class template that allows defining limiters via
 *        specialization of particular details.
 * \author P. Gomes
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


/*!
 * \brief A traits class for limiters, see notes for "computeLimiters_impl()".
 * \ingroup FvmAlgos
 * \note There is no default implementation (the code will compile but not
 *       link) specialization is mandatory.
 */
template<LIMITER LimiterKind>
struct CLimiterDetails
{
  /*!
   * \brief Compute any global value that may be needed by the other functions.
   * \note This function is called once by multiple threads.
   */
  template<class FieldType>
  inline void preprocess(CGeometry&, CConfig&, size_t varBegin,
                         size_t varEnd, const FieldType&);

  /*!
   * \brief Geometric modifier (e.g. increase limiting near sharp edges).
   * \note This function is called once per point inside an AD pre-
   *       -accumulation region, newly used variables should be registered.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry&) const;

  /*!
   * \brief Smooth (usually) function of the maximum/minimum (positive/negative)
   *        gradient projections onto the edges, and the deltas over direct neighbors.
   *        Both proj and delta may be 0.0, beware of divisions.
   * \note This function is called twice (min/max) per point per variable
   *       (also inside an AD pre-accumulation region).
   */
  inline su2double limiterFunction(size_t iVar, su2double proj, su2double delta) const;
};


/*!
 * \brief Common small functions used by limiters.
 * \ingroup FvmAlgos
 */
template<class Type = su2double>
struct LimiterHelpers
{
  FORCEINLINE static Type epsilon() {return std::numeric_limits<passivedouble>::epsilon();}

  FORCEINLINE static Type venkatFunction(const Type& proj, const Type& delta, const Type& eps2)
  {
    Type y = delta*(delta+proj) + eps2;
    return (y + delta*proj) / (y + 2*proj*proj);
  }

  FORCEINLINE static Type vanAlbadaFunction(const Type& proj, const Type& delta, const Type& eps)
  {
    return delta * (2*proj + delta) / (4*pow(proj, 2) + pow(delta, 2) + eps);
  }

  FORCEINLINE static Type raisedSine(const Type& dist)
  {
    Type factor = 0.5*(1.0+dist+sin(PI_NUMBER*dist)/PI_NUMBER);
    return max(0.0, min(factor, 1.0));
  }

  FORCEINLINE static Type r3Function(const Type& proj, const Type& delta, const Type& epsp)
  {
    Type Dp = fabs(delta);
    Type Dm = fabs(proj);
    if(Dp>(2.0*Dm)) return 1.0;
    Type y = pow(Dp, 3) + epsp;
    Type S3 = 4.0*Dm*Dm;
    return (y + Dp*S3) / (y + Dm*(delta*delta+S3));
  }

  FORCEINLINE static Type r4Function(const Type& proj, const Type& delta, const Type& epsp)
  {
    Type Dp = fabs(delta);
    Type Dm = fabs(proj);
    if(Dp>(2.0*Dm)) return 1.0;
    Type y = pow(Dp, 4) + epsp;
    Type S4 = 2.0*Dm*(Dp*Dp-2.0*Dm*(Dp-2.0*Dm));
    return (y + Dp*S4) / (y + Dm*(pow(delta,3)+S4));
  }

  FORCEINLINE static Type r5Function(const Type& proj, const Type& delta, const Type& epsp)
  {
    Type Dp = fabs(delta);
    Type Dm = fabs(proj);
    if(Dp>(2.0*Dm)) return 1.0;
    Type y = pow(Dp, 5) + epsp;
    Type S5 = 8.0*Dm*Dm*(Dp*Dp-2.0*Dm*(Dp-Dm));
    return (y + Dp*S5) / (y + Dm*(pow(delta,4)+S5));
  }

  FORCEINLINE static Type pipernoPhiFunction(const Type& R)
  {
    /*--- φ(R) = 0 if r <= 0 ---*/
    if (R <= 0.0) return 0.0;

    /*--- φ(R) = 1 + (3/2 r + 1)(r - 1)^3 if 0 <= r <= 1 ---*/
    if (R <= 1.0) {
      Type r_minus_1 = R - 1.0;
      return 1.0 + (1.5 * R + 1.0) * pow(r_minus_1, 3);
    }

    /*--- φ(R) = (3r^2 - 6r + 19) / (r^3 - 3r + 18) if 1 <= r ---*/
    Type r_squared = R * R;
    Type r_cubed = r_squared * R;
    Type numerator = 3.0 * r_squared - 6.0 * R + 19.0;
    Type denominator = r_cubed - 3.0 * R + 18.0;
    return numerator / max(denominator, std::numeric_limits<passivedouble>::epsilon());
  }

  FORCEINLINE static Type pipernoFunction(const Type& proj, const Type& delta, const Type& eps)
  {
    /*--- Avoid division by zero ---*/
    Type delta_safe = max(fabs(delta), eps);
    Type proj_safe = max(fabs(proj), eps);

    /*----------------------------------------------------------------*/
    /*--- In Piperno's notation:                                   ---*/
    /*---   proj = ∇u_i·Δx (gradient projection)                   ---*/
    /*---   delta = Δu_{i+1/2} (centered difference)               ---*/
    /*---   Δu_{i-1/2} = 2∇u_i·Δx - Δu_{i+1/2}                     ---*/
    /*---   R = Δu_{i+1/2} / Δu_{i-1/2} = delta / (2*proj - delta) ---*/
    /*----------------------------------------------------------------*/
    Type delta_upwind = 2.0 * proj_safe - delta_safe;
    Type R = delta_safe / max(fabs(delta_upwind), eps) *
             ((delta_upwind >= 0.0) ? 1.0 : -1.0);

    /*--- Compute inverse R for φ function ---*/
    Type inv_R = 1.0 / max(fabs(R), eps) * ((R >= 0.0) ? 1.0 : -1.0);

    /*--- Compute φ(1/R) ---*/
    Type phi_inv_R = pipernoPhiFunction(inv_R);

    /*--- Compute ψ(R) = (1/3 + 2/3 R) φ(1/R) ---*/
    /*--- When φ(R) = 1, this gives the beta scheme ---*/
    Type psi_R = (ONE3 + TWO3 * R) * phi_inv_R;

    /*--- Apply sign correction ---*/
    Type sign = (delta * proj >= 0.0) ? 1.0 : 0.0;

    return sign * psi_R;
  }
};


/*!
 * \brief Barth-Jespersen specialization.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::BARTH_JESPERSEN>
{
  su2double eps2;

  /*!
   * \brief Set a small epsilon to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(Ts&...) {eps2 = LimiterHelpers<>::epsilon();}

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Venkatakrishnan function with a numerical epsilon.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::venkatFunction(proj, delta, eps2);
  }
};


/*!
 * \brief Venkatakrishnan specialization.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::VENKATAKRISHNAN>
{
  su2double eps2;

  /*!
   * \brief Store the reference lenght based eps^2 parameter,
   *        limited to a small number to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, const CConfig& config, Ts&...)
  {
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    su2double eps1 = fabs(L*K);
    eps2 = max(eps1*eps1*eps1, LimiterHelpers<>::epsilon());
  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::venkatFunction(proj, delta, eps2);
  }
};


/*!
 * \brief Nishikawa's R3 limiter specialization.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::NISHIKAWA_R3>
{
  su2double epsp;

  /*!
   * \brief Store the reference lenght based eps^3 parameter,
   *        limited to a small number to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, const CConfig& config, Ts&...)
  {
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    su2double eps1 = fabs(L*K);
    epsp = max(pow(eps1, 4), LimiterHelpers<>::epsilon());
  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::r3Function(proj, delta, epsp);
  }
};


/*!
 * \brief Nishikawa's R4 limiter specialization.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::NISHIKAWA_R4>
{
  su2double epsp;

  /*!
   * \brief Store the reference lenght based eps^4 parameter,
   *        limited to a small number to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, const CConfig& config, Ts&...)
  {
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    su2double eps1 = fabs(L*K);
    epsp = max(pow(eps1, 5), LimiterHelpers<>::epsilon());
  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::r4Function(proj, delta, epsp);
  }
};


/*!
 * \brief Nishikawa's R5 limiter specialization.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::NISHIKAWA_R5>
{
  su2double epsp;

  /*!
   * \brief Store the reference lenght based eps^5 parameter,
   *        limited to a small number to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, const CConfig& config, Ts&...)
  {
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    su2double eps1 = fabs(L*K);
    epsp = max(pow(eps1, 6), LimiterHelpers<>::epsilon());
  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::r5Function(proj, delta, epsp);
  }
};

/*!
 * \brief Venkatakrishnan-Wang specialization.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::VENKATAKRISHNAN_WANG>
{
  static su2activevector sharedMin, sharedMax;
  su2activevector eps2;

  /*!
   * \brief Store the solution range based eps^2 parameter.
   */
  template<class FieldType>
  inline void preprocess(CGeometry& geometry, const CConfig& config, size_t varBegin,
                         size_t varEnd, const FieldType& field)
  {
    /*--- Determine the max and min global value for each variable. ---*/

    su2double largeNum = 0.1*std::numeric_limits<passivedouble>::max();

    /*--- Allocate the static members (shared between threads) to
     * perform the reduction across all threads in the rank. ---*/

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      sharedMin.resize(varEnd) = largeNum;
      sharedMax.resize(varEnd) =-largeNum;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    /*--- Per thread reduction. ---*/

    su2activevector localMin(varEnd), localMax(varEnd);
    localMin = largeNum;
    localMax =-largeNum;

    SU2_OMP_FOR_(schedule(static, 512) SU2_NOWAIT)
    for(size_t iPoint = 0; iPoint < geometry.GetnPointDomain(); ++iPoint)
    {
      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        localMin(iVar) = min(localMin(iVar), field(iPoint, iVar));
        localMax(iVar) = max(localMax(iVar), field(iPoint, iVar));
      }
    }
    END_SU2_OMP_FOR

    /*--- Per rank reduction. ---*/

    SU2_OMP_CRITICAL
    for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
    {
      sharedMin(iVar) = min(sharedMin(iVar), localMin(iVar));
      sharedMax(iVar) = max(sharedMax(iVar), localMax(iVar));
    }
    END_SU2_OMP_CRITICAL

    /*--- Global reduction. ---*/

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      localMin = sharedMin;
      SU2_MPI::Allreduce(localMin.data(), sharedMin.data(), varEnd, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());

      localMax = sharedMax;
      SU2_MPI::Allreduce(localMax.data(), sharedMax.data(), varEnd, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    /*--- Compute eps^2 (each thread has its own copy of it). ---*/

    eps2.resize(varEnd);
    su2double K = config.GetVenkat_LimiterCoeff();

    for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
    {
      su2double range = sharedMax(iVar) - sharedMin(iVar);
      eps2(iVar) = max(pow(K*range, 2), LimiterHelpers<>::epsilon());
    }
  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t iVar, su2double proj, su2double delta) const
  {
    AD::SetPreaccIn(eps2(iVar));
    return LimiterHelpers<>::venkatFunction(proj, delta, eps2(iVar));
  }
};


/*!
 * \brief Venkatakrishnan with sharp edge modification.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::SHARP_EDGES>
{
  su2double eps1, eps2, sharpCoeff;

  /*!
   * \brief Store the reference lenght based eps^2 parameter.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, const CConfig& config, Ts&...)
  {
    sharpCoeff = config.GetAdjSharp_LimiterCoeff();
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    eps1 = fabs(L*K);
    eps2 = max(eps1*eps1*eps1, LimiterHelpers<>::epsilon());
  }

  /*!
   * \brief Full limiting (1st order) near sharp edges.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry& geometry) const
  {
    AD::SetPreaccIn(geometry.nodes->GetSharpEdge_Distance(iPoint));
    su2double dist = geometry.nodes->GetSharpEdge_Distance(iPoint)/(sharpCoeff*eps1)-1.0;
    return LimiterHelpers<>::raisedSine(dist);
  }

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::venkatFunction(proj, delta, eps2);
  }
};


/*!
 * \brief Venkatakrishnan with wall distance modification.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::WALL_DISTANCE>
{
  su2double eps1, eps2, sharpCoeff;

  /*!
   * \brief Store the reference lenght based eps^2 parameter.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, const CConfig& config, Ts&...)
  {
    sharpCoeff = config.GetAdjSharp_LimiterCoeff();
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    eps1 = fabs(L*K);
    eps2 = max(eps1*eps1*eps1, LimiterHelpers<>::epsilon());
  }

  /*!
   * \brief Full limiting (1st order) near walls.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry& geometry) const
  {
    AD::SetPreaccIn(geometry.nodes->GetWall_Distance(iPoint));
    su2double dist = geometry.nodes->GetWall_Distance(iPoint)/(sharpCoeff*eps1)-1.0;
    return LimiterHelpers<>::raisedSine(dist);
  }

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::venkatFunction(proj, delta, eps2);
  }
};

/*!
 * \brief Piperno limiter specialization.
 * \ingroup FvmAlgos
 */
template<>
struct CLimiterDetails<LIMITER::PIPERNO>
{
  su2double eps;

  /*!
   * \brief Set a small epsilon to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, const CConfig&, Ts&...)
  {
    eps = LimiterHelpers<>::epsilon();
  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Piperno limiter function.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers<>::pipernoFunction(proj, delta, eps);
  }
};
