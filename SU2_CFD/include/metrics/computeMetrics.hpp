/*!
 * \file computeMetric.hpp
 * \brief Generic implementation of the metric tensor computation.
 * \note This allows the same implementation to be used for goal-oriented
 *       or feature-based mesh adaptation.
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

#include <vector>
#include <algorithm>
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"

namespace tensor {
  struct metric {
    template <class MatrixType, class MetricType>
    static void get(MetricType& metric_field, const unsigned long iPoint,
                    const unsigned short iSensor, MatrixType& mat, unsigned short nDim) {
      switch(nDim) {
        case 2: {
          mat[0][0] = metric_field(iPoint,0); mat[0][1] = metric_field(iPoint,1);
          mat[1][0] = metric_field(iPoint,1); mat[1][1] = metric_field(iPoint,2);
          break;
        }
        case 3: {
          mat[0][0] = metric_field(iPoint,0); mat[0][1] = metric_field(iPoint,1); mat[0][2] = metric_field(iPoint,2);
          mat[1][0] = metric_field(iPoint,1); mat[1][1] = metric_field(iPoint,3); mat[1][2] = metric_field(iPoint,4);
          mat[2][0] = metric_field(iPoint,2); mat[2][1] = metric_field(iPoint,4); mat[2][2] = metric_field(iPoint,5);
          break;
        }
      }
    }

    template <class ScalarType, class MatrixType, class MetricType>
    static void set(MetricType& metric_field, const unsigned long iPoint,
                    const unsigned short iSensor, MatrixType& mat, ScalarType scale,
                    unsigned short nDim) {
      switch(nDim) {
        case 2: {
          metric_field(iPoint,0) = mat[0][0] * scale;
          metric_field(iPoint,1) = mat[0][1] * scale;
          metric_field(iPoint,2) = mat[1][1] * scale;
          break;
        }
        case 3: {
          metric_field(iPoint,0) = mat[0][0] * scale;
          metric_field(iPoint,1) = mat[0][1] * scale;
          metric_field(iPoint,2) = mat[0][2] * scale;
          metric_field(iPoint,3) = mat[1][1] * scale;
          metric_field(iPoint,4) = mat[1][2] * scale;
          metric_field(iPoint,5) = mat[2][2] * scale;
          break;
        }
      }
    }
  };

  struct hessian {
    template <class MatrixType, class MetricType>
    static void get(MetricType& metric_field, const unsigned long iPoint,
                    const unsigned short iSensor, MatrixType& mat, unsigned short nDim) {
      switch(nDim) {
        case 2: {
          mat[0][0] = metric_field(iPoint,iSensor,0); mat[0][1] = metric_field(iPoint,iSensor,1);
          mat[1][0] = metric_field(iPoint,iSensor,1); mat[1][1] = metric_field(iPoint,iSensor,2);
          break;
        }
        case 3: {
          mat[0][0] = metric_field(iPoint,iSensor,0); mat[0][1] = metric_field(iPoint,iSensor,1); mat[0][2] = metric_field(iPoint,iSensor,2);
          mat[1][0] = metric_field(iPoint,iSensor,1); mat[1][1] = metric_field(iPoint,iSensor,3); mat[1][2] = metric_field(iPoint,iSensor,4);
          mat[2][0] = metric_field(iPoint,iSensor,2); mat[2][1] = metric_field(iPoint,iSensor,4); mat[2][2] = metric_field(iPoint,iSensor,5);
          break;
        }
      }
    }

    template <class ScalarType, class MatrixType, class MetricType>
    static void set(MetricType& metric_field, const unsigned long iPoint,
                    const unsigned short iSensor, MatrixType& mat, ScalarType scale,
                    unsigned short nDim) {
      switch(nDim) {
        case 2: {
          metric_field(iPoint,iSensor,0) = mat[0][0] * scale;
          metric_field(iPoint,iSensor,1) = mat[0][1] * scale;
          metric_field(iPoint,iSensor,2) = mat[1][1] * scale;
          break;
        }
        case 3: {
          metric_field(iPoint,iSensor,0) = mat[0][0] * scale;
          metric_field(iPoint,iSensor,1) = mat[0][1] * scale;
          metric_field(iPoint,iSensor,2) = mat[0][2] * scale;
          metric_field(iPoint,iSensor,3) = mat[1][1] * scale;
          metric_field(iPoint,iSensor,4) = mat[1][2] * scale;
          metric_field(iPoint,iSensor,5) = mat[2][2] * scale;
          break;
        }
      }
    }
  };
}

namespace detail {

/*!
 * \brief Make the eigenvalues of the metrics positive.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iSensor - Index of the sensor to work on.
 * \param[in,out] metric - Metric container.
 */
template<size_t nDim, class ScalarType, class Tensor, class MetricType>
void setPositiveDefiniteMetrics(CGeometry& geometry, const CConfig& config,
                               unsigned short iSensor, MetricType& metric) {

  const unsigned long nPointDomain = geometry.GetnPointDomain();

  static constexpr size_t MAXNDIM = 3;

  ScalarType A[MAXNDIM][MAXNDIM], EigVec[MAXNDIM][MAXNDIM], EigVal[MAXNDIM], work[MAXNDIM];

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    /*--- Get full metric tensor ---*/
    Tensor::get(metric, iPoint, iSensor, A, nDim);

    /*--- Compute eigenvalues and eigenvectors ---*/
    CBlasStructure::EigenDecomposition(A, EigVec, EigVal, nDim, work);

    /*--- If NaN detected, set values to zero ---*/
    /*--- Otherwise, store recombined matrix  ---*/
    bool check_hess = true;
    for (auto iDim = 0; iDim < nDim; iDim++) {
      if (EigVal[iDim] != EigVal[iDim] || fabs(EigVal[iDim]) < 1e-16) {
        EigVal[iDim] = 1e-16;
        check_hess = false;
      }
      EigVal[iDim] = fabs(EigVal[iDim]);
    }

    CBlasStructure::EigenRecomposition(A, EigVec, EigVal, nDim);

    /*--- Store upper half of metric tensor ---*/
    Tensor::set(metric, iPoint, iSensor, A, 1.0, nDim);
  }
}

/*!
 * \brief Integrate the Hessian field for the Lp-norm normalization of the metric.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iSensor - Index of the sensor to work on.
 * \param[in] metric - Metric container.
 * \return Integral of the metric tensor determinant.
*/
template<size_t nDim, class ScalarType, class Tensor, class MetricType>
ScalarType integrateMetrics(CGeometry& geometry, const CConfig& config,
                            unsigned short iSensor, MetricType& metric) {

  const unsigned long nPointDomain = geometry.GetnPointDomain();
  static constexpr size_t MAXNDIM = 3;

  /*--- Constants defining normalization ---*/
  const ScalarType p = config.GetMetric_Norm();
  const ScalarType normExp = p / (2.0 * p + nDim);

  ScalarType A[MAXNDIM][MAXNDIM], EigVec[MAXNDIM][MAXNDIM], EigVal[MAXNDIM], work[MAXNDIM];

  ScalarType localIntegral = 0.0;
  ScalarType globalIntegral = 0.0;
  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    auto nodes = geometry.nodes;

    /*--- Decompose metric ---*/
    Tensor::get(metric, iPoint, iSensor, A, nDim);
    CBlasStructure::EigenDecomposition(A, EigVec, EigVal, nDim, work);
    if (nDim == 2) EigVal[2] = 1.0;

    /*--- Integrate determinant ---*/
    const ScalarType det = EigVal[0] * EigVal[1] * EigVal[2];
    const ScalarType Vol = SU2_TYPE::GetValue(nodes->GetVolume(iPoint));
    localIntegral += pow(abs(det), normExp) * Vol;
  }

  CBaseMPIWrapper::Allreduce(&localIntegral, &globalIntegral, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  return globalIntegral;
}

/*!
 * \brief Perform an Lp-norm normalization of the metric.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iSensor - Index of the sensor to work on.
 * \param[in] integral - Integral of the metric tensor determinant.
 * \param[in,out] metric - Metric container.
 */
template<size_t nDim, class ScalarType, class Tensor, class MetricType>
void normalizeMetrics(CGeometry& geometry, const CConfig& config,
                      unsigned short iSensor, ScalarType integral,
                      MetricType& metric) {

  const unsigned long nPointDomain = geometry.GetnPointDomain();
  static constexpr size_t MAXNDIM = 3;

  const bool goal = (config.GetGoal_Oriented_Metric());

  /*--- Constants defining normalization ---*/
  const ScalarType p = config.GetMetric_Norm();
  const ScalarType N = ScalarType(config.GetMetric_Complexity() * config.GetnAdapt_Time_Subinterval());
  const ScalarType globalFactor = pow(N / integral, 2.0 / nDim);
  const ScalarType normExp = -1.0 / (2.0 * p + nDim);

  /*--- Size constraints ---*/
  const ScalarType hmin = SU2_TYPE::GetValue(config.GetMetric_Hmin());
  const ScalarType hmax = SU2_TYPE::GetValue(config.GetMetric_Hmax());
  const ScalarType eigmax = 1.0 / pow(hmin, 2.0);
  const ScalarType eigmin = 1.0 / pow(hmax, 2.0);
  const ScalarType armax2 = pow(SU2_TYPE::GetValue(config.GetMetric_ARmax()), 2.0);

  ScalarType A[MAXNDIM][MAXNDIM], EigVec[MAXNDIM][MAXNDIM], EigVal[MAXNDIM], work[MAXNDIM];

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    auto nodes = geometry.nodes;

    /*--- Decompose metric ---*/
    Tensor::get(metric, iPoint, iSensor, A, nDim);
    CBlasStructure::EigenDecomposition(A, EigVec, EigVal, nDim, work);
    if (nDim == 2) EigVal[2] = 1.0;

    /*--- Normalize eigenvalues ---*/
    const ScalarType det = EigVal[0] * EigVal[1] * EigVal[2];
    const ScalarType factor = globalFactor * pow(abs(det), normExp);
    for (auto iDim = 0u; iDim < nDim; ++iDim)
      EigVal[iDim] = factor * EigVal[iDim];

    /*--- Clip by user-specified size constraints ---*/
    for (auto iDim = 0u; iDim < nDim; ++iDim)
      EigVal[iDim] = min(max(abs(EigVal[iDim]), eigmin), eigmax);

    // /*--- Clip by user-specified aspect ratio ---*/
    // unsigned short iMax = 0;
    // for (auto iDim = 1; iDim < nDim; ++iDim)
    //   iMax = (EigVal[iDim] > EigVal[iMax])? iDim : iMax;

    // for (auto iDim = 0u; iDim < nDim; ++iDim)
    //   EigVal[iDim] = max(EigVal[iDim], EigVal[iMax]/armax2);

    /*--- Recompose and store metric ---*/
    CBlasStructure::EigenRecomposition(A, EigVec, EigVal, nDim);
    Tensor::set(metric, iPoint, iSensor, A, 1.0, nDim);
  }
}

} // end namespace detail

/*!
 * \brief Make the eigenvalues of the metrics positive.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iSensor - Index of the sensor to work on.
 * \param[in,out] metric - Metric container.
 */
template<class ScalarType, class Tensor, class MetricType>
void setPositiveDefiniteMetrics(CGeometry& geometry, const CConfig& config,
                               unsigned short iSensor, MetricType& metric) {
  switch (geometry.GetnDim()) {
    case 2:
      detail::setPositiveDefiniteMetrics<2, ScalarType, Tensor>(geometry, config, iSensor, metric);
      break;
    case 3:
      detail::setPositiveDefiniteMetrics<3, ScalarType, Tensor>(geometry, config, iSensor, metric);
      break;
    default:
      SU2_MPI::Error("Too many dimensions for metric computation.", CURRENT_FUNCTION);
      break;
  }
}

/*!
 * \brief Integrate the Hessian field for the Lp-norm normalization of the metric.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iSensor - Index of the sensor to work on.
 * \param[in] metric - Metric container.
 * \return Integral of the metric tensor determinant.
*/
template<class ScalarType, class Tensor, class MetricType>
ScalarType integrateMetrics(CGeometry& geometry, const CConfig& config,
                            unsigned short iSensor, MetricType& metric) {
  su2double integral;
  switch (geometry.GetnDim()) {
    case 2:
      integral = detail::integrateMetrics<2, ScalarType, Tensor>(geometry, config, iSensor, metric);
      break;
    case 3:
      integral = detail::integrateMetrics<3, ScalarType, Tensor>(geometry, config, iSensor, metric);
      break;
    default:
      SU2_MPI::Error("Too many dimensions for metric integration.", CURRENT_FUNCTION);
      break;
  }

  return integral;
}

/*!
 * \brief Perform an Lp-norm normalization of the metric.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iSensor - Index of the sensor to work on.
 * \param[in] integral - Integral of the metric tensor determinant.
 * \param[in,out] metric - Metric container.
 */
template<class ScalarType, class Tensor, class MetricType>
void normalizeMetrics(CGeometry& geometry, const CConfig& config,
                     unsigned short iSensor, ScalarType integral,
                     MetricType& metric) {
  switch (geometry.GetnDim()) {
    case 2:
      detail::normalizeMetrics<2, ScalarType, Tensor>(geometry, config, iSensor, integral, metric);
      break;
    case 3:
      detail::normalizeMetrics<3, ScalarType, Tensor>(geometry, config, iSensor, integral, metric);
      break;
    default:
      SU2_MPI::Error("Too many dimensions for metric normalization.", CURRENT_FUNCTION);
      break;
  }
}