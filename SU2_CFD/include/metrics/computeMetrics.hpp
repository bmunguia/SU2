/*!
 * \file computeMetrics.hpp
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
#include <limits>
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

namespace tensor {
  struct metric {
    template <class MatrixType, class MetricType>
    static void get(MetricType& metric_field, const unsigned long iPoint,
                    const unsigned short iSensor, MatrixType& mat, unsigned short nDim) {
      switch(nDim) {
        case 2: {
          mat[0][0] = metric_field(iPoint, 0); mat[0][1] = metric_field(iPoint, 1);
          mat[1][0] = metric_field(iPoint, 1); mat[1][1] = metric_field(iPoint, 2);
          break;
        }
        case 3: {
          mat[0][0] = metric_field(iPoint, 0); mat[0][1] = metric_field(iPoint, 1); mat[0][2] = metric_field(iPoint, 2);
          mat[1][0] = metric_field(iPoint, 1); mat[1][1] = metric_field(iPoint, 3); mat[1][2] = metric_field(iPoint, 4);
          mat[2][0] = metric_field(iPoint, 2); mat[2][1] = metric_field(iPoint, 4); mat[2][2] = metric_field(iPoint, 5);
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
          metric_field(iPoint, 0) = mat[0][0] * scale;
          metric_field(iPoint, 1) = mat[0][1] * scale;
          metric_field(iPoint, 2) = mat[1][1] * scale;
          break;
        }
        case 3: {
          metric_field(iPoint, 0) = mat[0][0] * scale;
          metric_field(iPoint, 1) = mat[0][1] * scale;
          metric_field(iPoint, 2) = mat[0][2] * scale;
          metric_field(iPoint, 3) = mat[1][1] * scale;
          metric_field(iPoint, 4) = mat[1][2] * scale;
          metric_field(iPoint, 5) = mat[2][2] * scale;
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
          mat[0][0] = metric_field(iPoint, iSensor, 0); mat[0][1] = metric_field(iPoint, iSensor, 1);
          mat[1][0] = metric_field(iPoint, iSensor, 1); mat[1][1] = metric_field(iPoint, iSensor, 2);
          break;
        }
        case 3: {
          mat[0][0] = metric_field(iPoint, iSensor, 0); mat[0][1] = metric_field(iPoint, iSensor, 1); mat[0][2] = metric_field(iPoint, iSensor, 2);
          mat[1][0] = metric_field(iPoint, iSensor, 1); mat[1][1] = metric_field(iPoint, iSensor, 3); mat[1][2] = metric_field(iPoint, iSensor, 4);
          mat[2][0] = metric_field(iPoint, iSensor, 2); mat[2][1] = metric_field(iPoint, iSensor, 4); mat[2][2] = metric_field(iPoint, iSensor, 5);
          break;
        }
      }
    }

    template <class ScalarType, class MatrixType, class MetricType>
    static void set(MetricType& metric_field, const unsigned long iPoint,
                    const unsigned short iSensor,  MatrixType& mat, ScalarType scale,
                    unsigned short nDim) {
      switch(nDim) {
        case 2: {
          metric_field(iPoint, iSensor, 0) = mat[0][0] * scale;
          metric_field(iPoint, iSensor, 1) = mat[0][1] * scale;
          metric_field(iPoint, iSensor, 2) = mat[1][1] * scale;
          break;
        }
        case 3: {
          metric_field(iPoint, iSensor, 0) = mat[0][0] * scale;
          metric_field(iPoint, iSensor, 1) = mat[0][1] * scale;
          metric_field(iPoint, iSensor, 2) = mat[0][2] * scale;
          metric_field(iPoint, iSensor, 3) = mat[1][1] * scale;
          metric_field(iPoint, iSensor, 4) = mat[1][2] * scale;
          metric_field(iPoint, iSensor, 5) = mat[2][2] * scale;
          break;
        }
      }
    }
  };
}

namespace detail {

/*!
 * \brief Compute determinant of eigenvalues for different dimensions.
 * \param[in] EigVal - Array of eigenvalues.
 * \return Determinant value.
 */
template<size_t nDim, class ScalarType>
ScalarType computeDeterminant(const ScalarType* EigVal) {
  if constexpr (nDim == 2) {
    return EigVal[0] * EigVal[1];
  } else if constexpr (nDim == 3) {
    return EigVal[0] * EigVal[1] * EigVal[2];
  } else {
    static_assert(nDim == 2 || nDim == 3, "Only 2D and 3D supported");
    return ScalarType(0.0);
  }
}

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

  ScalarType A[nDim][nDim], EigVec[nDim][nDim], EigVal[nDim], work[nDim];

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

  /*--- Constants defining normalization ---*/
  const ScalarType p = config.GetMetric_Norm();
  const ScalarType normExp = p / (2.0 * p + nDim);

  ScalarType A[nDim][nDim], EigVec[nDim][nDim], EigVal[nDim], work[nDim];

  ScalarType localIntegral = 0.0;
  ScalarType globalIntegral = 0.0;
  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    auto nodes = geometry.nodes;

    /*--- Decompose metric ---*/
    Tensor::get(metric, iPoint, iSensor, A, nDim);
    CBlasStructure::EigenDecomposition(A, EigVec, EigVal, nDim, work);

    /*--- Integrate determinant ---*/
    const ScalarType det = computeDeterminant<nDim>(EigVal);
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

  ScalarType A[nDim][nDim], EigVec[nDim][nDim], EigVal[nDim], work[nDim];

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    auto nodes = geometry.nodes;

    /*--- Decompose metric ---*/
    Tensor::get(metric, iPoint, iSensor, A, nDim);
    CBlasStructure::EigenDecomposition(A, EigVec, EigVal, nDim, work);

    /*--- Normalize eigenvalues ---*/
    const ScalarType det = computeDeterminant<nDim>(EigVal);
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

/*!
 * \brief Compute a metric tensor based on surface geometry.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[out] metric - The computed geometric metric tensor field.
 */
template<size_t nDim, class ScalarType, class Tensor, class MetricType>
void geometricSurfaceMetrics(CGeometry& geometry, const CConfig& config,
                             MetricType& metric) {

  const unsigned long nPointDomain = geometry.GetnPointDomain();

  /*--- Size constraints ---*/
  const ScalarType hmin = SU2_TYPE::GetValue(config.GetMetric_Hmin());
  const ScalarType hmax = SU2_TYPE::GetValue(config.GetMetric_Hmax());
  const ScalarType eigmax = 1.0 / pow(hmin, 2.0);
  const ScalarType eigmin = 1.0 / pow(hmax, 2.0);

  /*--- Constraint on deviation from tangent plane ---*/;
  const ScalarType deg2rad = M_PI / 180.0;

  /*--- Working arrays ---*/
  ScalarType M[nDim][nDim], R[nDim][nDim], EigVal[nDim], work[nDim];

  /*--- Initialize isotropic metric ---*/
  for(size_t i = 0; i < nDim; ++i) {
    for(size_t j = 0; j < nDim; ++j) {
      M[i][j] = (i == j) ? eigmin : 0.0;
    }
  }
  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    Tensor::set(metric, iPoint, 0, M, 1.0, nDim);
  }

  auto nodes = geometry.nodes;

  if constexpr (nDim == 2) {
    /*--- 2D: Use curvature computed by CGeometry::ComputeSurf_Curvature ---*/
    for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
      if (!nodes->GetPhysicalBoundary(iPoint)) continue;

      /*--- Find the minimum GeoDev value for this point across all applicable markers ---*/
      ScalarType alpha = std::numeric_limits<ScalarType>::max();
      ScalarType n[2] = {};
      bool foundMarker = false;

      for (unsigned short iMarkerGeoDev = 0; iMarkerGeoDev < config.GetnMarker_GeoDev(); ++iMarkerGeoDev) {
        const string& markerTag = config.GetMarker_GeoDev(iMarkerGeoDev);

        /*--- Find the geometry marker index corresponding to this tag ---*/
        for (unsigned short iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
          if (geometry.GetMarker_Tag(iMarker) == markerTag) {
            /*--- Check if this point belongs to this marker ---*/
            const auto iVertex = nodes->GetVertex(iPoint, iMarker);
            if (iVertex >= 0) {
              /*--- Point belongs to this marker, get the GeoDev value ---*/
              const ScalarType geodev_deg = SU2_TYPE::GetValue(config.GetMetric_GeoDev(iMarkerGeoDev));
              const ScalarType geodev_rad = geodev_deg * deg2rad;
              alpha = min(alpha, geodev_rad);
              foundMarker = true;

              /*--- Also add to the normal ---*/
              const auto* normal = geometry.vertex[iMarker][iVertex]->GetNormal();
              for (auto iDim = 0u; iDim < 2; ++iDim) n[iDim] += SU2_TYPE::GetValue(normal[iDim]);
            }
            break; /*--- Found the marker, no need to continue searching ---*/
          }
        }
      }

      /*--- If point doesn't belong to any GeoDev marker, skip it ---*/
      if (!foundMarker) continue;

      /*--- Get curvature directly from geometry ---*/
      const ScalarType curvature = SU2_TYPE::GetValue(nodes->GetCurvature(iPoint));

      /*--- Compute metric sizes ---*/
      const ScalarType h_curv = alpha / fabs(curvature);
      const ScalarType h_clipped = max(hmin, min(hmax, h_curv));

      /*--- Build metric eigenvalues ---*/
      EigVal[0] = 1.0 / (h_clipped * h_clipped);  // Tangent direction
      EigVal[1] = eigmin;                         // Normal direction

      /*--- Build rotation matrix ---*/
      const auto area = GeometryToolbox::Norm(2, n);
      for (auto iDim = 0u; iDim < 2; ++iDim) n[iDim] /= area;

      R[0][0] = -n[1]; R[0][1] = n[0];   // Tangent direction
      R[1][0] = n[0];  R[1][1] = n[1];   // Normal direction

      /*--- Recompose metric ---*/
      CBlasStructure::EigenRecomposition(M, R, EigVal, 2);
      Tensor::set(metric, iPoint, 0, M, 1.0, 2);
    }

  } else {
    /*--- 3D: Use methodology from Meyer, "Discrete Differential-Geometry Operators for ---*/
    /*--- Triangulated 2-Manifolds," 2003 for principal curvature computation ---*/
    /*--- Allocate arrays for curvature computation ---*/
    vector<ScalarType> angleDefect(geometry.GetnPoint(), 2.0 * M_PI);
    vector<ScalarType> areaVertex(geometry.GetnPoint(), 0.0);
    vector<ScalarType> angleAlpha(geometry.GetnEdge(), 0.0);
    vector<ScalarType> angleBeta(geometry.GetnEdge(), 0.0);
    vector<bool> checkEdge(geometry.GetnEdge(), true);

    /*--- Mean curvature vectors K(xi) = (1/2AMixed) * sum over 1-ring ---*/
    vector<vector<ScalarType>> meanCurvatureVector(geometry.GetnPoint(), vector<ScalarType>(3, 0.0));

    /*--- Data for least-squares curvature tensor fitting per vertex ---*/
    vector<vector<ScalarType>> tangentDirs_u(geometry.GetnPoint());  // u-components in tangent plane
    vector<vector<ScalarType>> tangentDirs_v(geometry.GetnPoint());  // v-components in tangent plane
    vector<vector<ScalarType>> normalCurvatures(geometry.GetnPoint()); // Kappaij_N values
    vector<vector<ScalarType>> weights(geometry.GetnPoint());         // cotangent weights

    /*--- Compute angle defects and vertex areas ---*/
    for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
      if (config.GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;

      for (size_t iElem = 0; iElem < geometry.GetnElem_Bound(iMarker); ++iElem) {
        if (geometry.bound[iMarker][iElem]->GetVTK_Type() != TRIANGLE) continue;

        for (size_t iNode = 0; iNode < 3; ++iNode) {
          const auto iPoint = geometry.bound[iMarker][iElem]->GetNode(iNode);

          /*--- Get triangle neighbors ---*/
          vector<unsigned long> trianglePoints;
          for (size_t jNode = 0; jNode < geometry.bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); ++jNode) {
            auto neighborNode = geometry.bound[iMarker][iElem]->GetNeighbor_Nodes(iNode, jNode);
            trianglePoints.push_back(geometry.bound[iMarker][iElem]->GetNode(neighborNode));
          }

          if (trianglePoints.size() == 2) {
            /*--- Compute vectors U, V and their cross product W ---*/
            ScalarType U[3], V[3], W[3];
            for (size_t iDim = 0; iDim < 3; ++iDim) {
              U[iDim] = nodes->GetCoord(trianglePoints[0], iDim) - nodes->GetCoord(iPoint, iDim);
              V[iDim] = nodes->GetCoord(trianglePoints[1], iDim) - nodes->GetCoord(iPoint, iDim);
            }

            W[0] = 0.5 * (U[1] * V[2] - U[2] * V[1]);
            W[1] = -0.5 * (U[0] * V[2] - U[2] * V[0]);
            W[2] = 0.5 * (U[0] * V[1] - U[1] * V[0]);

            /*--- Normalize and compute angle ---*/
            ScalarType lengthU = sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
            ScalarType lengthV = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
            ScalarType lengthW = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2]);

            if (lengthU > 1e-16 && lengthV > 1e-16) {
              for (size_t iDim = 0; iDim < 3; ++iDim) {
                U[iDim] /= lengthU;
                V[iDim] /= lengthV;
              }

              ScalarType cosValue = max(-1.0, min(1.0, U[0]*V[0] + U[1]*V[1] + U[2]*V[2]));
              ScalarType angleValue = acos(cosValue);

              areaVertex[iPoint] += lengthW;
              angleDefect[iPoint] -= angleValue;

              /*--- Store angles for edge processing ---*/
              auto iEdge = geometry.FindEdge(trianglePoints[0], trianglePoints[1]);
              if (angleAlpha[iEdge] == 0.0) {
                angleAlpha[iEdge] = angleValue;
              } else {
                angleBeta[iEdge] = angleValue;
              }
            }
          }
        }
      }
    }

    /*--- Compute mean curvature vectors AND normal curvatures for least squares ---*/
    for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
      if (config.GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;

      for (size_t iElem = 0; iElem < geometry.GetnElem_Bound(iMarker); ++iElem) {
        if (geometry.bound[iMarker][iElem]->GetVTK_Type() != TRIANGLE) continue;

        for (size_t iNode = 0; iNode < 3; ++iNode) {
          const auto iPoint = geometry.bound[iMarker][iElem]->GetNode(iNode);

          for (size_t jNode = 0; jNode < geometry.bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); ++jNode) {
            auto neighborNode = geometry.bound[iMarker][iElem]->GetNeighbor_Nodes(iNode, jNode);
            const auto jPoint = geometry.bound[iMarker][iElem]->GetNode(neighborNode);

            auto iEdge = geometry.FindEdge(iPoint, jPoint);

            if (checkEdge[iEdge]) {
              checkEdge[iEdge] = false;

              /*--- Compute cotangent weights ---*/
              ScalarType cotAlpha = (tan(angleAlpha[iEdge]) != 0.0) ? 1.0 / tan(angleAlpha[iEdge]) : 0.0;
              ScalarType cotBeta = (tan(angleBeta[iEdge]) != 0.0) ? 1.0 / tan(angleBeta[iEdge]) : 0.0;
              ScalarType cotSum = cotAlpha + cotBeta;

              /*--- Edge vector ---*/
              ScalarType edge[3];
              for (size_t iDim = 0; iDim < 3; ++iDim) {
                edge[iDim] = nodes->GetCoord(iPoint, iDim) - nodes->GetCoord(jPoint, iDim);
              }
              ScalarType edgeLength2 = edge[0]*edge[0] + edge[1]*edge[1] + edge[2]*edge[2];
              ScalarType edgeLength = sqrt(edgeLength2);

              /*--- Add contributions to mean curvature vector K(xi) ---*/
              /*--- K(xi) = (1/2AMixed) * sum((cot αij + cot βij) * (xi - xj)) ---*/
              for (size_t iDim = 0; iDim < 3; ++iDim) {
                if (areaVertex[iPoint] != 0.0) {
                  meanCurvatureVector[iPoint][iDim] += 0.5 * cotSum * edge[iDim] / areaVertex[iPoint];
                }
                if (areaVertex[jPoint] != 0.0) {
                  meanCurvatureVector[jPoint][iDim] += 0.5 * cotSum * (-edge[iDim]) / areaVertex[jPoint];
                }
              }
            }
          }
        }
      }
    }

    /*--- Build metric tensors from curvature information ---*/
    for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
      if (config.GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;

      /*--- Check if this marker is in the GeoDev list ---*/
      ScalarType alpha = std::numeric_limits<ScalarType>::max();
      bool foundGeoDevMarker = false;
      const string& markerTag = geometry.GetMarker_Tag(iMarker);

      for (unsigned short iMarkerGeoDev = 0; iMarkerGeoDev < config.GetnMarker_GeoDev(); ++iMarkerGeoDev) {
        if (config.GetMarker_GeoDev(iMarkerGeoDev) == markerTag) {
          /*--- This marker is in the GeoDev list ---*/
          const ScalarType geodev_deg = SU2_TYPE::GetValue(config.GetMetric_GeoDev(iMarkerGeoDev));
          const ScalarType geodev_rad = geodev_deg * deg2rad;
          alpha = min(alpha, geodev_rad);
          foundGeoDevMarker = true;
          break;
        }
      }

      /*--- If this marker is not in the GeoDev list, skip it ---*/
      if (!foundGeoDevMarker) continue;

      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
        const auto iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
        if (!nodes->GetDomain(iPoint)) continue;

        /*--- Get surface normal ---*/
        const auto* normal = geometry.vertex[iMarker][iVertex]->GetNormal();
        const auto area = GeometryToolbox::Norm(3, normal);
        ScalarType n[3];
        for (size_t i = 0; i < 3; ++i) {
          n[i] = SU2_TYPE::GetValue(normal[i] / area);
          R[i][2] = n[i];
        }

        /*--- Compute mean and Gaussian curvatures using Meyer et al. formulation ---*/
        /*--- Mean curvature: KappaH = (1/2) * ||K(xi)|| ---*/
        ScalarType meanK = 0.5 * GeometryToolbox::Norm(3, meanCurvatureVector[iPoint].data());

        /*--- Gaussian curvature: KappaG = (1/AMixed) * sum(angle defects) ---*/
        ScalarType gaussK = (areaVertex[iPoint] != 0.0) ? angleDefect[iPoint] / areaVertex[iPoint] : 0.0;

        /*--- Compute principal curvatures with safety check ---*/
        ScalarType delta = max(meanK * meanK - gaussK, 0.0);
        ScalarType kappa1 = meanK + sqrt(delta);  // Maximum principal curvature
        ScalarType kappa2 = meanK - sqrt(delta);  // Minimum principal curvature

        /*--- Compute mesh sizes ---*/
        ScalarType h1 = (fabs(kappa1) != 0.0) ? alpha / fabs(kappa1) : hmax;
        ScalarType h2 = (fabs(kappa2) != 0.0) ? alpha / fabs(kappa2) : hmax;

        h1 = max(hmin, min(hmax, h1));
        h2 = max(hmin, min(hmax, h2));

        /*----------------------------------------------------------------------------------*/
        /*--- Calculate the eigenvalues with Meyer's operators                           ---*/
        /*--- From Meyer: Although we could actually determine the principal curvatures  ---*/
        /*---     (and thus the mean and gaussian curvatures) using an unconstrained     ---*/
        /*---     least squares procedure, we use our operators to compute the           ---*/
        /*---     curvatures and only use the least squares for the principal directions ---*/
        /*---     as the curvature values computed from the least squares are often less ---*/
        /*---     accurate in practice while the directions are fairly robust.           ---*/
        /*----------------------------------------------------------------------------------*/
        EigVal[0] = 1.0 / pow(h1, 2.0);  // 1st principal direction
        EigVal[1] = 1.0 / pow(h2, 2.0);  // 2nd principal direction
        EigVal[2] = eigmin;              // Normal direction

        /*--- Now collect normal curvatures for least-squares tensor fitting ---*/
        /*--- We need to do this per vertex since we now have the surface normal ---*/
        tangentDirs_u[iPoint].clear();
        tangentDirs_v[iPoint].clear();
        normalCurvatures[iPoint].clear();
        weights[iPoint].clear();

        /*--- Build orthonormal tangent basis {u, v, n} ---*/
        ScalarType u[3], v[3];

        /*--- First tangent vector: remove normal component from coordinate axis ---*/
        size_t minComp = (fabs(n[0]) < fabs(n[1])) ? 0 : 1;
        minComp = (fabs(n[2]) < fabs(n[minComp])) ? 2 : minComp;

        for (size_t dim = 0; dim < 3; ++dim) u[dim] = 0.0;
        u[minComp] = 1.0;

        /*--- u = u - (u dot n)n ---*/
        ScalarType uDotN = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
        for (size_t dim = 0; dim < 3; ++dim) {
          u[dim] -= uDotN * n[dim];
        }

        /*--- Normalize u ---*/
        ScalarType uNorm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
        for (size_t dim = 0; dim < 3; ++dim) u[dim] /= uNorm;

        /*--- v = n cross u ---*/
        v[0] = n[1]*u[2] - n[2]*u[1];
        v[1] = n[2]*u[0] - n[0]*u[2];
        v[2] = n[0]*u[1] - n[1]*u[0];

        /*--- Find all neighboring vertices for normal curvature computation ---*/
        for (size_t jMarker = 0; jMarker < geometry.GetnMarker(); ++jMarker) {
          if (config.GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) continue;

          for (size_t jElem = 0; jElem < geometry.GetnElem_Bound(jMarker); ++jElem) {
            if (geometry.bound[jMarker][jElem]->GetVTK_Type() != TRIANGLE) continue;

            /*--- Check if iPoint is in this triangle ---*/
            bool pointInTriangle = false;
            size_t nodeIdx = 0;
            for (size_t kNode = 0; kNode < 3; ++kNode) {
              if (geometry.bound[jMarker][jElem]->GetNode(kNode) == iPoint) {
                pointInTriangle = true;
                nodeIdx = kNode;
                break;
              }
            }

            if (pointInTriangle) {
              /*--- Get the two neighboring points in this triangle ---*/
              auto jPoint1 = geometry.bound[jMarker][jElem]->GetNode((nodeIdx + 1) % 3);
              auto jPoint2 = geometry.bound[jMarker][jElem]->GetNode((nodeIdx + 2) % 3);

              /*--- Process both edges from iPoint ---*/
              for (auto jPoint : {jPoint1, jPoint2}) {
                /*--- Compute edge vector xi - xj ---*/
                ScalarType edge[3];
                for (size_t dim = 0; dim < 3; ++dim)
                  edge[dim] = nodes->GetCoord(iPoint, dim) - nodes->GetCoord(jPoint, dim);
                ScalarType edgeLength2 = GeometryToolbox::SquaredNorm(3, edge);
                if (edgeLength2 < 1e-24) continue; // Skip degenerate edges

                /*--- Normal curvature: Kappaij_N = 2(xi - xj) dot n / ||xi - xj||^2 ---*/
                ScalarType edgeDotN = GeometryToolbox::DotProduct(3, edge, n);
                ScalarType Kappaij_N = 2.0 * edgeDotN / edgeLength2;

                /*--- Project edge to tangent plane for least squares ---*/
                ScalarType tangentDir[3];
                for (size_t dim = 0; dim < 3; ++dim) {
                  tangentDir[dim] = edge[dim] - edgeDotN * n[dim];
                }

                /*--- Normalize tangent direction ---*/
                ScalarType tangentNorm = GeometryToolbox::Norm(3, tangentDir);
                if (tangentNorm < 1e-12) continue;

                for (size_t dim = 0; dim < 3; ++dim)
                  tangentDir[dim] /= tangentNorm;

                /*--- Project to 2D tangent coordinates ---*/
                ScalarType du = GeometryToolbox::DotProduct(3, tangentDir, u);
                ScalarType dv = GeometryToolbox::DotProduct(3, tangentDir, v);

                /*--- Compute cotangent weight ---*/
                auto iEdge = geometry.FindEdge(iPoint, jPoint);
                ScalarType cotAlpha = (tan(angleAlpha[iEdge]) != 0.0) ? 1.0 / tan(angleAlpha[iEdge]) : 0.0;
                ScalarType cotBeta = (tan(angleBeta[iEdge]) != 0.0) ? 1.0 / tan(angleBeta[iEdge]) : 0.0;
                ScalarType weight = (cotAlpha + cotBeta) * sqrt(edgeLength2) / (8.0 * areaVertex[iPoint]);

                /*--- Store data for least squares ---*/
                tangentDirs_u[iPoint].push_back(du);
                tangentDirs_v[iPoint].push_back(dv);
                normalCurvatures[iPoint].push_back(Kappaij_N);
                weights[iPoint].push_back(weight);
              }
            }
          }
        }

        /*--- Solve least squares problem for curvature tensor B = [[a,b],[b,c]] ---*/
        /*--- Minimize: sum_j w_j * (d_j^T B d_j - kappa_j)^2 ---*/
        /*--- where d_j^T B d_j = a*du_j^2 + 2*b*du_j*dv_j + c*dv_j^2 ---*/

        ScalarType B[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

        if (tangentDirs_u[iPoint].size() >= 3) {
          /*--- Set up normal equations: A^T A x = A^T b ---*/
          ScalarType ATA[3][3] = {{0.0}};  // For [a, b, c]
          ScalarType ATb[3] = {0.0};

          for (size_t j = 0; j < tangentDirs_u[iPoint].size(); ++j) {
            ScalarType du = tangentDirs_u[iPoint][j];
            ScalarType dv = tangentDirs_v[iPoint][j];
            ScalarType kappa = normalCurvatures[iPoint][j];
            ScalarType w = weights[iPoint][j];

            /*--- Row: [du^2, 2*du*dv, dv^2] ---*/
            ScalarType row[3] = {du*du, 2.0*du*dv, dv*dv};

            for (size_t i = 0; i < 3; ++i) {
              for (size_t k = 0; k < 3; ++k) {
                ATA[i][k] += w * row[i] * row[k];
              }
              ATb[i] += w * row[i] * kappa;
            }
          }

          /*--- Add trace constraint: a + c = 2*meanK using Lagrange multiplier approach ---*/
          /*--- Original system: ATA * x = ATb ---*/
          /*--- Augmented system: [ATA C^T; C 0] * [x; λ] = [ATb; 2*meanK] ---*/
          /*--- where C = [1 0 1] enforces a + c = 2*meanK ---*/

          ScalarType augmented[4][4] = {{0.0}};
          ScalarType rhsAugmented[4] = {0.0};

          /*--- Copy ATA into top-left 3x3 block ---*/
          for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
              augmented[i][j] = ATA[i][j];
            }
            rhsAugmented[i] = ATb[i];
          }

          /*--- Add constraint row: [1 0 1 0] ---*/
          augmented[3][0] = 1.0;  // a coefficient
          augmented[3][2] = 1.0;  // c coefficient
          rhsAugmented[3] = 2.0 * meanK;

          /*--- Add constraint column: transpose of constraint row ---*/
          augmented[0][3] = 1.0;
          augmented[2][3] = 1.0;

          /*--- Solve 4x4 augmented system ---*/
          ScalarType solution[4] = {0.0};

          /*--- Gaussian elimination with partial pivoting ---*/
          for (size_t i = 0; i < 4; ++i) {
            /*--- Find pivot ---*/
            size_t maxRow = i;
            for (size_t k = i + 1; k < 4; ++k) {
              if (fabs(augmented[k][i]) > fabs(augmented[maxRow][i])) maxRow = k;
            }

            /*--- Swap rows ---*/
            if (maxRow != i) {
              for (size_t k = 0; k < 4; ++k) {
                swap(augmented[i][k], augmented[maxRow][k]);
              }
              swap(rhsAugmented[i], rhsAugmented[maxRow]);
            }

            /*--- Eliminate column ---*/
            for (size_t k = i + 1; k < 4; ++k) {
              if (fabs(augmented[i][i]) > 1e-12) {
                ScalarType factor = augmented[k][i] / augmented[i][i];
                for (size_t j = i; j < 4; ++j) {
                  augmented[k][j] -= factor * augmented[i][j];
                }
                rhsAugmented[k] -= factor * rhsAugmented[i];
              }
            }
          }

          /*--- Back substitution ---*/
          for (int i = 3; i >= 0; --i) {
            solution[i] = rhsAugmented[i];
            for (size_t j = i + 1; j < 4; ++j) {
              solution[i] -= augmented[i][j] * solution[j];
            }
            if (fabs(augmented[i][i]) > 1e-12) {
              solution[i] /= augmented[i][i];
            }
          }

          /*--- Extract curvature tensor components (solution[3] is the Lagrange multiplier) ---*/
          B[0][0] = solution[0];  // a
          B[0][1] = B[1][0] = solution[1];  // b
          B[1][1] = solution[2];  // c
        } else {
          /*--- Fallback: use diagonal tensor scaled by principal curvatures ---*/
          B[0][0] = kappa1;
          B[1][1] = kappa2;
          B[0][1] = B[1][0] = 0.0;
        }

        /*--- Compute eigenvectors of 2x2 curvature tensor B ---*/
        ScalarType B2D[2][2] = {{B[0][0], B[0][1]}, {B[1][0], B[1][1]}};
        ScalarType eigVec2D[2][2], eigVal2D[2], work2D[2];

        CBlasStructure::EigenDecomposition(B2D, eigVec2D, eigVal2D, 2, work2D);

        /*--- Sort eigenvalues in descending order to match kappa1, kappa2 ---*/
        if (eigVal2D[1] > eigVal2D[0]) {
          swap(eigVal2D[0], eigVal2D[1]);
          /*--- Swap corresponding eigenvectors ---*/
          for (size_t i = 0; i < 2; ++i) {
            swap(eigVec2D[i][0], eigVec2D[i][1]);
          }
        }

        /*--- Build 3D principal directions from 2D eigenvectors ---*/
        /*--- R[:,0] = first principal direction (max curvature) ---*/
        for (size_t dim = 0; dim < 3; ++dim) {
          R[dim][0] = eigVec2D[0][0] * u[dim] + eigVec2D[1][0] * v[dim];
        }

        /*--- R[:,1] = second principal direction (min curvature) ---*/
        for (size_t dim = 0; dim < 3; ++dim) {
          R[dim][1] = eigVec2D[0][1] * u[dim] + eigVec2D[1][1] * v[dim];
        }

        /*--- R[:,2] = normal direction (already set above) ---*/

        /*--- Recompose metric tensor ---*/
        CBlasStructure::EigenRecomposition(M, R, EigVal, 3);
        Tensor::set(metric, iPoint, 0, M, 1.0, 3);
      }
    }
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

/*!
 * \brief Compute a metric tensor based on surface geometry.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[out] metric - The computed geometric metric tensor field.
 */
template<class ScalarType, class Tensor, class MetricType>
void geometricSurfaceMetrics(CGeometry& geometry, const CConfig& config,
                             MetricType& metric) {
  switch (geometry.GetnDim()) {
    case 2:
      detail::geometricSurfaceMetrics<2, ScalarType, Tensor>(geometry, config, metric);
      break;
    case 3:
      detail::geometricSurfaceMetrics<3, ScalarType, Tensor>(geometry, config, metric);
      break;
    default:
      SU2_MPI::Error("Too many dimensions for geometric metric computation.", CURRENT_FUNCTION);
      break;
  }
}