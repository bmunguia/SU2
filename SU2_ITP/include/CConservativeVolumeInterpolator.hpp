/*!
 * \file CConservativeVolumeInterpolator.hpp
 * \brief Headers of the main conservative solution interpolation subroutines.
 *        The subroutines and functions are in the <i>CConservativeVolumeInterpolator.cpp</i> file.
 * \author B. Munguía
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

#include <optional>
#include "CVolumeInterpolator.hpp"

struct IntersectionData {
    unsigned long srcElemID;   // ID of source element for the current intersection
    vector<su2double> coords;  // Coordinates of each intersection sub-element vertex
    vector<su2double> vols;    // Volume of each intersection sub-element
};

using IntersectionMesh = vector<IntersectionData>;
using IntersectionMeshMap = map<unsigned long, IntersectionMesh>;

/*!
 * \class CConservativeVolumeInterpolator
 * \brief Performs conservative solution interpolation between meshes using the Alauzet 2015 method.
 * \author B. Munguía
 */
class CConservativeVolumeInterpolator : public CVolumeInterpolator {
  public:
    /*!
     * \brief Constructor of the class.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CConservativeVolumeInterpolator(SU2_Comm MPICommunicator);

    /*!
     * \brief Default destructor.
     */
    ~CConservativeVolumeInterpolator() override = default;

    /*!
     * \brief Main interpolation routine using conservative interpolation.
     * \param[in] config - Configuration object
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_container_src - Source mesh solver
     * \param[in] solver_container_dst - Destination mesh solver
     */
    void Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                     CSolver** solver_container_src, CSolver** solver_container_dst) override;

  private:
    /*!
     * \brief Conservative interpolation (Alauzet 2015).
     * \param[in] config - Configuration object
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_src - Source mesh solver
     * \param[in] solver_dst - Destination mesh solver
     */
    void ConservativeInterpolation(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                   CSolver* solver_src, CSolver* solver_dst);

    /*!
     * \brief Containment search.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] coor_dst - Destination mesh coordinates (after curvature correction if applied)
     * \param[out] containingElems - Elements containing destination points
     * \param[out] containingElemRanks - Ranks of elements containing destination points
     * \param[out] pointsFailed - Nodes for which no containing element was found
     */
    void PointLocalization(CGeometry* geometry_src,
                           const vector<su2double>& coor_dst,
                           vector<optional<unsigned long>>& containingElems,
                           vector<int>& containingElemRanks,
                           vector<unsigned long>& pointsFailed);

    /*!
     * \brief Compute the mass and gradient of the solution variables.
     * \param[in] geometry - Mesh geometry
     * \param[in] solver - Solver definition
     * \param[out] elemMass - Mass at each element for each solution variable
     * \param[out] elemGrad - Gradient at each element for each solution variable
     */
    void ComputeSourceSolutionMass(CGeometry* geometry,
                                   CSolver* solver,
                                   vector<vector<su2double> >& elemMass,
                                   vector<vector<su2double> >& elemGrad);

    /*!
     * \brief Compute destination mesh mass and gradients using Gauss quadrature over intersection regions.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_src - Source mesh solver
     * \param[in] overlapMeshes - Map from dst element ID to vector of (src element ID, triangulated mesh) pairs
     * \param[in] incompleteOverlaps - List of destination element IDs with incomplete overlap
     * \param[in] srcElemMass - Source element masses
     * \param[in] srcElemGrad - Source element gradients
     * \param[out] dstElemMass - Destination element masses
     * \param[out] dstElemGrad - Destination element gradients
     */
    void ComputeDestinationMassAndGradient(CGeometry* geometry_src,
                                           CGeometry* geometry_dst,
                                           CSolver* solver_src,
                                           const IntersectionMeshMap& overlapMeshes,
                                           const vector<unsigned long>& incompleteOverlaps,
                                           const vector<vector<su2double>>& srcElemMass,
                                           const vector<vector<su2double>>& srcElemGrad,
                                           vector<vector<su2double>>& dstElemMass,
                                           vector<vector<su2double>>& dstElemGrad);

    /*!
     * \brief Compute overlapping elements between source and destination meshes, and mesh the intersection regions.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] coor_dst - Destination mesh coordinates (after curvature correction if applied)
     * \param[in] containingElems - Elements containing destination points (from point localization)
     * \param[out] overlapMeshes - Map from dst element ID to vector of (src element ID, triangulated mesh) pairs
     * \param[out] incompleteOverlaps - List of destination element IDs with incomplete overlap
     * \param[out] uncoveredNodes - Set of destination node IDs that are not covered by any source element
     */
    void CreateIntersectionMeshes(CGeometry* geometry_src,
                                  CGeometry* geometry_dst,
                                  const vector<su2double> &coor_dst,
                                  IntersectionMeshMap& overlapMeshes,
                                  vector<unsigned long>& incompleteOverlaps,
                                  set<unsigned long>& uncoveredNodes);

    /*!
     * \brief Triangle-triangle intersection using Alauzet method (signed distance functions).
     * \param[in] dstTri - Destination triangle vertex coordinates
     * \param[in] srcTri - Source triangle vertex coordinates
     * \param[out] intersectionPoints - Cloud of intersection points
     * \param[out] intersectionElemCoords - Flat array of meshed intersection coordinates
     * \param[out] intersectionElemVols - Volumes of meshed intersection region elements
     * \param[out] isCoveredVertexP - boolean array indicating which vertices of P are part of the intersection
     * \param[out] isCoveredVertexQ - boolean array indicating which vertices of Q are part of the intersection
     * \return True if triangles intersect
     */
    bool TriangleTriangleIntersection(const su2double dstTri[6],
                                      const su2double srcTri[6],
                                      vector<su2double>& intersectionPoints,
                                      vector<su2double>& intersectionElemCoords,
                                      vector<su2double>& intersectionElemVols,
                                      bool isCoveredVertexP[3],
                                      bool isCoveredVertexQ[3]);

    /*!
     * \brief Process degenerate edge-edge intersection cases following Alauzet's algorithm.
     * \param[in] P_edges - Destination triangle edges (3 edges, each with 2 vertices)
     * \param[in] Q_edges - Source triangle edges (3 edges, each with 2 vertices)
     * \param[in] power_P - Powers of P vertices w.r.t. Q edges [3][3]
     * \param[in] power_Q - Powers of Q vertices w.r.t. P edges [3][3]
     * \param[in] EPS - Tolerance for zero detection
     * \param[out] intersectionPoints - Cloud of intersection points
     * \param[out] isDegenerateVertexP - boolean array indicating which vertices of P are degenerate
     * \param[out] isDegenerateVertexQ - boolean array indicating which vertices of Q are degenerate
     * \param[out] isDegenerateEdgePair - 3x3 boolean array indicating which edge pairs are degenerate [iP][jQ]
     */
    void ProcessDegenerateEdgeIntersections(su2double* P_edges[3][2], su2double* Q_edges[3][2],
                                            const su2double power_P[3][3], const su2double power_Q[3][3],
                                            const su2double EPS,
                                            vector<su2double>& intersectionPoints,
                                            bool isDegenerateVertexP[3],
                                            bool isDegenerateVertexQ[3],
                                            bool isDegenerateEdgePair[3][3]);

    /*!
     * \brief Compute signed distance (power) of a point to a line.
     * \param[in] point - Point coordinates (x, y)
     * \param[in] lineStart - Line start point (x, y)
     * \param[in] lineEnd - Line end point (x, y)
     * \return Signed distance (positive if point is on left side of oriented line)
     */
    su2double ComputeSignedDistance(const su2double point[2], const su2double lineStart[2], const su2double lineEnd[2]);

    /*!
     * \brief Mesh a convex polygon into triangles following Alauzet's method.
     * \param[in] polygonPoints - Convex polygon vertices (x,y pairs in order)
     * \param[out] polygonElemCoords - Flat array of meshed intersection coordinates
     * \param[out] polygonElemVols - Volumes of meshed intersection region elements
     */
    void MeshConvexPolygon(const vector<su2double>& polygonPoints,
                           vector<su2double>& polygonElemCoords,
                           vector<su2double>& polygonElemVols);

    /*!
     * \brief Apply local maximum principle correction following Alauzet's method.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_src - Source mesh solver
     * \param[in] coor_dst - Destination mesh coordinates (after curvature correction if applied)
     * \param[in] overlapMeshes - Map from dst element ID to overlapping src element IDs
     * \param[in] srcElemMass - Source element masses
     * \param[in] srcElemGrad - Source element gradients
     * \param[in,out] dstElemMass - Destination element masses
     * \param[in,out] dstElemGrad - Destination element gradients (corrected)
     */
    void ApplyMaximumPrincipleCorrection(CGeometry* geometry_src,
                                         CGeometry* geometry_dst,
                                         CSolver* solver_src,
                                         const vector<su2double>& coor_dst,
                                         const IntersectionMeshMap& overlapMeshes,
                                         const vector<vector<su2double>>& srcElemMass,
                                         const vector<vector<su2double>>& srcElemGrad,
                                         vector<vector<su2double>>& dstElemMass,
                                         vector<vector<su2double>>& dstElemGrad);

    /*!
     * \brief Calculate the solution at destination nodes from the mass and gradient.
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_dst - Destination mesh solver
     * \param[in] coor_dst - Destination mesh coordinates (after curvature correction if applied)
     * \param[in] dstElemMass - Destination element masses
     * \param[in] dstElemGrad - Destination element gradients
     */
    void DistributeSolutionToNodes(CGeometry* geometry_dst,
                                   CSolver* solver_dst,
                                   const vector<su2double>& coor_dst,
                                   const vector<vector<su2double>>& dstElemMass,
                                   const vector<vector<su2double>>& dstElemGrad);

    /*!
     * \brief Handle solution distribution to destination nodes outside the source domain.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_src - Source mesh solver
     * \param[in] solver_dst - Destination mesh solver
     * \param[in] coor_dst - Destination mesh coordinates (after curvature correction if applied)
     * \param[in] overlapMeshes - Map from dst element ID to overlapping src element IDs
     * \param[in] uncoveredNodes - Set of destination node IDs that are not covered by any source element
     * \param[in] dstElemMass - Destination element masses
     * \param[in] dstElemGrad - Destination element gradients
     */
    void DistributeSolutionOutsideDomain(CGeometry* geometry_src,
                                         CGeometry* geometry_dst,
                                         CSolver* solver_src,
                                         CSolver* solver_dst,
                                         const vector<su2double>& coor_dst,
                                         const IntersectionMeshMap& overlapMeshes,
                                         const set<unsigned long>& uncoveredNodes,
                                         const vector<vector<su2double>>& dstElemMass,
                                         const vector<vector<su2double>>& dstElemGrad);

    /*!
     * \brief Perform surface interpolation for physical boundary nodes using nearest surface elements.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_src - Source mesh solver
     * \param[in,out] solver_dst - Destination mesh solver (solution will be updated)
     * \param[in] coor_dst - Destination mesh coordinates
     */
    void SurfaceInterpolation(CGeometry* geometry_src,
                              CGeometry* geometry_dst,
                              CSolver* solver_src,
                              CSolver* solver_dst,
                              const vector<su2double>& coor_dst);

  private:
    /*!
     * \brief Compute solution mass and gradient in a triangular element from vertex values.
     * \param[in] vertexCoor - Triangle vertex coordinates (6 values: x0,y0,x1,y1,x2,y2)
     * \param[in] vertexSol - Solution values at vertices (nVar values per vertex)
     * \param[in] elemVol - Volume of the element
     * \param[in] nVar - Number of variables
     * \param[out] mass - Computed mass integral for each variable
     * \param[out] grad - Computed gradient for each variable (nVar*nDim values)
     */
    void ComputeTriangleMassAndGradient(const su2double vertexCoor[6],
                                        const vector<vector<su2double>>& vertexSol,
                                        const su2double elemVol,
                                        const unsigned short nVar,
                                        vector<su2double>& mass,
                                        vector<su2double>& grad);

    /*!
     * \brief Compute solution mass and gradient in a triangular element from vertex values using CFEMStandardElement.
     * \param[in] vertexCoor - Triangle vertex coordinates (6 values: x0,y0,x1,y1,x2,y2)
     * \param[in] vertexSol - Solution values at vertices (nVar values per vertex)
     * \param[in] elemVol - Volume of the element
     * \param[in] nVar - Number of variables
     * \param[out] mass - Computed mass integral for each variable
     * \param[out] grad - Computed gradient for each variable (nVar*nDim values)
     */
    void ComputeTriangleMassAndGradientFEM(const su2double vertexCoor[6],
                                           const vector<vector<su2double>>& vertexSol,
                                           const su2double elemVol,
                                           const unsigned short nVar,
                                           vector<su2double>& mass,
                                           vector<su2double>& grad);
};
