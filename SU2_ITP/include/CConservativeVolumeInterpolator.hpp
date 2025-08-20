/*!
 * \file CConservativeVolumeInterpolator.hpp
 * \brief Headers of the main conservative solution interpolation subroutines.
 *        The subroutines and functions are in the <i>CConservativeVolumeInterpolator.cpp</i> file.
 * \author B. Munguía, E. van der Weide
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

#include "CVolumeInterpolator.hpp"

/*!
 * \enum IntersectionType
 * \brief Types of edge-edge intersections for triangle intersection detection.
 */
enum IntersectionType {
  NO_INTERSECTION,    /*!< \brief No intersection between edges. */
  VERTEX_ON_EDGE,     /*!< \brief One vertex lies on the other edge. */
  EDGE_CROSSING,      /*!< \brief Edges cross at interior points. */
  EDGE_OVERLAP        /*!< \brief Edges overlap (collinear case). */
};

/*!
 * \struct IntersectionResult
 * \brief Detailed information about edge-edge intersection.
 */
struct IntersectionResult {
  IntersectionType type;       /*!< \brief Type of intersection. */
  vector<su2double> points;    /*!< \brief Intersection points. */
  int vertexIndex;             /*!< \brief Which vertex lies on edge (for VERTEX_ON_EDGE). */
  bool isFirstEdge;            /*!< \brief Which edge has the vertex (for VERTEX_ON_EDGE). */
};

/*!
 * \class CConservativeVolumeInterpolator
 * \brief Performs conservative solution interpolation between meshes using the Alauzet 2015 method.
 * \author B. Munguía, E. van der Weide
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
     * \param[in] coor_corrected - Destination mesh coordinates after curvature correction
     * \param[out] containingElems - Elements containing destination points
     * \param[out] containingElemRanks - Ranks of elements containing destination points
     * \param[out] pointsFailed - Nodes for which no containing element was found
     */
    void PointLocalization(CGeometry* geometry_src,
                           const vector<su2double>& coor_corrected,
                           vector<unsigned long>& containingElems,
                           vector<int>& containingElemRanks,
                           vector<unsigned long>& pointsFailed);

    /*!
     * \brief Compute the mass and gradient of the solution variables.
     * \param[in] geometry - Mesh geometry
     * \param[in] solver - Solver definition
     * \param[out] elemMass - Mass at each element for each solution variable
     * \param[out] elemGrad - Gradient at each element for each solution variable
     */
    void ComputeSolutionMass(CGeometry* geometry,
                             CSolver* solver,
                             vector<vector<su2double> >& elemMass,
                             vector<vector<su2double> >& elemGrad);

    /*!
     * \brief Compute destination mesh mass and gradients using Gauss quadrature over intersection regions.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_src - Source mesh solver
     * \param[in] overlappingElements - Map from dst element ID to overlapping src element IDs
     * \param[in] intersectionMeshes - Map from dst element ID to triangulated intersection meshes
     * \param[in] srcElemMass - Source element masses
     * \param[in] srcElemGrad - Source element gradients
     * \param[out] dstElemMass - Destination element masses
     * \param[out] dstElemGrad - Destination element gradients
     */
    void ComputeDestinationMassAndGradient(CGeometry* geometry_src,
                                           CGeometry* geometry_dst,
                                           CSolver* solver_src,
                                           const map<unsigned long, vector<unsigned long>>& overlappingElements,
                                           const map<unsigned long, vector<vector<su2double>>>& intersectionMeshes,
                                           const vector<vector<su2double>>& srcElemMass,
                                           const vector<vector<su2double>>& srcElemGrad,
                                           vector<vector<su2double>>& dstElemMass,
                                           vector<vector<su2double>>& dstElemGrad);

    /*!
     * \brief Compute overlapping elements between source and destination meshes.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] containingElems - Elements containing destination points (from point localization)
     * \param[out] overlappingElements - Map from dst element ID to vector of overlapping src element IDs
     * \param[out] intersectionMeshes - Map from dst element ID to vector of triangulated intersection meshes
     */
    void ComputeOverlappingElements(CGeometry* geometry_src,
                                    CGeometry* geometry_dst,
                                    const vector<unsigned long>& containingElems,
                                    map<unsigned long, vector<unsigned long>>& overlappingElements,
                                    map<unsigned long, vector<vector<su2double>>>& intersectionMeshes);

    /*!
     * \brief Triangle-triangle intersection using Alauzet method (signed distance functions).
     * \param[in] geometry_src - Source mesh geometry (for neighbor detection)
     * \param[in] srcTri - Second triangle vertices (source triangle)
     * \param[in] dstTri - First triangle vertices (destination triangle)
     * \param[in] srcElemID - Source element ID
     * \param[out] intersectionPoints - Cloud of intersection points
     * \param[out] newCandidates - New candidate elements detected during intersection
     * \return True if triangles intersect
     */
    bool TriangleTriangleIntersection(CGeometry* geometry_src,
                                      const su2double dstTri[6],
                                      const su2double srcTri[6],
                                      unsigned long srcElemID,
                                      vector<su2double>& intersectionPoints,
                                      set<unsigned long>& newCandidates);

    /*!
     * \brief Add face (edge in 2D) neighbor to candidate list when face is intersected.
     * \param[in] geometry - Mesh geometry
     * \param[in] elemID - Element ID
     * \param[in] faceIndex - Local face index (0, 1, or 2 for triangles)
     * \param[out] newCandidates - Set to add new candidates to
     */
    void AddFaceNeighborToCandidates(CGeometry* geometry_src,
                                     unsigned long srcElemID,
                                     unsigned short edgeIndex,
                                     set<unsigned long>& newCandidates);

    /*!
     * \brief Add vertex ball to candidate list when vertex is inside triangle.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] srcElemID - Source element ID
     * \param[in] vertexIndex - Local vertex index (0, 1, or 2 for triangles)
     * \param[out] newCandidates - Set to add new candidates to
     */
    void AddVertexBallToCandidates(CGeometry* geometry_src,
                                   unsigned long srcElemID,
                                   unsigned short vertexIndex,
                                   set<unsigned long>& newCandidates);

    /*!
     * \brief Compute signed distance (power) of a point to a line.
     * \param[in] point - Point coordinates (x, y)
     * \param[in] lineStart - Line start point (x, y)
     * \param[in] lineEnd - Line end point (x, y)
     * \return Signed distance (positive if point is on left side of oriented line)
     */
    su2double ComputeSignedDistance(const su2double point[2], const su2double lineStart[2], const su2double lineEnd[2]);

    /*!
     * \brief Compute parametric coordinate of a point on a line segment.
     * \param[in] point - Point to find parametric coordinate for
     * \param[in] lineStart - Line segment start point
     * \param[in] lineEnd - Line segment end point
     * \return Parametric coordinate in range [-1, 1]
     */
    su2double ComputeLineParametricCoord(const su2double point[2], const su2double lineStart[2], const su2double lineEnd[2]);

     /*!
     * \brief Compute intersection of two line segments.
     * \param[in] P0 - First segment start point
     * \param[in] P1 - First segment end point
     * \param[in] Q0 - Second segment start point
     * \param[in] Q1 - Second segment end point
     * \return IntersectionResult with detailed intersection information
     */
    IntersectionResult LineSegmentIntersection(const su2double P0[2], const su2double P1[2],
                                               const su2double Q0[2], const su2double Q1[2]);

    /*!
     * \brief Compute convex hull of a set of 2D points using Graham scan.
     * \param[in,out] points - Input points, output convex hull points (in order)
     */
    void ComputeConvexHull(vector<su2double>& points);

    /*!
     * \brief Mesh a convex polygon into triangles following Alauzet's method.
     * \param[in] polygonPoints - Convex polygon vertices (x,y pairs in order)
     * \param[out] triangles - Output triangle vertices (6 coordinates per triangle: x0,y0,x1,y1,x2,y2)
     */
    void MeshConvexPolygon(const vector<su2double>& polygonPoints, vector<su2double>& triangles);

    /*!
     * \brief Apply local maximum principle correction following Alauzet's method.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_src - Source mesh solver
     * \param[in] overlappingElements - Map from dst element ID to overlapping src element IDs
     * \param[in] srcElemMass - Source element masses
     * \param[in] srcElemGrad - Source element gradients
     * \param[in,out] dstElemMass - Destination element masses
     * \param[in,out] dstElemGrad - Destination element gradients (corrected)
     */
    void ApplyMaximumPrincipleCorrection(CGeometry* geometry_src,
                                         CGeometry* geometry_dst,
                                         CSolver* solver_src,
                                         const map<unsigned long, vector<unsigned long>>& overlappingElements,
                                         const vector<vector<su2double>>& srcElemMass,
                                         const vector<vector<su2double>>& srcElemGrad,
                                         vector<vector<su2double>>& dstElemMass,
                                         vector<vector<su2double>>& dstElemGrad);

    /*!
     * \brief Calculate the solution at destination nodes from the mass and gradient.
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_dst - Destination mesh solver
     * \param[in] dstElemMass - Destination element masses
     * \param[in] dstElemGrad - Destination element gradients
     */
    void DistributeSolutionToNodes(CGeometry* geometry_dst,
                                   CSolver* solver_dst,
                                   const vector<vector<su2double>>& dstElemMass,
                                   const vector<vector<su2double>>& dstElemGrad);

  private:
    /*!
     * \brief Helper function to compute mass and gradient from vertex values using CFEMStandardElement.
     * \param[in] vertexCoords - Triangle vertex coordinates (6 values: x0,y0,x1,y1,x2,y2)
     * \param[in] vertexSolutions - Solution values at vertices (nVar values per vertex)
     * \param[in] nVar - Number of variables
     * \param[out] mass - Computed mass integral for each variable
     * \param[out] gradient - Computed gradient for each variable (nVar*nDim values)
     */
    void ComputeTriangleMassAndGradient(const su2double vertexCoords[6],
                                        const vector<vector<su2double>>& vertexSolutions,
                                        unsigned short nVar,
                                        vector<su2double>& mass,
                                        vector<su2double>& gradient);
};
