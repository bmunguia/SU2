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
    void Interpolate(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
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

    void PointLocalization(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                           const vector<su2double>& coor_corrected, vector<unsigned long>& containingElems,
                           vector<int>& containingElemRanks, vector<unsigned long>& pointsFailed);

    void ComputeSolutionMass(CGeometry* geometry, CSolver* solver, unsigned short nVar,
                             vector<vector<su2double> >& elemMass,
                             vector<vector<su2double> >& elemGrad);

    /*!
     * \brief Compute overlapping elements between source and destination meshes.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[out] overlappingElements - Map from dst element ID to vector of overlapping src element IDs
     * \param[in] containingElems - Elements containing destination points (from point localization)
     */
    void ComputeOverlappingElements(CGeometry* geometry_src, CGeometry* geometry_dst,
                                    map<unsigned long, vector<unsigned long>>& overlappingElements,
                                    const vector<unsigned long>& containingElems);

    /*!
     * \brief Robust triangle-triangle intersection using Alauzet method.
     * \param[in] dstTri - First triangle vertices (destination triangle)
     * \param[in] srcTri - Second triangle vertices (source triangle)
     * \param[out] intersectionPoints - Cloud of intersection points
     * \param[in] geometry_src - Source mesh geometry (for neighbor detection)
     * \param[in] srcElemID - Source element ID
     * \param[out] newCandidates - New candidate elements detected during intersection
     * \return True if triangles intersect
     */
    bool TriangleTriangleIntersection(const su2double dstTri[6], const su2double srcTri[6],
                                      vector<su2double>& intersectionPoints,
                                      CGeometry* geometry_src, unsigned long srcElemID,
                                      set<unsigned long>& newCandidates);

    /*!
     * \brief Add edge neighbor to candidate list when edge is intersected.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] srcElemID - Source element ID
     * \param[in] edgeIndex - Local edge index (0, 1, or 2 for triangles)
     * \param[out] newCandidates - Set to add new candidates to
     */
    void AddEdgeNeighborToCandidates(CGeometry* geometry_src, unsigned long srcElemID,
                                     unsigned short edgeIndex, set<unsigned long>& newCandidates);

    /*!
     * \brief Add vertex ball to candidate list when vertex is inside triangle.
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] srcElemID - Source element ID
     * \param[in] vertexIndex - Local vertex index (0, 1, or 2 for triangles)
     * \param[out] newCandidates - Set to add new candidates to
     */
    void AddVertexBallToCandidates(CGeometry* geometry_src, unsigned long srcElemID,
                                   unsigned short vertexIndex, set<unsigned long>& newCandidates);

    /*!
     * \brief Compute signed distance (power) of a point to a line.
     * \param[in] point - Point coordinates (x, y)
     * \param[in] lineStart - Line start point (x, y)
     * \param[in] lineEnd - Line end point (x, y)
     * \return Signed distance (positive if point is on left side of oriented line)
     */
    su2double ComputeSignedDistance(const su2double point[2], const su2double lineStart[2], const su2double lineEnd[2]);

     /*!
     * \brief Compute intersection of two line segments.
     * \param[in] P0 - First segment start point
     * \param[in] P1 - First segment end point
     * \param[in] Q0 - Second segment start point
     * \param[in] Q1 - Second segment end point
     * \param[out] intersections - Vector of intersection points (x,y pairs)
     * \return True if segments intersect
     */
    bool LineSegmentIntersection(const su2double P0[2], const su2double P1[2],
                                 const su2double Q0[2], const su2double Q1[2],
                                 vector<su2double>& intersections);

    /*!
     * \brief Compute convex hull of a set of 2D points using Graham scan.
     * \param[in,out] points - Input points, output convex hull points (in order)
     */
    void ComputeConvexHull(vector<su2double>& points);
};
