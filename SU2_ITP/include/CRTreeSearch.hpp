/*!
 * \file CRTreeSearch.hpp
 * \brief Declaration of the spatial search class using Boost R-tree.
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

#include "../../Common/include/parallelization/mpi_structure.hpp"
#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/geometry/CGeometry.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <vector>
#include <set>

/*!
 * \class CRTreeSearch
 * \brief Class for spatial search operations using Boost R-tree.
 *        Provides efficient bounding box queries for mesh interpolation.
 * \tparam nDim Number of spatial dimensions.
 */
template<unsigned short nDim>
class CRTreeSearch {
public:
  /*--- Boost.Geometry type definitions ---*/
  using Point = boost::geometry::model::point<su2double, nDim, boost::geometry::cs::cartesian>;
  using Box = boost::geometry::model::box<Point>;
  using NodeValue = std::pair<Point, unsigned long>;
  using RTree = boost::geometry::index::rtree<NodeValue, boost::geometry::index::quadratic<16>>;

private:
  int rank;
  RTree nodeTree;  /*!< \brief R-tree for spatial indexing of source mesh nodes. */
  bool treeBuilt;  /*!< \brief Whether the tree has been built. */

public:
  /*!
   * \brief Constructor.
   * \param[in] MPICommunicator - MPI communicator.
   */
  explicit CRTreeSearch(SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor.
   */
  ~CRTreeSearch() = default;

  /*!
   * \brief Build the R-tree from source mesh nodes.
   * \param[in] geometry_src - Source mesh geometry.
   */
  void BuildTree(CGeometry* geometry_src);
  
  /*!
   * \brief Clear the tree and reset the built flag.
   */
  void ClearTree();

  /*!
   * \brief Create a bounding box from element vertices.
   * \param[in] elemCoords - Array of element coordinates.
   * \return Bounding box encompassing the element.
   */
  Box CreateBoundingBox(const su2double* elemCoords) const;

  /*!
   * \brief Create a bounding box from element vertices with padding.
   * \param[in] elemCoords - Array of element coordinates.
   * \param[in] padding - Additional padding around the element.
   * \return Bounding box encompassing the element with padding.
   */
  Box CreateBoundingBox(const su2double* elemCoords, su2double padding) const;

private:
  /*!
   * \brief Create a bounding box from 2D triangle vertices with padding.
   * \param[in] tri - Flat array of triangle coordinates.
   * \param[in] padding - Additional padding around the triangle.
   * \return Bounding box encompassing the triangle with padding.
   */
  Box CreateBoundingBox2D(const su2double* tri, su2double padding) const;

  /*!
   * \brief Create a bounding box from 3D tetrahedron vertices with padding.
   * \param[in] tet - Flat array of tetrahedron coordinates.
   * \param[in] padding - Additional padding around the tetrahedron.
   * \return Bounding box encompassing the tetrahedron with padding.
   */
  Box CreateBoundingBox3D(const su2double* tet, su2double padding) const;

public:

  /*!
   * \brief Search for nodes within a bounding box.
   * \param[in] boundingBox - The bounding box to search within.
   * \param[out] containedNodes - Set of node IDs found within the bounding box.
   */
  void SearchNodesInBox(const Box& boundingBox, std::set<unsigned long>& containedNodes) const;

  /*!
   * \brief Check if the tree has been built.
   * \return True if tree is built, false otherwise.
   */
  bool IsTreeBuilt() const { return treeBuilt; }

  /*!
   * \brief Get the number of nodes in the tree.
   * \return Number of nodes stored in the tree.
   */
  size_t GetTreeSize() const { return nodeTree.size(); }
};
