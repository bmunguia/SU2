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

namespace bg = boost::geometry;

/*--- Define geometry types outside the class for concrete use ---*/
using Point2D = bg::model::point<su2double, 2, bg::cs::cartesian>;
using Point3D = bg::model::point<su2double, 3, bg::cs::cartesian>;
using Box2D = bg::model::box<Point2D>;
using Box3D = bg::model::box<Point3D>;
using NodeValue2D = std::pair<Point2D, unsigned long>;
using NodeValue3D = std::pair<Point3D, unsigned long>;
using RTree2D = bg::index::rtree<NodeValue2D, bg::index::quadratic<16>>;
using RTree3D = bg::index::rtree<NodeValue3D, bg::index::quadratic<16>>;

/*!
 * \class CRTreeSearchBase
 * \brief Base class for spatial search operations using R-tree.
 *        Provides a common interface for different dimensional R-trees.
 */
class CRTreeSearchBase {
public:
  /*!
   * \brief Virtual destructor.
   */
  virtual ~CRTreeSearchBase() = default;

  /*!
   * \brief Build the R-tree from mesh nodes.
   * \param[in] geometry - Mesh geometry.
   */
  virtual void BuildTree(CGeometry* geometry) = 0;

  /*!
   * \brief Clear the tree and reset the built flag.
   */
  virtual void ClearTree() = 0;

  /*!
   * \brief Check if the tree has been built.
   * \return True if tree is built, false otherwise.
   */
  virtual bool IsTreeBuilt() const = 0;

  /*!
   * \brief Get the number of nodes in the tree.
   * \return Number of nodes stored in the tree.
   */
  virtual size_t GetTreeSize() const = 0;

  /*!
   * \brief Search for nodes within a bounding box (element coordinates).
   * \param[in] elemCoor - Array of element vertex coordinates.
   * \param[out] containedNodes - Set of node IDs found within the bounding box.
   */
  virtual void SearchNodesInElement(const su2double* elemCoor, std::set<unsigned long>& containedNodes) const = 0;

  /*!
   * \brief Search for nodes within a bounding box (element coordinates with padding).
   * \param[in] elemCoor - Array of element vertex coordinates.
   * \param[out] containedNodes - Set of node IDs found within the bounding box.
   * \param[in] padding - Additional padding around the element.
   */
  virtual void SearchNodesInElement(const su2double* elemCoor, std::set<unsigned long>& containedNodes, su2double padding) const = 0;

protected:
  /*!
   * \brief Set coordinates for a 2D point.
   * \param[out] pt - Point to set coordinates for.
   * \param[in] coor - Array of coordinates [x, y].
   */
  void SetPoint2D(Point2D& pt, const su2double* coor) const;

  /*!
   * \brief Set coordinates for a 3D point.
   * \param[out] pt - Point to set coordinates for.
   * \param[in] coor - Array of coordinates [x, y, z].
   */
  void SetPoint3D(Point3D& pt, const su2double* coor) const;

  /*!
   * \brief Create a 2D bounding box from triangle vertices with padding.
   * \param[in] tri - Flat array of triangle coordinates.
   * \param[in] padding - Additional padding around the triangle.
   * \return Bounding box encompassing the triangle with padding.
   */
  Box2D CreateBoundingBox2D(const su2double* tri, su2double padding) const;

  /*!
   * \brief Create a 3D bounding box from tetrahedron vertices with padding.
   * \param[in] tet - Flat array of tetrahedron coordinates.
   * \param[in] padding - Additional padding around the tetrahedron.
   * \return Bounding box encompassing the tetrahedron with padding.
   */
  Box3D CreateBoundingBox3D(const su2double* tet, su2double padding) const;
};

/*!
 * \class CRTreeSearch
 * \brief Class for spatial search operations using Boost R-tree.
 *        Provides efficient bounding box queries for mesh interpolation.
 * \tparam nDim Number of spatial dimensions.
 */
template<unsigned short nDim>
class CRTreeSearch : public CRTreeSearchBase {
public:
  /*--- Dimension-specific type definitions ---*/
  using Point = typename std::conditional<nDim == 2, Point2D, Point3D>::type;
  using Box = typename std::conditional<nDim == 2, Box2D, Box3D>::type;
  using NodeValue = typename std::conditional<nDim == 2, NodeValue2D, NodeValue3D>::type;
  using RTree = typename std::conditional<nDim == 2, RTree2D, RTree3D>::type;

private:
  int rank;
  RTree nodeTree;  /*!< \brief R-tree for spatial indexing of mesh nodes. */
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
   * \brief Build the R-tree from mesh nodes.
   * \param[in] geometry - Mesh geometry.
   */
  void BuildTree(CGeometry* geometry) override;

  /*!
   * \brief Clear the tree and reset the built flag.
   */
  void ClearTree() override;

  /*!
   * \brief Check if the tree has been built.
   * \return True if tree is built, false otherwise.
   */
  bool IsTreeBuilt() const override { return treeBuilt; }

  /*!
   * \brief Get the number of nodes in the tree.
   * \return Number of nodes stored in the tree.
   */
  size_t GetTreeSize() const override { return nodeTree.size(); }

  /*!
   * \brief Set coordinates for a point.
   * \param[out] pt - Point to set coordinates for.
   * \param[in] coor - Array of coordinates.
   */
  void SetPoint(Point& pt, const su2double* coor) const;

  /*!
   * \brief Create a bounding box from element vertices.
   * \param[in] elemCoor - Array of element vertex coordinates.
   * \return Bounding box encompassing the element.
   */
  Box CreateBoundingBox(const su2double* elemCoor) const;

  /*!
   * \brief Create a bounding box from element vertices with padding.
   * \param[in] elemCoor - Array of element vertex coordinates.
   * \param[in] padding - Additional padding around the element.
   * \return Bounding box encompassing the element with padding.
   */
  Box CreateBoundingBox(const su2double* elemCoor, su2double padding) const;

  /*!
   * \brief Search for nodes within a bounding box.
   * \param[in] boundingBox - The bounding box to search within.
   * \param[out] containedNodes - Set of node IDs found within the bounding box.
   */
  void SearchNodesInBox(const Box& boundingBox, std::set<unsigned long>& containedNodes) const;

  /*!
   * \brief Search for nodes within a bounding box (element coordinates).
   * \param[in] elemCoor - Array of element vertex coordinates.
   * \param[out] containedNodes - Set of node IDs found within the bounding box.
   */
  void SearchNodesInElement(const su2double* elemCoor, std::set<unsigned long>& containedNodes) const override;

  /*!
   * \brief Search for nodes within a bounding box (element coordinates with padding).
   * \param[in] elemCoor - Array of element vertex coordinates.
   * \param[out] containedNodes - Set of node IDs found within the bounding box.
   * \param[in] padding - Additional padding around the element.
   */
  void SearchNodesInElement(const su2double* elemCoor, std::set<unsigned long>& containedNodes, su2double padding) const override;
};
