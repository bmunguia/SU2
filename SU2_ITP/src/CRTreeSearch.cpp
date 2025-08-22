/*!
 * \file CRTreeSearch.cpp
 * \brief Implementation of the spatial search class using Boost R-tree.
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

#include "../include/CRTreeSearch.hpp"
#include <algorithm>
#include <iostream>
#include <limits>

template<unsigned short nDim>
CRTreeSearch<nDim>::CRTreeSearch(SU2_Comm MPICommunicator) : treeBuilt(false) {
  rank = SU2_MPI::GetRank();
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::BuildTree(CGeometry* geometry_src) {
  /*--- Clear any existing tree ---*/
  ClearTree();

  if (rank == MASTER_NODE) {
    std::cout << "Building R-tree for source mesh nodes." << std::flush;
  }

  /*--- Get the number of points in the source mesh ---*/
  const unsigned long nPoint = geometry_src->GetnPoint();

  /*--- Vector to store all node-point pairs for bulk loading ---*/
  std::vector<NodeValue> nodeValues;
  nodeValues.reserve(nPoint);

  /*--- Loop over all source mesh nodes and add them to the tree ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Get node coordinates for all dimensions ---*/
    Point pt;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      boost::geometry::set<iDim>(pt, geometry_src->nodes->GetCoord(iPoint, iDim));
    }

    /*--- Add to vector ---*/
    nodeValues.emplace_back(pt, iPoint);
  }

  /*--- Build the R-tree using bulk loading for better performance ---*/
  nodeTree = RTree(nodeValues.begin(), nodeValues.end());
  treeBuilt = true;

  if (rank == MASTER_NODE) {
    std::cout << " Done. Added " << nodeTree.size() << " nodes to R-tree." << std::endl;
  }
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::ClearTree() {
  nodeTree.clear();
  treeBuilt = false;
}

template<unsigned short nDim>
typename CRTreeSearch<nDim>::Box CRTreeSearch<nDim>::CreateBoundingBox(const su2double* elemCoords) const {
  return CreateBoundingBox(elemCoords, 0.0);
}

template<unsigned short nDim>
typename CRTreeSearch<nDim>::Box CRTreeSearch<nDim>::CreateBoundingBox(const su2double* elemCoords, su2double padding) const {
  if (nDim == 2) {
    return CreateBoundingBox2D(elemCoords, padding);
  } else if (nDim == 3) {
    return CreateBoundingBox3D(elemCoords, padding);
  } else {
    std::cerr << "Error: Unsupported dimension " << nDim << " in CreateBoundingBox." << std::endl;
    return Box();
  }
}

template<unsigned short nDim>
typename CRTreeSearch<nDim>::Box CRTreeSearch<nDim>::CreateBoundingBox2D(const su2double* tri, su2double padding) const {
  /*--- Extract triangle vertices ---*/
  const su2double x0 = tri[0], y0 = tri[1];
  const su2double x1 = tri[2], y1 = tri[3];
  const su2double x2 = tri[4], y2 = tri[5];

  /*--- Find bounding box coordinates ---*/
  const su2double xmin = std::min({x0, x1, x2}) - padding;
  const su2double xmax = std::max({x0, x1, x2}) + padding;
  const su2double ymin = std::min({y0, y1, y2}) - padding;
  const su2double ymax = std::max({y0, y1, y2}) + padding;

  /*--- Create points for bounding box ---*/
  Point minCorner, maxCorner;
  boost::geometry::set<0>(minCorner, xmin);
  boost::geometry::set<1>(minCorner, ymin);
  boost::geometry::set<0>(maxCorner, xmax);
  boost::geometry::set<1>(maxCorner, ymax);

  return Box(minCorner, maxCorner);
}

template<unsigned short nDim>
typename CRTreeSearch<nDim>::Box CRTreeSearch<nDim>::CreateBoundingBox3D(const su2double* tet, su2double padding) const {
  /*--- Extract tetrahedron vertices ---*/
  const su2double x0 = tet[0], y0 = tet[1], z0 = tet[2];
  const su2double x1 = tet[3], y1 = tet[4], z1 = tet[5];
  const su2double x2 = tet[6], y2 = tet[7], z2 = tet[8];
  const su2double x3 = tet[9], y3 = tet[10], z3 = tet[11];

  /*--- Find bounding box coordinates ---*/
  const su2double xmin = std::min({x0, x1, x2, x3}) - padding;
  const su2double xmax = std::max({x0, x1, x2, x3}) + padding;
  const su2double ymin = std::min({y0, y1, y2, y3}) - padding;
  const su2double ymax = std::max({y0, y1, y2, y3}) + padding;
  const su2double zmin = std::min({z0, z1, z2, z3}) - padding;
  const su2double zmax = std::max({z0, z1, z2, z3}) + padding;

  /*--- Create points for bounding box ---*/
  Point minCorner, maxCorner;
  boost::geometry::set<0>(minCorner, xmin);
  boost::geometry::set<1>(minCorner, ymin);
  boost::geometry::set<2>(minCorner, zmin);
  boost::geometry::set<0>(maxCorner, xmax);
  boost::geometry::set<1>(maxCorner, ymax);
  boost::geometry::set<2>(maxCorner, zmax);

  return Box(minCorner, maxCorner);
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::SearchNodesInBox(const Box& boundingBox, std::set<unsigned long>& containedNodes) const {
  if (!treeBuilt) {
    std::cerr << "Error: R-tree has not been built. Call BuildTree() first." << std::endl;
    return;
  }

  /*--- Clear output set ---*/
  containedNodes.clear();

  /*--- Query the R-tree for all nodes within the bounding box ---*/
  std::vector<NodeValue> queryResults;
  nodeTree.query(boost::geometry::index::within(boundingBox), std::back_inserter(queryResults));

  /*--- Extract node IDs from query results ---*/
  for (const auto& nodeValue : queryResults) {
    containedNodes.insert(nodeValue.second);
  }
}

/*--- Explicit template instantiations for 2D and 3D ---*/
template class CRTreeSearch<2>;
template class CRTreeSearch<3>;