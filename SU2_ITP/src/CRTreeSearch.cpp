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
#include <limits>

/*--- Implementation of base class methods ---*/
void CRTreeSearchBase::SetPoint2D(Point2D& pt, const su2double* coor) const {
  bg::set<0>(pt, coor[0]);
  bg::set<1>(pt, coor[1]);
}

void CRTreeSearchBase::SetPoint3D(Point3D& pt, const su2double* coor) const {
  bg::set<0>(pt, coor[0]);
  bg::set<1>(pt, coor[1]);
  bg::set<2>(pt, coor[2]);
}

Box2D CRTreeSearchBase::CreateBoundingBox2D(const su2double* tri, su2double padding) const {
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
  Point2D minCorner(xmin, ymin);
  Point2D maxCorner(xmax, ymax);

  return Box2D(minCorner, maxCorner);
}

Box3D CRTreeSearchBase::CreateBoundingBox3D(const su2double* tet, su2double padding) const {
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
  Point3D minCorner(xmin, ymin, zmin);
  Point3D maxCorner(xmax, ymax, zmax);

  return Box3D(minCorner, maxCorner);
}

/*--- Template class implementation ---*/
template<unsigned short nDim>
CRTreeSearch<nDim>::CRTreeSearch(SU2_Comm MPICommunicator) : treeBuilt(false), surfaceTreeBuilt(false) {
  rank = SU2_MPI::GetRank();
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::BuildTree(CGeometry* geometry) {
  /*--- Clear any existing tree ---*/
  ClearTree();

  /*--- Get the number of points in the mesh ---*/
  const unsigned long nPoint = geometry->GetnPoint();

  /*--- Storage for point coordinates ---*/
  su2double coor[3];
  Point pt;

  /*--- Vector to store all node-point pairs ---*/
  std::vector<NodeValue> nodeValues;
  nodeValues.reserve(nPoint);

  /*--- Loop over all mesh nodes and add them to the tree ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Get node coordinates for all dimensions ---*/
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      coor[iDim] = geometry->nodes->GetCoord(iPoint, iDim);
    }

    /*--- Set point coordinates using template point type ---*/
    SetPoint(pt, coor);

    /*--- Add to vector ---*/
    nodeValues.emplace_back(pt, iPoint);
  }

  /*--- Build the R-tree ---*/
  nodeTree = RTree(nodeValues.begin(), nodeValues.end());
  treeBuilt = true;
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::ClearTree() {
  nodeTree.clear();
  treeBuilt = false;
  surfaceTreeBuilt = false;
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::SetPoint(Point& pt, const su2double* coor) const {
  if constexpr (nDim == 2) {
    SetPoint2D(pt, coor);
  } else if constexpr (nDim == 3) {
    SetPoint3D(pt, coor);
  }
}

template<unsigned short nDim>
typename CRTreeSearch<nDim>::Box CRTreeSearch<nDim>::CreateBoundingBox(const su2double* elemCoor) const {
  return CreateBoundingBox(elemCoor, 0.0);
}

template<unsigned short nDim>
typename CRTreeSearch<nDim>::Box CRTreeSearch<nDim>::CreateBoundingBox(const su2double* elemCoor, su2double padding) const {
  if constexpr (nDim == 2) {
    return CreateBoundingBox2D(elemCoor, padding);
  } else if constexpr (nDim == 3) {
    return CreateBoundingBox3D(elemCoor, padding);
  }
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
  nodeTree.query(bg::index::intersects(boundingBox), std::back_inserter(queryResults));

  /*--- Extract node IDs from query results ---*/
  for (const auto& nodeValue : queryResults) {
    containedNodes.insert(nodeValue.second);
  }
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::SearchNodesInElement(const su2double* elemCoor, std::set<unsigned long>& containedNodes) const {
  Box boundingBox = CreateBoundingBox(elemCoor);
  SearchNodesInBox(boundingBox, containedNodes);
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::SearchNodesInElement(const su2double* elemCoor, std::set<unsigned long>& containedNodes, su2double padding) const {
  Box boundingBox = CreateBoundingBox(elemCoor, padding);
  SearchNodesInBox(boundingBox, containedNodes);
}

template<unsigned short nDim>
void CRTreeSearch<nDim>::BuildSurfaceTree(CGeometry* geometry) {
  /*--- Clear any existing surface tree ---*/
  surfaceTree.clear();
  surfaceTreeBuilt = false;

  /*--- Storage for point coordinates ---*/
  su2double coor[3];
  Point pt;

  /*--- Map from nodeID to lists of markerIDs and elemIDs ---*/
  std::map<unsigned long, std::vector<unsigned short>> nodeToMarkers;
  std::map<unsigned long, std::vector<unsigned long>> nodeToElements;

  /*--- Loop over all boundary markers ---*/
  for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); ++iMarker) {
    /*--- Loop over all elements on this boundary marker ---*/
    for (auto iElem = 0ul; iElem < geometry->GetnElem_Bound(iMarker); ++iElem) {
      /*--- Loop over nodes of this boundary element ---*/
      for (auto iNode = 0u; iNode < geometry->bound[iMarker][iElem]->GetnNodes(); ++iNode) {
        const auto nodeID = geometry->bound[iMarker][iElem]->GetNode(iNode);

        /*--- Add this marker to the node's marker list ---*/
        nodeToMarkers[nodeID].push_back(static_cast<unsigned short>(iMarker));

        /*--- Add this element to the node's element list ---*/
        nodeToElements[nodeID].push_back(iElem);
      }
    }
  }

  /*--- Now build the vector of SurfaceNodeValue tuples ---*/
  std::vector<SurfaceNodeValue> surfaceNodeValues;
  surfaceNodeValues.reserve(nodeToMarkers.size());

  for (const auto& nodeEntry : nodeToMarkers) {
    const auto nodeID = nodeEntry.first;
    const auto& markerIDs = nodeEntry.second;
    const auto& elemIDs = nodeToElements[nodeID];

    /*--- Get node coordinates ---*/
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      coor[iDim] = geometry->nodes->GetCoord(nodeID, iDim);
    }

    /*--- Set point coordinates ---*/
    SetPoint(pt, coor);

    /*--- Create tuple: (point, nodeID, markerIDs, elemIDs) ---*/
    surfaceNodeValues.emplace_back(pt, nodeID, markerIDs, elemIDs);
  }

  /*--- Build the surface R-tree ---*/
  if (!surfaceNodeValues.empty()) {
    surfaceTree = SurfaceRTree(surfaceNodeValues.begin(), surfaceNodeValues.end());
    surfaceTreeBuilt = true;
  }
}

template<unsigned short nDim>
bool CRTreeSearch<nDim>::SearchNearestSurfaceNode(const su2double* coor, unsigned long& nearestNodeID,
                                                  std::vector<unsigned short>& markerIDs,
                                                  std::vector<unsigned long>& elemIDs) const {
  if (!surfaceTreeBuilt) {
    std::cerr << "Error: Surface R-tree has not been built. Call BuildSurfaceTree() first." << std::endl;
    return false;
  }

  /*--- Create query point ---*/
  Point queryPoint;
  SetPoint(queryPoint, coor);

  /*--- Search for nearest neighbor in surface tree ---*/
  std::vector<SurfaceNodeValue> result;
  surfaceTree.query(bg::index::nearest(queryPoint, 1), std::back_inserter(result));

  if (result.empty()) {
    return false;
  }

  /*--- Extract results from tuple: (point, nodeID, markerIDs, elemIDs) ---*/
  nearestNodeID = std::get<1>(result[0]);
  markerIDs = std::get<2>(result[0]);
  elemIDs = std::get<3>(result[0]);

  return true;
}

/*--- Explicit template instantiations for 2D and 3D ---*/
template class CRTreeSearch<2>;
template class CRTreeSearch<3>;