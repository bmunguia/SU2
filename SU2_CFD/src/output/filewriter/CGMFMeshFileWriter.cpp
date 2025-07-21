/*!
 * \file CGMFMeshFileWriter.hpp
 * \brief Filewriter class for GMF format mesh.
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

#include "../../../include/output/filewriter/CGMFMeshFileWriter.hpp"


const string CGMFMeshFileWriter::fileExt = ".meshb";

CGMFMeshFileWriter::CGMFMeshFileWriter(const CFVMDataSorter* valVolumeSorter,
                                       const CSurfaceFVMDataSorter* valSurfaceSorter,
                                       const std::vector<std::string>& markerList_,
                                       unsigned short valiZone, unsigned short valnZone)
  : CFileWriter(valVolumeSorter, fileExt), iZone(valiZone), nZone(valnZone),
    surfaceSorter(valSurfaceSorter), markerList(markerList_) {}

void CGMFMeshFileWriter::WriteData(string val_filename) {
#ifdef HAVE_GMF
  val_filename.append(fileExt);

  unsigned long nGlobalPoints = volumeSorter->GetnPointsGlobal();
  unsigned nDim = volumeSorter->GetnDim();
  const int ver = 2; // GMF file version

  if (rank == 0) {
    int64_t mesh_id = GmfOpenMesh(val_filename.c_str(), GmfWrite, ver, nDim);
    if (!mesh_id) SU2_MPI::Error("Could not open GMF file for writing.", CURRENT_FUNCTION);
    GmfSetKwd(mesh_id, GmfVertices, nGlobalPoints);
    /*--- Set element counts for boundaries and volume ---*/
    if (nDim == 2) {
      /*--- 2D: Triangles, quadrilaterals, edges ---*/
      unsigned long nTri = volumeSorter->GetnElemGlobal(TRIANGLE);
      unsigned long nQuad = volumeSorter->GetnElemGlobal(QUADRILATERAL);
      unsigned long nEdge = surfaceSorter->GetnElemGlobal(LINE);

      if (nTri > 0) GmfSetKwd(mesh_id, GmfTriangles, nTri);
      if (nQuad > 0) GmfSetKwd(mesh_id, GmfQuadrilaterals, nQuad);
      if (nEdge > 0) GmfSetKwd(mesh_id, GmfEdges, nEdge);
    } else {
      /*--- 3D: Tetrahedra, hexahedra, prisms, pyramids, triangles, and quadrilaterals ---*/
      unsigned long nTet = volumeSorter->GetnElemGlobal(TETRAHEDRON);
      unsigned long nHex = volumeSorter->GetnElemGlobal(HEXAHEDRON);
      unsigned long nPri = volumeSorter->GetnElemGlobal(PRISM);
      unsigned long nPyr = volumeSorter->GetnElemGlobal(PYRAMID);
      unsigned long nTri = surfaceSorter->GetnElemGlobal(TRIANGLE);
      unsigned long nQuad = surfaceSorter->GetnElemGlobal(QUADRILATERAL);

      if (nTet > 0) GmfSetKwd(mesh_id, GmfTetrahedra, nTet);
      if (nHex > 0) GmfSetKwd(mesh_id, GmfHexahedra, nHex);
      if (nPri > 0) GmfSetKwd(mesh_id, GmfPrisms, nPri);
      if (nPyr > 0) GmfSetKwd(mesh_id, GmfPyramids, nPyr);
      if (nTri > 0) GmfSetKwd(mesh_id, GmfTriangles, nTri);
      if (nQuad > 0) GmfSetKwd(mesh_id, GmfQuadrilaterals, nQuad);
    }
    GmfCloseMesh(mesh_id);
  }
  SU2_MPI::Barrier(SU2_MPI::GetComm());

  for (int iProc = 0; iProc < size; ++iProc) {
    if (rank == iProc) {
      int64_t mesh_id = GmfOpenMesh(val_filename.c_str(), GmfReadWrite, ver, nDim);
      if (!mesh_id) SU2_MPI::Error("Could not open GMF file for writing.", CURRENT_FUNCTION);
      WritePoints(mesh_id, nDim);
      WriteElements(mesh_id, nDim);
      WriteBoundaryElements(mesh_id, nDim);
      GmfCloseMesh(mesh_id);
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
#endif
}

#ifdef HAVE_GMF
void CGMFMeshFileWriter::WritePoints(int64_t mesh_id, unsigned short nDim) {
  unsigned long nLocalPoints = volumeSorter->GetnPoints();
  for (unsigned long i = 0; i < nLocalPoints; ++i) {
    double x = volumeSorter->GetData(0, i);
    double y = volumeSorter->GetData(1, i);
    double z = (nDim == 3) ? volumeSorter->GetData(2, i) : 0.0;
    int ref = 0;
    if (nDim == 2)
      GmfSetLin(mesh_id, GmfVertices, x, y, ref);
    else
      GmfSetLin(mesh_id, GmfVertices, x, y, z, ref);
  }
}

void CGMFMeshFileWriter::WriteElements(int64_t mesh_id, unsigned short nDim) {
  /*--- Volume elements with ref = 0 ---*/
  if (nDim == 2) {
    unsigned long nTri = volumeSorter->GetnElem(TRIANGLE);
    for (unsigned long i = 0; i < nTri; ++i) {
      int v1 = volumeSorter->GetElemConnectivity(TRIANGLE, i, 0);
      int v2 = volumeSorter->GetElemConnectivity(TRIANGLE, i, 1);
      int v3 = volumeSorter->GetElemConnectivity(TRIANGLE, i, 2);
      int ref = 0;
      GmfSetLin(mesh_id, GmfTriangles, v1, v2, v3, ref);
    }
    unsigned long nQuad = volumeSorter->GetnElem(QUADRILATERAL);
    for (unsigned long i = 0; i < nQuad; ++i) {
      int v1 = volumeSorter->GetElemConnectivity(QUADRILATERAL, i, 0);
      int v2 = volumeSorter->GetElemConnectivity(QUADRILATERAL, i, 1);
      int v3 = volumeSorter->GetElemConnectivity(QUADRILATERAL, i, 2);
      int v4 = volumeSorter->GetElemConnectivity(QUADRILATERAL, i, 3);
      int ref = 0;
      GmfSetLin(mesh_id, GmfQuadrilaterals, v1, v2, v3, v4, ref);
    }
  } else {
    unsigned long nTet = volumeSorter->GetnElem(TETRAHEDRON);
    for (unsigned long i = 0; i < nTet; ++i) {
      int v1 = volumeSorter->GetElemConnectivity(TETRAHEDRON, i, 0);
      int v2 = volumeSorter->GetElemConnectivity(TETRAHEDRON, i, 1);
      int v3 = volumeSorter->GetElemConnectivity(TETRAHEDRON, i, 2);
      int v4 = volumeSorter->GetElemConnectivity(TETRAHEDRON, i, 3);
      int ref = 0;
      GmfSetLin(mesh_id, GmfTetrahedra, v1, v2, v3, v4, ref);
    }
    unsigned long nHex = volumeSorter->GetnElem(HEXAHEDRON);
    for (unsigned long i = 0; i < nHex; ++i) {
      int v1 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 0);
      int v2 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 1);
      int v3 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 2);
      int v4 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 3);
      int v5 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 4);
      int v6 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 5);
      int v7 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 6);
      int v8 = volumeSorter->GetElemConnectivity(HEXAHEDRON, i, 7);
      int ref = 0;
      GmfSetLin(mesh_id, GmfHexahedra, v1, v2, v3, v4, v5, v6, v7, v8, ref);
    }
    unsigned long nPri = volumeSorter->GetnElem(PRISM);
    for (unsigned long i = 0; i < nPri; ++i) {
      int v1 = volumeSorter->GetElemConnectivity(PRISM, i, 0);
      int v2 = volumeSorter->GetElemConnectivity(PRISM, i, 1);
      int v3 = volumeSorter->GetElemConnectivity(PRISM, i, 2);
      int v4 = volumeSorter->GetElemConnectivity(PRISM, i, 3);
      int v5 = volumeSorter->GetElemConnectivity(PRISM, i, 4);
      int v6 = volumeSorter->GetElemConnectivity(PRISM, i, 5);
      int ref = 0;
      GmfSetLin(mesh_id, GmfPrisms, v1, v2, v3, v4, v5, v6, ref);
    }
    unsigned long nPyr = volumeSorter->GetnElem(PYRAMID);
    for (unsigned long i = 0; i < nPyr; ++i) {
      int v1 = volumeSorter->GetElemConnectivity(PYRAMID, i, 0);
      int v2 = volumeSorter->GetElemConnectivity(PYRAMID, i, 1);
      int v3 = volumeSorter->GetElemConnectivity(PYRAMID, i, 2);
      int v4 = volumeSorter->GetElemConnectivity(PYRAMID, i, 3);
      int v5 = volumeSorter->GetElemConnectivity(PYRAMID, i, 4);
      int ref = 0;
      GmfSetLin(mesh_id, GmfPyramids, v1, v2, v3, v4, v5, ref);
    }
  }
}

void CGMFMeshFileWriter::WriteBoundaryElements(int64_t mesh_id, unsigned short nDim) {
   /*--- Surface elements with ref = iMarker + 1 ---*/
  if (nDim == 2) {
    unsigned long nEdge = surfaceSorter->GetnElem(LINE);
    for (unsigned long i = 0; i < nEdge; ++i) {
      int v1 = surfaceSorter->GetElemConnectivity(LINE, i, 0);
      int v2 = surfaceSorter->GetElemConnectivity(LINE, i, 1);
      int iMarker = surfaceSorter->GetElemMarkerID(LINE, i);
      int ref = iMarker + 1;
      GmfSetLin(mesh_id, GmfEdges, v1, v2, ref);
    }
  } else {
    unsigned long nTri = surfaceSorter->GetnElem(TRIANGLE);
    for (unsigned long i = 0; i < nTri; ++i) {
      int v1 = surfaceSorter->GetElemConnectivity(TRIANGLE, i, 0);
      int v2 = surfaceSorter->GetElemConnectivity(TRIANGLE, i, 1);
      int v3 = surfaceSorter->GetElemConnectivity(TRIANGLE, i, 2);
      int iMarker = surfaceSorter->GetElemMarkerID(TRIANGLE, i);
      int ref = iMarker + 1;
      GmfSetLin(mesh_id, GmfTriangles, v1, v2, v3, ref);
    }
    unsigned long nQuad = surfaceSorter->GetnElem(QUADRILATERAL);
    for (unsigned long i = 0; i < nQuad; ++i) {
      int v1 = surfaceSorter->GetElemConnectivity(QUADRILATERAL, i, 0);
      int v2 = surfaceSorter->GetElemConnectivity(QUADRILATERAL, i, 1);
      int v3 = surfaceSorter->GetElemConnectivity(QUADRILATERAL, i, 2);
      int v4 = surfaceSorter->GetElemConnectivity(QUADRILATERAL, i, 3);
      int iMarker = surfaceSorter->GetElemMarkerID(QUADRILATERAL, i);
      int ref = iMarker + 1;
      GmfSetLin(mesh_id, GmfQuadrilaterals, v1, v2, v3, v4, ref);
    }
  }
}
#endif  // HAVE_GMF
