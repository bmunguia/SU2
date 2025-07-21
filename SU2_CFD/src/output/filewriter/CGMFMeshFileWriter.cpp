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


const string CGMFMeshFileWriter::fileExt = ".mesh";

CGMFMeshFileWriter::CGMFMeshFileWriter(CParallelDataSorter* valVolumeSorter,
                                       CParallelDataSorter* valSurfaceSorter,
                                       unsigned short valiZone, unsigned short valnZone)
  : CFileWriter(valVolumeSorter, fileExt), iZone(valiZone), nZone(valnZone),
    surfaceSorter(valSurfaceSorter) {}

void CGMFMeshFileWriter::WriteData(string val_filename) {
#ifdef HAVE_GMF
  /*--- TODO: currently doesn't work, since GmfOpenMesh seems to
              overwrite the existing mesh ---*/
  val_filename.append(fileExt);

  unsigned long nGlobalPoints = dataSorter->GetnPointsGlobal();
  unsigned nDim = dataSorter->GetnDim();
  const int ver = 2; // GMF file version

  if (rank == MASTER_NODE) {
    int64_t mesh_id = GmfOpenMesh(val_filename.c_str(), GmfWrite, ver, nDim);
    if (!mesh_id) SU2_MPI::Error("Could not open GMF file for writing.", CURRENT_FUNCTION);
    GmfSetKwd(mesh_id, GmfVertices, nGlobalPoints);
    /*--- Set element counts for boundaries and volume ---*/
    if (nDim == 2) {
      /*--- 2D: Triangles, quadrilaterals, edges ---*/
      unsigned long nTri = dataSorter->GetnElemGlobal(TRIANGLE);
      unsigned long nQuad = dataSorter->GetnElemGlobal(QUADRILATERAL);
      unsigned long nEdge = surfaceSorter->GetnElemGlobal(LINE);

      if (nTri > 0) GmfSetKwd(mesh_id, GmfTriangles, nTri);
      if (nQuad > 0) GmfSetKwd(mesh_id, GmfQuadrilaterals, nQuad);
      if (nEdge > 0) GmfSetKwd(mesh_id, GmfEdges, nEdge);
    } else {
      /*--- 3D: Tetrahedra, hexahedra, prisms, pyramids, triangles, and quadrilaterals ---*/
      unsigned long nTet = dataSorter->GetnElemGlobal(TETRAHEDRON);
      unsigned long nHex = dataSorter->GetnElemGlobal(HEXAHEDRON);
      unsigned long nPri = dataSorter->GetnElemGlobal(PRISM);
      unsigned long nPyr = dataSorter->GetnElemGlobal(PYRAMID);
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
      int64_t mesh_id = GmfOpenMesh(val_filename.c_str(), GmfWrite, ver, nDim);
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
  unsigned long nLocalPoints = dataSorter->GetnPoints();
  for (unsigned long i = 0; i < nLocalPoints; ++i) {
    double x = dataSorter->GetData(0, i);
    double y = dataSorter->GetData(1, i);
    double z = (nDim == 3) ? dataSorter->GetData(2, i) : 0.0;
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
    unsigned long nTri = dataSorter->GetnElem(TRIANGLE);
    for (unsigned long i = 0; i < nTri; ++i) {
      int v1 = dataSorter->GetElemConnectivity(TRIANGLE, i, 0);
      int v2 = dataSorter->GetElemConnectivity(TRIANGLE, i, 1);
      int v3 = dataSorter->GetElemConnectivity(TRIANGLE, i, 2);
      int ref = 0;
      GmfSetLin(mesh_id, GmfTriangles, v1, v2, v3, ref);
    }
    unsigned long nQuad = dataSorter->GetnElem(QUADRILATERAL);
    for (unsigned long i = 0; i < nQuad; ++i) {
      int v1 = dataSorter->GetElemConnectivity(QUADRILATERAL, i, 0);
      int v2 = dataSorter->GetElemConnectivity(QUADRILATERAL, i, 1);
      int v3 = dataSorter->GetElemConnectivity(QUADRILATERAL, i, 2);
      int v4 = dataSorter->GetElemConnectivity(QUADRILATERAL, i, 3);
      int ref = 0;
      GmfSetLin(mesh_id, GmfQuadrilaterals, v1, v2, v3, v4, ref);
    }
  } else {
    unsigned long nTet = dataSorter->GetnElem(TETRAHEDRON);
    for (unsigned long i = 0; i < nTet; ++i) {
      int v1 = dataSorter->GetElemConnectivity(TETRAHEDRON, i, 0);
      int v2 = dataSorter->GetElemConnectivity(TETRAHEDRON, i, 1);
      int v3 = dataSorter->GetElemConnectivity(TETRAHEDRON, i, 2);
      int v4 = dataSorter->GetElemConnectivity(TETRAHEDRON, i, 3);
      int ref = 0;
      GmfSetLin(mesh_id, GmfTetrahedra, v1, v2, v3, v4, ref);
    }
    unsigned long nHex = dataSorter->GetnElem(HEXAHEDRON);
    for (unsigned long i = 0; i < nHex; ++i) {
      int v1 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 0);
      int v2 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 1);
      int v3 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 2);
      int v4 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 3);
      int v5 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 4);
      int v6 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 5);
      int v7 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 6);
      int v8 = dataSorter->GetElemConnectivity(HEXAHEDRON, i, 7);
      int ref = 0;
      GmfSetLin(mesh_id, GmfHexahedra, v1, v2, v3, v4, v5, v6, v7, v8, ref);
    }
    unsigned long nPri = dataSorter->GetnElem(PRISM);
    for (unsigned long i = 0; i < nPri; ++i) {
      int v1 = dataSorter->GetElemConnectivity(PRISM, i, 0);
      int v2 = dataSorter->GetElemConnectivity(PRISM, i, 1);
      int v3 = dataSorter->GetElemConnectivity(PRISM, i, 2);
      int v4 = dataSorter->GetElemConnectivity(PRISM, i, 3);
      int v5 = dataSorter->GetElemConnectivity(PRISM, i, 4);
      int v6 = dataSorter->GetElemConnectivity(PRISM, i, 5);
      int ref = 0;
      GmfSetLin(mesh_id, GmfPrisms, v1, v2, v3, v4, v5, v6, ref);
    }
    unsigned long nPyr = dataSorter->GetnElem(PYRAMID);
    for (unsigned long i = 0; i < nPyr; ++i) {
      int v1 = dataSorter->GetElemConnectivity(PYRAMID, i, 0);
      int v2 = dataSorter->GetElemConnectivity(PYRAMID, i, 1);
      int v3 = dataSorter->GetElemConnectivity(PYRAMID, i, 2);
      int v4 = dataSorter->GetElemConnectivity(PYRAMID, i, 3);
      int v5 = dataSorter->GetElemConnectivity(PYRAMID, i, 4);
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
