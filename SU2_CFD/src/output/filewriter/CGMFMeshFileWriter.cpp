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

  /*--- Open the mesh file for writing. */
  int64_t mesh_id = GmfOpenMesh(val_filename.c_str(), GmfWrite, ver, nDim);
  if (!mesh_id) SU2_MPI::Error("Could not open GMF file for writing.", CURRENT_FUNCTION);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

  /*--- Write global number of points and coordinates. ---*/
  WritePoints(mesh_id, nDim);
  /*--- Write volume and boundary elements. */
  if (nDim == 2) {
    /*--- 2D: Triangles, quadrilaterals, edges ---*/
    WriteElements(mesh_id, TRIANGLE, nDim);
    WriteElements(mesh_id, QUADRILATERAL, nDim);
    WriteElements(mesh_id, LINE, nDim, true);
  } else {
    /*--- 3D: Tetrahedra, hexahedra, prisms, pyramids, triangles, and quadrilaterals ---*/
    WriteElements(mesh_id, TETRAHEDRON, nDim);
    WriteElements(mesh_id, HEXAHEDRON, nDim);
    WriteElements(mesh_id, PRISM, nDim);
    WriteElements(mesh_id, PYRAMID, nDim);
    WriteElements(mesh_id, TRIANGLE, nDim, true);
    WriteElements(mesh_id, QUADRILATERAL, nDim, true);
  }

  /*--- Close the mesh file. ---*/
  GmfCloseMesh(mesh_id);
#endif
}

#ifdef HAVE_GMF
void CGMFMeshFileWriter::WritePoints(int64_t mesh_id, unsigned short nDim) {
  if (rank == MASTER_NODE) {
    unsigned long nGlobalPoints = dataSorter->GetnPointsGlobal();
    GmfSetKwd(mesh_id, GmfVertices, nGlobalPoints);
  }
  unsigned long nLocalPoints = dataSorter->GetnPoints();
  for (int iProc = 0; iProc < size; iProc++) {
    if (rank == iProc) {
      for (unsigned long i = 0; i < nLocalPoints; i++) {
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
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
}

void CGMFMeshFileWriter::WriteElements(int64_t mesh_id, GEO_TYPE type, unsigned short nDim, bool isSurf) {
  const CParallelDataSorter* sorter = isSurf ? surfaceSorter : dataSorter;
  unsigned long nGlobalElems = sorter->GetnElemGlobal(type);

  if (nGlobalElems == 0) return;

  /*--- Set keyword for element type. ---*/
  auto GmfKwd = GetElementKwd(type);
  if (rank == MASTER_NODE) GmfSetKwd(mesh_id, GmfKwd, nGlobalElems);

  /*--- Default to ref=0 for volume elements. ---*/
  int ref = 0;
  int v[8]; // max vertices for hexahedron
  int nNodes = nPointsOfElementType(type);
  unsigned long nLocalElems = sorter->GetnElem(type);
  for (int iProc = 0; iProc < size; iProc++) {
    if (rank == iProc) {
      for (auto i = 0u; i < nLocalElems; i++) {
        /*--- Get element vertices. ---*/
        for (auto j = 0u; j < nNodes; j++) {
          v[j] = sorter->GetElemConnectivity(type, i, j);
        }
        if (isSurf) {
          /*--- Update ref for surface elements. We store iMarker+1 as the MarkerID
                in the surface datasorter. ---*/
          ref = sorter->GetElemMarkerID(type, i);
        }
        switch(type) {
          case LINE: {
            GmfSetLin(mesh_id, GmfKwd, v[0], v[1], ref);
            break;
          }
          case TRIANGLE: {
            GmfSetLin(mesh_id, GmfKwd, v[0], v[1], v[2], ref);
            break;
          }
          case QUADRILATERAL: {
            GmfSetLin(mesh_id, GmfKwd, v[0], v[1], v[2], v[3], ref);
            break;
          }
          case TETRAHEDRON: {
            GmfSetLin(mesh_id, GmfKwd, v[0], v[1], v[2], v[3], ref);
            break;
          }
          case HEXAHEDRON: {
            GmfSetLin(mesh_id, GmfKwd, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], ref);
            break;
          }
          case PRISM: {
            GmfSetLin(mesh_id, GmfKwd, v[0], v[1], v[2], v[3], v[4], v[5], ref);
            break;
          }
          case PYRAMID: {
            GmfSetLin(mesh_id, GmfKwd, v[0], v[1], v[2], v[3], v[4], ref);
            break;
          }
          default:
            break;
        }
      }
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
}

int CGMFMeshFileWriter::GetElementKwd(GEO_TYPE type) {
  switch (type) {
    case LINE:
      return GmfEdges;
      break;
    case TRIANGLE:
      return GmfTriangles;
      break;
    case QUADRILATERAL:
      return GmfQuadrilaterals;
      break;
    case TETRAHEDRON:
      return GmfTetrahedra;
      break;
    case HEXAHEDRON:
      return GmfHexahedra;
      break;
    case PRISM:
      return GmfPrisms;
      break;
    case PYRAMID:
      return GmfPyramids;
      break;
    default:
      return -1;
      break;
  }
}
#endif  // HAVE_GMF
