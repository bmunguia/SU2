/*!
 * \file CGMFFileWriter.hpp
 * \brief Filewriter class for GMF format solution.
 * \author B. Munguía
 * \version 8.2.0 "Harrier"
 * \bug
 *   Only supports serial mesh writing; parallel output is not
 *   supported due to GMF/libMeshb limitations.
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

#include "../../../include/output/filewriter/CGMFFileWriter.hpp"

const string CGMFFileWriter::fileExt = ".solb";

CGMFFileWriter::CGMFFileWriter(CParallelDataSorter *valDataSorter)  :
  CFileWriter(valDataSorter, fileExt){}


CGMFFileWriter::~CGMFFileWriter()= default;

void CGMFFileWriter::WriteData(string val_filename){
#ifdef HAVE_GMF
  /*--- TODO: currently doesn't work, since GmfOpenMesh seems to
              overwrite the existing mesh; hopefully libmeshb v8
              can resolve this ---*/
  val_filename.append(fileExt);

  const int nDim = dataSorter->GetnDim();
  const int ver = 2; // GMF file version

  /*--- Open the solution file for writing. ---*/
  int64_t sol_id = GmfOpenMesh(val_filename.c_str(), GmfWrite, ver, nDim);
  if (!sol_id) SU2_MPI::Error("Could not open GMF solution file for writing.", CURRENT_FUNCTION);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

  /*--- Prepare solution types and buffer. ---*/
  std::vector<int> sol_types;
  const auto& fieldNames = dataSorter->GetFieldNames();
  const unsigned short nFields = fieldNames.size();
  for (unsigned short iField = 0u; iField < nFields; ) {
    auto GmfFieldKwd = GetFieldKwd(fieldNames[iField], nDim);
    sol_types.push_back(GmfFieldKwd);
    auto nComp = GetFieldSize(fieldNames[iField], nDim);
    iField += nComp;
  }

  /*--- Set the keyword for the solution. ---*/
  if (rank == MASTER_NODE) {
    const unsigned long nGlobalPoints = dataSorter->GetnPointsGlobal();
    GmfSetKwd(sol_id, GmfSolAtVertices, nGlobalPoints, sol_types.size(), sol_types.data());
  }

  /*--- Write solution data for each point. ---*/
  std::vector<double> bufDbl;
  const unsigned long nLocalPoints = dataSorter->GetnPoints();
  for (int iProc = 0u; iProc < size; iProc++) {
    if (rank == iProc) {
      for (auto i = 0u; i < nLocalPoints; ++i) {
        bufDbl.clear();
        /*--- No special processing required. Just add all the field data. ---*/
        for (auto j = 0u; j < nFields; ++j) {
          bufDbl.push_back(dataSorter->GetData(j, i));
        }
        GmfSetLin(sol_id, GmfSolAtVertices, bufDbl.data());
      }
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  /*--- Close the solution file. ---*/
  GmfCloseMesh(sol_id);
#endif
}

#ifdef HAVE_GMF
int CGMFFileWriter::GetFieldKwd(const std::string& fieldname, int nDim) {
  /*--- Check for tensor components. ---*/
  if (fieldname.find("_xx") != std::string::npos) return GmfSymMat;

  /*--- Check for vector components. ---*/
  if (fieldname.find("_x") != std::string::npos && fieldname.find("_xx") == std::string::npos) return GmfVec;

  /*--- Otherwise, it's scalar. ---*/
  return GmfSca;
}

int CGMFFileWriter::GetFieldSize(const std::string& fieldname, int nDim) {
  /*--- Symmetric tensor: 3 components for 2D, 6 for 3D. ---*/
  if (fieldname.find("_xx") != std::string::npos) return 3 * (nDim - 1);

  /*--- Vector: 2 components for 2D, 3 for 3D. ---*/
  if (fieldname.find("_x") != std::string::npos && fieldname.find("_xx") == std::string::npos) return nDim;

  /*--- Scalar: 1 component. ---*/
  return 1;
}
#endif