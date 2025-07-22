/*!
 * \file CGMFMeshFileWriter.hpp
 * \brief Headers for GMF mesh file writer class.
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

#ifdef HAVE_GMF
extern "C" {
#include "libmeshb7.h"
}
#endif

#include "CFileWriter.hpp"

class CGMFMeshFileWriter final: public CFileWriter{

private:
  unsigned short iZone, nZone;

#ifdef HAVE_GMF
  const CParallelDataSorter* surfaceSorter;
  std::vector<std::string> markerList;
#endif
public:
  const static string fileExt;

  /*!
   * \brief Construct a file writer using both volume and surface data sorters and marker list.
   * \param[in] valVolumeSorter - The parallel sorted volume data
   * \param[in] valSurfaceSorter - The parallel sorted surface data
   * \param[in] valiZone - The index of the current zone
   * \param[in] valnZone - The total number of zones
   */
  CGMFMeshFileWriter(CParallelDataSorter* valVolumeSorter,
                     CParallelDataSorter* valSurfaceSorter,
                     unsigned short valiZone, unsigned short valnZone);

  /*!
   * \brief Write sorted data to file
   * \param[in] val_filename - The name of the file
   */
  void WriteData(string val_filename) override;

#ifdef HAVE_GMF
  /*!
   * \brief Write mesh vertex coordinates to GMF file.
   * \param[in] mesh_id - GMF mesh file handle
   * \param[in] nDim - Mesh spatial dimension
   */
  void WritePoints(int64_t mesh_id, unsigned short nDim);

  /*!
   * \brief Write mesh elements of a given type to GMF file.
   * \param[in] mesh_id - GMF mesh file handle
   * \param[in] type - Element type
   * \param[in] nDim - Mesh spatial dimension
   * \param[in] isSurf - True if writing surface/boundary elements
   */
  void WriteElements(int64_t mesh_id, GEO_TYPE type, unsigned short nDim, bool isSurf = false);

  /*!
   * \brief Get GMF keyword for a given element type.
   * \param[in] type - Element type
   * \return GMF keyword integer for the element type
   */
  int GetElementKwd(GEO_TYPE type);
#endif
};