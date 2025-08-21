/*!
 * \file CLinearVolumeInterpolator.cpp
 * \brief Implementation of the main linear solution interpolation subroutines.
 *        This file contains all core interpolation logic and utility functions.
 * \author B. Munguía, E. van der Weide
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

#include "../include/CLinearVolumeInterpolator.hpp"

CLinearVolumeInterpolator::CLinearVolumeInterpolator(SU2_Comm MPICommunicator)
    : CVolumeInterpolator(MPICommunicator) { }

void CLinearVolumeInterpolator::Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                            CSolver** solver_container_src, CSolver** solver_container_dst) {
  if (rank == MASTER_NODE) {
    cout << endl << "----------------------------- Interpolation -----------------------------" << endl;
    cout << "Linear solution interpolation from source mesh to destination mesh." << endl;
    cout << "Source mesh: " << geometry_src->GetGlobal_nPointDomain() << " points, ";
    cout << geometry_src->GetGlobal_nElemDomain() << " elements." << endl;
    cout << "Destination mesh: " << geometry_dst->GetGlobal_nPointDomain() << " points,";
    cout << geometry_dst->GetGlobal_nElemDomain() << " elements." << endl;
  }

  /*--- Build the ADTs ---*/
  InitializeADTs(config, geometry_src, geometry_dst);

  /*--- Call the internal interpolation method ---*/
  LinearInterpolation(config, geometry_src, geometry_dst, solver_container_src[FLOW_SOL], solver_container_dst[FLOW_SOL]);
}

void CLinearVolumeInterpolator::LinearInterpolation(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                                    CSolver* solver_src, CSolver* solver_dst) {
  /*--------------------------------------------------------------------------*/
  /*--- Step 0: Initialize destination coordinate vector. If applying a    ---*/
  /*---         curvature correction, this vector will be modified.        ---*/
  /*---         Otherwise, it will just contain the original coordinates.  ---*/
  /*--------------------------------------------------------------------------*/
  nVar = solver_src->GetnVar();
  vector<su2double> coorDst;
  InitializeCoords(geometry_dst, coorDst);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Apply the curvature correction to the destination nodes.   ---*/
  /*--------------------------------------------------------------------------*/
  // if (rank == MASTER_NODE) cout << "Applying curvature correction." << endl;
  // ApplyCurvatureCorrection(config, geometry_src, geometry_dst, nDim, coorDst);

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Volume interpolation, via a containment search.            ---*/
  /*--------------------------------------------------------------------------*/

  if (rank == MASTER_NODE) cout << "Performing volume interpolation." << endl;
  vector<unsigned long> pointsFailed;
  VolumeInterpolation(geometry_src, solver_src, solver_dst, coorDst, pointsFailed);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Carry out a surface interpolation, via a minimum distance  ---*/
  /*---         search, for the points that could not be interpolated via  ---*/
  /*---         the regular volume interpolation. Print a warning.         ---*/
  /*--------------------------------------------------------------------------*/
  if (pointsFailed.size()) {
    if (rank == MASTER_NODE) {
      cout << pointsFailed.size() << " DOFs for which the containment search failed." << endl;
      cout << "A minimum distance search to the boundary of the domain is used for these points. " << endl;
    }
    unsigned long nPointsBeforeSurface = pointsFailed.size();
    SurfaceInterpolation(geometry_src, geometry_dst, solver_src, solver_dst, pointsFailed);
  }
}

void CLinearVolumeInterpolator::VolumeInterpolation(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                                                    const vector<su2double>& coor_dst, vector<unsigned long>& pointsFailed) {
  /*--- Search for donor elements for the given coordinates ---*/
  CADTElemClass& volumeADT = GetSourceVolumeADT();
  const unsigned long nDOFsDst = coor_dst.size() / nDim;

  /*--- Loop over the DOFs to be interpolated ---*/
  pointsFailed.clear();
  for (auto l = 0u; l < nDOFsDst; ++l) {
    /*--- Set a pointer to the coordinates to be searched ---*/
    const su2double* coor = coor_dst.data() + l * nDim;

    /*--- Carry out the containment search and check if it was successful ---*/
    unsigned short subElemID;
    unsigned long elemID;
    int rankID;
    su2double parCoor[3], weightsInterpol[8];

    if (volumeADT.DetermineContainingElement(coor, subElemID, elemID, rankID, parCoor,
                                             weightsInterpol)) {
      /*--- Get element information ---*/
      unsigned short nNodes = geometry_src->elem[elemID]->GetnNodes();

      /*--- Initialize interpolated solution to zero ---*/
      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        solver_dst->GetNodes()->SetSolution(l, iVar, 0.0);
      }

      /*--- Interpolate using shape function weights ---*/
      for (auto iNode = 0u; iNode < nNodes; ++iNode) {
        unsigned long nodeID = geometry_src->elem[elemID]->GetNode(iNode);

        for (auto iVar = 0u; iVar < nVar; ++iVar) {
          su2double val = solver_src->GetNodes()->GetSolution(nodeID, iVar);
          solver_dst->GetNodes()->Add_DeltaSolution(l, iVar, weightsInterpol[iNode] * val);
        }
      }
    } else {
      /*--- Containment search failed - store the index ---*/
      pointsFailed.push_back(l);
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Volume search finished. " << pointsFailed.size() << " points failed." << endl << flush;
  }
}

void CLinearVolumeInterpolator::SurfaceInterpolation(CGeometry* geometry_src, CGeometry* geometry_dst, CSolver* solver_src,
                                                     CSolver* solver_dst, vector<unsigned long> &pointsFailed) {
  if (pointsFailed.empty()) return;

  /*--- Check if surface ADT was built successfully ---*/
  CADTElemClass& surfaceADT = GetSourceSurfaceADT();
  if (!surfaceADT.IsEmpty()) {

    if (rank == MASTER_NODE) cout << " Done." << endl;

    /*--- Search for donor elements for failed points ---*/
    if (rank == MASTER_NODE) {
      cout << "Performing minimum distance search for " << pointsFailed.size()
           << " failed points." << endl << flush;
    }

    /*--- Loop over failed points for minimum distance search ---*/
    unsigned long nExtrapolated = 0;
    for (auto l = 0u; l < pointsFailed.size(); ++l) {
      /*--- Get coordinates of failed point ---*/
      const unsigned long pointID = pointsFailed[l];
      const su2double* coor = geometry_dst->nodes->GetCoord(pointID);

      /*--- Find nearest surface element ---*/
      unsigned short markerID;
      unsigned long elemID;
      int rankID;
      su2double dist;
      su2double surfCoor[3];

      /*--- Find the closest point on the source surface mesh ---*/
      surfaceADT.DetermineNearestElement(coor, dist, markerID, elemID, rankID);
      NearestPointOnElement(geometry_src, markerID, elemID, coor, surfCoor,
                            dist, nDim);

      /*--- Get surface element information ---*/
      unsigned short nNodes = geometry_src->bound[markerID][elemID]->GetnNodes();

      /*--- Use nearest surface element nodes for interpolation ---*/
      su2double weightsInterpol[4];
      if (geometry_src->GetnDim() == 3) {
        /*--- Use surface ADT to get interpolation weights*/
        su2double parCoor[3];
        surfaceADT.DetermineContainingElement(surfCoor, markerID, elemID, rankID, parCoor, weightsInterpol);
      } else {
        /*--- For 2D case (LINE elements), use inverse distance weighting ---*/
        su2double totalWeight = 0.0;

        for (auto iNode = 0u; iNode < nNodes; ++iNode) {
          unsigned long nodeID = geometry_src->bound[markerID][elemID]->GetNode(iNode);

          /*--- Compute distance from interpolation point to node ---*/
          su2double dist2 = 0.0;
          for (auto k = 0u; k < nDim; ++k) {
            su2double diff = coor[k] - geometry_src->nodes->GetCoord(nodeID, k);
            dist2 += diff * diff;
          }

          /*--- Inverse distance weighting (with small epsilon to avoid division by zero) ---*/
          weightsInterpol[iNode] = 1.0 / (sqrt(dist2) + 1e-12);
          totalWeight += weightsInterpol[iNode];
        }

        /*--- Normalize weights ---*/
        for (auto iNode = 0u; iNode < nNodes; ++iNode) {
          weightsInterpol[iNode] /= totalWeight;
        }
      }

      /*--- Initialize interpolated solution to zero ---*/
      for (auto iVar = 0u; iVar < nVar; ++iVar) {
        solver_dst->GetNodes()->SetSolution(pointID, iVar, 0.0);
      }

      /*--- Interpolate using shape function weights ---*/
      for (auto iNode = 0u; iNode < nNodes; ++iNode) {
        unsigned long nodeID = geometry_src->bound[markerID][elemID]->GetNode(iNode);

        for (auto iVar = 0u; iVar < nVar; ++iVar) {
          su2double val = solver_src->GetNodes()->GetSolution(nodeID, iVar);
          solver_dst->GetNodes()->Add_DeltaSolution(pointID, iVar, weightsInterpol[iNode] * val);
        }
      }

      nExtrapolated++;
    }

    if (rank == MASTER_NODE) {
      cout << "Surface search finished. " << nExtrapolated << " points extrapolated."
           << endl << flush;
    }

  } else {
    /*--- No surface elements found ---*/
    if (rank == MASTER_NODE) {
      cout << " No surface elements found for minimum distance search." << endl;
    }
  }
}