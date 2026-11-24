/*!
 * \file CMetricBaselineOutput.cpp
 * \brief Main subroutines for metric field output
 * \author B. Munguía
 * \version 8.4.0 "Harrier"
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



#include "../../include/output/CMetricBaselineOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CMetricBaselineOutput::CMetricBaselineOutput(CConfig *config, unsigned short nDim, CSolver* solver) : CBaselineOutput(config, nDim, solver) {

  /*--- Set the requested volume fields to all fields in the solver ---*/

  requestedVolumeFields.clear();

  requestedVolumeFields.emplace_back("COORDINATES");
  requestedVolumeFields.emplace_back("MESH_ADAPT");

  nRequestedVolumeFields = requestedVolumeFields.size();

  /*--- Keep the coordinates and the adaptive metric of every sensor (AdaptiveMetric_<Sensor>_*),
   *    and store the solution column of each metric field ---*/
  std::vector<string> kept(fields.begin(), fields.begin() + nDim);
  for (auto iField = static_cast<unsigned long>(nDim); iField < fields.size(); ++iField) {
    if (fields[iField].rfind("AdaptiveMetric_", 0) == 0) {
      kept.push_back(fields[iField]);
      metricColumns.push_back(iField);
    }
  }
  fields = kept;

}

CMetricBaselineOutput::~CMetricBaselineOutput() = default;

void CMetricBaselineOutput::SetVolumeOutputFields(CConfig *config){

  unsigned short iField = 0;

  // Grid coordinates
  AddVolumeOutput(fields[0], fields[0], "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput(fields[1], fields[1], "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput(fields[2], fields[2], "COORDINATES", "z-component of the coordinate vector");

  // Add all the remaining fields

  for (iField = nDim; iField < fields.size(); iField++){
    AddVolumeOutput(fields[iField], fields[iField], "MESH_ADAPT", "");
  }

}

void CMetricBaselineOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  /*--- Take the solver at index 0 --- */

  CVariable* Node_Sol  = solver[0]->GetNodes();
  auto* Node_Geo = geometry->nodes;

  // Grid coordinates
  SetVolumeOutputValue(fields[0], iPoint, Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue(fields[1], iPoint, Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue(fields[2], iPoint, Node_Geo->GetCoord(iPoint, 2));

  for (auto iMetric = 0ul; iMetric < metricColumns.size(); ++iMetric)
    SetVolumeOutputValue(fields[nDim + iMetric], iPoint, Node_Sol->GetSolution(iPoint, metricColumns[iMetric]));

}
