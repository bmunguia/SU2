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

  /*--- Remove fields that aren't the coordinates or metric ---*/
  auto it = fields.begin() + nDim;
  while (it != fields.end()) {
    if (it->rfind("Metric_", 0) != 0) {
      it = fields.erase(it);
    } else {
      ++it;
    }
  }

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

  unsigned short iField = 0;

  /*--- Take the solver at index 0 --- */

  CVariable* Node_Sol  = solver[0]->GetNodes();
  auto* Node_Geo = geometry->nodes;

  // Grid coordinates
  SetVolumeOutputValue(fields[0], iPoint, Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue(fields[1], iPoint, Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue(fields[2], iPoint, Node_Geo->GetCoord(iPoint, 2));

  if(config->GetCompute_Metric()) {
    vector<string> solver_fields = solver[0]->GetSolutionFields();
    solver_fields.erase(solver_fields.begin()); // remove Point_ID
    /*--- Get the first index that has Metric_ in the name ---*/
    int first_ind = -1;
    for (size_t i = 0; i < solver_fields.size(); ++i) {
      string field_name = solver_fields[i].substr(1, solver_fields[i].size() - 2);
      if (field_name.rfind("Metric_", 0) == 0) {
        first_ind = i;
        break;
      }
    }

    if (first_ind >= 0) {
      // Common metric components for both 2D and 3D
      SetVolumeOutputValue("Metric_xx", iPoint, Node_Sol->GetSolution(iPoint, first_ind + 0));
      SetVolumeOutputValue("Metric_xy", iPoint, Node_Sol->GetSolution(iPoint, first_ind + 1));
      SetVolumeOutputValue("Metric_yy", iPoint, Node_Sol->GetSolution(iPoint, first_ind + 2));

      // Additional components for 3D
      if (nDim == 3) {
        SetVolumeOutputValue("Metric_xz", iPoint, Node_Sol->GetSolution(iPoint, first_ind + 3));
        SetVolumeOutputValue("Metric_yz", iPoint, Node_Sol->GetSolution(iPoint, first_ind + 4));
        SetVolumeOutputValue("Metric_zz", iPoint, Node_Sol->GetSolution(iPoint, first_ind + 5));
      }
    }
  }

}
