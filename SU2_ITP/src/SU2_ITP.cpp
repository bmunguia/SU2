/*!
 * \file SU2_ITP.cpp
 * \brief Main file for the solution interpolation code (SU2_ITP).
 *        All interpolation logic is implemented in the <i>interpolation.cpp</i> file.
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

#include "../include/CInterpolatorDriver.hpp"

int main(int argc, char* argv[]) {
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;

  char config_file_name[MAX_STRING_SIZE];

  /*--- MPI initialization ---*/

  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Comm MPICommunicator = SU2_MPI::GetComm();

  const int rank = SU2_MPI::GetRank();
  const int size = SU2_MPI::GetSize();

  if (argc == 2 || argc == 3) {
    strcpy(config_file_name, argv[1]);
  } else {
    strcpy(config_file_name, "default.cfg");
  }

  /*--- Set up a timer for performance benchmarking ---*/

  StartTime = SU2_MPI::Wtime();

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Solution Postprocessing -----------------------" << endl;

  /*--- Initialize the driver: reads config, builds geometry, lazy-inits solvers ---*/

  CInterpolatorDriver driver(config_file_name, MPICommunicator);

  /*--- Run the interpolation ---*/

  if (driver.IsUnsteady()) {
    unsigned long TimeIter = driver.GetStartTimeIter();
    while (!driver.StopCalc(TimeIter)) {
      if (driver.ShouldOutput(TimeIter)) {
        driver.Preprocess(TimeIter);
        driver.Postprocess();
        driver.Output(TimeIter);
      }
      TimeIter++;
    }
  } else {
    driver.Preprocess(0);
    driver.Postprocess();
    driver.Output(0);
  }

  driver.Finalize();

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Finalize Solver -------------------------" << endl;

  /*--- Compute/print the total time for performance benchmarking ---*/

  StopTime = SU2_MPI::Wtime();
  UsedTime = StopTime - StartTime;

  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on " << size;
    if (size == 1)
      cout << " core." << endl;
    else
      cout << " cores." << endl;
  }

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Exit Success (SU2_ITP) ------------------------" << endl << endl;

  SU2_MPI::Finalize();

  return EXIT_SUCCESS;
}
