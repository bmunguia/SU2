#!/usr/bin/env python

## \file metrics.py
#  \brief Helper utilities for mesh adaptation metric computation with generalized sensors
#  \author B. Munguía
#  \version 8.2.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

from typing import Optional


def resolve_sensors(
    driver, sensor_list: list[str], verbose: bool = True
) -> list[tuple[int, int, str]]:
    """
    Resolve sensor name strings to (solver_idx, var_idx, name) tuples.

    This function queries the driver for available primitive variables in all
    solvers, then matches the requested sensor names to their indices.

    Args:
        driver: SU2 driver object with GetSolverVariables() method
        sensor_list: list of primitive variable names (e.g., "DENSITY", "PRESSURE")
        verbose: Print resolution information

    Returns:
        list of (solver_idx, var_idx, name) tuples for resolved sensors

    Raises:
        ValueError: If a sensor cannot be resolved

    Note:
        Only primitive variables are currently supported.
    """
    # Get available variables from all solvers
    solver_vars = driver.GetSolverVariables()

    if verbose:
        print("\n--- Python Sensor Resolution ---")
        print("Available solver variables:")
        for solver_name, var_names in solver_vars.items():
            print(
                f"  {solver_name}: {', '.join(var_names) if var_names else '(empty)'}"
            )

    resolved = []

    # Resolve each sensor by searching all solvers
    for sensor in sensor_list:
        found = False
        solver_idx = 0

        for solver_name, var_list in solver_vars.items():
            if sensor in var_list:
                var_idx = var_list.index(sensor)
                resolved.append((solver_idx, var_idx, sensor))
                found = True

                if verbose:
                    print(
                        f"Resolved '{sensor}' -> Solver {solver_idx} ({solver_name}), Variable {var_idx}"
                    )
                break
            solver_idx += 1

        if not found:
            raise ValueError(
                f"Could not resolve sensor '{sensor}' in any solver.\n"
                f"Available solvers: {list(solver_vars.keys())}"
            )

    return resolved


def setup_metric_computation(
    driver, verbose: bool = True
) -> list[tuple[int, int, str]]:
    """
    Set up metric computation by resolving sensors and configuring the driver.

    This is the main entry point for setting up mesh adaptation metrics.
    It reads the METRIC_SENSOR config option, resolves sensor names,
    and passes the resolved indices back to the C++ driver.

    Args:
        driver: SU2 driver object
        verbose: Print detailed resolution information

    Returns:
        list of (solver_idx, var_idx, name) tuples for resolved sensors

    Note:
        Only primitive variables are supported in METRIC_SENSOR.
        Use variables like DENSITY, PRESSURE, VELOCITY_X, etc.
    """
    # Get sensor list from config
    sensor_list = driver.GetMetricSensorList()

    if not sensor_list:
        if verbose:
            print("No metric sensors specified in config (METRIC_SENSOR is empty)")
        return []

    if verbose:
        print(f"Requested sensors: {', '.join(sensor_list)}")

    # Resolve sensors to (solver_idx, var_idx, name) tuples
    resolved = resolve_sensors(driver, sensor_list, verbose=verbose)

    # Pass resolved sensors back to C++ driver
    driver.SetResolvedSensors(resolved)

    if verbose:
        print(f"Successfully resolved {len(resolved)} sensor(s)")
        print()

    return resolved


def print_sensor_info(driver, resolved_sensors: Optional[list] = None):
    """
    Print information about available and resolved sensors.

    Args:
        driver: SU2 driver object
        resolved_sensors: Optional list of resolved sensors from setup_metric_computation()
    """
    print("\n=== Metric Sensor Information ===")

    # Available variables
    solver_vars = driver.GetSolverVariables()
    print("\nAvailable primitive variables by solver:")
    for solver_name, var_names in solver_vars.items():
        print(f"  {solver_name}:")
        if var_names:
            for var in var_names:
                print(f"    - {var}")
        else:
            print("    (none)")

    # Resolved sensors
    if resolved_sensors:
        print(f"\nResolved sensors ({len(resolved_sensors)}):")
        for solver_idx, var_idx, name in resolved_sensors:
            print(f"  {name} -> Solver {solver_idx}, Index {var_idx}")

    print("=" * 35 + "\n")
