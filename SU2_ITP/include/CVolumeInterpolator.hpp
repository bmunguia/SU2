/*!
 * \file CVolumeInterpolator.hpp
 * \brief Headers of the main solution interpolation subroutines.
 *        The subroutines and functions are in the <i>CVolumeInterpolator.cpp</i> file.
 * \author B. Munguía, E. van der Weide
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

#include <cmath>
#include <memory>
#include <optional>

#include "../../SU2_CFD/include/solvers/CBaselineSolver.hpp"
#include "../../SU2_CFD/include/solvers/CBaselineSolver_FEM.hpp"
#include "../../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../../SU2_CFD/include/output/CBaselineOutput.hpp"
#include "../../SU2_CFD/include/output/COutputFactory.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/adt/CADTElemClass.hpp"

/*!
 * \class CVolumeInterpolator
 * \brief Base class for performing solution interpolation between meshes.
 * \author B. Munguía, E. van der Weide
 */
class CVolumeInterpolator {
  protected:
    int rank;  /*!< \brief MPI Rank. */
    int size;  /*!< \brief MPI Size. */

    unsigned short nDim = 0;       /*!< \brief Problem dimension. */
    unsigned short nVar = 0;       /*!< \brief Number of solution variables. */
    unsigned long nPoint_src = 0;  /*!< \brief Number of points on the source mesh. */
    unsigned long nPoint_dst = 0;  /*!< \brief Number of points on the destination mesh. */
    unsigned long nElem_src = 0;   /*!< \brief Number of elements on the source mesh. */
    unsigned long nElem_dst = 0;   /*!< \brief Number of elements on the destination mesh. */

    vector<su2double> coorDst;     /*!< \brief Coordinates of destination mesh DOFs. */

    vector<optional<unsigned long>> containingElems;  /*!< \brief Map of node index to source element containing the node. */
    vector<int> containingElemRanks;                  /*!< \brief Map of node index to rank of source element containing the node. */
    vector<unsigned long> uncontainedNodes;             /*!< \brief List of points for which there is no containing source element. */
    vector<optional<unsigned long>> nearestNodes;    /*!< \brief List of the nearest contained node to each uncontained node. */
    vector<unsigned long> nearestElemRanks;           /*!< \brief List of the rank of the nearest destination element to each uncontained node. */

    std::unique_ptr<CFEMStandardElement> stdElement_ptr; /*!< \brief Standard element object used for Gauss quadrature. */

    /*--- Volume ADT data structure for volume interpolation ---*/
    std::unique_ptr<CADTElemClass> srcVolumeADT_ptr; /*!< \brief ADT for source surface mesh. */

    /*--- Surface ADT data structures for curvature correction and surface interpolation ---*/
    std::unique_ptr<CADTElemClass> srcSurfaceADT_ptr; /*!< \brief ADT for source surface mesh. */
    std::unique_ptr<CADTElemClass> dstSurfaceADT_ptr; /*!< \brief ADT for destination surface mesh. */

    /*!
     * \brief Constructor of the class.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CVolumeInterpolator(SU2_Comm MPICommunicator);

  public:
    /*!
     * \brief Virtual destructor.
     */
    virtual ~CVolumeInterpolator() = default;

    void InitializeConfig(CConfig* driver_config, CConfig** config_container, char* zone_file_name,
                          char* config_file_name, SU2_COMPONENT val_software, int iZone, int nZone,
                          SU2_MPI::Comm MPICommunicator, bool isSource);

    void InitializeGeometry(CConfig* config, CGeometry*& geometry, int iZone, int iInst,
                            int nZone, bool isSource);

    void InitializeSolver(CConfig* config, CGeometry* geometry, CSolver**& solver_container, int iZone,
                          int iInst, int nZone);

    void InitializeOutput(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput*& output,
                          int iZone, int iInst, int nZone);

    void LoadRestarts(CConfig* config, CGeometry** geometry_container, CSolver*** solver_container, int iZone,
                      int iInst, int TimeIter, bool UpdateGeo);

  protected:
    /*!
     * \brief Initialize volume and surface ADTs for both source and destination meshes.
     * \param[in] config - Configuration object
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] update - <code>TRUE</code> means to re-initialize the ADTs
     */
    void InitializeADTs(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst, bool update = false);

    std::unique_ptr<CADTElemClass> BuildVolumeADT(CGeometry* geometry);

    std::unique_ptr<CADTElemClass> BuildSurfaceADT(const CConfig* config, CGeometry* geometry);

    std::unique_ptr<CFEMStandardElement> InitializeFEMStandardElement(CConfig* config) {
      unsigned short nPoly = 1;
      bool constJac = true;
      unsigned short orderExact = 0; // Let CFEMStandardElementBase determine the order
      unsigned short VTK_Type = TRIANGLE;
      if (nDim == 3) {
        VTK_Type = TETRAHEDRON;
      }

      return std::unique_ptr<CFEMStandardElement>(
          new CFEMStandardElement(VTK_Type, nPoly, constJac, config, orderExact));
    }

    /*!
     * \brief Get reference to source volume ADT.
     * \return Reference to source volume ADT
     */
    CADTElemClass& GetSourceVolumeADT() {
      if (!srcVolumeADT_ptr) {
        SU2_MPI::Error("Source volume ADT not initialized. Call InitializeADT first.", CURRENT_FUNCTION);
      }
      return *srcVolumeADT_ptr;
    }

    /*!
     * \brief Get reference to source surface ADT.
     * \return Reference to source surface ADT
     */
    CADTElemClass& GetSourceSurfaceADT() {
      if (!srcSurfaceADT_ptr) {
        SU2_MPI::Error("Source surface ADT not initialized. Call InitializeADT first.", CURRENT_FUNCTION);
      }
      return *srcSurfaceADT_ptr;
    }

    /*!
     * \brief Get reference to destination surface ADT.
     * \return Reference to destination surface ADT
     */
    CADTElemClass& GetDestinationSurfaceADT() {
      if (!dstSurfaceADT_ptr) {
        SU2_MPI::Error("Destination surface ADT not initialized. Call InitializeADT first.", CURRENT_FUNCTION);
      }
      return *dstSurfaceADT_ptr;
    }

    /*!
     * \brief Get reference to destination surface ADT.
     * \return Reference to destination surface ADT
     */
    CFEMStandardElement& GetFEMStandardElement() {
      if (!stdElement_ptr) {
        SU2_MPI::Error("FEM standard element not initialized. Call InitializeFEMStandardElement first.", CURRENT_FUNCTION);
      }
      return *stdElement_ptr;
    }

    void InitializeCoords(CGeometry* geometry, vector<su2double>& coords) {
      const unsigned long nPoint = geometry->GetnPoint();
      coords.resize(nPoint * nDim, 0.0);
      for (auto l = 0u; l < nPoint; ++l) {
        for (auto k = 0u; k < nDim; ++k) {
          coords[l * nDim + k] = geometry->nodes->GetCoord(l, k);
        }
      }
    }

  public:
    /*!
     * \brief Main interpolation routine - pure virtual function.
     * \param[in] config - Configuration object
     * \param[in] geometry_src - Source mesh geometry
     * \param[in] geometry_dst - Destination mesh geometry
     * \param[in] solver_container_src - Source mesh solver
     * \param[in] solver_container_dst - Destination mesh solver
     * \param[in] initial_interp - <code>TRUE</code> means this is the first interpolation for the zone
     */
    virtual void Interpolate(CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                             CSolver** solver_container_src, CSolver** solver_container_dst,
                             bool initial_interp) = 0;

  protected:
    virtual void LinearInterpolation(const CConfig *config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                     CSolver* solver_src, CSolver* solver_dst) { }

    virtual void VolumeInterpolation(CGeometry* geometry_src, CSolver* solver_src, CSolver* solver_dst,
                                     const vector<su2double> &coor_corrected, vector<unsigned long> &uncontainedNodes) { }

    void SurfaceInterpolation(CGeometry* geometry_src, CGeometry* geometry_dst, CSolver* solver_src,
                              CSolver* solver_dst);

    void NearestPointOnElement(CGeometry* geometry, unsigned short markerID, unsigned long elemID,
                              const su2double* coor, su2double* surfCoor, su2double& dist2Elem,
                              const unsigned short nDim);

    void NearestPointOnLine(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                            const su2double* coor, su2double* surfCoor, su2double& dist2Line,
                            const unsigned short nDim);

    bool NearestPointOnTriangle(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                                const unsigned long i2, const su2double* coor, su2double* surfCoor,
                                su2double& dist2Tria, const unsigned short nDim);

    bool NearestPointOnQuadrilateral(CGeometry* geometry, const unsigned long i0, const unsigned long i1,
                                    const unsigned long i2, const unsigned long i3, const su2double* coor,
                                    su2double* surfCoor, su2double& dist2Quad, const unsigned short nDim);

    void ApplyCurvatureCorrection(const CConfig* config, CGeometry* geometry_src, CGeometry* geometry_dst,
                                  const unsigned short nDim, vector<su2double> &coor_dst);

  public:
    void WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container,
                    COutput* output, unsigned long TimeIter);
};
