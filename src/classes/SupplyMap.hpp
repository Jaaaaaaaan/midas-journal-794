/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: SupplyMap.h,v $
    Language: C++
    Date: $Date: 2011/02/08 10:43:00 $
    Version: $Revision: 1.0 $

   Copyright (c) 2011 Medical Imaging Analysis Lab, Simon Fraser University,
   British Columbia, Canada.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

* The name of the Insight Consortium, nor the names of any consortium members,
   nor of any contributors, may be used to endorse or promote products derived
   from this software without specific prior written permission.

* Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

* Free for non-commercial use only.  For commercial use, explicit approval
   must be requested by contacting the Authors.

* If you use the code in your work, you must acknowledge it

* Modifications of the source code must also be released as open source

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   =========================================================================*/

/*=========================================================================
   This source code was modified as part of the thesis
   "Generating Synthetic Vasculature in Organ-Like 3D Meshes" at NCT TSO Dresden.
   =========================================================================*/

#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <string>
#include <vector>

/**
 * \brief Describes how the oxygen demand changes when a new
 * terminal node is added to a vascular tree.
 *
 * For each new terminal node, the oxygenation map calls
 * \see SupplyMap::reduction to update the oxygen demand in the vicinity
 * of the terminal node.
 *
 * The supply map is described by a supply vector
 * (w0, w1, ..., wn) and a treshold distance d.
 * Consider a newly added terminal node at position p_t.
 * Then the reduction in oxygen demand for a voxel at position p_v is
 * computed as
 * 0                                            (if ||p_t - p_v|| > d)
 * 1 / (sum_{i = 0}^n w_i * ||p_t - p_v||^i)    (otherwise).
 */
class SupplyMap
{
private:
    std::vector<double> m_supplyVector;
    double m_tresholdDistance;

public:

    /**
     * \brief Loads the supply map parameters from a file.
     *
     * The layout of the supply map file is assumed to be
     * <line 1>:
     * <line 2>:
     * <line 3>: w1 w2 ... wn d
     *
     * The contents of the first two lines are ignored for compatibility.
     */
    void loadMap(const std::string &filepath);

    /**
     * \brief Calculates the reduction in oxygen demand at
     * "target" if "supplied" becomes a terminal node.
     *
     * \param supplied The position of a new terminal node.
     * \param target The position at which to compute the reduction.
     *
     * \return The remaining oxygen demand at "target"
     * as a perecentage of its current oxygen demand.
     */
    double reduction(Eigen::Vector3i supplied, Eigen::Vector3i target) const;

    /**
     * \brief Getter function for the treshold distance.
     *
     * Outside of the treshold distance, the oxygen demand around
     * a new terminal node remains unchanged.
     */
    int getTresholdDistance() const
    {
        return std::ceil(m_tresholdDistance);
    }
};
