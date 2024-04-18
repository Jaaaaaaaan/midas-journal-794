/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: OxygenationMap.h,v $
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

#include <Eigen/Dense>
#include <memory>
#include <random>
#include <string>

#include "SupplyMap.hpp"

/**
 * \brief Describes the oxygen demand in the volume to be vascularized.
 *
 * The oxygen demand map (ODM) is a 3D array that stores an
 * oxygen demand between 0.0 and 1.0 for each voxel.
 * The oxygen demand must be updated as nodes are being
 * added to the vascular tree.
 * The updated map is stored in the m_effectiveMap array,
 * whereas the original map is stored in the m_map array.
 *
 * No vessel branches will be generated in areas with an oxygen demand of 0.
 */
class OxygenationMap
{
private:
    std::mt19937 m_rand;
    std::uniform_real_distribution<double> m_randDist;

    Eigen::Vector3i m_dim;

    std::shared_ptr<SupplyMap> m_supplyMap;
    std::shared_ptr<std::vector<double>> m_map;
    std::shared_ptr<std::vector<double>> m_effectiveMap;

public:

    /**
     * \brief Constructs a new map instance.
     *
     * \param supplyMap The supply map used to compute reductions to the
     * oxygen demand.
     * \param mapArray A pointer to the array describing the oxygen demand.
     * \param dim A vector describing the extents of the map array.
     * \param randomSeed The random seed based on which new candidate
     * voxels are selected. If randomSeed is negative, a random seed
     * is selected pseudo-randomly instead.
     */
    OxygenationMap(std::shared_ptr<SupplyMap> supplyMap,
                   std::shared_ptr<std::vector<double>> mapArray,
                   Eigen::Vector3i dim, int randomSeed);

    /**
     * \brief Getter function for the extents of the map array.
     */
    Eigen::Vector3i getDim() const
    {
        return m_dim;
    }

    /**
     * \brief Computes the sum of the m_map array.
     *
     * \return The computed sum.
     */
    double sum() const;

    /**
     * \brief Randomly selects a voxel from the map array as
     * a new candidate terminal node.
     *
     * The probability of any given voxel being selected
     * is proportional to its oxygen demand.
     *
     * \return The selected voxel position.
     */
    Eigen::Vector3i candidate();

    /**
     * Updates the m_effectiveMap array upon selection of
     * "candidate" as a new terminal node position.
     *
     * \param candidate The position of a new terminal node.
     */
    void applyCandidate(Eigen::Vector3i candidate);

    /**
     * Determines whether the given positions can be connected by
     * a vessel branch without crossing regions of
     * zero oxygen demand.
     * This computation is based on the m_map array,
     * not on the m_effectiveMap array.
     *
     * \param source The source position.
     * \param target The target position.
     *
     * \return True if a connection is possible and false otherwise.
     */
    bool visible(Eigen::Vector3d source, Eigen::Vector3d target) const;
};