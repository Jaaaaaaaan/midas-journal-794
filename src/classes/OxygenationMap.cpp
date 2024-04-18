/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: OxygenationMap.cpp,v $
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

#include "OxygenationMap.hpp"
#include "SupplyMap.hpp"
#include "UtilsOxygenationMap.hpp"

// A parametric macro to allow the use of dynamic multi-dimensional arrays
#define arr(arr, x, y, z, dim) arr->at(z + dim[2] * (y + dim[1] * (x)))

OxygenationMap::OxygenationMap(std::shared_ptr<SupplyMap> supplyMap,
                               std::shared_ptr<std::vector<double>> mapArray,
                               Eigen::Vector3i dim, int randomSeed) :
                               m_dim(dim), m_supplyMap(supplyMap)
{
    m_randDist = std::uniform_real_distribution<double>(0.0, 1.0);

    // Use shallow copy for m_map member because it is never written to
    m_map = mapArray;
    m_effectiveMap = Utils::OxygenationMap::copyMap(mapArray, dim);
    m_rand = std::mt19937(randomSeed);
}

double OxygenationMap::sum() const
{
    double acc = 0.0;

    for (int i = 0; i < m_dim[0]; i++)
    {
        for (int j = 0; j < m_dim[1]; j++)
        {
            for (int k = 0; k < m_dim[2]; k++)
            {
                acc += arr(m_effectiveMap, i, j, k, m_dim);
            }
        }
    }

    return acc;
}

Eigen::Vector3i OxygenationMap::candidate()
{
    const double r = m_randDist(m_rand) * sum();
    double acc = 0.0;

    for (int i = 0; i < m_dim[0]; i++)
    {
        for (int j = 0; j < m_dim[1]; j++)
        {
            for (int k = 0; k < m_dim[2]; k++)
            {
                acc += arr(m_effectiveMap, i, j, k, m_dim);
                if (acc >= r)
                {
                    return Eigen::Vector3i(i, j, k);
                }
            }
        }
    }

    throw std::runtime_error("Unreachable");
}

void OxygenationMap::applyCandidate(Eigen::Vector3i candidate)
{
    const int tresholdDist = m_supplyMap->getTresholdDistance();

    auto isValid = [this](Eigen::Vector3i pos)
    {
        return (pos[0] > 0 && pos[0] < m_dim[0]) &&
               (pos[1] > 0 && pos[1] < m_dim[1]) &&
               (pos[2] > 0 && pos[2] < m_dim[2]);
    };

    for (int rx = -tresholdDist; rx <= tresholdDist; ++rx)
    {
        for (int ry = -tresholdDist; ry <= tresholdDist; ++ry)
        {
            for (int rz = -tresholdDist; rz <= tresholdDist; ++rz)
            {
                Eigen::Vector3i pos = candidate + Eigen::Vector3i(rx, ry, rz);

                if (!isValid(pos))
                {
                    continue;
                }

                arr(m_effectiveMap, pos[0], pos[1], pos[2], m_dim) *=
                    m_supplyMap->reduction(candidate, pos);
            }
        }
    }
}

bool OxygenationMap::visible(Eigen::Vector3d source, Eigen::Vector3d target) const
{
    const Eigen::Vector3d vect = target - source;
    Eigen::Vector3d pos = source;

    Eigen::Vector3i voxel(
        (int)(source[0] + 0.5),
        (int)(source[1] + 0.5),
        (int)(source[2] + 0.5));

    const Eigen::Vector3i targetVoxel(
        (int)(target[0] + 0.5),
        (int)(target[1] + 0.5),
        (int)(target[2] + 0.5));

    const double epsilon = 1e-6;

    while (!(voxel == targetVoxel) &&
           !((std::abs(pos[0] - target[0]) < epsilon) &&
             (std::abs(pos[1] - target[1]) < epsilon) &&
             (std::abs(pos[2] - target[2]) < epsilon)))
    {
        double minMult = 1e200;

        for (int i = 0; i < 3; ++i)
        {
            const double dir = vect[i] < 0 ? 0.5 + epsilon : -0.5 - epsilon;
            const double mult = std::abs((voxel[i] - (pos[i] + dir)) / vect[i]);

            if (mult < minMult)
            {
                minMult = mult;
            }
        }

        pos += minMult * vect;
        voxel = Eigen::Vector3i(
            (int)(pos[0] + 0.5),
            (int)(pos[1] + 0.5),
            (int)(pos[2] + 0.5));

        if (arr(m_map, voxel[0], voxel[1], voxel[2], m_dim) == 0)
        {
            return false;
        }
    }

    return true;
}