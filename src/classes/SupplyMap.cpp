/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: SupplyMap.cpp,v $
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

#include <fstream>
#include <sstream>

#include "SupplyMap.hpp"

void SupplyMap::loadMap(const std::string &filepath)
{
    std::ifstream filestream(filepath);
    std::string line;

    // Ignore first two lines of supply map file
    // These are only present for compatability
    std::getline(filestream, line);
    std::getline(filestream, line);
    std::getline(filestream, line);

    double val;
    std::istringstream linestream(line);
    linestream >> val;

    // The first N entries of the supply map file specify the supply vector
    while (linestream.peek() != std::istream::traits_type::eof())
    {
        m_supplyVector.push_back(val);
        linestream >> val;
    }

    // The final entry of the supply map file is the treshold distance
    m_tresholdDistance = val;
}

double SupplyMap::reduction(Eigen::Vector3i supplied, Eigen::Vector3i target) const
{
    if (supplied == target)
    {
        return 0.0;
    }

    const double dist = (supplied - target).cast<double>().norm();

    if (m_tresholdDistance <= dist)
    {
        return 1.0;
    }

    double acc = 0.0;

    for (std::size_t i = 0; i < m_supplyVector.size(); ++i)
    {
        acc += m_supplyVector[i] * pow(dist, i);
    }

    return std::max(1.0 - (1.0 / acc), 0.0);
}
