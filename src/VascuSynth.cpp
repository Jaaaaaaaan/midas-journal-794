/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: VascuSynth.cpp,v $
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

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <random>
#include <vector>

#include "OBB.hpp"
#include "OxygenationMap.hpp"
#include "SupplyMap.hpp"
#include "UtilsOxygenationMap.hpp"
#include "VascularTree.hpp"

/**
 * \brief Utility function to read a text file and store the lines into a vector
 * that is returned.
 *
 * \param filepath Path to the file to read.
 *
 * \return A vector of lines read from the file.
 *
 * \throws Runtime error if the file cannot be read
 */
std::vector<std::string> readFileLines(const std::string &filepath)
{
    std::ifstream oFile;

    oFile.open(filepath, std::ios::in);
    std::vector<std::string> lines;
    std::string line;

    if (!oFile.is_open())
    {
        throw std::runtime_error("Could not open file " + (filepath));
    }

    while (!oFile.eof())
    {
        std::getline(oFile, line);
        lines.push_back(line);
    }

    oFile.close();
    return lines;
}

int parseRandomSeed(const std::string &seedString)
{
    int randomSeed = std::stoi(seedString);

    if (randomSeed < 0)
    {
        std::random_device device;
        randomSeed = device();

        std::cerr << "Warning: The program parameters contain a negative "
                        "random seed. A pseudo-randomly selected "
                        "seed will be used instead." << std::endl;
    }

    return randomSeed;
}

/**
 *  \brief Reads parameters for a vascular tree from a parameter file and
 * generates a corresponding \see VascularTree instance.
 *
 * Parameter File Entries:
 *
 * PERF_POINT: Perforation point of the tree.
 * PERF_PRESSURE: Perfussion pressure.
 * TERM_PRESSURE: Terminal pressure.
 * PERF_FLOW: Perfusion flow.
 * RHO: Rho.
 * GAMMA: Gamma.
 * LAMBDA: Lambda.
 * MU: Mu.
 * MIN_DISTANCE: Minimum distance between two vessel branches.
 * NUM_NODES: Number of terminal nodes.
 * VOXEL_WIDTH: Voxel width of the oxygen demand map in SI units.
 * CLOSEST_NEIGHBOURS: Number of neighbours to test during optimization.
 * CHECK_SELF_INTERSECTION: If set to zero, self-intersection tests are
 * omitted.
 * OUTPUT_FILENAME: Filename for the XML representation of the vascular tree.
 *
 * For more information on these parameters, see:
 * Preet Jassi and Ghassan Hamarneh
 * VascuSynth: Vascular Tree Synthesis Software.
 * Insight Journal (2011).
 *
 * \param filepath Path to the parameter file.
 * \param randomSeed Variable in which to store the random seed for this tree.
 *
 * \return A pointer to the initialized vascular tree.
 */
std::shared_ptr<VascularTree> initializeTree(const std::string &filepath)
{
    const std::vector<std::string> lines = readFileLines(filepath);
    std::map<std::string, std::string> args;

    for (const auto &line : lines)
    {
        if (line.compare("") == 0)
        {
            continue;
        }

        const int colonPosition = line.find(":");
        const std::string name = line.substr(0, colonPosition);
        const std::string value = line.substr(colonPosition + 2);
        args[name] = value;
    }

    auto checkDefined = [args](const std::string &name)
    {
        if (args.find(name) != args.end())
        {
            return;
        }

        throw std::invalid_argument("Not all parameters have been defined");
    };

    auto getInt = [args, checkDefined](const std::string &name)
    {
        checkDefined(name);
        return std::stoi(args.at(name));
    };

    auto getDouble = [args, checkDefined](const std::string &name)
    {
        checkDefined(name);
        return std::stof(args.at(name));
    };

    // Parse perforation point
    checkDefined("PERF_POINT");
    std::istringstream perfStream(args.at("PERF_POINT"));
    std::vector<double> perfVec;
    double val;

    while (perfStream >> val)
    {
        perfVec.push_back(val);
    }

    if (perfVec.size() != 3)
    {
        throw std::invalid_argument("Received invalid perforation point");
    }

    const Eigen::Vector3d perforationPoint(perfVec[0], perfVec[1], perfVec[2]);

    // Parse remaining arguments
    const double Pperf = getDouble("PERF_PRESSURE");
    const double Pterm = getDouble("TERM_PRESSURE");
    const double Qperf = getDouble("PERF_FLOW");
    const double rho = getDouble("RHO");
    const double gamma = getDouble("GAMMA");
    const double lambda = getDouble("LAMBDA");
    const double mu = getDouble("MU");
    const double minDistance = getDouble("MIN_DISTANCE");
    const double voxelWidth = getDouble("VOXEL_WIDTH");
    const int numTerminalNodes = getInt("NUM_NODES");
    const int closestNeighbours = getInt("CLOSEST_NEIGHBOURS");
    bool checkSelfIntersection = false;
    std::string outputFilename = "";

    // Parse optional arguments
    if (args.find("CHECK_SELF_INTERSECTION") != args.end())
    {
        checkSelfIntersection = getInt("CHECK_SELF_INTERSECTION") != 0;
    }

    if (args.find("OUTPUT_FILENAME") != args.end())
    {
        outputFilename = args.at("OUTPUT_FILENAME");
    }

    return std::make_shared<VascularTree>(perforationPoint, outputFilename,
                                          Pperf, Pterm, Qperf,
                                          rho, gamma, lambda, mu,
                                          minDistance, voxelWidth,
                                          numTerminalNodes, closestNeighbours,
                                          checkSelfIntersection);
}

/**
 * \brief Iteratively builds non-intersecting vascular trees ´
 * with the desired parameters.
 *
 * \param vts A vector of pointers to the initialized
 * vascular tree instances.
 */
void buildTrees(std::vector<std::shared_ptr<VascularTree>> vts)
{
    for (const auto &vt : vts)
    {
        vt->setVascularTrees(vts);
    }
    std::deque<std::shared_ptr<VascularTree>> queue(vts.begin(), vts.end());

    // Schedule iteration steps for the vascular trees
    // using round-robin scheduling
    while (queue.size() > 0)
    {
        std::shared_ptr<VascularTree> vt = std::move(queue.front());
        queue.pop_front();
        vt->buildTreeStep();

        if (vt->getNumberOfNodes() < vt->getTargetNumberOfNodes())
        {
            queue.push_back(std::move(vt));
        }
    }

    for (const auto &vt : vts)
    {
        vt->calculateRadius();
    }

    // Check whether interesections have been introduced due to increases
    // of radius during tree growth
    for (const auto &vt : vts)
    {
        for (int i = 1; i < vt->getNumberOfNodes(); ++i)
        {
            const Node *node = vt->getNode(i);
            const Eigen::Vector3d target = node->pos;
            const Eigen::Vector3d source = node->parent->pos;
            const double radius = node->radius / vt->getMapVoxelWidth();

            if (vt->intersectsTrees(target, source, radius, 1.0))
            {
                throw std::runtime_error("Intersection detected");
            }
        }

        if (vt->intersectsSelfTree())
        {
            throw std::runtime_error("Self-intersection detected");
        }
    }
}


/**
 * \brief Assigns a direction of curvature and a curve height
 * to each branch.
 * Based on this information, it is possbile to generate
 * non-intersecting, curved vessel branches.
 *
 * \param vts The vascular trees for which to compute curvature parameters.
 * \param dim The extents of the volume that is being vascularized
 * in units of voxels.
 * \param randomSeed A random seed for the random generator used during the
 * computation. If randomSeed is negative, a random seed is selected
 * pseudo-randomly instead.
 */
std::deque<std::shared_ptr<CurvedBranchBoundingBox>>
generateCurveBoundingBoxes(
    std::vector<std::shared_ptr<VascularTree>> vts,
    Eigen::Vector3i dim, int randomSeed)
{
    std::shared_ptr<CurvedBranchBoundingBox> current;
    std::deque<std::shared_ptr<CurvedBranchBoundingBox>> boxes;
    std::deque<std::shared_ptr<CurvedBranchBoundingBox>> queue;

    int iter = 1;
    const int maxIter = 1000;

    const double diameter = pow(dim[0] * dim[0] + dim[1] * dim[1] + dim[2] * dim[2], 0.5);
    const double startWidth = 7.0 * diameter / vts[0]->getTargetNumberOfNodes();
    const double increment = 0.25;

    for (const auto &vt : vts)
    {
        for (int id = 1; id < vt->getNumberOfNodes(); ++id)
        {
            const Node *node = vt->getNode(id);
            const Eigen::Vector3d target = node->pos;
            const Eigen::Vector3d source = node->parent->pos;
            const double radius = node->radius / vt->getMapVoxelWidth();

            current = std::make_shared<CurvedBranchBoundingBox>(target, source, radius,
                                                                startWidth, id);
            boxes.push_back(std::move(current));
        }
    }

    int offsetLeft = 0;
    int offsetRight = 0;
    int index = 0;

    for (const auto &vt : vts)
    {
        offsetRight += vt->getNumberOfNodes() - 1;
        for (int id = 1; id < vt->getNumberOfNodes(); ++id)
        {
            current = boxes[index];
            index++;

            current->initializeIntersections({ boxes.begin(), boxes.begin() + offsetLeft });
            current->initializeIntersections({ boxes.begin() + offsetRight, boxes.end() });
            current->initializeSelfIntersections({ boxes.begin() + offsetLeft, boxes.begin() + offsetRight });

            if (current->hasIntersections())
            {
                queue.push_back(std::move(current));
            }
        }

        offsetLeft = offsetRight;
    }

    std::mt19937 rand(randomSeed);
    shuffle(queue.begin(), queue.end(), rand);

    while (queue.size() > 0)
    {
        iter++;
        current = std::move(queue.front());
        queue.pop_front();
        current->contract(increment * current->getHeight());
        current->updateIntersections();

        if (current->hasIntersections())
        {
            queue.push_back(std::move(current));
        }

        if (iter > maxIter)
        {
            for (auto box : queue)
            {
                box->contract(startWidth);
            }
            break;
        }
    }

    for (auto box : boxes)
    {
        box->updateSelfIntersections();
        if (box->hasSelfIntersections())
        {
            box->contract(startWidth);
        }
    }

    return boxes;
}

/**
 * \brief Prints a node into XML/GXL format
 *
 * \param node The node to print.
 * \param os The stream to which the node is printed.
 */
void printNode(const Node *node, std::ofstream &os)
{
    const NodeType type = node->type;
    std::string typeStr;

    switch (type)
    {
    case NodeType::ROOT:
        typeStr = "root node";
        break;

    case NodeType::TERM:
        typeStr = "terminal node";
        break;

    case NodeType::BIF:
        typeStr = "bifurcation";
        break;

    default:
        typeStr = "unknown type";
    }

    os << "  <node id=\"n" << node->index << "\">" << std::endl;
    os << "    <attr name=\" nodeType\">" << std::endl;
    os << "      <string> " << typeStr << " </string>" << std::endl;
    os << "    </attr>" << std::endl;

    os << "    <attr name=\" position\">" << std::endl;
    os << "      <tup>" << std::endl;
    os << "        <float>" << node->pos[0] << "</float>" << std::endl;
    os << "        <float>" << node->pos[1] << "</float>" << std::endl;
    os << "        <float>" << node->pos[2] << "</float>" << std::endl;
    os << "      </tup>" << std::endl;
    os << "    </attr>" << std::endl;
    os << "  </node>" << std::endl;

    if (type != NodeType::TERM)
    {
        printNode(node->leftChild, os);
        if (type != NodeType::ROOT)
        {
            printNode(node->rightChild, os);
        }
    }
}

/**
 * \brief Prints an edge into XML/GXL format
 *
 * \param node The target node of the edge to print.
 * \param os The stream to which the edge is printed.
 */
void printEdge(const Node *node, std::ofstream &os,
               const std::deque<std::shared_ptr<CurvedBranchBoundingBox>> &curves)
{
    const NodeType type = node->type;
    const int idx = node->index;

    if (type != NodeType::ROOT)
    {
        os << "  <edge id=\"e" << idx << "\" to=\"n" << idx;
        os << "\" from=\"n" << node->parent->index << "\">" << std::endl;
        os << "    <attr name=\" flow\">" << std::endl;
        os << "      <float>" << node->flow << "</float>" << std::endl;
        os << "    </attr>" << std::endl;

        os << "    <attr name=\" radius\">" << std::endl;
        os << "      <float>" << node->radius << "</float>" << std::endl;
        os << "    </attr>" << std::endl;

        if (!curves.empty())
        {
            double curveHeight = curves[idx - 1]->getCurveHeight();
            os << "    <attr name=\" curveHeight\">" << std::endl;
            os << "      <float>" << curveHeight << "</float>" << std::endl;
            os << "    </attr>" << std::endl;

            Eigen::Vector3d dir(curves[idx - 1]->getCurveDirection());

            os << "    <attr name=\" curveDirection\">" << std::endl;
            os << "      <tup>" << std::endl;
            os << "        <float>" << dir[0] << "</float>" << std::endl;
            os << "        <float>" << dir[1] << "</float>" << std::endl;
            os << "        <float>" << dir[2] << "</float>" << std::endl;
            os << "      </tup>" << std::endl;
            os << "    </attr>" << std::endl;
        }

        os << "  </edge>" << std::endl;
    }

    if (type != NodeType::TERM)
    {
        printEdge(node->leftChild, os, curves);
        if (type != NodeType::ROOT)
        {
            printEdge(node->rightChild, os, curves);
        }
    }
}

/**
 * \brief Writes a tree structure to a GXL file.
 *
 * \param vt A pointer to the vascular tree that is being printed.
 * \param filepath Path to the output file.
 * \param curves Contains a curve height and direction of curvature for
 * each branch, which can be used to generate curved branches.
 */
void printTreeStructure(std::shared_ptr<VascularTree> vt,
                        const std::string &filepath,
                        const std::deque<std::shared_ptr<CurvedBranchBoundingBox>> &curves)
{
    std::ofstream output;

    output.open(filepath);
    output << "<gxl><graph id=\"" << filepath
           << "\" edgeids=\" true\" edgemode=\" directed\" hypergraph=\" false\">"
           << std::endl;

    const Node *root = vt->getNode(0);
    printNode(root, output);
    printEdge(root, output, curves);
    output << "</graph></gxl>" << std::endl;
    output.close();
}

/**
 * \brief VascuSynth: takes a series of parameter files,
 * and generates non-intersecting vascular structures
 * based on the parameters.
 * Information about the vascular structure is saved as
 * a GXL file that can be visualized using software such as GraphViz.
 *
 * Arguments: oxygenationMap.txt supplyMap.txt randomSeed paramFile1.txt ... paramFileN.txt
 */
int main(int argc, const char **argv)
{
    std::vector<std::shared_ptr<VascularTree>> vts;

    if (argc < 5)
    {
        // not enough parameters specified
        std::cout << "An error has occured: incorrect number of arguments" << std::endl;
        std::cout << "Usage: VascuSynth [oxygenDemandMap] [supplyMap] [randomSeed] [paramFile 1] ... [paramFile N]" << std::endl;
        return 0;
    }

    const std::string omPath = argv[1];
    const std::string smPath = argv[2];
    int randomSeed = parseRandomSeed(argv[3]);
    Eigen::Vector3i dim;

    std::cout << "Reading parameters and building the trees..." << std::endl;

    std::shared_ptr<SupplyMap> supplyMap = std::make_shared<SupplyMap>();
    supplyMap->loadMap(smPath);
    std::shared_ptr<std::vector<double>> mapArray = Utils::OxygenationMap::loadMap(omPath, dim);

    for (int i = 4; i < argc; ++i)
    {
        std::shared_ptr<VascularTree> vt = initializeTree(argv[i]);
        const int treeSeed = randomSeed + i;

        auto oxygenationMap = std::make_unique<OxygenationMap>(supplyMap, mapArray, dim, treeSeed);
        vt->setOxygenationMap(std::move(oxygenationMap));
        vts.push_back(std::move(vt));
    }

    CurvedBranchBoundingBox::initRandom(randomSeed);

    try
    {
        buildTrees(vts);
    }
    catch (const std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        std::cout << "Exiting VascuSynth" << std::endl;
        return 1;
    }

    std::cout << "The vascular trees have been built sucessfully..." << std::endl;

    std::deque<std::shared_ptr<CurvedBranchBoundingBox>> curves;
    // Stores curvature information for an individual vascular tree
    curves = generateCurveBoundingBoxes(vts, dim, randomSeed + 1);

    // create the GXL files
    #pragma omp parallel for
    for (std::size_t i = 0; i < vts.size(); ++i)
    {
        // Compute start and end indices for curves vector
        int start = 0;
        for (std::size_t j = 0; j < i; ++j)
        {
            start += vts[i]->getNumberOfNodes() - 1;
        }

        int end = start + (vts[i]->getNumberOfNodes() - 1);

        const std::deque<std::shared_ptr<CurvedBranchBoundingBox>> curvesVt =
            { curves.begin() + start, curves.begin() + end };

        std::string filename = vts[i]->getOutputFilename();
        if (filename.empty())
        {
            filename = "tree_structure" + std::to_string(i + 1);
        }

        filename += ".xml";
        printTreeStructure(vts[i], filename.c_str(), curvesVt);
    }

    std::cout << "Information about the vascular structures has been saved..." << std::endl;
    return 0;
}
