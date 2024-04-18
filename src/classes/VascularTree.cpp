/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: VascularTree.cpp,v $
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

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <queue>

#include "VascularTree.hpp"
#include "OxygenationMap.hpp"
#include "NodeTable.hpp"
#include "OBB.hpp"

VascularTree::VascularTree(
    Eigen::Vector3d perforationPoint, std::string outputFilename,
    double Pperf, double Pterm, double Qperf,
    double rho, double gamma, double lambda, double mu,
    double minDistance, double mapVoxelWidth, int numTerminalNodes,
    int closestNeighbours, bool checkSelfIntersection) :
    m_nodeTable(numTerminalNodes), m_outputFilename(outputFilename),
    m_Pperf(Pperf), m_Pterm(Pterm), m_Qperf(Qperf), m_Qterm(Qperf / numTerminalNodes),
    m_rho(rho), m_gamma(gamma), m_lambda(lambda), m_mu(mu),
    m_minDistance(minDistance), m_mapVoxelWidth(mapVoxelWidth),
    m_numTerminalNodes(numTerminalNodes), m_closestNeighbours(closestNeighbours),
    m_checkSelfIntersection(checkSelfIntersection)
{
    m_nodeTable.addNode(perforationPoint, 0.0, 1.0, 1.0, NULL, NodeType::ROOT);
}

void VascularTree::setVascularTrees(std::vector<std::shared_ptr<VascularTree>> vascularTrees)
{
    for (const auto &vt : vascularTrees)
    {
        m_vascularTrees.push_back(vt);
    }
}

double VascularTree::getRadiusHeuristic(Eigen::Vector3d target, Eigen::Vector3d source,
                                        double flow) const
{
    const double distance = (target - source).norm() * m_mapVoxelWidth;
    const double reducedResistance = 8.0 * m_rho * distance / M_PI;

    double radiusHeuristic = flow * reducedResistance / (m_Pperf - m_Pterm);
    radiusHeuristic = pow(radiusHeuristic, 0.25);
    return radiusHeuristic;
}

double VascularTree::distanceBetweenTwoNodes(const Node *source, const Node *target) const
{
    return (source->pos - target->pos).norm() * m_mapVoxelWidth;
}

void VascularTree::calculateReducedResistance(Node *node)
{
    double acc = 0.0;

    if (node->type != NodeType::TERM)
    {
        acc += pow(node->leftRatio, 4) / (node->leftChild->reducedResistance);
        acc += pow(node->rightRatio, 4) / (node->rightChild->reducedResistance);
        acc = 1.0 / acc;
    }

    acc += (8.0 * m_rho * distanceBetweenTwoNodes(node, node->parent) / M_PI);
    node->reducedResistance = acc;
}

void VascularTree::calculateRatios(Node *node)
{
    Node *left = node->leftChild;
    Node *right = node->rightChild;

    double leftOverRight = (left->flow * left->reducedResistance) /
                           ((right->flow * right->reducedResistance));

    leftOverRight = pow(leftOverRight, 0.25);

    node->leftRatio = pow((1 + pow(leftOverRight, -m_gamma)), (-1.0) / m_gamma);
    node->rightRatio = pow((1 + pow(leftOverRight, m_gamma)), (-1.0) / m_gamma);
}

void VascularTree::updateAtBifurcation(Node *node)
{
    calculateReducedResistance(node);
    if (node->parent->type != NodeType::ROOT)
    {
        calculateRatios(node->parent);
        updateAtBifurcation(node->parent);
    }
}

void VascularTree::calculateRadius(NodeTable &nt)
{
    Node *rootChild = nt.getNode(0)->leftChild;

    if (rootChild == NULL)
    {
        return;
    }

    double radius = (rootChild->flow * rootChild->reducedResistance) / (m_Pperf - m_Pterm);
    radius = pow(radius, 0.25);
    rootChild->radius = radius;
    calculateRadius(rootChild);
}

void VascularTree::calculateRadius(Node *node)
{
    if (node->type == NodeType::TERM)
    {
        return;
    }

    Node *left = node->leftChild;
    Node *right = node->rightChild;
    left->radius = node->radius * node->leftRatio;
    right->radius = node->radius * node->rightRatio;

    calculateRadius(left);
    calculateRadius(right);
}

double VascularTree::calculateFitness(const NodeTable &nt) const
{
    double acc = 0;

    for (int i = 1; i < nt.numNodes(); ++i)
    {
        const Node *node = nt.getNode(i);
        acc += pow(distanceBetweenTwoNodes(node, node->parent), m_mu) *
               pow(node->radius, m_lambda);
    }

    return acc;
}

bool VascularTree::isCandidateValid(Eigen::Vector3d candidate,
                                     const Node *ignored) const
{
    for (int i = 1; i < m_nodeTable.numNodes(); ++i)
    {
        const Node *node = m_nodeTable.getNode(i);

        if (pointSegmentDistance(candidate, node) >= m_minDistance)
        {
            continue;
        }

        if (node != ignored)
        {
            return false;
        }
    }

    return true;
}

void VascularTree::connectPoint(Eigen::Vector3d point,
                                Eigen::Vector3d bifurcationPos,
                                Node *node, NodeTable &nt)
{
    if (node->type == NodeType::ROOT)
    {
        Node *child = nt.addNode(point, m_Qterm, 1.0, 1.0, node, NodeType::TERM);

        // The above call to nt.addNode sets "child" as the left child of "node".
        // The call to node->addChild then also sets the right child
        // of the root node to "child".
        // This is done because the root node is supposed to have
        // only one child node.
        node->addChild(child);
        calculateReducedResistance(child);
        return;
    }

    /**
     *      parent * ----------- * node
     *
     *          becomes
     *
     *	    parent * --------- * bifurcationNode ---------- * node
     *                         \
     *                          \
     *                           \
     *                            \
     *                             \
     *                              * terminal Node
     */
    Node *parent = node->parent;
    Node *left = parent->leftChild;
    Node *right = parent->rightChild;
    parent->unlinkChildren();

    Node *bifurcationNode = nt.addNode(bifurcationPos, node->flow + m_Qterm, 1.0, 1.0, 
                                       parent, NodeType::BIF);
    bifurcationNode->addChild(node);

    Node *terminalNode = nt.addNode(point, m_Qterm, 1.0, 1.0, 
                                    bifurcationNode, NodeType::TERM);

    if (left == node)
    {
        parent->leftChild = bifurcationNode;
        parent->rightChild = right;
    }
    else
    {
        parent->rightChild = bifurcationNode;
        parent->leftChild = left;
    }

    node->parent = bifurcationNode;

    incrementFlow(parent);
    calculateReducedResistance(node);
    updateAtBifurcation(terminalNode);
}

void VascularTree::incrementFlow(Node *node)
{
    node->flow += m_Qterm;
    if (node->parent != NULL)
    {
        incrementFlow(node->parent);
    }
}

double VascularTree::pointSegmentDistance(Eigen::Vector3d point,
                                          const Node *endNode) const
{
    const Eigen::Vector3d pos1 = endNode->pos;
    const Eigen::Vector3d pos2 = endNode->parent->pos;

    double t = (pos2 - pos1).dot(point - pos1) / (pos2 - pos1).dot(pos2 - pos1);

    t = std::clamp(t, 0.0, 1.0);
    return (pos1 + t * (pos2 - pos1) - point).norm();
}

Eigen::Vector3d VascularTree::localOptimization(Eigen::Vector3d point,
                                                const Node *node, double &bestFitness)
{
    const Eigen::Vector3d nodePos = node->pos;
    const Eigen::Vector3d parentPos = node->parent->pos;
    // Initial value for the position that is being optimized
    Eigen::Vector3d bifurcationPos = (nodePos + parentPos) / 2.0;

    bestFitness = 1e200;
    const int numSteps = 20;
    const double stepSize = ((parentPos + nodePos + point) / 3.0 - bifurcationPos).
                            dot(Eigen::Vector3d(1.0, 1.0, 1.0)) * 2.0 / (double)numSteps;

    // Make sure point is visible from bifurcationPos - otherwise return null
    if (!m_oxygenationMap->visible(bifurcationPos, point))
    {
        bestFitness = -1.0;
        return Eigen::Vector3d();
    }

    for (int i = 0; i < numSteps; ++i)
    {
        Eigen::Vector3d localBest = bifurcationPos;

        /**
         * Test all 6 neighbouring positions (+x, -x, +y, -y, +z, -z)
         * of the current bifurcationPos,
         * choose the optimal neighbouring position as the new bifurcationPos.
         *
         * A position can become a new bifurcationPos only if:
         *  - parentPos is visbible from bifurcationPos
         *  - nodePos is visible from bifurcationPos
         *  - point is visible from bifurcationPos
         *  - bifurcationPos satifies the minimum distance criterion
         *    (see VascularTree::isCandidateValid)
         *
         * If all neighbouring positions yield a higher cost than
         * the current bifurcationPos, a local minimum has been found.
         */
        #pragma omp parallel for collapse(2) shared(bestFitness) shared(localBest)
        for (int index = 0; index < 3; ++index)
        {
            for (int sgn = -1; sgn <= 1; sgn += 2)
            {
                double fitness;
                Eigen::Vector3d bifurcationNeighbour = bifurcationPos;
                bifurcationNeighbour[index] += sgn * stepSize;

                if (!(inVolume(bifurcationNeighbour) &&
                      m_oxygenationMap->visible(parentPos, bifurcationNeighbour) &&
                      m_oxygenationMap->visible(nodePos, bifurcationNeighbour) &&
                      m_oxygenationMap->visible(point, bifurcationNeighbour) &&
                      isCandidateValid(bifurcationNeighbour, node)))
                {
                    continue;
                }

                // Copy the node table, connect the candidate point
                // in the copied node table, and evaluate the cost function.
                // Since a throwaway copy of the node table is used,
                // no changes to the actual node table are being made.
                const std::unique_ptr<NodeTable> ntCopy = m_nodeTable.getCopy();
                Node *nodeCopy = ntCopy->getNode(node->index);

                connectPoint(point, bifurcationNeighbour, nodeCopy, *ntCopy);
                calculateRadius(*ntCopy);
                fitness = calculateFitness(*ntCopy);

                #pragma omp critical
                if (fitness < bestFitness)
                {
                    localBest = bifurcationNeighbour;
                    bestFitness = fitness;
                }
            }
        }

        if (localBest != bifurcationPos)
        {
            bifurcationPos = localBest;
        }
        else
        {
            break;
        }
    }

    return bifurcationPos;
}

bool VascularTree::inVolume(Eigen::Vector3d point) const
{
    for (int idx = 0; idx < 3; ++idx)
    {
        if (point[idx] < 0 || point[idx] >= m_oxygenationMap->getDim()[idx])
        {
            return false;
        }
    }

    return true;
}

bool VascularTree::connectCandidate(Eigen::Vector3d point)
{
    // Candiate is too close to an existing segment
    if (!isCandidateValid(point, NULL))
    {
        return false;
    }
    
    if (m_nodeTable.numNodes() == 1)
    {
        Node *root = m_nodeTable.getNode(0);

        // Calculate heuristic upper bound for radius of root segment
        double radiusHeuristic = 2.0 * getRadiusHeuristic(point, root->pos, m_Qperf);
        if (!m_oxygenationMap->visible(root->pos, point) ||
            intersectsTrees(point, root->pos, radiusHeuristic))
        {
            return false;
        }

        // The "bifurcationPos" argument of connectPoint is ignore for
        // the root node, so just set it to (0.0, 0.0, 0.0)
        connectPoint(point, Eigen::Vector3d(0.0, 0.0, 0.0), root, m_nodeTable);
        return true;
    }

    std::vector<double> distances(m_nodeTable.numNodes());
    std::vector<int> indices(m_nodeTable.numNodes() - 1);
    std::iota(indices.begin(), indices.end(), 1); // indices = 1, 2, 3, ...

    // Determine the distance between the point and each segment
    for (int i = 1; i < m_nodeTable.numNodes(); ++i)
    {
        distances[i] = pointSegmentDistance(point, m_nodeTable.getNode(i));
    }

    // Sort segments by distance
    std::sort(indices.begin(),
              indices.end(),
              [&distances](int i, int j)
    {
        return distances[i] < distances[j];
    });

    int count = 0;

    double fitness;
    double bestFitness = 1e200;

    Eigen::Vector3d bifurcationPos;
    Eigen::Vector3d bestBifurcationPos;

    // Pointers to target nodes of bifurcated segments
    Node *bifurcationSegment;
    Node *bestBifurcationSegment = NULL;

    // Try the first "m_closestNeighbours" segments
    for (std::size_t j = 0; j < indices.size() && count < m_closestNeighbours; ++j)
    {
        bifurcationSegment = m_nodeTable.getNode(indices[j]);
        bifurcationPos = localOptimization(point, bifurcationSegment, fitness);

        if (fitness < 0)
        {
            continue;
        }

        count++;

        if (bifurcationHasIntersections(point, bifurcationPos, bifurcationSegment))
        {
            continue;
        }

        if (fitness < bestFitness)
        {
            bestFitness = fitness;
            bestBifurcationPos = bifurcationPos;
            bestBifurcationSegment = bifurcationSegment;
        }
    }

    // Could not connect candidate to ANY segment
    if (bestBifurcationSegment == NULL)
    {
        return false;
    }

    connectPoint(point, bestBifurcationPos, bestBifurcationSegment, m_nodeTable);
    return true;
}

void VascularTree::buildTreeStep()
{
    const int maxIter = 50;

    for (int i = 0; i < maxIter; ++i)
    {
        const Eigen::Vector3i term = m_oxygenationMap->candidate();

        if (connectCandidate(term.cast<double>()))
        {
            m_oxygenationMap->applyCandidate(term);
            return;
        }
    }

    throw std::runtime_error("Could not find suitable terminal node.");
}

double VascularTree::minDistance(Eigen::Vector3d target1, Eigen::Vector3d source1,
                                 Eigen::Vector3d target2, Eigen::Vector3d source2)
{
    const Eigen::Vector3d v1 = source1 - target1;
    const Eigen::Vector3d v2 = source2 - target2;

    // The line segments are given by p1(t1) = target1 + t1 * v1 and
    // p2(t2) = target2 + t2 * v2, with t1 and t2 in the interval [0, 1].
    // Phrase the seach for values (t1, t2) that yield the minimum
    // distance as a linear system of equations A * (t1, t2) = b.
    Eigen::Matrix2d A;
    A << v1.dot(v1), -v1.dot(v2),
         v1.dot(v2), -v2.dot(v2);
    const Eigen::Vector2d b(v1.dot(target2 - target1), v2.dot(target2 - target1));

    // Solve linear system A * t = b using the Eigen library
    Eigen::Vector2d t = A.colPivHouseholderQr().solve(b);
    t[0] = std::clamp(t[0], 0.0, 1.0);
    t[1] = std::clamp(t[1], 0.0, 1.0);

    const Eigen::Vector3d p1t = target1 + t[0] * v1;
    const Eigen::Vector3d p2t = target2 + t[1] * v2;
    return (p2t - p1t).norm();
}

std::vector<Node *> VascularTree::getNeighbours(Node *node)
{
    std::vector<Node *> neighbours = { node, node->leftChild, node->rightChild };

    if (node->parent != NULL)
    {
        neighbours.push_back(node->parent);

        Node *sibling = node->parent->leftChild == node ?
                        node->parent->rightChild : node->parent->leftChild;
        neighbours.push_back(sibling);
    }

    return neighbours;
}

bool VascularTree::bifurcationHasIntersections(Eigen::Vector3d point,
                                               Eigen::Vector3d bifurcationPos,
                                               Node *node)
{
    const double scale = 3.0;
    // Radius of the bifurcated segment
    double radius = node->radius / m_mapVoxelWidth;
    
    const Eigen::Vector3d oldChildPos = node->pos;
    const Eigen::Vector3d oldParentPos = node->parent->pos;
    std::vector<Node *> neighbours;

    const auto intersectionTests = [this, scale]
    (auto source, auto target, auto radius, auto neighbours)
    {
        return intersectsTrees(source, target, radius, scale) ||
               intersectsSelfSegment(source, target, radius, scale, neighbours);
    };

    neighbours = { node, node->leftChild, node->rightChild };
    if (intersectionTests(oldChildPos, bifurcationPos, radius, neighbours))
    {
        return true;
    }

    neighbours = { node->parent, node->parent->leftChild, node->parent->rightChild };
    if (intersectionTests(bifurcationPos, oldParentPos, radius, neighbours))
    {
        return true;
    }

    neighbours = { node };
    // Approximation for the radius of the new terminal segment
    radius = getRadiusHeuristic(point, bifurcationPos, m_Qterm) / m_mapVoxelWidth;
    if (intersectionTests(point, bifurcationPos, radius, neighbours))
    {
        return true;
    }

    return false;
}

bool VascularTree::intersectsTree(VascularTree &vt,
                                  Eigen::Vector3d target, Eigen::Vector3d source,
                                  double radius, double scale,
                                  const std::vector<Node *> &ignored) const
{
    const int size = vt.getNumberOfNodes();

    if (size == 1)
    {
        return false;
    }

    vt.calculateRadius();

    volatile bool intersection = false;
    // The "ignored" argument is used only when testing for self-intersections
    const bool selfIntersectionTest = ignored.size() > 0;

    #pragma omp parallel for shared(intersection)
    for (int i = 1; i < size; ++i)
    {
        // Simulate early return from the loop by skipping
        // all iterations once an intersection has been detected.
        // This way, the loop is suitable for use with omp.
        if (intersection)
        {
            continue;
        }

        const Node *node = vt.getNode(i);
        if (std::find(ignored.begin(), ignored.end(), node) != ignored.end())
        {
            continue;
        }

        const Eigen::Vector3d target2 = node->pos;
        const Eigen::Vector3d source2 = node->parent->pos;
        const double radius2 = node->radius / vt.getMapVoxelWidth();

        if (branchOBBIntersect(target, source, target2, source2, scale * radius, scale * radius2))
        {
            if (minDistance(target, source, target2, source2) > scale * (radius + radius2))
            {
                continue;
            }

            const Eigen::Vector3d targetRotated = source + (source - target);
            const Eigen::Vector3d target2Rotated = source2 + (source2 - target2);

            if (!selfIntersectionTest)
            {
                #pragma omp critical
                intersection = true;
            }
            else if (!branchOBBIntersect(target, source, source2, target2Rotated, radius, radius2) &&
                     !branchOBBIntersect(source, targetRotated, target2, source2, radius, radius2))
            {
                #pragma omp critical
                intersection = true;
            }
        }
    }

    return intersection;
}

bool VascularTree::intersectsTrees(Eigen::Vector3d target, Eigen::Vector3d source,
                                   double radius, double scale) const
{
    for (const auto &vt : m_vascularTrees)
    {
        auto vtLocked = vt.lock();

        if (vtLocked && vtLocked.get() != this &&
            intersectsTree(*vtLocked, target, source, radius, scale))
        {
            return true;
        }
    }

    return false;
}

bool VascularTree::intersectsSelfSegment(Eigen::Vector3d target, Eigen::Vector3d source,
                                  double radius, double scale,
                                  const std::vector<Node *> &ignored)
{
    if (!m_checkSelfIntersection)
    {
        return false;
    }

    return intersectsTree(*this, target, source, radius, scale, ignored);
}

bool VascularTree::intersectsSelfTree()
{
    if (!m_checkSelfIntersection)
    {
        return false;
    }

    for (int i = 1; i < getNumberOfNodes(); ++i)
    {
        Node *node = m_nodeTable.getNode(i);
        const Eigen::Vector3d target = node->pos;
        const Eigen::Vector3d source = node->parent->pos;
        const double radius = node->radius / m_mapVoxelWidth;
        
        std::vector<Node *> ignored = getNeighbours(node);
        if (intersectsSelfSegment(target, source, radius, 1.0, ignored))
        {
            return true;
        }
    }

    return false;
}
