/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: VascularTree.h,v $
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

#include <memory>
#include <vector>

#include "NodeTable.hpp"
#include "OxygenationMap.hpp"

/**
 * \brief Class that iteratively builds the vascular tree.
 *
 * Iteratively builds the vascular structure by selecting a new candidate
 * node and connecting it to an existing edge,
 * creating a bifurcation location.
 *
 * Since vascular trees are assumed to be binary trees here,
 * each child node has exactly one parent node.
 * Thus, the edge connecting the two nodes can be uniquely identified using
 * only its target node.
 * 
 * As a result, information about edges in the tree can be stored
 * together with the respective target nodes, and there is no need for a separate
 * edge data structure.
 * Similarly, iteration over the edges of the tree is implemented using
 * a loop such as "for (int i = 1; i < m_nodeTable.numNodes(); ++i)",
 * where one iterates over all child nodes in the tree, i.e. all nodes
 * except for the root node.
 */
class VascularTree
{
private:
    NodeTable m_nodeTable;
    std::string m_outputFilename;
    std::vector<std::weak_ptr<VascularTree>> m_vascularTrees;

    double m_Pperf;
    double m_Pterm;
    double m_Qperf;
    double m_Qterm;

    double m_rho;
    double m_gamma;
    double m_lambda;
    double m_mu;

    double m_minDistance;
    double m_mapVoxelWidth; // in SI units

    std::unique_ptr<OxygenationMap> m_oxygenationMap;

    int m_numTerminalNodes;
    int m_closestNeighbours;
    bool m_checkSelfIntersection;

public:

    /**
     * \brief Constructs a new vascular tree instance.
     *
     * For more information on the parameters describing a vascular
     * tree, see:
     * Preet Jassi and Ghassan Hamarneh
     * VascuSynth: Vascular Tree Synthesis Software.
     * Insight Journal (2011).
     */
    VascularTree(Eigen::Vector3d perforationPoint, std::string outputFilename,
                 double Pperf, double Pterm, double Qperf,
                 double rho, double gamma, double lambda, double mu,
                 double minDistance, double mapVoxelWidth, int numTerminalNodes,
                 int closestNeighbours, bool checkSelfIntersection);

    /**
     * \brief Setter for the m_vascularTrees vector that stores pointers
     * to all the vascular tree instances that are
     * being generated simultaneously.
     */
    void setVascularTrees(std::vector<std::shared_ptr<VascularTree>> vascularTrees);

    /**
     * \brief Setter for the oxygenation map of the vascular tree.
     */
    void setOxygenationMap(std::unique_ptr<OxygenationMap> oxygenationMap)
    {
        m_oxygenationMap = std::move(oxygenationMap);
    }

    /**
     * \brief Getter for the number of nodes currently in the tree.
     */
    int getNumberOfNodes() const
    {
        return (int)m_nodeTable.numNodes();
    }

    /**
     * \brief Getter for the number of nodes in the final tree.
     */
    int getTargetNumberOfNodes() const
    {
        return 2 * m_numTerminalNodes;
    }

    /**
     * \brief Getter for the map voxel width.
     */
    double getMapVoxelWidth() const
    {
        return m_mapVoxelWidth;
    }

    /**
     * \brief Getter for individual nodes in the node table of the tree.
     *
     * \param index The index in the node table.
     */
    Node *getNode(int index)
    {
        return m_nodeTable.getNode(index);
    }

    const Node *getNode(int index) const
    {
        return m_nodeTable.getNode(index);
    }

    /**
     * \brief Getter for the output filename of the vascular tree's
     * XML representation.
     */
    std::string getOutputFilename() const
    {
        return m_outputFilename;
    }

    /**
     * \brief Calculates a heuristic upper bound for the radius of a segment.
     *
     * \param target The target position of the segment.
     * \param source The source position of the segment.
     * \param flow The blood flow through the segment.
     *
     * \return The computed upper bound.
     */
    double getRadiusHeuristic(Eigen::Vector3d target, Eigen::Vector3d source, double flow) const;

    /**
     * \brief Calculates the distance between two nodes.
     *
     * \param from A pointer to the source node.
     * \param to A pointer to the target node.
     *
     * \return The computed distance.
     */
    double distanceBetweenTwoNodes(const Node *source, const Node *target) const;

    /**
     * \brief Calculates the reduced resistance of the vessel segment
     * with target node "node" and stores it in "node".
     *
     * \param node: A pointer to the target node of the segment
     * under consideration.
     */
    void calculateReducedResistance(Node *node);

    /**
     * \brief Calculates the radius ratios for the vessel segment
     * with target node "node" and stores it in "node".
     *
     * \param node: A pointer to the target node of the segment
     * under consideration.
     */
    void calculateRatios(Node *node);

    /** \brief Recusively updates the reduced resistance and radius
     * ratios throughout the vascular tree.
     *
     * \param node: Pointer to the node that is to be updated.
     * The parent of this node will be updated recursively.
     */
    void updateAtBifurcation(Node *node);

    /**
     * \brief Recursively computes segment radii throughout the tree
     */
    void calculateRadius()
    {
        calculateRadius(m_nodeTable);
    }

    /**
     * \brief Recursively computes segment radii in a given node table.
     * The node table does not necessarily have to be
     * the tree's own node table.
     *
     * \param nt The node table for which to compute segment radii.
     */
    void calculateRadius(NodeTable &nt);

    /**
     * \brief Recursively computes segment radii for the subtree starting
     * at the given node.
     *
     * \param node The starting node for the radius computation.
     */
    void calculateRadius(Node *node);

    /**
     * \brief Evaluates the cost function for a vascular tree.
     *
     * \param nt The node table describing the vascular tree.
     *
     * \return The value of the cost function.
     */
    double calculateFitness(const NodeTable &nt) const;

    /**
     * \brief Checks whether the distance from a candidate terminal node
     * to all vessel segments is greater than m_minDistance.
     *
     * \param candidate The position of the candidate terminal node.
     * \param ignored Specifies the target node of a segment for which
     * the distance criterion should not be tested.
     * If ignored is set to NULL, no segment is ignored.
     *
     * \return True if the distance criterion is satisfied for all segments,
     * false otherwise.
     */
    bool isCandidateValid(Eigen::Vector3d candidate, const Node *ignored) const;

    /**
     * \brief Given a segment, a bifurcation location, and a candidate
     * terminal node position, this function connects the candidate node
     * to the segment at the bifurcation location.
     *
     * \param point The position of the candidate terminal node.
     * \param bifurcationPos The position at which to place the bifurcation.
     * \param node The target node of the segment to which the candidate
     * node is being connected,
     * \param nt The node table on which to apply the connection.
     */
    void connectPoint(Eigen::Vector3d point, Eigen::Vector3d bifurcationPos,
                      Node *node, NodeTable &nt);

    /**
     * \brief Recursively updates the flow throughout the vascular tree.
     *
     * \param node Pointer to the node that is to be updated.
     * The parent of this node will be updated recursively.
     */
    void incrementFlow(Node *node);

    /**
     * \brief Computes the distance between a point and a vessel segment.
     *
     * \param point The point for which to compute the distance.
     * \param node A pointer to the terminal node of the vessel
     * segment for which to compute the distance.
     *
     * \return The computed distance.
     */
    double pointSegmentDistance(Eigen::Vector3d point, const Node *node) const;

    /**
     * \brief Attempts to find the optimal bifurcation location
     * for connecting a new candidate node to a given segment.
     *
     * \param point The position of the candidate terminal node to
     * be connected to the vascular tree.
     * \param node The target node of the segment to which "point"
     * is to be connected.
     * \param bestFitness Used to store the computed value of the cost
     * function for the optimized bifurcation location.
     *
     * \return The optimized bifurcation location.
     */
    Eigen::Vector3d localOptimization(Eigen::Vector3d point, const Node *node,
                                      double &bestFitness);

    /**
     * \brief Determines whether a point is positioned within the bounds
     * of the oxygen demand map of this tree.
     *
     * \param point The point to be tested.
     *
     * \return True if the point is within bounds, and false otherwise.
     */
    bool inVolume(Eigen::Vector3d point) const;

    /**
     * \brief Attempts to connect a candidate terminal node
     * to an existing segment of the vascular tree
     * such that the cost function is minimized.
     *
     * Connections to the m_closestNeighbours closest vessel segments
     * to point are evaluated.
     * The function VascularTree::localOptimization is used to optimized
     * the bifurcation location.
     *
     * \param point The position of the candidate terminal node.
     *
     * \return True if the candidate node was successfully connected to
     * the vascular tree, false otherwise.
     */
    bool connectCandidate(Eigen::Vector3d point);

    /**
     * \brief Executes one iteration step of vessel tree generation.
     */
    void buildTreeStep();

    /**
     * \brief Given two line segments with endpoints p00, p01 and p10, p11,
     * respectively, this function determines the minimum distance between
     * the two line segments.
     * 
     * \param p00 The starting point of the first line segment.
     * \param p01 The end point of the first line segment.
     * \param p10 The starting point of the second line segment.
     * \param p11 The end point of the second line segment.
     * 
     * \return The determined minimum distance.
    */
    static double minDistance(Eigen::Vector3d p00, Eigen::Vector3d p01,
                              Eigen::Vector3d p10, Eigen::Vector3d p11);

    /**
     * \brief Finds all adjacent segments of the segment with target node
     * "node".
     * 
     * \param node A pointer to the target node of the segment
     * for which to find neighbours.
     * 
     * \return Target nodes of all neighbouring segments, including
     * the starting segment itself.
    */
    static std::vector<Node *> getNeighbours(Node *node);

    /**
     * \brief When a bifurcation is added to the vascular tree,
     * the bifurcated segment is replaced by three new segments.
     * This method tests whether any of these three segments
     * intersects existing segments in the other vascular trees.
     * In addition, the function tests for intersections with existing
     * segments in this vascular tree if the m_checkSelfIntersections flag
     * is set.
     * 
     * \param point The position of a candidate terminal node that is to be
     * connected to the tree.
     * \param bifurcationPos The position of the bifurcation.
     * \param node A pointer to the terminal node of the segment that
     * is being bifurcated.
     * 
     * \return True if any intersections were found, and false otherwise.
    */
    bool bifurcationHasIntersections(Eigen::Vector3d point,
                                     Eigen::Vector3d bifurcationPos,
                                     Node *node);

    /**
     * \brief Checks whether a candidate segment intersects
     * any existing segment in a given vascular tree.
     *
     * \param vt The vascular tree which to test for intersections.
     * \param target The target position of the candidate segment.
     * \param source The source position of the candidate segment.
     * \param radius The radius of the candidate segment.
     * \param scale As the tree grows, segment radii in the tree become larger
     * because the blood flow through segments increases.
     * To account for this, all radii are scaled up by this scale factor
     * when testing for intersections  during tree growth.
     * \param ignored A vector of nodes in the vascular tree "vt"
     * for which no intersection tests should be perfomed.
     * This is used when checking for self-intersections.
     *
     * \return True if an intersection is detected, and false otherwise.
     */
    bool intersectsTree(VascularTree &vt, Eigen::Vector3d target, Eigen::Vector3d source,
                        double radius, double scale,
                        const std::vector<Node *> &ignored = std::vector<Node *>()) const;

    /**
     * \brief Wrapper function for \see VascularTree::intersectsTree that
     * checks for intersections with all of the other vascular trees.
     */
    bool intersectsTrees(Eigen::Vector3d target, Eigen::Vector3d source,
                         double radius, double scale = 1.0) const;

    /**
     * \brief Checks whether a candidate segment intersects
     * any existing segment in this vascular tree.
     *
     * \param target The target position of the candidate segment.
     * \param source The source position of the candidate segment.
     * \param radius The radius of the candidate segment.
     * \param scale As the tree grows, segment radii in the tree become larger
     * because the blood flow through segments increases.
     * To account for this, all radii are scaled up by this scale factor
     * when testing for intersections  during tree growth.
     * \param ignored A vector of nodes in the vascular tree
     * for which no intersection tests should be perfomed.
     *
     * \return True if an intersection is detected, and false otherwise.
     */
    bool intersectsSelfSegment(Eigen::Vector3d target, Eigen::Vector3d source,
                        double radius, double scale,
                        const std::vector<Node *> &ignored);

    /**
     * \brief Wrapper function for \see VascularTree::intersectsSelfSegment
     * that checks the entire tree for self intersections.
     */
    bool intersectsSelfTree();
};
