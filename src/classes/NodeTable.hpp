/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: NodeTable.h,v $
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
#include <vector>

/**
 * \brief Describes the type of a node.
 *
 * All nodes are either terminal nodes or bifurcation nodes or the root node.
 */
enum class NodeType
{
    /** A terminal node. */
    TERM,
    /** The root node. */
    ROOT,
    /** A bifurcation node. */
    BIF,
};

/**
 * \brief Represents a node in a vascular tree.
 *
 * Nodes are described using a number of parameters such as their radius
 * and the blood flow through them.
 * These parameters are updated as the tree grows and the node
 * gains more children.
 *
 * This struct also contains connectivity information,
 * i.e. pointers to a parent node and direct child nodes.
 * Since we are modelling vascular trees as binary trees,
 * each node can have up to two children.
 * The children are labelled as a "left child" and "right child",
 * but these labels are interchangable.
 * They are used only to distinguish the two children.
 */
struct Node
{
    /**
     * \brief The position of the node.
     */
    Eigen::Vector3d pos;

    /**
     * \brief The blood flow through the segment ending on this node.
     */
    double flow;

    /**
     * \brief The radius ratio for the left child segment.
     *
     * Let r_parent be the radius of the segment ending at this node.
     * Let r_left be the radius of the segment from this node to its
     * left child. Then leftRatio = r_left / r_parent.
     */
    double leftRatio;

    /**
     * \brief The radius ratio for the right child segment.
     *
     * Let r_parent be the radius of the segment ending at this node.
     * Let r_right be the radius of the segment from this node to its
     * right child. Then rightRatio = r_right / r_parent.
     */
    double rightRatio;

    /**
     * \brief The radius of the segment ending on this node.
     */
    double radius = 0.0;

    /**
     * \brief The reduced resistance of the segment ending on this node.
     */
    double reducedResistance = 0.0;

    /**
     * \brief A pointer to the parent node.
     */
    Node *parent;

    /**
     * \brief A pointer to the left child node.
     */
    Node *leftChild = NULL;

    /**
     * \brief A pointer to the right child node.
     */
    Node *rightChild = NULL;

    /**
     * \brief The index of the node in its node table.
     *
     * This is used when generating an XML file representing a vascular tree.
     */
    int index;

    /**
     * \brief The type of the node.
     */
    NodeType type;

    Node() {}

    Node(Eigen::Vector3d pos, double flow, double leftRatio, double rightRatio,
         Node *parent, int index, NodeType type) :
         pos(pos), flow(flow), leftRatio(leftRatio), rightRatio(rightRatio),
         parent(parent), index(index), type(type)
         {}

    /**
     * \brief If this node has neither a leftChild nor a rightChild,
     * this method sets the \see Node::leftChild variable.
     * If there already is a left child, this method sets
     * the \see Node::rightChild variable.
     *
     * \throws Invalid argument exception if the node already has
     * both a left and right child.
     *
     * \param child A pointer to the new child node.
     */
    void addChild(Node *child)
    {
        if (leftChild == NULL)
        {
            leftChild = child;
            return;
        }
        if (rightChild == NULL)
        {
            rightChild = child;
            return;
        }

        throw std::invalid_argument("Added too many children to node!");
    }

    /**
     * \brief Resets the \see Node::leftChild and
     * \see Node::rightChild members of this class.
     */
    void unlinkChildren()
    {
        leftChild = NULL;
        rightChild = NULL;
    }
};

/**
 * \brief Stores nodes and edge information that makes up a vascular tree.
 *
 * For each node, this class stores an instance of the \see Node struct.
 */
class NodeTable
{
private:
    std::vector<Node> m_nodes;
    int m_index = 0;
    int m_numTerminalNodes;

public:

    /**
     * \brief Constructs a new node table.
     *
     * The total number of nodes that can be stored in the table is
     * 2 * numTerminalNodes + 1:
     * The number of nodes in a binary tree is twice the number of terminal
     * nodes. An additional slot is needed during construction of the tree.
     *
     * \param numTerminalNodes The number of terminal nodes in the vascular
     * tree represented by this table.
     */
    NodeTable(int numTerminalNodes);

    /**
     * \brief Getter for the number of nodes currently stored in the table.
     */
    int numNodes() const
    {
        return m_index;
    }

    /**
     * \brief Fetches a pointer to a node from the table.
     *
     * \param index The index of the desired node in the table.
     *
     * \return A pointer to the node at the given index.
     */
    Node *getNode(int index);
    const Node *getNode(int index) const;

    /**
     * \brief Adds a new node to the table.
     *
     * Note that this method does not check whether there is enough space
     * in the table for a new node.
     *
     * \param pos The position of the node.
     * \param flow The blood flow through the node.
     * \param leftRatio The radius ratio for the left child.
     * See \see Node::leftRatio for details.
     * \param rightRatio The radius ratio for the right child.
     * See \see Node::rightRatio for details.
     * \param parent A pointer to the parent node. This function automatically
     * calls the \see Node::addChild method of the parent node.
     * Hence, a pointer to the newly generated node is automataically stored
     * in the parent node.
     * \param type The node type.
     *
     * \return A pointer to the newly added node.
     */
    Node *addNode(Eigen::Vector3d pos, double flow, double leftRatio,
                  double rightRatio, Node *parent, NodeType type);

    /**
     * \brief Generates a deep copy of the table.
     *
     * \return A pointer to a copy of the node table.
     */
    std::unique_ptr<NodeTable> getCopy() const;
};
