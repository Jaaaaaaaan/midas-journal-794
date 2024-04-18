/*=========================================================================

    Program: VascuSynth
    Module: $RCSfile: NodeTable.cpp,v $
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

#include "NodeTable.hpp"

NodeTable::NodeTable(int numTerminalNodes) : m_numTerminalNodes(numTerminalNodes)
{
    // Note that it is crucial that this vector is initialized
    // with a fixed size and never resized!
    // Otherwise, producing pointers to elements of this vectors
    // becomes unsafe.
    m_nodes = std::vector<Node>(2 * numTerminalNodes + 1);
}

Node *NodeTable::getNode(int index)
{
    return &m_nodes[index];
}

const Node *NodeTable::getNode(int index) const
{
    return &m_nodes[index];
}

Node *NodeTable::addNode(Eigen::Vector3d pos, double flow, double leftRatio,
                         double rightRatio, Node *parent, NodeType type)
{
    m_nodes[m_index] = Node(pos, flow, leftRatio, rightRatio, parent,
                            numNodes(), type);

    auto *node = getNode(m_index);
    m_index++;

    if (parent != NULL)
    {
        parent->addChild(node);
    }

    return node;
}

std::unique_ptr<NodeTable> NodeTable::getCopy() const
{
    auto res = std::make_unique<NodeTable>(m_numTerminalNodes);

    res->m_index = m_index;

    // Copy nodes into new table
    for (int i = 0; i < numNodes(); ++i)
    {
        res->m_nodes[i] = Node(m_nodes[i]);
    }

    // Update pointers
    for (int i = 0; i < numNodes(); ++i)
    {
        Node *node = res->getNode(i);

        if (node->parent != NULL)
        {
            node->parent = res->getNode(node->parent->index);
        }

        if (node->leftChild != NULL)
        {
            node->leftChild = res->getNode(node->leftChild->index);
        }

        if (node->rightChild != NULL)
        {
            node->rightChild = res->getNode(node->rightChild->index);
        }
    }

    return res;
}
