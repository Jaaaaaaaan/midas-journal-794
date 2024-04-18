#pragma once

#include <deque>
#include <Eigen/Dense>
#include <memory>
#include <random>
#include <vector>

class CurvedBranchBoundingBox;

/**
 * \brief Bounding box for vessel branches, based on the BoundingBox class
 * in the MathGeoLib library (https://github.com/juj/MathGeoLib).
 */
class BranchBoundingBox
{
protected:

    /**
     * \brief Center of the bounding box.
     */
    Eigen::Vector3d m_pos;

    /**
     * \brief Half lengths of the bounding box along its axes.
     *
     * Describes the extents of the bounding box.
     * This information is stored as half lengths by convention.
     */
    Eigen::Vector3d m_halfLengths;

    /**
     * \brief Normalized axes of the bounding box.
     */
    Eigen::Vector3d m_axes[3];

public:
    /**
     * \brief Computes center, half lengths and axes
     * of a bounding box for the vessel branch described by the arguments.
     *
     * \param target The target position of the branch.
     * \param source The source position of the branch.
     * \param radius The radius of the branch.
     */
    BranchBoundingBox(Eigen::Vector3d target, Eigen::Vector3d source, double radius);

    /**
     * \brief Getter for the bounding box height, i.e. its extent along the
     * direction orthogonal to (source - target).
     */
    double getHeight() const
    {
        return m_halfLengths[2];
    }

    /**
     * \brief Getter for the bounding box center.
     */
    Eigen::Vector3d getPos() const
    {
        return m_pos;
    }

    /**
     * \brief Getter for the bounding box half lengths.
     */
    Eigen::Vector3d getHalfLengths() const
    {
        return m_halfLengths;
    }

    /**
     * \brief Getter for the bounding box axes.
     * 
     * \param index Specifies which axis vector is to be retrieved.
     */
    Eigen::Vector3d getAxis(std::size_t index) const
    {
        return m_axes[index];
    }
};

/**
 * \brief Bounding box for curved vessel branches.
 *
 * To generate non-intersecting curved branches, a direction
 * along which curvature is allowed is assigned randomly to each branch.
 * This is called the "direction of curvature".
 *
 * In addition, a maximum height of the curved branch along this direction
 * is determined such that there are no intersections with
 * other curved branches. This is refered to as the "curve height".
 *
 * To compute the curve height, bounding boxes are incrementally updated
 * when generating curved branches. During this process,
 * the bounding box keeps track of its remaining intersections with other
 * bounding boxes in the m_intersectsWith vector.
 */
class CurvedBranchBoundingBox : public BranchBoundingBox
{
private:
    std::vector<std::weak_ptr<CurvedBranchBoundingBox>> m_intersectsWith;
    std::vector<std::weak_ptr<CurvedBranchBoundingBox>> m_selfIntersectsWith;

    static std::mt19937 s_rand;
    static std::uniform_real_distribution<double> s_randDist;

    double m_radius;
    double m_initialPosition;
    int m_id;

    /**
     * \brief Executes updates to the m_intersectsWith and
     * m_selfIntersectsWith vectors.
     * 
     * \param intersectsWith The vector to be updated.
     */
    void updateIntersectionsDispatcher(
        std::vector<std::weak_ptr<CurvedBranchBoundingBox>> &intersectsWith);

public:

    /**
     * \brief Constructs a bounding box instance for the computation of branch
     * curvature parameters for a given branch.
     *
     * \param target The target position of the branch.
     * \param source The source position of the branch.
     * \param radius The radius of the branch.
     * \param startHeight The initial curve height.
     */
    CurvedBranchBoundingBox(Eigen::Vector3d target, Eigen::Vector3d source,
                            double radius, double startHeight, int id);

    /**
     * \brief Getter for the target node id of the branch bounded by this
     * bounding box.
     */
    int getId() const
    {
        return m_id;
    }

    /**
     * \brief Getter for the direction of curvature determined by
     * the \see CurvedBranchBoundingBox::randomizeAxes function.
     */
    Eigen::Vector3d getCurveDirection() const;

    /**
     * \brief Getter for the computed curve height.
     */
    double getCurveHeight() const;

    /**
     * \brief Creates a bounding box for the branch that is obtained by
     * "reflecting" this branch, i.e. by inverting the direction of the
     * branch curvature in the local coordinate system.
     */
    CurvedBranchBoundingBox getReflected() const
    {
        CurvedBranchBoundingBox res(*this);

        res.m_axes[1] *= -1;
        res.m_axes[2] *= -1;
        res.m_pos += (2.0 * m_halfLengths[1] - m_radius) * m_axes[1];

        return res;
    }

    /**
     * \brief Initializes the random generator for
     * \see CurvedBranchBoundingBox::randomizeAxes.
     *
     * \param randomSeed The random seed for the random generator.
     * If this is a negative number, a random seed is selected
     * pseudo-randomly instead.
     */
    static void initRandom(int randomSeed);

    /**
     * \brief Rotates the bounding box coordinate system by a
     * random angle around the vector source - target, where target and source
     * refer to the arguments of
     * \see CurvedBranchBoundingBox::CurvedBranchBoundingBox.
     *
     * This function is used to randomly determine a direction of curvature.
     */
    void randomizeAxes();

    /**
     * \brief Increases the curve height.
     *
     * \param increment The increment by which to expand the height.
     */
    void expand(double increment);

    /**
     * \brief Decreases the curve height.
     * This function has no effect if the contracted height
     * of the bounding box would be less than the branch diameter.
     *
     * \param increment The increment by which to contract the height.
     */
    void contract(double increment);

    /**
     * \brief Tests all bounding boxes of other vascular trees for
     * intersections with this bounding box.
     * 
     * Let T be the vascular tree that contains this branch.
     * This function tests branches in all vascular trees T' != T
     * for intersections with this branch. Intersections that have been found
     * are kept track of in the m_intersectsWith vector.
     *
     * \param branches Contains pointers to curved branch bounding boxes
     * for all branches in other vascular trees.
     */
    void initializeIntersections(
        std::vector<std::weak_ptr<CurvedBranchBoundingBox>> branches);

    /**
     * \brief Tests all bounding boxes of this vascular tree for
     * intersections with this bounding box.
     * 
     * Let T be the vascular tree that contains this branch.
     * This function tests branches in T for intersections with this branch.
     * Intersections that have been found are kept track of in the
     * m_selfIntersectsWith vector.
     *
     * \param branches Contains pointers to curved branch bounding boxes
     * for all branches in this vascular tree.
     */
    void initializeSelfIntersections(
        std::vector<std::weak_ptr<CurvedBranchBoundingBox>> branches);

    /**
     * \brief Updates the m_intersectsWith with vector.
     */
    void updateIntersections();

    /**
     * \brief Updates the m_selfIntersectsWith with vector.
     */
    void updateSelfIntersections();

    /**
     * \brief Checks whether the m_intersectsWith vector is empty.
     *
     * \return True if there are intersections with bounding boxes of other
     * vascular trees, and false otherwise.
     */
    bool hasIntersections() const;

    /**
     * \brief Checks whether the m_selfIntersectsWith vector is empty.
     *
     * \return True if there are intersections with bounding boxes
     * of this vascular tree, and false otherwise.
     */
    bool hasSelfIntersections() const;
};

/**
 * \brief Implementation of oriented bounding box (OBB) intersection test.
 *
 * Adapted from the MathGeoLib library
 * (https://github.com/juj/MathGeoLib) under the Apache license.
 *
 * \param box1 A bounding box to test for intersections
 * \param box2 A bounding box to test for intersections
 *
 * \return True if an intersection was detected, and false otherwise.
 */
bool OBBIntersect(const BranchBoundingBox &box1, const BranchBoundingBox &box2);

/**
 * \brief Wrapper function for \see OBBIntersect that
 * applies a bounding box intersection test to two vessel branches.
 *
 * \param target1 The target position of the first vessel branch.
 * \param source1 The source position of the first vessel branch.
 * \param target2 The target position of the second vessel branch.
 * \param source2 The source position of the second vessel branch.
 * \param radius1 The radius of the first vessel branch.
 * \param radius2 The radius of the second vessel branch.
 *
 * \return True if an intersection was detected, and false otherwise.
 */
bool branchOBBIntersect(Eigen::Vector3d target1, Eigen::Vector3d source1,
                        Eigen::Vector3d target2, Eigen::Vector3d source2,
                        double radius1, double radius2);
