#define _USE_MATH_DEFINES
#include <cmath>

#include "OBB.hpp"

BranchBoundingBox::BranchBoundingBox(Eigen::Vector3d target, Eigen::Vector3d source,
                                     double radius)
{
    const Eigen::Vector3d v = source - target;
    const double h = v.norm();

    m_halfLengths = Eigen::Vector3d(0.5 * h, radius, radius);
    m_pos = target + v / 2.0;

    // Create helper vector v0 with v0[argmin_i abs(v[i])] = 1
    // and all other entries set to 0
    Eigen::Vector3d v0(0.0, 0.0, 0.0);
    if (std::abs(v[0]) <= std::abs(v[1]) &&
        std::abs(v[0]) <= std::abs(v[2]))
    {
        v0[0] = 1.0;
    }
    else if (std::abs(v[1]) <= std::abs(v[0]) &&
             std::abs(v[1]) <= std::abs(v[2]))
    {
        v0[1] = 1.0;
    }
    else
    {
        v0[2] = 1.0;
    }

    // Set up OBB axes
    m_axes[0] = v.normalized();
    m_axes[1] = v0.cross(m_axes[0]).normalized();
    m_axes[2] = m_axes[0].cross(m_axes[1]);
}

CurvedBranchBoundingBox::CurvedBranchBoundingBox(
    Eigen::Vector3d target, Eigen::Vector3d source,
    double radius, double startHeight, int id) :
    BranchBoundingBox(target, source, radius), m_radius(radius),
    m_initialPosition(m_pos[2]), m_id(id)
{
    randomizeAxes();
    expand(startHeight - 2.0 * radius);
}

Eigen::Vector3d CurvedBranchBoundingBox::getCurveDirection() const
{
    return Eigen::Vector3d(m_axes[1]);
}

double CurvedBranchBoundingBox::getCurveHeight() const
{
    return 2.0 * m_halfLengths[2] - 2.0 * m_radius;
}

void CurvedBranchBoundingBox::initRandom(int randomSeed)
{
    s_rand = std::mt19937(randomSeed);
}

void CurvedBranchBoundingBox::randomizeAxes()
{
    // Randomly choose rotation angle
    const double theta = s_randDist(s_rand) * 2.0 * M_PI;

    // Implementation of Rodrigues' rotation formula for axes[1]
    m_axes[1] = cos(theta) * m_axes[1] + sin(theta) * m_axes[2];
    // Update axes[2]
    m_axes[2] = m_axes[0].cross(m_axes[1]);
}

void CurvedBranchBoundingBox::expand(double increment)
{
    m_pos[2] += increment / 2.0;
    m_halfLengths[2] += increment / 2.0;
}

void CurvedBranchBoundingBox::contract(double increment)
{
    if (m_halfLengths[2] - increment / 2.0 < m_radius)
    {
        m_pos[2] = m_initialPosition;
        m_halfLengths[2] = m_radius;
        return;
    }

    m_pos[2] -= increment / 2.0;
    m_halfLengths[2] -= increment / 2.0;
}

void CurvedBranchBoundingBox::initializeIntersections(
    std::vector<std::weak_ptr<CurvedBranchBoundingBox>> branches)
{
    for (const auto &box : branches)
    {
        if (auto boxLocked = box.lock(); boxLocked && OBBIntersect(*this, *boxLocked))
        {
            m_intersectsWith.push_back(box);
        }
    }
}

void CurvedBranchBoundingBox::initializeSelfIntersections(
    std::vector<std::weak_ptr<CurvedBranchBoundingBox>> branches)
{
    for (const auto &box : branches)
    {
        if (auto boxLocked = box.lock())
        {
            CurvedBranchBoundingBox thisReflected = this->getReflected();
            CurvedBranchBoundingBox boxReflected = boxLocked->getReflected();

            if (OBBIntersect(*this, *boxLocked) &&
                !OBBIntersect(*this, boxReflected) &&
                !OBBIntersect(thisReflected, *boxLocked))
            {
                m_selfIntersectsWith.push_back(box);
            }

        }
    }
}

void CurvedBranchBoundingBox::updateIntersections()
{
    updateIntersectionsDispatcher(m_intersectsWith);
}

void CurvedBranchBoundingBox::updateSelfIntersections()
{
    updateIntersectionsDispatcher(m_selfIntersectsWith);
}

void CurvedBranchBoundingBox::updateIntersectionsDispatcher(
    std::vector<std::weak_ptr<CurvedBranchBoundingBox>> &intersectsWith)
{
    std::vector<std::weak_ptr<CurvedBranchBoundingBox>> updated;

    for (const auto &box : intersectsWith)
    {
        if (auto boxLocked = box.lock(); boxLocked && OBBIntersect(*this, *boxLocked))
        {
            updated.push_back(box);
        }
    }

    intersectsWith = updated;
}

bool CurvedBranchBoundingBox::hasIntersections() const
{
    return m_intersectsWith.size() > 0;
}

bool CurvedBranchBoundingBox::hasSelfIntersections() const
{
    return m_selfIntersectsWith.size() > 0;
}

bool OBBIntersect(const BranchBoundingBox &box1, const BranchBoundingBox &box2)
{
    const double epsilon = std::numeric_limits<double>::epsilon();
    const Eigen::Vector3d r1 = box1.getHalfLengths();
    const Eigen::Vector3d r2 = box2.getHalfLengths();
    const Eigen::Vector3d axes1[3] = { box1.getAxis(0), box1.getAxis(1), box1.getAxis(2) };
    const Eigen::Vector3d axes2[3] = { box2.getAxis(0), box2.getAxis(1), box2.getAxis(2) };

    // Generate a rotation matrix that transforms
    // from world space to this OBB's coordinate space.
    double R[3][3];

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            R[i][j] = (axes1[i]).dot(axes2[j]);
        }
    }

    // Compute translation vector
    const Eigen::Vector3d t0 = box2.getPos() - box1.getPos();
    // Express the translation vector in a's coordinate frame.
    const Eigen::Vector3d t(t0.dot(axes1[0]), t0.dot(axes1[1]), t0.dot(axes1[2]));

    double AbsR[3][3];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            AbsR[i][j] = std::abs(R[i][j]) + epsilon;
        }
    }

    // Test the three major axes of this OBB.
    for (int i = 0; i < 3; ++i)
    {
        const double ra = r1[i];
        const double ti = t[i];
        const double rb = r2[0] * AbsR[i][0] + r2[1] * AbsR[i][1] + r2[2] * AbsR[i][2];

        if (std::abs(ti) > ra + rb)
        {
            return false;
        }
    }

    // Test the three major axes of the OBB b.
    for (int i = 0; i < 3; ++i)
    {
        const double ra = r1[0] * AbsR[0][i] + r1[1] * AbsR[1][i] + r1[2] * AbsR[2][i];
        const double rb = r2[i];

        if (std::abs(t[0] * R[0][i] + t[1] * R[1][i] + t[2] * R[2][i]) > ra + rb)
        {
            return false;
        }
    }

    // Test the 9 different cross-axes.

    // A.x <cross> B.x
    double ra = r1[1] * AbsR[2][0] + r1[2] * AbsR[1][0];
    double rb = r2[1] * AbsR[0][2] + r2[2] * AbsR[0][1];
    if (std::abs(t[2] * R[1][0] - t[1] * R[2][0]) > ra + rb)
    {
        return false;
    }

    // A.x < cross> B.y
    ra = r1[1] * AbsR[2][1] + r1[2] * AbsR[1][1];
    rb = r2[0] * AbsR[0][2] + r2[2] * AbsR[0][0];
    if (std::abs(t[2] * R[1][1] - t[1] * R[2][1]) > ra + rb)
    {
        return false;
    }

    // A.x <cross> B.z
    ra = r1[1] * AbsR[2][2] + r1[2] * AbsR[1][2];
    rb = r2[0] * AbsR[0][1] + r2[1] * AbsR[0][0];
    if (std::abs(t[2] * R[1][2] - t[1] * R[2][2]) > ra + rb)
    {
        return false;
    }

    // A.y <cross> B.x
    ra = r1[0] * AbsR[2][0] + r1[2] * AbsR[0][0];
    rb = r2[1] * AbsR[1][2] + r2[2] * AbsR[1][1];
    if (std::abs(t[0] * R[2][0] - t[2] * R[0][0]) > ra + rb)
    {
        return false;
    }

    // A.y <cross> B.y
    ra = r1[0] * AbsR[2][1] + r1[2] * AbsR[0][1];
    rb = r2[0] * AbsR[1][2] + r2[2] * AbsR[1][0];
    if (std::abs(t[0] * R[2][1] - t[2] * R[0][1]) > ra + rb)
    {
        return false;
    }

    // A.y <cross> B.z
    ra = r1[0] * AbsR[2][2] + r1[2] * AbsR[0][2];
    rb = r2[0] * AbsR[1][1] + r2[1] * AbsR[1][0];
    if (std::abs(t[0] * R[2][2] - t[2] * R[0][2]) > ra + rb)
    {
        return false;
    }

    // A.z <cross> B.x
    ra = r1[0] * AbsR[1][0] + r1[1] * AbsR[0][0];
    rb = r2[1] * AbsR[2][2] + r2[2] * AbsR[2][1];
    if (std::abs(t[1] * R[0][0] - t[0] * R[1][0]) > ra + rb)
    {
        return false;
    }

    // A.z <cross> B.y
    ra = r1[0] * AbsR[1][1] + r1[1] * AbsR[0][1];
    rb = r2[0] * AbsR[2][2] + r2[2] * AbsR[2][0];
    if (std::abs(t[1] * R[0][1] - t[0] * R[1][1]) > ra + rb)
    {
        return false;
    }

    // A.z <cross> B.z
    ra = r1[0] * AbsR[1][2] + r1[1] * AbsR[0][2];
    rb = r2[0] * AbsR[2][1] + r2[1] * AbsR[2][0];
    if (std::abs(t[1] * R[0][2] - t[0] * R[1][2]) > ra + rb)
    {
        return false;
    }

    // No separating axis exists, so the two OBB don't intersect.
    return true;
}

bool branchOBBIntersect(Eigen::Vector3d target1, Eigen::Vector3d source1,
                        Eigen::Vector3d target2, Eigen::Vector3d source2,
                        double radius1, double radius2)
{
    const BranchBoundingBox box1(target1, source1, radius1);
    const BranchBoundingBox box2(target2, source2, radius2);

    return OBBIntersect(box1, box2);
}

// Initialize static member variables.
// Initialize rand with dummy seed because it is re-initialized later by
// a call to CurvedBranchBoundingBox::initRandom
std::mt19937 CurvedBranchBoundingBox::s_rand(0);
std::uniform_real_distribution<double>
CurvedBranchBoundingBox::s_randDist(0.0, 1.0);
