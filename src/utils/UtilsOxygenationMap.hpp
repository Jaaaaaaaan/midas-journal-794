#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Utils::OxygenationMap
{
    /**
     * \brief Loads an ODM from a file.
     *
     * This function returns a pointer to the loaded map array
     * so that multiple vascular trees can share the same map array.
     *
     * For more information on the layout of ODM files, see:
     * Preet Jassi and Ghassan Hamarneh
     * VascuSynth: Vascular Tree Synthesis Software.
     * Insight Journal (2011).
     *
     * \param filepath Path to the ODM file.
     * \param dim Empty vector to which the extents of the map are written.
     *
     * \return A pointer to the flattened three dimensional array
     * describing the oxygen demand.
     */
    std::shared_ptr<std::vector<double>>
    loadMap(const std::string &filepath, Eigen::Vector3i &dim);

    /**
     * \brief Deep copies a map array
     *
     * \param map A pointer to the array that is being copied.
     * \param dim The extents of the array to copy.
     *
     * \return A pointer to the copied array.
     */
    std::shared_ptr<std::vector<double>>
    copyMap(std::shared_ptr<std::vector<double>> map, Eigen::Vector3i dim);
}