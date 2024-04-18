#include <fstream>
#include <iostream>

#include "UtilsOxygenationMap.hpp"

// A parametric macro to allow the use of dynamic multi-dimensional arrays
#define arr(arr, x, y, z, dim) arr->at(z + dim[2] * (y + dim[1] * (x)))

std::shared_ptr<std::vector<double>>
Utils::OxygenationMap::loadMap(const std::string &filepath, Eigen::Vector3i &dim)
{
    std::ifstream mapFile;

    mapFile.open(filepath, std::ios::in);
    std::string line;

    if (!mapFile.is_open())
    {
        throw std::runtime_error("Could not read the OxygenationMap file");
    }

    std::getline(mapFile, line);
    char *start = const_cast<char *>(line.c_str());
    int index = 0;

    while (*start != '\0' && *start != '\n')
    {
        dim[index] = strtol(start, &start, 10);
        index++;
    }

    auto map = std::make_shared<std::vector<double>>(
        std::vector<double>(dim[0] * dim[1] * dim[2]));

    for (int i = 0; i < dim[0]; i++)
    {
        for (int j = 0; j < dim[1]; j++)
        {
            for (int k = 0; k < dim[2]; k++)
            {
                arr(map, i, j, k, dim) = 0.0;
            }
        }
    }

    std::getline(mapFile, line);

    while (!mapFile.eof())
    {
        char *start = const_cast<char *>(line.c_str());
        int index = 0;
        int region[6];

        while (*start != '\0' && *start != '\n')
        {
            region[index] = strtol(start, &start, 10);
            index++;
        }

        std::getline(mapFile, line);
        const double value = atof(line.c_str());

        for (int i = region[0]; i < region[3]; i++)
        {
            for (int j = region[1]; j < region[4]; j++)
            {
                for (int k = region[2]; k < region[5]; k++)
                {
                    arr(map, i, j, k, dim) = value;
                }
            }
        }

        std::getline(mapFile, line);
    }

    mapFile.close();
    return map;
}

std::shared_ptr<std::vector<double>>
Utils::OxygenationMap::copyMap(std::shared_ptr<std::vector<double>> map, Eigen::Vector3i dim)
{
    auto copy = std::make_shared<std::vector<double>>(
        std::vector<double>(dim[0] * dim[1] * dim[2]));

    for (int i = 0; i < dim[0]; ++i)
    {
        for (int j = 0; j < dim[1]; ++j)
        {
            for (int k = 0; k < dim[2]; ++k)
            {
                arr(copy, i, j, k, dim) = arr(map, i, j, k, dim);
            }
        }
    }

    return copy;
}