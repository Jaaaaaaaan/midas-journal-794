# Description

This repository contains modified source files of the VascuSynth software [1][2] that I used during my undergraduate thesis at NCT TSO Dresden.
Additional modifications were made as part of a student assistant position at NCT TSO Dresden.

Differences to the original VascuSynth source code include support for generating multiple non-intersecting vascular trees in a shared volume, simple parallelization, and basic documentation of all functions.
In addition, the modified software exports information on maximal bounding boxes for vessel branches.
Based on this information, I was able to generate 3D models of non-intersecting, curved vessel trees with the help of a Python wrapper.
The `TreeDrawer` class from the original project was removed in order to allow for compilation without the ITK library.

The syntax used to call the program is similar to the one in the original VascuSynth project:
```
VascuSynth [oxygenDemandMap] [supplyMap] [randomSeed] [paramFile 1] ... [paramFile N]
```
For more information on the contents of the map and parameter files, see Ref. [1] as well as the `VascuSynth.cpp` file.

# Installation

On Linux, the project can be compiled using the CMake files provided in this repository.
Note that the source code depends on the Eigen3 library.

# Third Party Code

A routine for intersection tests on oriented bounding boxes was adapted from the open source MathGeoLib library [3].
The code is contained in the `OBB.cpp` and `OBB.hpp` files and explicitly marked as third party code.

# References

[1] P. Jassi and G. Hamarneh, “VascuSynth: Vascular Tree Synthesis Software,” The Insight Journal, Apr. 2011.

[2] P. Jassi and G. Hamarneh, "VascuSynth", 2011, available at https://github.com/midas-journal/midas-journal-794, visited on 18/04/2024

[3] Jukka Jylänki, "MathGeoLib", 2023, available at https://github.com/juj/MathGeoLib, visited on 18/04/2024
