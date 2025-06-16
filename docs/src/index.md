
# ThermovisorImages.jl

## General description

**ThermovisorImages.jl** is designed to process static thermal images stored as matrices in CSV files or as image files. It treats each matrix element as a temperature value. ThermovisorImages.jl provides functions to calculate temperature distributions and perform statistical analyses of temperatures within Regions of Interest (ROIs), such as circles, squares, rectangles, or along lines. ROI objects can be fitted to image patterns (regions that stand out from the background). It is also possible to evaluate statistics across multiple ROIs, including distributions of side length, area, and perimeter.

**ThermovisorImages.jl** also provides functions to recalculate the temperature distribution of the entire image (or its part within the ROI or labeled pattern), taking into account the emissivity of the surface and the spectral range of the infrared camera.


## Contact

Contacts [GitHub repository](https://github.com/Manarom/ThermovisorImages.jl).

## License

Copyright (c) 2025 Roman Mironov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.

