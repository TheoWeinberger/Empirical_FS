# README for empirical_fs.py

Copyright 2023 T I Weinberger

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Code to empirically fit QO data to a Fermi surface defined by a geometric shape as opposed to band structure calculations

Note the code currently doesn't perform any symmetry operations/checking
and so this must be imposed by hand.

The code is currently set up for the calculation of QOs from a minimally corrugated
UTe2 fermi surface. This means that the symmetries used are for the case of Immm UTe2
and any other user must check the code carefully to remove undesired symmetry
operations.

Furthermore, current the code to extract extremal frequencies will only extract the 
maximum and minimum frequencies. This is because, in the case of the UTe2 cylinders, there
is no z-dependent warping of the cross-sectional area meaning that there is a maximum 
of two frequencies per sheet. This code must be adapted for more complex FS geometries

## Installation

This script is best run in a conda environment in which all the listed dependencies have been installed.

Typicall this should take < 10 minutes to set up.

## Usage

In order to calculate the frequencies for an empirical fermi surface, the user must define functions which describe the shape of the fermi surface while obyeing the symmetry operations of the system. 

For UTe2 this consists of two cylinders with different radial extents and warping parameters.

For each surface a list of parameters and surface functions must be defined. This can then be iterated over by the create_surface function to create a list of surfaces for plotting with pyvista. 

Note currently create_surface is set up for usage with only warped cylinders but could be adapted to take varying parameters through the use of **kwargs.

Once the surfaces have been defined. The predicted FS can be plotted within the 1st BZ. The BZ is calculated using Voronoi analysis and requires an input of the basis vectors of reciprocal space in a 3x3 matrix.

Once a satisfactory FS has been determined, the maximal and minimal frequencies for each sheet at each field angle in the c-a and c-b plane can be calculated using calculate_extremal_frequencies.

This function takes the coordinates of the surface as well as plane_res and rotations_points which determine the number of field intersection samples as well as number of angles that will be sampled as the field rotates through 90 degrees.

## Dependencies

This code has been developed and tested with the following packages on a Linux 22.04 LTS system:
    - python 3.9.10
    - numpy 1.20.3
    - matplotlib 3.7.1
    - scipy 1.7.3
    - pyvista 0.37.0
    - triangle 1.6
    - pandas 1.5.3 
    - seaborn 0.11.2
    - copy from python 3.9.10

## Instruction for UTe2

Current the code has the user parameters set to recreate the data including in our manuscript. Running the script as is should recreate the data in UTe2-band{0/1}.{c-a/c-b}.out.

The degree of warping can be changed by altering r_h/e

The radial extent of the super ellipses can be altered by changing Ra_h/e or Rb_h/e.

The extent of the sheets in the z direction is controlled by the t_h parameter.

Finally the z-resolution of the surfaces is controlled by res wheras the radial resolution is controlled by theta_res.

The expected run time for these parametes is < 10 minutes.