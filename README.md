# CRUST
CRUST: Climatic Research Unit Standardisation of Tree-ring data

Copyright &copy; 2013, Thomas M. Melvin and Keith R. Briffa, see the GNU General Public License.
Climatic Research Unit, School of Environmental Sciences, University of East Anglia, Norwich, NR4 7TJ, U.K.

The program CRUST program is designed to perform the specific processing necessary to construct tree-ring chronologies and undertake  various analyses. Executable versions of the program are provided for LINUX and Windows. The program uses specific files and needs to be run within the specific directory structure provided here. The "data" directory contains sample tree-ring measurement files. The "f90" directory contains the source code and makefile.

The program can be run by clicking on the appropriate executable. If errors occur (e.g. program disappears from the screen) then run the program from the "Command" line or "DOS", both of which should provide some error diagnostics. The makefile, amended to select for machine type from LIN or WIN can be used to compile the program. For LINUX, gfortran and Dislin need to have been installed. For Windows, MGW32 (GCC) with gfortran and Dislin need to have been installed. The program has also been compiled on MacOS (thanks to Paul Krusic) but getting the correct components installed for the appropriate version of software has proved awkward.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version, see <http://www.gnu.org/licenses/>.

We are not able to offer any technical support for this program. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Acknowledgements:
-----------------

We wish to thank Ed Cook and Paul Krusic for their contribution of computer code, algorithms and helpful advice.

The DISLIN graphics software, by Helmut Michels, was obtained from the Max Planck Institute for Solar System Research, Lindau at <http://www.mps.mpg.de/dislin/examples.html>

Attribution:
------------
Whenever you publish results obtained using this software you must acknowledge our work by including the following references:

* Melvin T. M. and Briffa K. R. (2014) CRUST: Software for the implementation of Regional Chronology Standardisation: Part 1, Signal-Free RCS. Dendrochronologia, 32, 7-20, doi: <http://dx.doi.org/10.1016/j.dendro.2013.06.002>

* Melvin T. M. and Briffa K. R. (2014) CRUST: Software for the implementation of Regional Chronology Standardisation: Part 2, Further RCS options and recommendations. Dendrochronologia, 32, 343-356, doi: <http://dx.doi.org/10.1016/j.dendro.2014.07.008>
