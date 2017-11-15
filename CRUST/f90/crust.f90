! Copyright (C) 2013, Thomas M. Melvin and Keith R. Briffa, see the GNU 
! General Public License. Climatic Research Unit, School of Environmental
! Sciences, University of East Anglia, Norwich, NR4 7TJ, U.K.
!
! CRUST - A tree-ring standardisations program to generate the chronologies
! and figures see:

! CRUST: Software for the implementation of Regional Chronology 
! Standardisation: Part 1. Signal-Free RCS Melvin, T. M. and Briffa, K. R.
! Dendrochronologia 32 (2014) 7-20 doi: 
! http://dx.doi.org/10.1016/j.dendro.2013.06.002
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version, see <http://www.gnu.org/licenses/>.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

      PROGRAM CRUST
      USE crustprocs
      IMPLICIT NONE
      CALL crust_setup()     ! CRUST buttons and graphics
      CALL start_menu()      ! Main program
      CALL DISFIN()
      STOP
      END PROGRAM CRUST
