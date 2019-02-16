
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! BFM ECOLOGICAL MODEL COMMON TRANSPORTED VARIABLES
!
! revision log :
!
! 31.03.2008    acc     arrays for BFM TRANSPORT
!
!*****************************************************************************
!
!-------------------------------- BFM MODEL-----------------------------------
!
!---------------------------- TRASNPORTED ARRAYS------------------------------
!
!-----------------------------------------------------------------------------
!        fO2o                                            Oxgen      mmol O2/m3
!        fN1p                                        Phosphate       mmol P/m3
!        fN3n                                          Nitrate       mmol N/m3
!        fN4n                                         Ammonium       mmol N/m3
!        fO4n                                     NitrogenSink       mmol N/m3
!        fN5s                                         Silicate      mmol Si/m3
!        fN6r                            Reduction Equivalents     mmol S--/m3
!        fB1c                                 Pelagic Bacteria         mg C/m3
!        fB1n                                 Pelagic Bacteria       mmol N/m3
!        fB1p                                 Pelagic Bacteria       mmol P/m3
!        fP1c                                          Diatoms         mg C/m3
!        fP1n                                          Diatoms       mmol N/m3
!        fP1p                                          Diatoms       mmol P/m3
!        fP1l                                          Diatoms       mg Chl/m3
!        fP1s                                          Diatoms     mmmol Si/m3
!        fP2c                                      Flagellates         mg C/m3
!        fP2n                                      Flagellates       mmol N/m3
!        fP2p                                      Flagellates       mmol P/m3
!        fP2l                                      Flagellates       mg Chl/m3
!        fP3c                                PicoPhytoPlankton         mg C/m3
!        fP3n                                PicoPhytoPlankton       mmol N/m3
!        fP3p                                PicoPhytoPlankton       mmol P/m3
!        fP3l                                PicoPhytoPlankton       mg Chl/m3
!        fP4c                                  Dinoflagellates         mg C/m3
!        fP4n                                  Dinoflagellates       mmol N/m3
!        fP4p                                  Dinoflagellates       mmol P/m3
!        fP4l                                  Dinoflagellates       mg Chl/m3
!        fZ3c                      Carnivorous mesozooplankton         mg C/m3
!        fZ3n                      Carnivorous mesozooplankton       mmol N/m3
!        fZ3p                      Carnivorous mesozooplankton       mmol P/m3
!        fZ4c                       Omnivorous mesozooplankton         mg C/m3
!        fZ4n                       Omnivorous mesozooplankton       mmol N/m3
!        fZ4p                       Omnivorous mesozooplankton       mmol P/m3
!        fZ5c                                 Microzooplankton         mg C/m3
!        fZ5n                                 Microzooplankton       mmol N/m3
!        fZ5p                                 Microzooplankton       mmol P/m3
!        fZ6c             Heterotrophic nanoflagellates (HNAN)         mg C/m3
!        fZ6n             Heterotrophic nanoflagellates (HNAN)       mmol N/m3
!        fZ6p             Heterotrophic nanoflagellates (HNAN)       mmol P/m3
!        fR1c                      Labile Organic Carbon (LOC)         mg C/m3
!        fR1n                      Labile Organic Carbon (LOC)       mmol N/m3
!        fR1p                      Labile Organic Carbon (LOC)       mmol P/m3
!        fR2c                           CarboHydrates (sugars)         mg C/m3
!        fR6c                 Particulate Organic Carbon (POC)         mg C/m3
!        fR6n                 Particulate Organic Carbon (POC)       mmol N/m3
!        fR6p                 Particulate Organic Carbon (POC)       mmol P/m3
!        fR6s                 Particulate Organic Carbon (POC)     mmmol Si/m3
!        fR7c                                   Refractory DOC         mg C/m3
!-----------------------------------------------------------------------------

	real fO2o(nkndim)
        common /fO2o/fO2o
        real fN1p(nkndim)
	common /fN1p/fN1p
	real fN3n(nkndim)
        common /fN3n/fN3n
	real  fN4n(nkndim)
        common /fN4n/fN4n
	real fO4n(nkndim)
        common /fO4n/fO4n
	real fN5s(nkndim)
        common /fN5s/fN5s
	real fN6r(nkndim)
        common /fN6r/fN6r
	real fB1c(nkndim)
        common /fB1c/fB1c
	real  fB1n(nkndim)
        common /fB1n/fB1n
        real fB1p(nkndim)
	common /fB1p/fB1p
	real fP1c(nkndim)
        common /fP1c/fP1c
	real fP1n(nkndim)
        common /fP1n/fP1n
	real fP1p(nkndim)
        common /fP1p/fP1p
	real fP1l(nkndim)
        common /fP1l/fP1l
	real fP1s(nkndim)
        common /fP1s/fP1s
	real fP2c(nkndim)
        common /fP2c/fP2c
	real fP2n(nkndim)
        common /fP2n/fP2n
        real fP2p(nkndim)
	common /fP2p/fP2p
	real fP2l(nkndim)
        common /fP2l/fP2l
	real  fP3c(nkndim)
        common /fP3c/fP3c
	real fP3n(nkndim)
        common /fP3n/fP3n
	real fP3p(nkndim)
        common /fP3p/fP3p
	real fP3l(nkndim)
        common /fP3l/fP3l
	real fP4c(nkndim)
        common /fP4c/fP4c
	real fP4n(nkndim)
        common /fP4n/fP4n
	real  fP4p(nkndim)
        common /fP4p/fP4p
        real fP4l(nkndim)
	common /fP4l/fP4l
	real fZ3c(nkndim)
        common /fZ3c/fZ3c
	real fZ3n(nkndim)
        common /fZ3n/fZ3n
	real fZ3p(nkndim)
        common /fZ3p/fZ3p
	real fZ4c(nkndim)
        common /fZ4c/fZ4c
	real fZ4n(nkndim)
        common /fZ4n/fZ4n
	real fZ4p(nkndim)
        common /fZ4p/fZ4p
	real fZ5c(nkndim)
        common /fZ5c/fZ5c
	real fZ5n(nkndim)
        common /fZ5n/fZ5n
	real fZ5p(nkndim)
        common /fZ5p/fZ5p
	real fZ6c(nkndim)
        common /fZ6c/fZ6c
	real fZ6n(nkndim)
        common /fZ6n/fZ6n
	real fZ6p(nkndim)
        common /fZ6p/fZ6p
	real fR1c(nkndim)
        common /fR1c/fR1c
	real fR1n(nkndim)
        common /fR1n/fR1n
	real fR1p(nkndim)
        common /fR1p/fR1p
	real fR2c(nkndim)
        common /fR2c/fR2c
	real  fR6c(nkndim)
        common /fR6c/fR6c
        real fR6n(nkndim)
	common /fR6n/fR6n
	real fR6p(nkndim)
        common /fR6p/fR6p
	real fR6s(nkndim)
        common /fR6s/fR6s
	real fR7c(nkndim)
        common /fR7c/fR7c

!*****************************************************************************

