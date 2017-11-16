
Program name: GenEx Light

Set of GenEx-based MC generators


Description:
	
	Set of Monte Carlo tools for generation phase space:
	1) Spherical Decay - TDecayFoam
	2) Decay with two leading particles - TEventMakerTFoam
	3) Proof of concept of some generators that should not be used in production - Helper


Authors: 

- Radoslaw Kycia (kycia.radoslaw@gmail.com)
Tadeusz Kościuszko Cracow University of Technology,
Warszawska 24, 31-155 Kraków, Poland.
AND
Faculty of Science, Masaryk University,
Kotlářská 2, 602 00 Brno, Czech Republic.

- Jacek Turnau

- Janusz Chwastowski (Janusz.Chwastowski@ifj.edu.pl)
Institute of Nuclear Physics Polish Academy of Sciences,
Radzikowskiego 152, 31-342 Kraków, Poland

- Rafal Staszewski (Rafal.Staszewski@ifj.edu.pl)
Institute of Nuclear Physics Polish Academy of Sciences,
Radzikowskiego 152, 31-342 Kraków, Poland

- Maciej Trzebinski (Maciej.Trzebinski@ifj.edu.pl)
Institute of Nuclear Physics Polish Academy of Sciences,
Radzikowskiego 152, 31-342 Kraków, Poland



License:

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>



Detailed description of the directory:

-- TDecayTFoam - spherically symmetric decay using TDecay class using adaptive integration by FOAM.

-- TEventMakerTFoam - adaptive simulation using two leading particles. One can select bouncing or nonbouncing solution setting isol parameter.

- Helper:

- TEventMakerTFoam2Solutions - adaptive simulation using two leading particles. Bouncing and nonbouncing solution is selected randomly.


- TDeacyTRandomSimplifiedRejectionSamling - TRandom with simplified Rejection Sampling instead of FOAM; This implementation is only for comparison with TDecayTFoam and it is messy.

- TDeacyTRandom - decay of particles using TDecay without adaptative integration; It is not meant to be used in applications, only for testing and learning - the simplest version of all generators in this directory.

- TDecayTRandomRemoveUnordered  - Spherical decay with rejecting uordered random numbers for generating intermediate masses.



Commands:

Every generator is supplied with GNU Makefile that allows to perform basic operations:

- make run - (re)compile and run generator;

- make clean - remove executables from directory;

- make cleanest - make clean + removing figures(eps) and root files(events);



