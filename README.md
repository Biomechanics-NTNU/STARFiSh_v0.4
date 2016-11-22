# STARFiSh: STochastic ARterial Flow Simulator

<http://www.ntnu.no/starfish>

Copyright 2012-2016 Vinzenz Eck

## License

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.



## Contributors:

* Leif Rune Hellevik
* Vinzenz Gregor Eck
* Jacob Sturdy
* Knut Petter Maraak
* Paul Roger Leinan 
* Fredrik Eikeland Fossan 
* Einar Nyberg Karlsen
* Yvan Gugler
* Yapi Donatien Achou
* Hallvard Moian Nydal 


# STARFiSh: 
is a python2, shell-based, scientific simulation program for blood flow in mammals!

## Components

* vascular network creator - vnc - is a shell-based tool to create arterial networks for the simulation tool STARFiSh.
* VascularPolynomialChaos - is a module allowing users to simulate blood flow subject to uncertainties in model parameters
* Visualisation - is a module relying on GTK to allow visualisation of the results from simulations both as 2D and 3D plots of pressures, flows, areas, etc. throughout the simulated network.

## Installation

* To get started with starfish download the code to your desired location
* There are two scripts `ubuntu_dependencies.sh` and `fedora_dependencies.sh` that should install the required dependencies. (Note these likely require administrator privileges
* Check the status of the install with `python systemCheck.py`
* You should be able to start with `python starfish.py` in the same director as you have downloaded the files
* When you first run `starfish.py` it will ask you to create a working directory. This working directory is where all output from simulations and network creation will be saved.

## Coding
For programmers going to work on this code:

The standard convention for writing Python code:
<https://www.python.org/dev/peps/pep-0008/>

Using Sphinx with some extensions, docstrings in the code 
will be used to autogenerate documentation. To generate it, 
enter the AutoDocumentation folder and run "make html". If 
you've altered modules, run "make clean" first. If you've 
altered module names, folders, or top level modules, edit 
AutoDocumentation/sources/index.rst accordingly, and then 
run "make fullupdate".

We will use a slightly modified Google standard 
for writing docstrings. (return is different)

```
def foo(input1, input2):
	"""
	describe function here
	
	Args:
		input1 (type): description of input1
		input2 (type): description of input2
	
	Returns:   ### This has different syntax depending ###
               ### on how many outputs it has.         ###

        type of soleOutput: description of soleOutput

		type of output1
			description of output1
		type of output2
			description of output2
	
	Raises:
		IOError: An error occured loading myClass.myStuff
	"""
``` 

Example can be found in UtilityLib/moduleCSV.py
in the function readBCFromCSV



