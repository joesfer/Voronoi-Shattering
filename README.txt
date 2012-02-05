Description:
	
	Maya plugin implementing geometry fracturing based on 
	Voronoi Diagrams (Voronoi Shattering).
	
	Visit http://www.joesfer.com/?p=60 for further information.

License:

	This software is released under the LGPL-3.0 license: http://www.opensource.org/licenses/lgpl-3.0.html	

	Copyright (c) 2012, Jose Esteve. http://www.joesfer.com
	
	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 3.0 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  
	USA

Compilation:

	- Clone the git repository into a local folder:

		mkdir <shatter_folder>
		cd <shatter_folder>
		git clone git://github.com/joesfer/Shatter.git 

	- The project depends on the CoreLib and RenderLib shared libraries, 
	  included as submodules.
	  
		in <shatter_folder>
		git submodule init
		git submodule update
		
		it should create <shatter_folder>/RenderLib and 
		<shatter_folder>/CoreLib otherwise, add them manually:
		
		git submodule add git://github.com/joesfer/CoreLib.git CoreLib
		git submodule add git://github.com/joesfer/RenderLib.git RenderLib		
		
	- Build the Maya plugin using CMake:

		cd <shatter_folder>
		mkdir .build
		cd .build
		cmake ..
		
		On Windows: cmake will generate a Visual studio solution on .build
		On Linux: cmake will generate a GCC makefile

		Build the plugin using visual studio or make.

		Note the solution will contain 3 projects: CoreLib, RenderLib 
		and Shatter, which need to be built in that order. If the process 
		succeeds, the resulting .mll plugin file will be located under 
		<shatter_folder>/bin

	- Load the .mll file in Maya's plugin manager.
	
	- Run the provided MEL script for an example on how to use the plugin