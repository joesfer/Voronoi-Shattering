Description:
	
	Maya plugin implementing the Space Colonization algorithm to grow
	vein-like structures around objects. 
	
	Visit http://www.joesfer.com/?p=46 for further information.

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
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

Compilation:

	- Clone the git repository into a local folder:

		mkdir <shatter_folder>
		cd <shatter_folder>
		git clone git://github.com/joesfer/Shatter.git 

	- The project depends on the RenderLib shared library, included as a submodule.
		in <shatter_folder>
		git submodules init
		git submodules update
		
		git should say Cloning into 'src/RenderLib'...
		
	- Build RenderLib using CMake:

		cd <shatter_folder>/src/RenderLib
		mkdir .build
		cd .build
		cmake ..

		Under windows: cmake will generate a Visual studio solution on .build
		Under linux: cmake will generate a GCC makefile

		Build the library. If everything went well, a new folder structure 
		<shatter_folder>/src/RenderLib/lib containing the static library 
		should have been generated.
		
	- Build the Grower plugin using CMake:

		cd <shatter_folder>
		mkdir .build
		cd .build
		cmake ..
		
		Under windows: cmake will generate a Visual studio solution on .build
		Under linux: cmake will generate a GCC makefile

		Build the plugin using visual studio or make.

		This will find the precompiled renderLib and build the .mll plugin
		under <shatter_folder>/bin

	- Load the .mll file in Maya's plugin manager.