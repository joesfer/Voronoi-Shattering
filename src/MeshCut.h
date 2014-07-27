/*
	================================================================================
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
	================================================================================
*/


#pragma once

class MPointArray;
class MIntArray;
class MPlane;
class MObject;

#include <maya/MTypes.h>
#include "typedef.h"
#include <coreLib.h>

enum MeshCutType {
	MESHCUT_DISCARD_FRONT,
	MESHCUT_DISCARD_BACK,
	MESHCUT_DISCARD_NONE
};

bool meshCut( MPointArray& vertices, MIntArray& indices, LIST( bool )& triangleIsInterior, const MPlane& plane, MeshCutType cutType );
