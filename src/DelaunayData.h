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

#ifndef DelaunayData_h__
#define DelaunayData_h__

#include <maya/MPxGeometryData.h>
#include <maya/MTypeId.h>
#include <maya/MString.h>
#include <maya/MPointArray.h>
#include <maya/MPointArray.h>

#include <CoreLib.h>
#include <RenderLib.h>

/////////////////////////////////////////////////////////////////////
//
// class DelaunayData
//
/////////////////////////////////////////////////////////////////////

class ShatterNodeData : public MPxGeometryData {

public:
	//////////////////////////////////////////////////////////////////
	//
	// Overrides from MPxData
	//
	//////////////////////////////////////////////////////////////////
	ShatterNodeData();
	virtual					~ShatterNodeData();

	virtual	void			copy ( const MPxData& );

	virtual MTypeId         typeId() const;
	virtual MString         name() const;

	//////////////////////////////////////////////////////////////////
	//
	// Helper methods
	//
	//////////////////////////////////////////////////////////////////

	static void *	creator();

public:
	static const MString typeName;
	static const MTypeId id;

	CoreLib::List< RenderLib::Geometry::Delaunay::tetrahedron_t > tetrahedra;
};
#endif // DelaunayData_h__
