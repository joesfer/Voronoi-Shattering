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

#include "DelaunayData.h"

const MTypeId ShatterNodeData::id( 0x80777 );
const MString ShatterNodeData::typeName( "ShatterData" );

//////////////////////////////////////////////////////////////////////////
// DelaunayData::DelaunayData()
//////////////////////////////////////////////////////////////////////////

ShatterNodeData::ShatterNodeData() {
}

//////////////////////////////////////////////////////////////////////////
// DelaunayData::~DelaunayData()
//////////////////////////////////////////////////////////////////////////

ShatterNodeData::~ShatterNodeData() {
}

//////////////////////////////////////////////////////////////////////////
// DelaunayData::copy (override)
//////////////////////////////////////////////////////////////////////////

void ShatterNodeData::copy ( const MPxData& other ) {
	if ( &other != this ) {
		const ShatterNodeData& _other = (const ShatterNodeData &)other;
		tetrahedra = _other.tetrahedra;
	}
}

//////////////////////////////////////////////////////////////////////////
// DelaunayData::typeId (override)
//
//	Binary tag used to identify this kind of data
//////////////////////////////////////////////////////////////////////////

MTypeId ShatterNodeData::typeId() const {
	return ShatterNodeData::id;
}

//////////////////////////////////////////////////////////////////////////
// DelaunayData::typeId (override)
//
//	String name used to identify this kind of data
//////////////////////////////////////////////////////////////////////////
MString ShatterNodeData::name() const {
	return ShatterNodeData::typeName;
}

//////////////////////////////////////////////////////////////////////////
// DelaunayData::creator
//
//	This method exists to give Maya a way to create new objects
//	of this type. 
//////////////////////////////////////////////////////////////////////////
void * ShatterNodeData::creator() {
	return new ShatterNodeData;
}
