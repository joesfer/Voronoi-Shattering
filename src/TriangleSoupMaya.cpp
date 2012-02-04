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


#include <maya/MPointArray.h>
#include <maya/MIntArray.h>

#include "TriangleSoupMaya.h"
#include <assert.h>

TriangleSoupMaya::TriangleSoupMaya ( const MFnMesh& mesh, MSpace::Space space ) {
	verts.resize( mesh.numVertices() );
	
	MPointArray positions;
	mesh.getPoints( positions, space );
	for( int i = 0; i < mesh.numVertices(); i++ ) {
		const MPoint& p = positions[ i ];
		verts[ i ].position = RenderLib::Math::Point3f( (float)p.x, (float)p.y, (float)p.z );
	}

	MIntArray triangleCounts, triangleVertices;
	mesh.getTriangles( triangleCounts, triangleVertices );
	indices.resize( triangleVertices.length() );
	for( int i = 0; i < (int)triangleVertices.length(); i++ ) {
		const int index = triangleVertices[ i ];
		assert( index < (int)verts.size() );
		indices[ i ] = index;
	}

}

RenderLib::Geometry::BoundingBox TriangleSoupMaya::calculateBounds() const {
	RenderLib::Geometry::BoundingBox bounds;
	for( unsigned int i = 0; i < numVertices(); i++ ) {
		bounds.expand( getVertices()[ i ].position );
	}
	return bounds;
}
