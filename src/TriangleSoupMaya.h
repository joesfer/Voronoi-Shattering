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


#ifndef TriangleSoupMaya_h__
#define TriangleSoupMaya_h__

#include <RenderLib.h>
#include <maya/MFnMesh.h>
#include <vector>

typedef struct {
	RenderLib::Math::Point3f position;
} VertexPosition;

class TriangleSoupMaya : public RenderLib::DataStructures::ITriangleSoup<VertexPosition> {
public:
	TriangleSoupMaya ( const MFnMesh& mesh, MSpace::Space space = MSpace::kObject );

	virtual size_t		numIndices() const { return indices.size(); }
	virtual const int*	getIndices() const { return &indices[ 0 ]; }

	virtual size_t					numVertices() const { return verts.size(); }
	virtual const VertexPosition*	getVertices() const { return &verts[ 0 ]; }

	RenderLib::Geometry::BoundingBox calculateBounds() const;

private:
	std::vector< int > indices;
	std::vector< VertexPosition > verts;
};

#endif // TriangleSoupMaya_h__