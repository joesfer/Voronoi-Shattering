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

#ifndef DelaunayNode_h__
#define DelaunayNode_h__

#include <maya/MPxNode.h>
#include <maya/MTypeId.h> 
#include <maya/MFnMesh.h>
#include <maya/MPointArray.h>
#include <maya/MVectorArray.h>

#include <coreLib.h>
#include <renderLib.h>

class ShatterNode : public MPxNode
{
public:
	ShatterNode();
	virtual				~ShatterNode(); 

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static  void*		creator();
	static  MStatus		initialize();

public:

	// There needs to be a MObject handle declared for each attribute that
	// the node will have.  These handles are needed for getting and setting
	// the values later.
	//
	static	MObject		inputSamples;	// input vector array
	static	MObject		inputPoints;	
	static	MObject		inputNormals;
	static	MObject		inputMesh;
	static  MObject		shrink; // units to shrink the voronoi cells by to avoid perfect cell overlapping

	static	MObject		outData;		// DelaunayData


	// The typeid is a unique 32bit identifier that describes this node.
	// It is used to save and retrieve nodes of this type from the binary
	// file format.  If it is not unique, it will cause file IO problems.
	//
	static const MTypeId	id;
	static const MString	typeName;

private:
	void tetrahedralize( CoreLib::List< RenderLib::Math::Point3d > &pointArray, const MVectorArray& normals, CoreLib::List< RenderLib::Geometry::Delaunay::tetrahedron_t >& tetrahedra );
	void voronoiCutting( CoreLib::List< RenderLib::Math::Point3d > &pointArray, CoreLib::List< RenderLib::Geometry::Delaunay::tetrahedron_t > &tetrahedra, const float shrinkUnits, MFnMesh& mesh );
};

#endif // DelaunayNode_h__
