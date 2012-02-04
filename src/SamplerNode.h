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


#ifndef _SamplerNode
#define _SamplerNode

#include <maya/MPxNode.h>
#include <maya/MTypeId.h> 
#include <maya/MFnMesh.h>
#include <maya/MPointArray.h>
#include <maya/MVectorArray.h>
#include <renderlib.h>
 
class TriangleSoupMaya;

class Sampler : public MPxNode
{
public:
						 Sampler( MSpace::Space space );
	virtual				~Sampler(); 

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static  void*		creator();
	static  MStatus		initialize();

public:

	// There needs to be a MObject handle declared for each attribute that
	// the node will have.  These handles are needed for getting and setting
	// the values later.
	//
	static  MObject		nSamples;		// number of desired samples
	static	MObject		inputMesh;		// input mesh to sample
	static  MObject		outputSamples;	// output samples array
	static	MObject		outputPoints;	// child of outputSamples
	static	MObject		outputNormals;	// child of outputSamples
	static	MObject		worldToLocal;	// output copy of mesh transform


	// The typeid is a unique 32bit identifier that describes this node.
	// It is used to save and retrieve nodes of this type from the binary
	// file format.  If it is not unique, it will cause file IO problems.
	//
	static	MTypeId		id;

private:

	void sampleMesh( const MFnMesh& mesh, int numSamples, MPointArray& points, MVectorArray& normals );

	void sampleShell( const MFnMesh& mesh, int numSamples, MPointArray& points, MVectorArray& normals );
	
	void sampleInterior( const MFnMesh& mesh, int numSamples, MPointArray& points, MVectorArray& normals );
	void buildAccelerator( const MFnMesh& mesh, MSpace::Space space );
	int traceRay( const TriangleSoupMaya* geometry, const RenderLib::Raytracing::Ray ray, float depthStep, MPointArray& points, MVectorArray& normals );

private:
	TriangleSoupMaya*					geometry;
	RenderLib::DataStructures::KdTree*	accelerator;
	MSpace::Space						sampleSpace;
};

#endif
