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


#include "SamplerNode.h"

#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MVectorArray.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MIntArray.h>
#include <maya/MPointArray.h>
#include <maya/MVector.h>
#include <maya/MFloatVector.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MDagPath.h>
#include <maya/MFnTransform.h>
#include <maya/MMatrix.h>
#include <maya/MFnMatrixData.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MGlobal.h>

#include <vector>
#include <math.h>

#include "TriangleSoupMaya.h"

// You MUST change this to a unique value!!!  The id is a 32bit value used
// to identify this type of node in the binary file format.  
//
MTypeId     Sampler::id( 0x80099 );

// Attributes
MObject		Sampler::nSamples;
MObject     Sampler::inputMesh;        
MObject     Sampler::outputSamples;
MObject		Sampler::outputPoints;
MObject		Sampler::outputNormals;
MObject		Sampler::worldToLocal;

Sampler::Sampler( MSpace::Space space ) : sampleSpace( space ) {
	accelerator = NULL;
	geometry = NULL;
}
Sampler::~Sampler() {
	delete( accelerator );
	delete( geometry );
}

MStatus Sampler::compute( const MPlug& plug, MDataBlock& data )
//
//	Description:
//		This method computes the value of the given output plug based
//		on the values of the input attributes.
//
//	Arguments:
//		plug - the plug to compute
//		data - object that provides access to the attributes for this node
//
{
	MStatus returnStatus;
 
	// Check which output attribute we have been asked to compute.  If this 
	// node doesn't know how to compute it, we must return 
	// MS::kUnknownParameter.
	// 
	if( plug == outputSamples ) {
		// Get a handle to the input attribute that we will need for the
		// computation.  If the value is being supplied via a connection 
		// in the dependency graph, then this call will cause all upstream  
		// connections to be evaluated so that the correct value is supplied.
		// 
		MDataHandle inputMeshHandle = data.inputValue( inputMesh, &returnStatus );
		if ( returnStatus != MS::kSuccess ) return MStatus::kInvalidParameter;
		int numSamples = data.inputValue( nSamples, &returnStatus ).asInt();
		if ( returnStatus != MS::kSuccess ) return MStatus::kInvalidParameter;

		MFnMesh mesh( inputMeshHandle.asMesh() );
		{					
			MDataHandle outputHandle = data.outputValue( Sampler::outputSamples );
			MFnPointArrayData pointsHandle( outputHandle.child( Sampler::outputPoints ).data() );
			MPointArray points = pointsHandle.array();
			MFnVectorArrayData normalsHandle( outputHandle.child( Sampler::outputNormals ).data() );
			MVectorArray normals = normalsHandle.array();
						
			MFnDependencyNode thisNode( thisMObject() );
			MPlug meshPlug = thisNode.findPlug( Sampler::inputMesh );
			MPlugArray connected;
			meshPlug.connectedTo( connected, true, false );
			if ( connected.length() > 0 ) {
				MPlug meshShapePlug = connected[ 0 ];
				MFnMesh meshShapeNode( meshShapePlug.node() );
				MDagPath meshPath;
				meshShapeNode.getPath( meshPath );
				MFnTransform meshTransformNode( meshPath.transform() );
				MMatrix meshTransform = meshPath.inclusiveMatrix();

				MDataHandle world2LocalHandle = data.outputValue( Sampler::worldToLocal );
				MMatrix worldToLocal = meshTransform.inverse();
				MFnMatrixData matrixData;
				MObject matrixDataObject = matrixData.create();
				matrixData.set( worldToLocal );
				world2LocalHandle.set( matrixDataObject );	
			}
		
			sampleMesh( mesh, numSamples, points, normals );

			// Mark the destination plug as being clean.  This will prevent the
			// dependency graph from repeating this calculation until an input 
			// of this node changes.
			// 
			data.setClean(plug);
		}
	} else if ( plug == inputMesh ) {
		// rebuild accelerator structure

		MDataHandle inputMeshHandle = data.inputValue( inputMesh, &returnStatus );
		if ( returnStatus != MS::kSuccess ) return MStatus::kInvalidParameter;
		int numSamples = data.inputValue( nSamples, &returnStatus ).asInt();
		if ( returnStatus != MS::kSuccess ) return MStatus::kInvalidParameter;

		MFnMesh mesh( inputMeshHandle.asMesh() );

		buildAccelerator( mesh, sampleSpace );

	} else {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}

struct triSampling_t {
	int triangle;
	float importance;
};

int ImportanceSort( const void* a, const void* b ) {
	const float importanceA = static_cast< const triSampling_t* >(a)->importance;
	const float importanceB = static_cast< const triSampling_t* >(b)->importance;
	if ( importanceA < importanceB ) {
		return 1;
	} else if( importanceA > importanceB ) {
		return -1;
	} else {
		return 0;
	}
}

void Sampler::sampleMesh( const MFnMesh& mesh, int numSamples, MPointArray& points, MVectorArray& normals ) {
	points.clear();
	normals.clear();
	
	//int shellSamples = numSamples / 3;
	//int interiorSamples = 2 * numSamples / 3;
	//sampleShell( mesh, shellSamples, points, normals );
	//sampleInterior( mesh, interiorSamples, points, normals );
	sampleInterior( mesh, numSamples, points, normals );
}

void Sampler::sampleShell( const MFnMesh& mesh, int numSamples, MPointArray& points, MVectorArray& normals ) {
	
#if 1 // Sample at vertices
	MPointArray positions;
	mesh.getPoints( positions, sampleSpace );
	for( int i = 0; i < mesh.numVertices(); i++ ) {
		points.append( positions[ i ] );
		normals.append( MVector(0,0,0) );
	}

#else // Sample shell
	MIntArray triangleCounts, triangleVertices;
	mesh.getTriangles( triangleCounts, triangleVertices );
	unsigned int numTriangles = triangleVertices.length() / 3;
	triSampling_t* triSampling = (triSampling_t*)malloc( numTriangles * sizeof( triSampling_t ) );

	MFloatVectorArray vNormals;
#if ( MAYA_API_VERSION == 200806 )
	// maya 2008
	MVector n;
	for( int i = 0; i < mesh.numVertices(); i++ ) {
		mesh.getVertexNormal( i, n );
		vNormals.append( n );
	}
#else 
	// maya 2011
	mesh.getVertexNormals( true, vNormals );
#endif

	MPointArray verts;
	mesh.getPoints( verts, sampleSpace );
	for( unsigned int i = 0; i < numTriangles; i ++ ) {
		const int iA = triangleVertices[ 3 * i + 0 ];
		const int iB = triangleVertices[ 3 * i + 1 ];
		const int iC = triangleVertices[ 3 * i + 2 ];
		const MVector AB = verts[ iB ] - verts[ iA ];
		const MVector AC = verts[ iC ] - verts[ iA ];
		const float area = 0.5f * (float)( AB ^ AC ).length();

		const float importance = area;	

		triSampling[ i ].triangle = i;
		triSampling[ i ].importance = importance;
	}

	qsort( triSampling, numTriangles, sizeof(triSampling_t), ImportanceSort );
	while( numTriangles > 0 && triSampling[ numTriangles - 1 ].importance < 1e-5f ) { numTriangles--; }
	if ( numTriangles == 0 ) {
		free( triSampling );
		return;
	}

	// cumulative probability distribution for faces
	std::vector< int > triangleId;	
	// normalize sorted areas against the smaller triangle (so smaller importance is 1)
	const float commonDenominator = triSampling[ numTriangles - 1 ].importance;
	for( unsigned int i = 0; i < numTriangles; i++ ) {
		int area = (int)ceilf( triSampling[ i ].importance / commonDenominator );
		for( int j = 0; j < area; j++ ) {
			triangleId.push_back( triSampling[ i ].triangle );
		}
	}

	free( triSampling );

	// Sample triangles

	for( int i = 0; i < numSamples; i++ ) {
		float r = (float)rand() / RAND_MAX;
		int triId = triangleId[ (int)( r * ( (int)triangleId.size() - 1 ) ) ];

		// sample using barycentric coordinates
		float u, v;
		do {
			u = (float)rand() / RAND_MAX;
			v = (float)rand() / RAND_MAX;
		} while( u + v > 1 );

		const int iA = triangleVertices[ 3 * triId + 0 ];
		const int iB = triangleVertices[ 3 * triId + 1 ];
		const int iC = triangleVertices[ 3 * triId + 2 ];
		const MPoint A = verts[ iA ];
		const MVector AB = verts[ iB ] - A;
		const MVector AC = verts[ iC ] - A;

		points.append( A + AB * u + AC * v );
		MVector n = vNormals[ iA ] * (1.0f - u - v ) + vNormals[ iB ] * u + vNormals[ iC ] * v; 
		n.normalize();
		normals.append( n );
	}
#endif
}

void Sampler::sampleInterior( const MFnMesh& mesh, int numSamples, MPointArray& points, MVectorArray& normals ) {
	using namespace RenderLib::Math;
	using namespace RenderLib::Geometry;
	using namespace RenderLib::Raytracing;	

	if ( accelerator == NULL ) {
		buildAccelerator( mesh, sampleSpace );
	}
	
	BoundingBox bounds = accelerator->bounds();
	Vector3f extents = bounds.extents();
	int samples = 0;
	int face = 0;

	const int raysPerFace = std::min( numSamples, std::max( 10, numSamples / 6 ) );

	// ideally we'll need fewer than numSamples well-distributed points to fill the requested samples
	// although if this is not the case (a lot of empty space in the bounding box causing very few
	// successful hits?), we'll just switch to random numbers
	const int numOffsetSamples = numSamples * 100; // FIXME
	float* originOffsets = (float*)malloc( numOffsetSamples * 2 * sizeof( float ) ); 
	const int PRIME_NUMBER = 7;
	RenderLib::Math::planeHalton( originOffsets, numOffsetSamples, PRIME_NUMBER );
	int originOffset = 0;


	Point3f o;
	Vector3f u, v, w;
	float depthStep = bounds.extents()[ bounds.shortestAxis() ] / sqrtf( (float)raysPerFace );
	if ( depthStep < 1.0e-4f ) {
		MGlobal::displayError( "Sampler:: input mesh is degenerate" );
		return;
	}

	while ( samples < numSamples ) {
		face = ( face + 1 ) % 6;
		switch( face ) {
			case 0: // bottom face
				o = bounds.min();
				u = Vector3f( extents.x, 0, 0 );
				v = Vector3f( 0, 0, extents.z );
				w = Vector3f( 0, 1, 0 );				
				break;
			case 1: // top face
				o.x = bounds.min().x; o.y = bounds.max().y; o.z = bounds.min().z;
				u = Vector3f( extents.x,  0, 0 );
				v = Vector3f( 0,  0, extents.z );
				w = Vector3f( 0, -1, 0 );
				break;		
			case 2: // left face
				o = bounds.min();
				u = Vector3f( 0, extents.y, 0 );
				v = Vector3f( 0, 0, extents.z );
				w = Vector3f( 1, 0, 0 );
				break;
			case 3: // right face
				o.x = bounds.max().x; o.y = bounds.min().y; o.z = bounds.min().z;
				u = Vector3f( 0, extents.y, 0 );
				v = Vector3f( 0, 0, extents.z );
				w = Vector3f( -1, 0, 0 );
				break;
			case 4: // front face
				o = bounds.min();
				u = Vector3f( extents.x, 0, 0 );
				v = Vector3f( 0, extents.y, 0 );
				w = Vector3f( 0, 0, 1 );
				break;
			case 5: // back face
				o.x = bounds.min().x; o.y = bounds.min().y; o.z = bounds.max().z;
				u = Vector3f( extents.x, 0, 0 );
				v = Vector3f( 0, extents.y, 0 );
				w = Vector3f( 0, 0, -1 );
				break;
			default: break;
		}

		const float safetyOffset = -1.0f; // prevents the ray from starting ON a face adjacent to the bounding box
		Ray ray;
		for( int i = 0; i < raysPerFace; i++ ) {
			if ( originOffset < numOffsetSamples ) {
				ray.origin = o + w * safetyOffset + u * originOffsets[ 2 * originOffset ] + v * originOffsets[ 2 * originOffset + 1 ];
				originOffset++;
			} else {
				ray.origin = o + w * safetyOffset + u * randomFloat01() + v * randomFloat01();
			}
			ray.direction = w;
			ray.tMax = 99999;
			ray.tMin = 0;

			samples += traceRay( geometry, ray, depthStep, points, normals );
			if ( samples >= numSamples ) break;
		}
	}

	free( originOffsets );
}

int Sampler::traceRay( const TriangleSoupMaya* geometry, RenderLib::Raytracing::Ray ray, float depthStep, MPointArray& points, MVectorArray& normals ) {
	using namespace RenderLib::Math;
	RenderLib::DataStructures::TraceDesc traceDesc;
	RenderLib::DataStructures::TraceIsectDesc isect;
	
	traceDesc.doubleSided = true;
	traceDesc.startPoint = ray.origin;
	traceDesc.testOnly = false;
	Point3f endpoint = ray.at( 999999 );
	traceDesc.endPoint = endpoint;

	int samples = 0;

	bool inside = false; // this code assumes the ray starts OUTSIDE the object
	Point3f lastHit = ray.origin;
	float lastHitT = 0;
	
	while( accelerator->traceClosest( traceDesc, geometry, isect ) ) {
		inside = !inside;
		const float t = isect.t + lastHitT; // isect.t is relative to the new origin, make t relative to ray.origin
		Point3f hit = traceDesc.startPoint + ray.direction * isect.t;

		if ( !inside ) {
			const int segments = (int)ceilf( ( t - lastHitT ) / depthStep );
			if ( segments > 0 ) {
				float segmentStartT = lastHitT;
				float segmentEndT = std::min( t, lastHitT + depthStep );
				for( int i = 0; i < segments; i++ ) {
					// avoid laying the samples regularly in the segments by stratifying them
					float innerT = randomFloat01() * ( segmentEndT - segmentStartT ) + segmentStartT;

					const Point3f p = ray.at( innerT );
					points.append( p.x, p.y, p.z );
					normals.append( MVector( 0, 0, 0 ) );

					segmentStartT = segmentEndT;
					segmentEndT = std::min( t, segmentEndT + depthStep );

					samples++;
				}
			}
		}

		lastHit = hit;
		lastHitT = t;

		traceDesc.startPoint = hit + ray.direction * 1.0e-1f;
		traceDesc.endPoint = endpoint;	
	}

	return samples;
}

void Sampler::buildAccelerator( const MFnMesh& mesh, MSpace::Space space ) {
	delete( accelerator );
	delete( geometry );

	accelerator = new RenderLib::DataStructures::KdTree();
	geometry = new TriangleSoupMaya( mesh, space );
	
	if ( !accelerator->init<VertexPosition>( geometry, 24, 16 ) ) {
		MGlobal::displayError( "Sampler: Unable to build accelerator for the provided mesh" );
		delete(accelerator);
		return;
	}
}

void* Sampler::creator()
//
//	Description:
//		this method exists to give Maya a way to create new objects
//      of this type. 
//
//	Return Value:
//		a new object of this type
//
{
	return new Sampler( MSpace::kWorld );
}

MStatus Sampler::initialize()
//
//	Description:
//		This method is called to create and initialize all of the attributes
//      and attribute dependencies for this node type.  This is only called 
//		once when the node type is registered with Maya.
//
//	Return Values:
//		MS::kSuccess
//		MS::kFailure
//		
{
	// This sample creates a single input float attribute and a single
	// output float attribute.
	//
	MFnTypedAttribute	tAttr;
	MFnNumericAttribute nAttr;
	MFnCompoundAttribute cAttr;
	MStatus				stat;

	nSamples = nAttr.create( "sampleCount", "s", MFnNumericData::kInt, 1000, &stat );
	if ( !stat ) return stat;
	nAttr.setMin( 1 );
	nAttr.setWritable( true );
	nAttr.setStorable( true );

	
	inputMesh = tAttr.create( "inputMesh", "m", MFnData::kMesh, MObject::kNullObj, &stat );
	if ( !stat ) return stat;
	tAttr.setWritable( true );
	tAttr.setStorable( false );
	tAttr.setHidden( true );


	MFnPointArrayData pCreator;
	MObject pa = pCreator.create();	
	outputPoints = tAttr.create( "outPoints", "osp", MFnData::kPointArray, pa, &stat );
	if ( !stat ) return stat;
	// Attribute is read-only because it is an output attribute
	tAttr.setWritable(false);
	// Attribute will not be written to files when this type of node is stored
	tAttr.setStorable(false);
	tAttr.setHidden( true );

	MFnVectorArrayData nCreator;
	MObject va = nCreator.create();	
	outputNormals = tAttr.create( "outNormals", "osn", MFnData::kVectorArray, va, &stat );
	if ( !stat ) return stat;
	// Attribute is read-only because it is an output attribute
	tAttr.setWritable(false);
	// Attribute will not be written to files when this type of node is stored
	tAttr.setStorable(false);
	//tAttr.setHidden( true );

	outputSamples = cAttr.create( "outSamples", "os", &stat );
	if ( !stat ) return stat;

	MFnCompoundAttribute os(outputSamples );
	os.addChild( outputPoints );
	os.addChild( outputNormals );
	os.setHidden( true );

	worldToLocal = tAttr.create( "worldToLocal", "wtl", MFnData::kMatrix, MObject::kNullObj, &stat );
	if ( !stat ) return stat;
	tAttr.setWritable( false );
	tAttr.setStorable( false );
	tAttr.setHidden( true );

	// Add the attributes we have created to the node
	//
	stat = addAttribute( nSamples );
	if (!stat) { stat.perror("addAttribute"); return stat;}
	stat = addAttribute( inputMesh );
	if (!stat) { stat.perror("addAttribute"); return stat;}
	stat = addAttribute( outputSamples );
	if (!stat) { stat.perror("addAttribute"); return stat;}
	stat = addAttribute( worldToLocal );
	if (!stat) { stat.perror("addAttribute"); return stat;}

	// Set up a dependency between the input and the output.  This will cause
	// the output to be marked dirty when the input changes.  The output will
	// then be recomputed the next time the value of the output is requested.
	//
	stat = attributeAffects( nSamples, outputSamples );
	if (!stat) { stat.perror("attributeAffects"); return stat;}
	stat = attributeAffects( inputMesh, outputSamples );
	if (!stat) { stat.perror("attributeAffects"); return stat;}
	stat = attributeAffects( inputMesh, worldToLocal );
	if (!stat) { stat.perror("attributeAffects"); return stat;}

	return MS::kSuccess;

}

