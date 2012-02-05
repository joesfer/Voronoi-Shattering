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

#include "DelaunayNode.h"
#include "DelaunayData.h"
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MFnPluginData.h>
#include <maya/MFnMeshData.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MFloatVectorArray.h>
#include "MPlane.h"

#include <RenderLib.h>
#include "MeshCut.h"
#include "typedef.h"

#include <ostream>
#include <iostream>
#include <sstream>

#define DEBUG_STEPS		0
#define DUMP_DT			0
#define DUMP_VD_CELLS	0
#define DUMP_INTERMEDIATE_SPLITS 0
int globalPrefix = 0;

//////////////////////////////////////////////////////////////////////
//
// Error checking
//
//    MCHECKERROR       - check the status and print the given error message
//    MCHECKERRORNORET  - same as above but does not return
//
//////////////////////////////////////////////////////////////////////

#define MCHECKERROR(STAT,MSG)       \
	if ( MS::kSuccess != STAT ) {   \
	cerr << MSG << endl;        \
	return MS::kFailure;    \
		}

#define MCHECKERRORNORET(STAT,MSG)  \
	if ( MS::kSuccess != STAT ) {   \
	cerr << MSG << endl;        \
		}

// You MUST change this to a unique value!!!  The typeId is a 32bit value used
// to identify this type of node in the binary file format.  
//
const MTypeId   ShatterNode::id( 0x80097 );
const MString	ShatterNode::typeName( "DelaunayNode" );

// Attributes
MObject		ShatterNode::inputSamples;
MObject		ShatterNode::inputPoints;
MObject		ShatterNode::inputNormals;
MObject		ShatterNode::inputMesh;
MObject		ShatterNode::shrink;
MObject		ShatterNode::outData;

ShatterNode::ShatterNode() {}
ShatterNode::~ShatterNode() {}

MStatus ShatterNode::compute( const MPlug& plug, MDataBlock& data )
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
	MStatus stat;
	if ( plug == outData ) {

		MDataHandle inputPointsHandle = data.inputValue( inputSamples, &stat );

		float shrinkUnits = data.inputValue( shrink, &stat ).asFloat();
		if ( stat != MS::kSuccess ) return MStatus::kInvalidParameter;

		MObject pointArrayObj = inputPointsHandle.child( inputPoints ).data();
		MFnPointArrayData pointVecData;
		pointVecData.setObject(pointArrayObj);
		MPointArray pointVec = pointVecData.array();

		MObject normalArrayObj = inputPointsHandle.child( inputNormals ).data();
		MFnVectorArrayData normalVecData;
		normalVecData.setObject(normalArrayObj);
		MVectorArray normalVec = normalVecData.array();

		MFnPluginData fnDataCreator;
		MTypeId tmpid( ShatterNodeData::id );
		ShatterNodeData * newData = NULL;

		MDataHandle outHandle = data.outputValue( outData );	
		newData = (ShatterNodeData*)outHandle.asPluginData();

		if ( newData == NULL ) {
			// Create some output data
			fnDataCreator.create( tmpid, &stat );
			MCHECKERROR( stat, "compute : error creating DelaunayData")
				newData = (ShatterNodeData*)fnDataCreator.data( &stat );
			MCHECKERROR( stat, "compute : error gettin at proxy DelaunayData object")
		}

		MDataHandle inputMeshHandle = data.inputValue( inputMesh, &stat );
		if ( stat != MS::kSuccess ) return MStatus::kInvalidParameter;
		MFnMesh mesh( inputMeshHandle.asMesh() );		

		// compute the output values			

		newData->tetrahedra.resize( 0, false );
		
		if ( pointVec.length() > 0 ) { 
			CoreLib::List< RenderLib::Math::Point3d > pointArray;
			pointArray.resize( (int)pointVec.length() );
			for( unsigned int i = 0; i < pointVec.length(); i++ ) {
				pointArray[ i ] = RenderLib::Math::Point3d(	pointVec[ i ].x, 
														pointVec[ i ].y,
														pointVec[ i ].z );
			}
		
			CoreLib::Memory::StaticMemoryPoolBase::init( 16 * 1024 * 1024 );

			tetrahedralize( pointArray, normalVec, newData->tetrahedra );
			voronoiCutting( pointArray, newData->tetrahedra, shrinkUnits, mesh );

			CoreLib::Memory::StaticMemoryPoolBase::destroy();
		}

		// Assign the new data to the outputSurface handle

		if ( newData != outHandle.asPluginData() ) {
			outHandle.set( newData );
		}

		data.setClean(plug);
		return MS::kSuccess;

	}

	return MS::kInvalidParameter;
}


void* ShatterNode::creator()
//
//	Description:
//		this method exists to give Maya a way to create new objects
//      of this type. 
//
//	Return Value:
//		a new object of this type
//
{
	return new ShatterNode();
}

MStatus ShatterNode::initialize()
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
	MFnNumericAttribute nFn;
	MFnTypedAttribute	typedFn;	
	MFnCompoundAttribute cFn;
	MStatus				stat;

	inputPoints = typedFn.create( "samplesPoints", "sp", MFnData::kPointArray );
	typedFn.setStorable( false );
	typedFn.setWritable( true );

	inputNormals = typedFn.create( "samplesNormals", "sn", MFnData::kVectorArray );
	typedFn.setStorable( false );
	typedFn.setWritable( true );

	inputSamples = cFn.create( "samples", "s" );
	cFn.setWritable( true );
	cFn.addChild( inputPoints );
	cFn.addChild( inputNormals );
	cFn.setHidden( true );

	inputMesh = typedFn.create( "inputMesh", "m", MFnData::kMesh, MObject::kNullObj, &stat );
	if ( !stat ) return stat;
	typedFn.setWritable( true );
	typedFn.setStorable( false );
	typedFn.setHidden( true );

	shrink = nFn.create( "shrinking", "sh", MFnNumericData::kFloat, 1 );
	if ( !stat ) return stat;
	typedFn.setWritable( true );
	typedFn.setStorable( true );
	typedFn.setHidden( false );


	outData = typedFn.create( "output", "out", ShatterNodeData::id );
	typedFn.setWritable( false );
	typedFn.setStorable(false);
	typedFn.setHidden( true );

	// Add the attributes we have created to the node
	//
	stat = addAttribute( inputSamples );
	if (!stat) { stat.perror("addAttribute"); return stat;}
	stat = addAttribute( inputMesh );
	if (!stat) { stat.perror("addAttribute"); return stat;}
	stat = addAttribute( shrink );
	if (!stat) { stat.perror("addAttribute"); return stat;}
	stat = addAttribute( outData );
	if (!stat) { stat.perror("addAttribute"); return stat;}

	attributeAffects( inputSamples, outData );
	attributeAffects( inputMesh, outData );
	attributeAffects( shrink, outData );

	return MS::kSuccess;

}

void ShatterNode::tetrahedralize( CoreLib::List< RenderLib::Math::Point3d >& pointArray, const MVectorArray& normals, CoreLib::List<  RenderLib::Geometry::Delaunay::tetrahedron_t >& tetrahedra ) {
	
	using namespace RenderLib::Geometry::Delaunay;
	
	Delaunay3D delaunay;

	if ( !delaunay.tetrahedralize( pointArray, tetrahedra ) ) {
		return;
	}
#if _DEBUG

#if 0 // dump sampling points
	for( int i = 0; i < pointArray.Num(); i++ ) {
		MString cmd = CoreLib::va( "spaceLocator -p %f %f %f;", pointArray[ i ].x, pointArray[ i ].y, pointArray[ i ].z );
		MGlobal::executeCommand( cmd );
	}
#endif
#if DUMP_DT 
	{
		MFnMeshData md;
		MPointArray points( 4 );
		MIntArray polygonCounts( 4, 3 );
		MIntArray indices( 4 * 3 );
		indices[0] = 0; indices[1] = 2; indices[2] = 1;
		indices[3] = 0; indices[4] = 1; indices[5] = 3;
		indices[6] = 1; indices[7] = 2; indices[8] = 3;
		indices[9] = 2; indices[10] = 0; indices[11] = 3;

		for( size_t i = 0; i < tetrahedra.size(); i++ ) {
			const RenderLib::Geometry::Delaunay::tetrahedron_t& t = tetrahedra[ i ];
			if ( !t.isValid() ) {
				continue;
			}
			const RenderLib::Math::Point3d& v0 = pointArray[ t.v[ 0 ] ];
			const RenderLib::Math::Point3d& v1 = pointArray[ t.v[ 1 ] ];
			const RenderLib::Math::Point3d& v2 = pointArray[ t.v[ 2 ] ];
			const RenderLib::Math::Point3d& v3 = pointArray[ t.v[ 3 ] ];

			MObject meshObject = md.create();
			MFnMesh mesh( meshObject );

			points[ 0 ] = MPoint( v0.x, v0.y, v0.z );
			points[ 1 ] = MPoint( v1.x, v1.y, v1.z );
			points[ 2 ] = MPoint( v2.x, v2.y, v2.z );
			points[ 3 ] = MPoint( v3.x, v3.y, v3.z );

			mesh.create( 4, 4, points, polygonCounts, indices );	
			std::stringstream name; 
			name << globalPrefix << "_TD_" << i;
			mesh.setName( name.str().c_str() );
			MGlobal::addToModel( meshObject );

		}
	}
#endif
#endif
}

#include <random>
float Random01() {
	return (float)rand() / RAND_MAX;
}

void ShatterNode::voronoiCutting( CoreLib::List< RenderLib::Math::Point3d > &pointArray, CoreLib::List< RenderLib::Geometry::Delaunay::tetrahedron_t > &tetrahedra, const float shrinkUnits, MFnMesh& mesh ) {
	using namespace RenderLib::Math;

	// Extract Voronoi Diagram from Delaunay Triangulation
	// ---------------------------------------------------
	
	CoreLib::List< CoreLib::List< int > > pointTetrahedraAdjacencyInfo; // tetrahedra adjacent to each source point
	pointTetrahedraAdjacencyInfo.resize( pointArray.size() ); 

	for( size_t i = 0; i < tetrahedra.size(); i++ ) {
		const RenderLib::Geometry::Delaunay::tetrahedron_t& t = tetrahedra[ i ];
		if ( !t.isValid() ) {
			continue;
		}
		for( int j = 0; j < 4; j++ ) {
			pointTetrahedraAdjacencyInfo[ t.v[ j ] ].append( i );
		}
	}

	CoreLib::List< int > edgesToP( 32 );
	int a,b,c;
	MFnMeshData meshData;

#if _DEBUG
#if DUMP_VD_CELLS
	CoreLib::List< Vector3d > planeNormals;
	CoreLib::List< double > planeDists;
#endif
#endif

	for( size_t i = 0; i < pointArray.size(); i++ ) {

		CoreLib::Memory::StaticMemoryPoolBase::clearMemory();

#if _DEBUG
#if DUMP_VD_CELLS
		planeNormals.resize( 0, false );
		planeDists.resize( 0, false );
#endif
#endif
		// Clone the source mesh; we'll trim it with the VD faces
		MPointArray vertices;	
		MIntArray triangleCounts, triangleVertices;
		mesh.getPoints( vertices, MSpace::kWorld );
		mesh.getTriangles( triangleCounts, triangleVertices );
		
		LIST( bool ) triangleIsInterior;
		triangleIsInterior.resize( triangleVertices.length() / 3 );
		triangleIsInterior.setGranularity( triangleIsInterior.size() );
		memset( &triangleIsInterior[ 0 ], false, triangleIsInterior.size() * sizeof( bool ) );

		int cuts = 0;

		// 1. extract all DT edges arriving to the point (we just need to store the other end)
		edgesToP.resize( 0, false );
		const CoreLib::List< int >& adjacentTetrahedra = pointTetrahedraAdjacencyInfo[ i ];

		for( size_t j = 0; j < adjacentTetrahedra.size(); j++ ) {
			const RenderLib::Geometry::Delaunay::tetrahedron_t& adjacentT = tetrahedra[ adjacentTetrahedra[ j ] ];
			assert( adjacentT.isValid() );

			for( int f = 0; f < 4; f++ ) {
				adjacentT.getFaceVertices( f, a, b, c );
				// if face contains our point, extract the two adjacent edges
				if ( a == i ) {
					edgesToP.addUnique( b );
					edgesToP.addUnique( c );
				} else if ( b == i ) {
					edgesToP.addUnique( a );
					edgesToP.addUnique( c );
				} else if ( c == i ) {
					edgesToP.addUnique( a );
					edgesToP.addUnique( b );
				}
			}
		}
		
		// 2. for each adjacent edge, build the dual VD face by "rotating"
		//	  around the edge, collecting the shared DT tetrahedra.

		for( size_t j = 0; j < edgesToP.size(); j++ ) {
#if _DEBUG
			//if ( j > 5 ) continue;
#if DUMP_VD_CELLS && DEBUG_STEPS && 0
			{
				MString cmd = CoreLib::varArgsStr<1024>( "$c = `curve -d 1 -p %f %f %f -p %f %f %f -k 0 -k 1`;\
														  rename $c P%dE%d;", 
														  pointArray[ i ].x, pointArray[ i ].y, pointArray[ i ].z,
														  pointArray[ edgesToP[ j ] ].x, pointArray[ edgesToP[ j ] ].y, pointArray[ edgesToP[ j ] ].z, 
														  i, j );
				MGlobal::executeCommand( cmd );
			}
#endif
#endif
			int candidates = 0;
			double smallerSphereRadius = DBL_MAX;
			int betterCandidate;
			RenderLib::Math::Point3d facePoint;	

			// Optimization: we don't need to calculate the whole VD face
			// since we're cutting the geometry with a plane: We know the
			// plane normal will be parallel to the edge, and just need
			// ONE of the face points to calculate the plane distance

			// nevertheless, instead of just picking the first adjacent tetrahedron
			// we'll try to choose "best" (less-distorted) one for that point
			// since flat tetrahedra tend to generate huge spheres which situate
			// the face far-off from the actual VD cell

			for( size_t t = 0; t < adjacentTetrahedra.size(); t++ ) {
				// Pick the less-distorted adjacent tetrahedron 
				// By less-distorted we're looking for a tetrahedron which doesn't tend to be flat
				// that is, the points are not coplanar and the resulting sphere is smaller
				const RenderLib::Geometry::Delaunay::tetrahedron_t& adjacentT = tetrahedra[ adjacentTetrahedra[ t ] ];
				if( adjacentT.containsVertex( edgesToP[ j ] ) ) {
					candidates++;
					RenderLib::Math::Point3d center;
					double radius;

					RenderLib::Geometry::SphereFrom4Points(	pointArray[ adjacentT.v[0] ], 
															pointArray[ adjacentT.v[1] ], 
															pointArray[ adjacentT.v[2] ], 
															pointArray[ adjacentT.v[3] ], 
															center, 
															radius );

					if ( radius < smallerSphereRadius) {
						smallerSphereRadius = radius;
						facePoint = center;
						betterCandidate = adjacentTetrahedra[ t ];
					}
				}
			}
			
			if ( candidates == 0 ) {
				continue;
			}
#if _DEBUG
#if DUMP_VD_CELLS && 0		
			{				
				MString cmd = CoreLib::va( "$s = `sphere -r %f -pivot %f %f %f`;\
											rename $s P%dE%dT%d_sphere;", smallerSphereRadius, facePoint.x, facePoint.y, facePoint.z, i, j, betterCandidate );
				MGlobal::executeCommand( cmd );
				cmd = CoreLib::va( "$sl = `spaceLocator -p %f %f %f`; rename $sl P%dE%dT%d_center;", facePoint.x, facePoint.y, facePoint.z, i, j, betterCandidate );
				MGlobal::executeCommand( cmd );
			}
#endif
#endif
	
			Vector3d planeNormal = pointArray[ edgesToP[ j ] ] - pointArray[ i ];
			planeNormal.normalize();
			double planeDist = Vector3d::dot( planeNormal, facePoint.fromOrigin() );

#if 0
			MString cmd = CoreLib::va("$p = `polyPlane -ax %f %f %f -w 100 -h 100 -sx 1 -sy 1`;\
									  move -absolute %f %f %f $p;",
									  planeNormal.x, planeNormal.y, planeNormal.z,
									  planeNormal.x * planeDist, planeNormal.y * planeDist, planeNormal.z * planeDist );
			MGlobal::executeCommand( cmd );				
#endif
#if 0
			if ( j == 0 ) {
				srand( 1234 );
				planeNormal.x = Random01(); planeNormal.y = Random01(); planeNormal.z = Random01();
				planeDist = 0;
			} else if ( j < 10 ) {
				planeNormal.x = Random01(); planeNormal.y = Random01(); planeNormal.z = Random01();
				planeDist = 0;
			} else {
				break;
			}
			planeNormal.Normalize();
#endif
#if _DEBUG
#if DUMP_VD_CELLS
			planeNormals.Append( planeNormal );
			planeDists.Append( planeDist );
#endif 
#endif
			MPlane plane;
			plane.setPlane( planeNormal.x, planeNormal.y, planeNormal.z, -planeDist + shrinkUnits );
		
			if ( !meshCut( vertices, triangleVertices, triangleIsInterior, plane, MESHCUT_DISCARD_FRONT ) ) {
				cuts = 0;
				break;
			} else {
				cuts++;
			}
#if DUMP_INTERMEDIATE_SPLITS
			const int resultingTris = triangleVertices.length() / 3;
			if ( resultingTris > 0 ) {
				MObject res = meshData.create();
				MIntArray polygonCounts( resultingTris, 3 );
				MFnMesh resMesh( res );
				resMesh.create( vertices.length(), resultingTris, vertices, polygonCounts, triangleVertices );
				MGlobal::addToModel( res );
			}
#endif

		}
#if _DEBUG
#if DUMP_VD_CELLS && 1
		for( int i = 0; i < planeNormals.Num(); i++ ) {
			const Vector3d& planeNormal = planeNormals[ i ];
			double planeDist = planeDists[ i ];
			static int vdc = 0;
			vdc++;
			MString cmd = CoreLib::va("$p = `polyPlane -ax %f %f %f -w 1000 -h 1000 -sx 1 -sy 1`;\
									  move -absolute %f %f %f $p;\
									  select $p;\
									  rename VDP_%d;",
									  planeNormal.x, planeNormal.y, planeNormal.z,
									  planeNormal.x * planeDist, planeNormal.y * planeDist, planeNormal.z * planeDist,
									  vdc );
			MGlobal::executeCommand( cmd );				
			for( int j = 0; j < planeNormals.Num(); j++ ) {
				if ( j == i ) {
					continue;
				}
				const Vector3d& cutPlaneNormal = -planeNormals[ j ];
				const double Xrot = -asin( cutPlaneNormal.y ) * 180.0 / 3.1415926535898;
				const double Yrot = atan2( cutPlaneNormal.x, cutPlaneNormal.z ) * 180.0 / 3.1415926535898;
				double cutPlaneDist = -planeDists[ j ];
				const Vector3d cutPlaneCenter = cutPlaneNormal * cutPlaneDist;

				MString cmd = CoreLib::va( "polyCut -df 1 -pc %f %f %f -rx %f -ry %f;",
											cutPlaneCenter.x, cutPlaneCenter.y, cutPlaneCenter.z,
											Xrot, Yrot );
				MGlobal::executeCommand( cmd );				
			}
		}
#endif

#endif
		// After all the trimming, is there any piece of the mesh left?
		// if so, create a model for it

		const int resultingTris = triangleVertices.length() / 3;
		if ( cuts > 0 && resultingTris > 0 ) {
			MObject res = meshData.create();
			MIntArray polygonCounts( resultingTris, 3 );
			MFnMesh resMesh( res );			
			resMesh.create( vertices.length(), resultingTris, vertices, polygonCounts, triangleVertices );		
			
			// tweak interior face normals to make them planar
			MIntArray faces;
			for( size_t i = 0; i < triangleIsInterior.size(); i++ ) {
				if ( triangleIsInterior[ i ] ) {
					faces.append( i );
				}
			}
			resMesh.extractFaces( faces, NULL );

			//MVectorArray normals;
			for( size_t i = 0; i < triangleIsInterior.size(); i++ ) {
				if ( triangleIsInterior[ i ] ) {
					MVector normal;
					resMesh.getPolygonNormal( i, normal );
					
					int triangleVertices[3];
					resMesh.getPolygonTriangleVertices( i, 0, triangleVertices );

					for( int j = 0; j < 3; j++ ) {
						//resMesh.setFaceColor( MColor(  1, 0, 0 ), i );
						resMesh.setFaceVertexNormal( normal, i, triangleVertices[ j ] );
					}
				}
			}

			MGlobal::addToModel( res );
		}
	}

}