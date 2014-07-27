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


#include <maya/MGlobal.h>
#include "FillHoles.h"
#include <stdlib.h>
#include <math.h>

// temp //
#include <maya/MIntArray.h>
#include <maya/MPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
//
#include <renderLib.h>

#define DEBUG_STEPS	0

typedef std::pair< int, int > IntPair_t;
typedef IntPair_t edge_t;

//////////////////////////////////////////////////////////////////////////


/*
================
 fillPlanarHoles
================
*/

void fillPlanarHoles( const LIST( RenderLib::Math::Vector2f )& vertices, LIST( int )& edges, const float planeNormalZSign, LIST( int )& triangles ) {
	using namespace RenderLib::Geometry;
	using namespace RenderLib::Math;

	if ( vertices.size() < 3 || edges.size() < 6 ) { 
		return;
	}
	
	LIST(int) tris; // resulting tris from the CDT

	Delaunay::AdjacencyInfo adjacency;

	Delaunay::constrainedDelaunay2D( vertices, 
									 edges,
									 &tris,
									 &adjacency );

	// The output from the CDT is a complete triangulation of the input vertices
	// which includes the provided segments as well. Nevertheless this means that
	// both outlines and holes are triangulated, and therefore we now need to
	// delete all the triangles within the holes.
	//
	// What we'll do is tracing an "infinite" (rather finite) line and count the
	// intersections with all the segments (both from the outlines and the holes)
	// if the number is even, we know he triangle is outside and we can discard it.

	Bounds2D bounds; // use the bounds of the vertices to get a sense for what's the smaller "infinite" line we can use (better numerical precision) 
	for( size_t i = 0; i < adjacency.vertices.size(); i++ ) {
		bounds.expand( adjacency.vertices[ i ] );
	}

	LIST( Vector2f )& verts = adjacency.vertices;
	for( size_t i = 0; i < adjacency.triangles.size(); i++ ) {
		Delaunay::AdjacencyInfo::Triangle_t& triangle = adjacency.triangles[ i ];
		if ( !triangle.valid ) {
			continue;
		}

		const float diagonal = bounds.max().distanceTo( bounds.min() );
		Vector2f centroid = ( verts[ triangle.vertices[ 0 ] ] + verts[ triangle.vertices[ 1 ] ] + verts[ triangle.vertices[ 2 ] ] ) / 3;

		if ( Delaunay::isTriOutside(centroid, diagonal, edges, verts, i) ) {
			 adjacency.removeTriangle( i );
		}
	}

#if _DEBUG && DEBUG_STEPS
		// dump result to maya
		{				
			MPointArray vertices;
			for( size_t i = 0; i < adjacency.vertices.size(); i++ ) {
				MPoint v( adjacency.vertices[ i ].x, adjacency.vertices[ i ].y, 0 );
				vertices.append( v );
			}

			MIntArray triangleVertices;
			for( size_t i = 0; i < adjacency.triangles.size(); i++ ) {
				const Delaunay::AdjacencyInfo::Triangle_t& triangle = adjacency.triangles[ i ];
				if ( !triangle.valid ) {
					continue;
				}
				triangleVertices.append( triangle.vertices[ 0 ] );
				triangleVertices.append( triangle.vertices[ 1 ] );
				triangleVertices.append( triangle.vertices[ 2 ] );
			}
			const int resultingTris = triangleVertices.length() / 3;
			MFnMeshData meshData;
			MObject meshObj = meshData.create();
			MFnMesh mesh( meshObj );
			MIntArray polygonCounts( resultingTris, 3 );
			MStatus status;
			mesh.create( vertices.length(), resultingTris, vertices, polygonCounts, triangleVertices, MObject::kNullObj, &status );
			MGlobal::addToModel( meshObj );			
		}
#endif

	// write the output triangles
	for( size_t i = 0; i < adjacency.triangles.size(); i++ ) {
		if ( !adjacency.triangles[ i ].valid ) {
			continue;
		}
		triangles.append( adjacency.triangles[ i ].vertices[ 0 ] );
		triangles.append( adjacency.triangles[ i ].vertices[ 1 ] );
		triangles.append( adjacency.triangles[ i ].vertices[ 2 ] );
	}

}
