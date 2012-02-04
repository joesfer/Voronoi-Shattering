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


#include "MeshCut.h"
#include "FillHoles.h"

#include <maya/MFnMesh.h>
#include <maya/MBoundingBox.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include "MPlane.h"
#include <maya/MGlobal.h>

#include <CoreLib.h>

#define DEBUG_STEPS	0


#define CUT_EPSILON 1e-4
#define ON_PLANE_EPSILON 2 * CUT_EPSILON
#define MERGE_EPSILON 1e-2
#define COLLINEAR_EPSILON 1e-5

typedef std::pair< int, int > IntPair_t;

struct edge_t {
	int v[2]; // vertices
	int f[2]; // faces: f[0] has the edge v[0]->v[1]; f[1] has the edge v[1]->v[0]
};

struct adjacencyInfo_t {
	LIST( edge_t ) edges;
	LIST( LIST( int ) ) vertEdges; // edges adjacent to each vertex
	adjacencyInfo_t() {
		edges.setGranularity( 512 );
		vertEdges.setGranularity( 512 );
	}
};

bool BoundsPlaneIntersect( const MBoundingBox& bounds, const MPlane& plane ) {
	MPoint p[8] = { MPoint( bounds.min().x, bounds.min().y, bounds.min().z ),
					MPoint( bounds.max().x, bounds.min().y, bounds.min().z ),
					MPoint( bounds.min().x, bounds.max().y, bounds.min().z ),
					MPoint( bounds.max().x, bounds.max().y, bounds.min().z ),
					MPoint( bounds.min().x, bounds.min().y, bounds.max().z ),
					MPoint( bounds.max().x, bounds.min().y, bounds.max().z ),
					MPoint( bounds.min().x, bounds.max().y, bounds.max().z ),
					MPoint( bounds.max().x, bounds.max().y, bounds.max().z ) };

	bool front = false, back = false;

	for( int i = 0; i < 8; i++ ) {
		const double dist = plane.directedDistance( p[ i ] );
		if ( dist >= CUT_EPSILON ) {
			front = true;
		} else {
			back = true;
		}
		if ( front && back ) {
			return true;
		}
	}
	return false;
}

void splitEdge( edge_t& edge, const int iEdge, MIntArray &triangles, LIST( bool )& triangleIsInterior, int iCut, adjacencyInfo_t &adjacency ) {
	
	/*
		 	  vi          edge: vi -> vj
		    / |\          edge.f[0] = tf
	      /   |  \        edge.f[1] = tb
	    / tb' | tf \      C = edge cut point
	 vkb------C------vkf  
		\  tb | tf' /     tf is split into tf, tf'
		  \   |   /       tb is split into tb, tb'
			\ |  /        edge is split into vi->C, C->vj
			 \|/
			  vj
	*/

	const int vi = edge.v[0]; 
	const int vj = edge.v[1];
	const int tf = edge.f[0];
	assert( tf != -1 );
	const int tb = edge.f[1];

	// fix edge : vi --> iCut
	edge.v[1] = iCut;
	adjacency.vertEdges[ vj ].removeFast( iEdge );
	adjacency.vertEdges[ iCut ].append( iEdge );

	{ // split front tri

		int edgeTriVertOffset = 0;

		for( int j = 0; j < 3; j++ ) {
			if ( triangles[ 3 * tf + edgeTriVertOffset ] == vi &&
				triangles[ 3 * tf + ( edgeTriVertOffset + 1 ) % 3 ] == vj ) {
					break;
			}
			edgeTriVertOffset++;
		}

		if ( edgeTriVertOffset == 3 ) {
			assert( false );
			return;
		}
	
		const int ii = edgeTriVertOffset;
		const int ij = ( edgeTriVertOffset + 1 ) % 3;
		const int ik = ( edgeTriVertOffset + 2 ) % 3;
		assert( triangles[ 3 * tf + ii ] == vi );
		assert( triangles[ 3 * tf + ij ] == vj );

		const int vkf = triangles[ 3 * tf + ik ];

		triangles[ 3 * tf + ij ] = iCut;

		// append new face: fT'

		int newTri = triangles.length() / 3;

		triangles.append( vj );
		triangles.append( vkf );
		triangles.append( iCut );

		const bool tfInterior = triangleIsInterior[ tf ];
		triangleIsInterior.append( tfInterior );


		// fix edge vj --> vkf			
		for( size_t i = 0; i < adjacency.vertEdges[ vj ].size(); i++ ) {
			edge_t& e = adjacency.edges[ adjacency.vertEdges[ vj ][ i ] ];
			if ( e.v[ 0 ] == vj && e.v[ 1 ] == vkf && e.f[ 0 ] == tf ) {
				e.f[ 0 ] = newTri;
				break;
			} else if ( e.v[ 1 ] == vj && e.v[ 0 ] == vkf && e.f[ 1 ] == tf ) {
				e.f[ 1 ] = newTri;
				break;
			}
		}

		// newEdge1 : iCut --> vj
		const int iNewEdge1 = (int)adjacency.edges.size();
		edge_t& newEdge1 = adjacency.edges.append();

		newEdge1.v[ 0 ] = iCut;	
		newEdge1.v[ 1 ] = vj;
		newEdge1.f[ 0 ] = newTri;
		newEdge1.f[ 1 ] = edge.f[1];				
		adjacency.vertEdges[ iCut ].append( iNewEdge1 );
		adjacency.vertEdges[ vj ].append( iNewEdge1 );

		// newEdge 2 : vkf --> iCut
		const int iNewEdge2 = (int)adjacency.edges.size();
		edge_t& newEdge2 = adjacency.edges.append();

		newEdge2.v[ 0 ] = vkf;	
		newEdge2.v[ 1 ] = iCut;
		newEdge2.f[ 0 ] = newTri;
		newEdge2.f[ 1 ] = tf;				
		adjacency.vertEdges[ vkf ].append( iNewEdge2 );
		adjacency.vertEdges[ iCut ].append( iNewEdge2 );
	}

	{ // split back tri
		if ( tb == -1 ) {
			return;
		}

		int edgeTriVertOffset = 0;

		for( int j = 0; j < 3; j++ ) {
			if ( triangles[ 3 * tb + edgeTriVertOffset ] == vj &&
				triangles[ 3 * tb + ( edgeTriVertOffset + 1 ) % 3 ] == vi ) {
					break;
			}
			edgeTriVertOffset++;
		}

		if ( edgeTriVertOffset == 3 ) {
			assert( false );
			return;
		}

		const int ii = ( edgeTriVertOffset + 1 ) % 3;
		const int ij = edgeTriVertOffset;
		const int ik = ( edgeTriVertOffset + 2 ) % 3;
		assert( triangles[ 3 * tb + ii ] == vi );
		assert( triangles[ 3 * tb + ij ] == vj );

		const int vkb = triangles[ 3 * tb + ik ];

		triangles[ 3 * tb + ii ] = iCut;
		
		// append new face: fB'

		int newTri = triangles.length() / 3;

		triangles.append( iCut );
		triangles.append( vi );
		triangles.append( vkb );

		const bool tbInterior = triangleIsInterior[ tb ];
		triangleIsInterior.append( tbInterior );

		// fix edge vi --> vkf			
		for( size_t i = 0; i < adjacency.vertEdges[ vi ].size(); i++ ) {
			edge_t& e = adjacency.edges[ adjacency.vertEdges[ vi ][ i ] ];
			if ( e.v[ 0 ] == vi && e.v[ 1 ] == vkb && e.f[ 0 ] == tb ) {
				e.f[ 0 ] = newTri;
				break;
			} else if ( e.v[ 1 ] == vi && e.v[ 0 ] == vkb && e.f[ 1 ] == tb ) {
				e.f[ 1 ] = newTri;
				break;
			}
		}

		// fix edge vi --> iCut			
		edge.f[ 1 ] = newTri;

		// newEdge 2 : vkf --> iCut
		const int iNewEdge2 = (int)adjacency.edges.size();
		edge_t& newEdge2 = adjacency.edges.append();

		newEdge2.v[ 0 ] = vkb;	
		newEdge2.v[ 1 ] = iCut;
		newEdge2.f[ 0 ] = newTri;
		newEdge2.f[ 1 ] = tb;				
		adjacency.vertEdges[ vkb ].append( iNewEdge2 );
		adjacency.vertEdges[ iCut ].append( iNewEdge2 );
	}
}

void splitEdges( const int vi, const int vj, MPointArray& vertices, MIntArray& triangles, LIST( bool )& triangleIsInterior, double distVi, double distVj, adjacencyInfo_t& adjacency ) {
	
	const double edgeLength = abs( distVi ) + abs( distVj );
	const double bias = abs( distVi ) / edgeLength;
	
	MPoint cut = vertices[ ((unsigned int)vi) ] + ( vertices[ (unsigned int)vj ] - vertices[ (unsigned int)vi ] ) * bias;
	int iCut = vertices.length();
	vertices.append( cut );
	LIST( int )& ve = adjacency.vertEdges.append();
	ve.setGranularity( 32 );

	// split every edge linking vi with vj
	for( size_t i = 0; i < adjacency.vertEdges[ vi ].size(); i++ ) {
		const int iEdge = adjacency.vertEdges[ vi ][ i ];
		edge_t& edge = adjacency.edges[ iEdge ];
		if( !( ( edge.v[ 0 ] == vi && edge.v[ 1 ] == vj ) || ( edge.v[ 0 ] == vj && edge.v[ 1 ] == vi ) ) ) {
			continue;
		}
					
		// split edge
		splitEdge(edge, iEdge, triangles, triangleIsInterior, iCut, adjacency);				

	}
}

bool splitTriangle( MIntArray& triangles, const int t, MPointArray& vertices, LIST( bool )& triangleIsInterior, const MPlane& Plane, adjacencyInfo_t& adjacency ) {
	bool front = false, back = false;
	int isectOffset = 0;
	double dist[3];
	int iv[ 3 ];
	for( int i = 0; i < 3; i++ ) {
		iv[ i ] = triangles[ 3 * t + i ]; // copy the triangle before we start splitting it
		const MPoint v = vertices[ iv[ i ] ];
		dist[ i ] = Plane.directedDistance( v );
	}
	
	// split all the affected edges
	for( int i = 0; i < 3; i++ ) {
		const int j = ( i + 1 ) % 3;
		if ( ( dist[ j ] < -CUT_EPSILON && dist[ i ] >= CUT_EPSILON ) || 
			 ( dist[ i ] < -CUT_EPSILON && dist[ j ] >= CUT_EPSILON ) ) {
			splitEdges( iv[ i ], iv[ j ], vertices, triangles, triangleIsInterior, dist[ i ], dist[ j ], adjacency );
		}
	}

	return true;
}

void PurgeTris( MIntArray& triangleVertices, MPointArray& vertices, LIST( bool )& triangleIsInterior, const MPlane& plane, const MeshCutType cutType ) {
	
	if ( cutType == MESHCUT_DISCARD_NONE ) {
		return;
	}

	MIntArray vertexMapping( vertices.length(), -1 );

	{	
		MPointArray srcVertices = vertices;
		vertices.clear();
		for( unsigned int i = 0; i < srcVertices.length(); i++ ) {
			const MPoint& vertex = srcVertices[ i ];
			const double dist = plane.directedDistance( vertex );
			bool discarded = ( cutType == MESHCUT_DISCARD_FRONT && dist > CUT_EPSILON ) ||
							 ( cutType == MESHCUT_DISCARD_BACK && dist < -CUT_EPSILON );

			if ( !discarded ) {
				vertexMapping[ i ] = vertices.length();
				vertices.append( vertex );
			}
		}
	}

	{
		MIntArray srcTriangles = triangleVertices;
		LIST( bool ) triangleWasInterior;
		triangleWasInterior.swap( triangleIsInterior );
		triangleIsInterior.setGranularity( triangleWasInterior.getGranularity() );
		const unsigned int numTris = triangleVertices.length() / 3;
		triangleVertices.clear();
	
		int a,b,c;
		for( unsigned int i = 0; i < numTris; i++ ) {
			a = srcTriangles[ 3 * i + 0 ];
			b = srcTriangles[ 3 * i + 1 ];
			c = srcTriangles[ 3 * i + 2 ];
			if ( vertexMapping[ a ] == -1 || vertexMapping[ b ] == -1 || vertexMapping[ c ] == -1 ) {
				// discarded triangle
				continue;
			}
			triangleVertices.append( vertexMapping[ a ] );
			triangleVertices.append( vertexMapping[ b ] );
			triangleVertices.append( vertexMapping[ c ] );
			triangleIsInterior.append( triangleWasInterior[ i ] );
		}
	}
}

void MergeVertices( MIntArray& triangleVertices, MPointArray& vertices, LIST( bool )& triangleIsInterior ) {
	
	LIST( IntPair_t ) vertexRemap; // we'll merge close vertices

	for( unsigned int i = 0; i < vertices.length() - 1; i++ ) {
		for( unsigned int j = i + 1; j < vertices.length(); j++ ) {
			if ( vertices[ i ].distanceTo( vertices[ j ] ) <= 2 * MERGE_EPSILON ) {
				IntPair_t& remap = vertexRemap.append();
				remap.first = j;
				remap.second = i;
			}
		}
	}

	bool merges;
	do {
		merges = false;
		for( int i = 0; i < (int)vertexRemap.size() - 1; i++ ) {
			for( int j = i + 1; j < (int)vertexRemap.size(); j++ ) {
				if ( vertexRemap[ j ].second == vertexRemap[ i ].first ) {
					vertexRemap[ j ].second = vertexRemap[ i ].second;
					merges = true;
				}
				if ( vertexRemap[ j ].first == vertexRemap[ i ].first && 
					vertexRemap[ j ].second == vertexRemap[ i ].second ) {
						vertexRemap.removeIndexFast( j );
						j--;
						merges = true;
				}
			}
		}
	} while( merges );

	for( size_t r = 0; r < vertexRemap.size(); r++ ) {
		IntPair_t& remap = vertexRemap[ r ];

		for( unsigned int t = 0; t < triangleVertices.length(); t += 3 ) {
			
			for( int i = 0; i < 3; i++ ) {
				if ( triangleVertices[ t + i ] == remap.first ) {
					triangleVertices[ t + i ] = remap.second;
				}
			}

			if ( triangleVertices[ t + 0 ] == triangleVertices[ t + 1 ] ||
				triangleVertices[ t + 0 ] == triangleVertices[ t + 2 ] ) {
					// degenerate triangle. Remove
				triangleVertices[ t + 0 ] = triangleVertices[ triangleVertices.length() - 3 ];
				triangleVertices[ t + 1 ] = triangleVertices[ triangleVertices.length() - 2 ];
				triangleVertices[ t + 2 ] = triangleVertices[ triangleVertices.length() - 1 ];
				triangleVertices.setLength( triangleVertices.length() - 3 );
				triangleIsInterior.removeIndexFast( t / 3 );
				t -= 3;
			}
		}
	}
}

int EdgePairSort( const void* a, const void* b ) {
	const int* edge1 = static_cast< const int* >(a);
	const int* edge2 = static_cast< const int* >(b);
	if ( edge1[ 0 ] != edge2[ 0 ] ) {
		return edge1[ 0 ] - edge2[ 0 ]; 
	} else {
		return edge1[ 1 ] - edge2[ 1 ]; 
	}
}

void MergeCollinearEdges( LIST( RenderLib::Math::Vector2f )& verts, LIST( int )& edges ) {
	
	using namespace RenderLib::Math;

	bool touched = false;
	bool changed;
	do {
		const size_t orgEdgeNum = edges.size();

		for( size_t i = 0; i < edges.size() - 2 ; i += 2 ) {
			int from = edges[ i ];
			int to = edges[ i + 1 ];
			const Vector2f normalEdge = Vector2f::normalize( verts[ to ] - verts[ from ] );
			size_t j = i + 2;
			for( ; j < edges.size(); j += 2 ) {
				const int fromJ = edges[ j ];
				const int toJ = edges[ j + 1 ];

				if ( !( from == fromJ || from == toJ || to == fromJ || to == toJ ) ) {
					continue;
				} 
				if ( ( from == fromJ && to == toJ ) || ( from == toJ && to == fromJ ) ) {
					continue;
				}

				const int commonPoint = ( from == fromJ || from == toJ ) ? from : to;
				const int src1 = from != commonPoint ? from : to;
				const int src2 = fromJ != commonPoint ? fromJ : toJ;

				assert( src1 != src2 );
				assert( src1 != commonPoint );

				// src1 ------ commonPoint | commonPoint -------- src2 ------ next

				const float cosine = Vector2f::dot(normalEdge, Vector2f::normalize( verts[ commonPoint ] - verts[ src2 ] ));
				const bool collinear = fabsf( cosine ) > 1.f - COLLINEAR_EPSILON;
				if ( collinear ) {

					// ensure src2 only has adjacency 2 (that is, it's only linked to commonPoint and 
					// another vertex "next", before we remove commonPoint--src2

					int adjacencySrc2 = 0;
					int adjacencyCommon = 0;
					for( size_t k = 0; k < edges.size(); k++ ) {
						if ( edges[ k ] == src2 ) {
							adjacencySrc2++;
							if ( adjacencySrc2 > 2 ) {
								break;
							}
						} else if ( edges[ k ] == commonPoint ) {
							adjacencyCommon++;
							if ( adjacencyCommon > 2 ) {
								break;
							}
						}

					}

					if ( adjacencySrc2 > 2 || adjacencyCommon > 2 ) {
						continue;
					}

					if ( commonPoint == from ) {
						edges[ i ] = src2;
						from = src2;
					} else {
						edges[ i + 1 ] = src2;
						to = src2;
					}
					edges.removeIndexFast( j + 1 );
					edges.removeIndexFast( j );
					j -= 2;
				}

			}			
		}

		changed = edges.size() != orgEdgeNum;
		touched |= changed;
	} while( changed );

	if ( touched ) {
		// need to sort the edges array again
		qsort( &edges[ 0 ], edges.size() / 2, 2 * sizeof(int), EdgePairSort );		
	}
}

bool LongestCycles( const LIST( RenderLib::Math::Vector2f )& verts, LIST( int )& edges ) {

	int correctVertices = 0;
	int singleLinked = 0;
	int multipleLinked = 0;
	LIST( LIST( int ) ) adjacency;
	adjacency.resize( verts.size() );
	LIST( int ) uniqueVertices( verts.size() );
	for( size_t i = 0; i < edges.size(); i += 2 ) {
		const int from = edges[ i ];
		const int to = edges[ i + 1 ];
		uniqueVertices.addUnique( from );
		uniqueVertices.addUnique( to );
		const int numAdjacents_i = adjacency[ from ].append( to ) + 1;
		const int numAdjacents_i1 = adjacency[ to ].append( from ) + 1;
		if ( numAdjacents_i == 1 ) {
			singleLinked++;
		} else if ( numAdjacents_i == 2 ) {
			correctVertices++;
			singleLinked--;
		} else if ( numAdjacents_i == 3 ) {
			// it was correct, but now it's not
			correctVertices--;
			multipleLinked++;
		}

		if ( numAdjacents_i1 == 1 ) {
			singleLinked++;
		} else if ( numAdjacents_i1 == 2 ) {
			correctVertices++;
			singleLinked--;
		} else if ( numAdjacents_i1 == 3 ) {
			// it was correct, but now it's not
			correctVertices--;
			multipleLinked++;
		}		
	}
	
	if ( correctVertices == (int)uniqueVertices.size() ) {
		// further processing is not necessary
		return true;
	}

	if ( singleLinked > 0 && multipleLinked == 0 ) {
		// if there are single linked segments (dead ends), but no
		// segments with more than 2 links (multiple linked) we know
		// for sure the segments describe disjoint lines with no loops
		return false;
	}

	LIST( int ) candidateEdges;
	candidateEdges.swap( edges );

	LIST( int ) path;
	path.preAllocate( verts.size() );
	LIST( int ) longestCycle;
	longestCycle.preAllocate( verts.size() );

	LIST( int ) bifurcations( 16 );
	while( !candidateEdges.empty() ) {

		int source = candidateEdges[ 0 ];
		longestCycle.resize( 0, false );

		do {

			bifurcations.resize( 0, false );
			path.resize( 0, false );
			
			int v = source;
			int prev = source;		
		
			do {
				path.append( v );
				int next = v;

				if ( prev != v ) {
					// delete the edge we just traversed
					for( size_t i = 0; i < candidateEdges.size(); i += 2 ) {
						if ( ( candidateEdges[ i ] == prev && candidateEdges[ i + 1 ] == v ) ||
							( candidateEdges[ i ] == v && candidateEdges[ i + 1 ] == prev ) ) {
								candidateEdges.removeIndexFast( i + 1 );
								candidateEdges.removeIndexFast( i );
								break;
						}
					}					
				}
				
				for( size_t i = 0; i < adjacency[ v ].size(); i++ ) {
					if( adjacency[ v ][ i ] != prev ) {
						next = adjacency[ v ][ i ];
						if ( adjacency[ v ].size() > 2 ) {
							adjacency[ v ].removeIndexFast( i );
							adjacency[ next ].removeFast( v );
							bifurcations.append( path.size() - 1 );
						}
						break;
					}
				}				
				if ( next == v ) { 	
					if ( bifurcations.size() > 0 ) {
						const int bifurcationIdx = bifurcations[ bifurcations.size() - 1 ];
						// delete the latest part of the path
						for( int j = bifurcationIdx; j < (int)path.size() - 1; j ++ ) {
							for( size_t i = 0; i < candidateEdges.size(); i += 2 ) {
								// this is a 1-adjacent path, we're fine just checking for 1 of the endpoints
								if ( candidateEdges[ i ] == path[ j ] || candidateEdges[ i + 1 ] == path[ j ] ) {
									candidateEdges.removeIndexFast( i + 1 );
									candidateEdges.removeIndexFast( i );
									i -= 2;
								}
							}
						}

						// reset the path to the latest bifurcation
						bifurcations.resize( bifurcations.size() - 1, false );
						v = path[ std::max( 0, bifurcationIdx - 1  ) ];
						prev = v;
						next = path[ bifurcationIdx ]; 
						path.resize( bifurcationIdx, false ); // pop the bifurcation element out as well, as it will be reinserted when we loop over						
					} else {
						// dead end path, delete it
						for( size_t i = 0; i < candidateEdges.size(); i += 2 ) {
							// this is a 1-adjacent path, we're fine just checking for 1 of the endpoints
							if ( candidateEdges[ i ] == source || candidateEdges[ i + 1 ] == source ) {
								candidateEdges.removeIndexFast( i + 1 );
								candidateEdges.removeIndexFast( i );
								i -= 2;
							}
						}
						path.resize( 0, false );						
					}					
				}

				prev = v;
				v = next;

				if ( v != source ) {
					const int idx = path.findIndex( v );
					if ( idx >= 0 ) {
						// we found a loop which doesn't start on our source vertex (like a "6" shape)
						// this is a symmetrical case to when we were deleting the latest part of a path
						// (that part is the beginning of the path now).
						for( int j = 0; j < bifurcations[ 0 ]; j++ ) {
							for( size_t i = 0; i < candidateEdges.size(); i += 2 ) {
								// this is a 1-adjacent path, we're fine just checking for 1 of the endpoints
								if ( candidateEdges[ i ] == path[ j ] || candidateEdges[ i + 1 ] == path[ j ] ) {
									candidateEdges.removeIndexFast( i + 1 );
									candidateEdges.removeIndexFast( i );
									i -= 2;
								}
							}
						}
						// FIXME: reuse this loop?
						path.resize( 0, false );
						break;
					}
				}

			} while( v != source && !path.empty() );			

			if ( path.size() > 0 && ( longestCycle.empty() || path.size() > longestCycle.size() ) ) {
				longestCycle.swap( path );
			}		
		} while( !bifurcations.empty() );

		for( size_t i = 0; i < longestCycle.size(); i++ ) {
			for( size_t j = 0; j < candidateEdges.size(); j += 2 ) {
				if ( ( candidateEdges[ j ] == longestCycle[ i ] && candidateEdges[ j + 1 ] == longestCycle[ ( i + 1 ) % longestCycle.size() ] ) ||
					( candidateEdges[ j ] == longestCycle[ ( i + 1 ) % longestCycle.size() ] && candidateEdges[ j + 1 ] == longestCycle[ i ] ) ) {
						candidateEdges.removeIndexFast( j + 1 );
						candidateEdges.removeIndexFast( j );
						break;
				}
			}			
			edges.append( longestCycle[ i ] );
			edges.append( longestCycle[ ( i + 1 ) % longestCycle.size() ] );
			assert( edges.size() < uniqueVertices.size() * uniqueVertices.size() ); // check for nasty loops going on here
		}
	}

	return true;
}

bool meshCut( MPointArray& vertices, MIntArray& triangleVertices, LIST( bool )& triangleIsInterior, const MPlane& plane, MeshCutType cutType ) {

	using namespace RenderLib::Math;

	const int numTriangles = triangleVertices.length() / 3;

	if ( numTriangles == 0 ) { 
		return false;
	}

	// trivial discard: do nothing if the plane doesn't intersect the mesh bounds
	MBoundingBox bounds;	

	{ // mesh.boundingbox() doesn't seem to be calculated 
		bounds.clear();
		for( unsigned int i = 0; i < vertices.length(); i++ ) {
			bounds.expand( vertices[ i ] );
		}
	}

	if ( !BoundsPlaneIntersect( bounds, plane ) ) { // skip the geometry splitting if the plane does not even intersect the bounds (we might purge away all the triangles, though)

		switch( cutType ) {
			case MESHCUT_DISCARD_NONE:
				return triangleVertices.length() > 0;
			case MESHCUT_DISCARD_FRONT:
				{
					const double dist = plane.directedDistance( bounds.center() );
					return dist > CUT_EPSILON ? false : true;
				}
			case MESHCUT_DISCARD_BACK:
				{
					const double dist = plane.directedDistance( bounds.center() );
					return dist < CUT_EPSILON ? false : true;
				}
		}

	}
	
	// Build adjacency caches
	
	adjacencyInfo_t adjacency;
	adjacency.edges.preAllocate( 3 * numTriangles ); // in theory 2 * numTriangles should suffice given the shared edges
	adjacency.vertEdges.resize( triangleVertices.length() );

	for( int t = 0; t < numTriangles; t++ ) {

		const int v[3] = {  triangleVertices[ 3 * t + 0 ],
							triangleVertices[ 3 * t + 1 ], 
							triangleVertices[ 3 * t + 2 ] };
		
		int ie[3] = { -1, -1, -1 };
		for( int i = 0; i < 3; i++ ) {
			const int vi = v[ i ];
			const int vj = v[ ( i + 1 ) % 3 ];

			// check if the reversed edge does exist already
			for( size_t j = 0; j < adjacency.vertEdges[ vj ].size(); j++ ) {
				const int iedge = adjacency.vertEdges[ vj ][ j ];
				edge_t& edge = adjacency.edges[ iedge ];				
				if ( edge.v[ 0 ] == vj && edge.v[ 1 ] == vi && edge.f[ 1 ] == -1 ) {
					ie[ i ] = iedge;
					edge.f[ 1 ] = t;
					break;
				} 
			}
			if ( ie[ i ] == -1 ) { // the edge does not exist - create it
				ie[ i ] = (int)adjacency.edges.size();
				edge_t& e = adjacency.edges.append();
				e.v[ 0 ] = vi;
				e.v[ 1 ] = vj;
				e.f[ 0 ] = t;
				e.f[ 1 ] = -1;
				adjacency.vertEdges[ e.v[ 0 ] ].append( ie[ i ] );
				adjacency.vertEdges[ e.v[ 1 ] ].append( ie[ i ] );						
			} 
		}
	}	

	// split mesh
	for( unsigned int t = 0; t < triangleVertices.length() / 3 ; t++ ) { // triangleVertices will grow while we split each triangle
		splitTriangle( triangleVertices, t, vertices, triangleIsInterior, plane, adjacency );
	}

	// identify the vertices which lie on the plane. It is not enough to
	// assume that they'll be the new vertices created by the triangle splitting
	// as they may be previously existing vertices that were already on the plane.
	// Example degenerate case: a plane splitting right along the middle of a symmetric model
	// (along a ring of edges) will not generate any new vertices.
	LIST( int ) verticesOnPlane( triangleVertices.length() / 3 * 2 );

	LIST( IntPair_t ) vertexRemap; // we'll merge close vertices
	
	LIST( int ) edgeVertices( triangleVertices.length() / 3 * 2 );

	for( unsigned int v = 0; v < triangleVertices.length(); v++ ) { 
		const int candidate = triangleVertices[ v ];

		const double planeDist = plane.distance( vertices[ candidate ] );
		if ( planeDist <= ON_PLANE_EPSILON ) {

			edgeVertices.addUnique( candidate );

			bool accepted = true;
			for( size_t v2 = 0; v2 < verticesOnPlane.size(); v2 ++ ) {
				const int vop = verticesOnPlane[ v2 ];
				
				if ( vop == candidate ) {
					accepted = false;
					break;
				} else if ( vertices[ vop ].distanceTo( vertices[ candidate ] ) < MERGE_EPSILON ) {				

					IntPair_t remap( candidate, vop );
					vertexRemap.addUnique( remap );
					accepted = false;
					break;
				}
			}
			if ( accepted ) {
				verticesOnPlane.append( candidate );
			}
		}
	}

#if _DEBUG && DEBUG_STEPS
	static int cutEdge = 0;
	cutEdge++;
	for( size_t i = 0; i < edgeVertices.size(); i ++ ) {
		const int v0 = edgeVertices[ i ];
		for( size_t j = 0; j < adjacency.vertEdges[ v0 ].size(); j++ ) {
			const edge_t& edge = adjacency.edges[ adjacency.vertEdges[ v0 ][ j ] ];
			const int v1 = edge.v[ 0 ] == v0 ? edge.v[ 1 ] : edge.v[ 0 ];
			if ( edgeVertices.FindIndex( v1 ) == -1 ) {
				continue;
			}
			MPoint& a = vertices[ v0 ];
			MPoint& b = vertices[ v1 ];
			MString cmd = CoreLib::varArgsStr<1024>( "$c = `curve -d 1 -p %f %f %f -p %f %f %f -k 0 -k 1`; rename $c cutEdge%d_%d_%d;", a.x, a.y, a.z, b.x, b.y, b.z, cutEdge, v0, v1 );
			MGlobal::executeCommand( cmd );

		}
	}
#endif
	if ( verticesOnPlane.size() > 2 ) {
		// We know for sure the new vertices are coplanar since they've been generated from
		// splitting the triangles by the plane. This reduces our hole filling problem to 2D

		LIST( Vector2f ) vertices2D;
		vertices2D.preAllocate( verticesOnPlane.size() );
		LIST( int ) edges;
		edges.setGranularity( verticesOnPlane.size() );
		LIST( int ) fillingTriangles;
		fillingTriangles.setGranularity( verticesOnPlane.size() );

		// Obtain the projected points in a 2D basis
		MVector planeU = vertices[ verticesOnPlane[ 1 ] ] - vertices[ verticesOnPlane[ 0 ] ]; // these points are coplanar, so the resulting vector should be perpendicular to the normal
		planeU.normalize();
		MVector planeV = plane.normal() ^ planeU;
		planeU = planeV ^ plane.normal(); // ensure the base is orthonormal

		for( size_t i = 0; i < verticesOnPlane.size(); i++ ) {
			const MVector& v3D = vertices[ verticesOnPlane[ i ] ];
			assert( plane.distance( v3D ) <= ON_PLANE_EPSILON );
			const float u = (float)(v3D * planeU);
			const float v = (float)(v3D * planeV);
			vertices2D.append( Vector2f( u, v ) );
		}

		for( size_t i = 0; i < edgeVertices.size(); i++ ) {
			int vIdx = edgeVertices[ i ];
			int remappedI = verticesOnPlane.findIndex( vIdx );
			if ( remappedI == -1 ) { // the candidate vertex doesn't exist on verticesOnPlane, therefore it must have been merged with another vertex
				for( size_t r = 0; r < vertexRemap.size(); r++ ) {
					if ( vertexRemap[ r ].first == vIdx ) {							
						const int remappedVIdx = vertexRemap[ r ].second;
						remappedI = verticesOnPlane.findIndex( remappedVIdx );
						break;
					}
				}
			}
			assert( remappedI >= 0 );

			const size_t nEdges = adjacency.vertEdges[ vIdx ].size();
			for( size_t j = 0; j < nEdges; j++ ) {
				const edge_t& edge = adjacency.edges[ adjacency.vertEdges[ vIdx ][ j ] ];
				int otherVertex = edge.v[ 0 ] != vIdx ? edge.v[ 0 ] : edge.v[ 1 ];
				
				for( size_t r = 0; r < vertexRemap.size(); r++ ) {
					if ( vertexRemap[ r ].first == otherVertex ) {							
						otherVertex = vertexRemap[ r ].second;
						break;
					}
				}
				
				const int vIdxOther = verticesOnPlane.findIndex( otherVertex );
				if ( vIdxOther >= 0 ) {
					const int from = std::min( remappedI, vIdxOther );					
					const int to = std::max( remappedI, vIdxOther );					
					edges.append( from );
					edges.append( to );					
				}
			}
		}
		
		{ // remove duplicate and degenerate edges
			LIST( int ) uniqueEdges( edges.size() );

			for( size_t i = 0; i < edges.size(); i += 2 ) {
				if ( edges[ i ] == edges[ i + 1 ] ) {
					// degenerate edge
					edges.removeIndexFast( i + 1 );
					edges.removeIndexFast( i );
					i -= 2;
				}

			}
			qsort( &edges[ 0 ], edges.size() / 2, 2 * sizeof(int), EdgePairSort );
			for( size_t i = 0; i < edges.size(); i += 2 ) {
				size_t j = i + 2;
				for( ; j < edges.size(); j += 2 ) {
					if ( ( edges[ j ] != edges[ i ] ) || ( edges[ j + 1 ] != edges[ i + 1 ] ) ) {
						break;
					}
				}
				uniqueEdges.append( edges[ i ] );
				uniqueEdges.append( edges[ i + 1 ] );
				i = j - 2;				
			}
			edges.swap( uniqueEdges );
		}		

#if _DEBUG
		LIST( int ) edgesCopy = edges; // LongestCycles is a destructive operation
#endif

#if _DEBUG && DEBUG_STEPS
		static int holeId = 0;
		holeId++;
		for( size_t i = 0; i < edges.size(); i += 2 ) {
			Vector2f a = vertices2D[ edges[ i ] ];
			Vector2f b = vertices2D[ edges[ i + 1 ] ];
			MString cmd = CoreLib::varArgsStr<1024>( "$c = `curve -d 1 -p %f %f 0 -p %f %f 0 -k 0 -k 1`; rename $c inputEdge%d_%d_%d;", a.x, a.y, b.x, b.y, holeId, edges[ i ], edges[ i + 1 ] );
			MGlobal::executeCommand( cmd );
		}
#endif
		MergeCollinearEdges( vertices2D, edges );
#if _DEBUG && DEBUG_STEPS
		for( size_t i = 0; i < edges.size(); i += 2 ) {
			Vector2f a = vertices2D[ edges[ i ] ];
			Vector2f b = vertices2D[ edges[ i + 1 ] ];
			MString cmd = CoreLib::varArgsStr<1024>( "$c = `curve -d 1 -p %f %f 0 -p %f %f 0 -k 0 -k 1`; rename $c mergedEdge%d_%d_%d;", a.x, a.y, b.x, b.y, holeId, edges[ i ], edges[ i + 1 ] );
			MGlobal::executeCommand( cmd );
		}
#endif

		// merge close vertices

		LIST( IntPair_t ) vertexRemap( vertices2D.size() );

		for( size_t i = 0; i < vertices2D.size(); i++ ) {
			for( size_t j = i + 1; j < vertices2D.size(); j++ ) {
				if ( ( vertices2D[ i ] - vertices2D[ j ] ).length() <= 2 * MERGE_EPSILON ) {

					// remap j to i
					IntPair_t& remap1 = vertexRemap.append();
					remap1.first = j;
					remap1.second = i;										
				}
			}
		}		

		bool merges;
		do {
			merges = false;
			for( int i = 0; i < (int)vertexRemap.size() - 1; i++ ) {
				for( int j = i + 1; j < (int)vertexRemap.size(); j++ ) {
					if ( vertexRemap[ j ].second == vertexRemap[ i ].first ) {
						vertexRemap[ j ].second = vertexRemap[ i ].second;
						merges = true;
					}
					if ( vertexRemap[ j ].first == vertexRemap[ i ].first && 
						vertexRemap[ j ].second == vertexRemap[ i ].second ) {
							vertexRemap.removeIndexFast( j );
							j--;
							merges = true;
					}
				}
			}
		} while( merges );

		LIST( int ) activeVerticesRemap( vertices2D.size() );
		LIST( Vector2f ) activeVertices( vertices2D.size() );
		for( size_t i = 0; i < vertices2D.size(); i++ ) {
			bool include = true;
			for( size_t j = 0; j < vertexRemap.size(); j++ ) {
				if ( vertexRemap[ j ].first == i ) {
					include = false;
					break;
				}
			}
			if ( include ) {
				activeVertices.append( vertices2D[ i ] );
				activeVerticesRemap.append( i );
			}
		}
		for( size_t i = 0; i < edges.size(); i++ ) {
			for( size_t j = 0; j < vertexRemap.size(); j++ ) {
				if ( edges[ i ] == vertexRemap[ j ].first ) {
					edges[ i ] = vertexRemap[ j ].second;
					break;
				}
			}
			edges[ i ] = activeVerticesRemap.findIndex( edges[ i ] );
		}

		fillPlanarHoles( activeVertices, edges, 1.0f, fillingTriangles );

		// add hole filling triangles, remapping the vertex indices
		for( size_t i = 0; i < fillingTriangles.size(); i++ ) {
			int v = fillingTriangles[ i ];
			// Undo remapping to get the original triangle index
			// FIXME: when several vertices i,j,k have been merged to a single
			// vertex m, there'll be several entries (i, m) (j, m) (k, m) in the
			// vertex remap list. Since our search retrieves the first entry
			// which second == m, vertices will be restored as either i, j, or k
			// however this shouldn't matter as they're very close (that's why they
			// were merged!)
			v = activeVerticesRemap[ v ];
			/*for( int r = 0; r < vertexRemap.Num(); r++ ) {
				const IntPair_t& remap = vertexRemap[ r ];
				if ( remap.second == v ) {
					v = remap.first;
					break;
				}
			}*/
			assert( verticesOnPlane[ v ] < (int)vertices.length() );
			triangleVertices.append( verticesOnPlane[ v ] );
			triangleIsInterior.append( true );
		}
	}	

	PurgeTris( triangleVertices, vertices, triangleIsInterior, plane, cutType );

	MergeVertices( triangleVertices, vertices, triangleIsInterior );

	const int newTriCount = triangleVertices.length() / 3 ;
	
	return newTriCount > 0;
}
