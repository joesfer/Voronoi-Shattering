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


#define WIN32_LEAN_AND_MEAN

#include "SamplerNode.h"
#include "DelaunayData.h"
#include "DelaunayNode.h"
#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is loaded into Maya.  It 
//		registers all of the services that this plug-in provides with 
//		Maya.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{ 
	MStatus   status;
	MFnPlugin plugin( obj, "Jose Esteve. www.joesfer.com", "2011", "Any");

	// Add plug-in feature registration here
	//
	status = plugin.registerNode( "SamplerNode", Sampler::id, Sampler::creator, Sampler::initialize );
	if (!status) {
		status.perror("registerNode");
		return status;
	}
	
	status = plugin.registerData( ShatterNodeData::typeName, ShatterNodeData::id, ShatterNodeData::creator, MPxData::kGeometryData );
	if (!status) {
		status.perror("registerData");
		return status;
	}

	status = plugin.registerNode( "ShatterNode", ShatterNode::id, ShatterNode::creator, ShatterNode::initialize );
	if (!status) {
		status.perror("registerNode");
		return status;
	}


	return status;
}

MStatus uninitializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is unloaded from Maya. It 
//		deregisters all of the services that it was providing.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{
	MStatus   status;
	MFnPlugin plugin( obj );

	// Add plug-in feature deregistration here
	//

	status = plugin.deregisterNode( Sampler::id );
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

	status = plugin.deregisterData( ShatterNodeData::id );
	if (!status) {
		status.perror("deregisterData");
		return status;
	}

	status = plugin.deregisterNode( ShatterNode::id );
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

	return status;
}
