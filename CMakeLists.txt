cmake_minimum_required(VERSION 2.6)

project(Shatter)

set ( MAYA_PLUGIN_NAME "Shatter" )

# Need to supply MAYA_ROOT
set( MAYA_ROOT "/usr/autodesk/maya/" CACHE STRING "path to maya installation folder")

set( CORELIB_INCLUDE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/CoreLib/include CACHE STRING "path to CoreLib include folder" )
set( CORELIB_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/CoreLib/lib CACHE STRING "path to CoreLib libray folder" )

set( RENDERLIB_INCLUDE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/RenderLib/include CACHE STRING "path to RenderLib include folder" )
set( RENDERLIB_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/RenderLib/lib CACHE STRING "path to RenderLib libray folder" )

# Platform specific definitions
if(WIN32)
    if(MSVC)
        # To enable /MP, parralel build on MSVC
        option(WIN32_USE_MP "Set to OFF to diable /MP." ON)
        MARK_AS_ADVANCED(WIN32_USE_MP)
        if(WIN32_USE_MP)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
        endif(WIN32_USE_MP)

        # Other MSVC compilation flags
        add_definitions(-D_SCL_SECURE_NO_WARNINGS)
        add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
		add_definitions(-D NOMINMAX) # prevent windef.h min/max macros from interfering with stl's
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")

    endif(MSVC)
endif (WIN32)


if(WIN32)
	set( MAYA_USER_DIR $ENV{USERPROFILE}/Documents/maya CACHE PATH "Maya user home directory" )
	set( MAYA_DEFINITIONS "-D_AFXDLL -D_MBCS -DNT_PLUGIN -DREQUIRE_IOSTREAM" )
else()
	set( MAYA_USER_DIR $ENV{HOME}/maya CACHE PATH "Maya user home directory" )
	add_definitions( -DREQUIRE_IOSTREAM )
	add_definitions( -D_BOOL )
	add_definitions( -fPIC )
	set( MAYA_DEFINITIONS "" )
endif()
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	message( "DEBUG BUILD" )
	add_definitions( -D_DEBUG )
endif( ${CMAKE_BUILD_TYPE} EQUAL "DEBUG" )

if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()

find_package(OpenGL)


set( MAYA_HEADERS_DIR ${MAYA_ROOT}/include )
set( MAYA_LIBRARY_DIR ${MAYA_ROOT}/lib )

set ( LOCAL_WARNING_FLAGS /W3 )
set ( LOCAL_RTTI_FLAGS /GR )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# add base to include directories
include_directories ( ${PUBLIC_INCLUDE_DIRS} ${MAYA_HEADERS_DIR} ${CORELIB_INCLUDE_DIR} ${RENDERLIB_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR})
link_directories ( ${PUBLIC_LIBRARY_DIRS} ${MAYA_LIBRARY_DIR} ${CORELIB_INCLUDE_DIR} ${RENDERLIB_LIBRARY_DIR})

set( CORE_LIB CoreLib )
set( RENDER_LIB RenderLib )
if (WIN32)
set( OPENGL_LIB opengl32 )
else()
set( OPENGL_LIB  GL )
endif()
set( MAYASDK_LIBRARIES "Foundation" "OpenMaya" "OpenMayaUI" "OpenMayaAnim" "OpenMayaFX" "OpenMayaRender" )


#specify app sources

function(create_source_group sourceGroupName relativeSourcePath)
	foreach(currentSourceFile ${ARGN})
		file(RELATIVE_PATH folder ${relativeSourcePath} ${currentSourceFile})
		get_filename_component(filename ${folder} NAME)
		string(REPLACE ${filename} "" folder ${folder})
		if(NOT folder STREQUAL "")
			string(REGEX REPLACE "/+$" "" folderlast ${folder})
			string(REPLACE "/" "\\" folderlast ${folderlast})
			source_group("${sourceGroupName}\\${folderlast}" FILES ${currentSourceFile})
		endif(NOT folder STREQUAL "")
	endforeach(currentSourceFile ${ARGN})
endfunction(create_source_group)

file(GLOB_RECURSE RENDERLIB_SOURCES 
	RenderLib/source/*.cpp 
	RenderLib/include/*.h 
	RenderLib/include/*.inl)
file(GLOB_RECURSE CORELIB_SOURCES 
	CoreLib/source/*.cpp 
	CoreLib/include/*.h 
	CoreLib/include/*.inl)
file(GLOB_RECURSE MAYA_PLUGIN_SOURCES
	src/*.cpp
	src/*.h)

# mimic disk folder structure on the project
#create_source_group( "coreLib" "${CMAKE_CURRENT_SOURCE_DIR}/CoreLib/" ${CORELIB_SOURCES} )
#create_source_group( "renderLib" "${CMAKE_CURRENT_SOURCE_DIR}/RenderLib/" ${RENDERLIB_SOURCES} )

# Core Lib
add_library(${CORE_LIB} STATIC ${CORELIB_SOURCES})
set_target_properties( ${CORE_LIB} PROPERTIES PREFIX "" )
set_target_properties( ${CORE_LIB} PROPERTIES OUTPUT_NAME ${CORE_LIB} )
set_target_properties( ${CORE_LIB} PROPERTIES LINKER_LANGUAGE CXX )
set_target_properties( ${CORE_LIB} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CORELIB_LIBRARY_DIR} )

# Render Lib

add_library(${RENDER_LIB} STATIC ${RENDERLIB_SOURCES})
target_link_libraries( ${RENDER_LIB} ${CORE_LIB})

set_target_properties( ${RENDER_LIB} PROPERTIES PREFIX "" )
set_target_properties( ${RENDER_LIB} PROPERTIES OUTPUT_NAME ${RENDER_LIB} )
set_target_properties( ${RENDER_LIB} PROPERTIES LINKER_LANGUAGE CXX )
set_target_properties( ${RENDER_LIB} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${RENDERLIB_LIBRARY_DIR} )


# Shatter Plugin
add_library( ${MAYA_PLUGIN_NAME} SHARED ${MAYA_PLUGIN_SOURCES} )
target_link_libraries( ${MAYA_PLUGIN_NAME} ${MAYASDK_LIBRARIES} ${OPENGL_LIB} ${CORE_LIB} ${RENDER_LIB})
#set_target_properties( ${MAYA_PLUGIN_NAME} PROPERTIES COMPILE_DEFINITIONS ${MAYA_DEFINITIONS} )
set_target_properties( ${MAYA_PLUGIN_NAME} PROPERTIES OUTPUT_NAME ${MAYA_PLUGIN_NAME} )
set_target_properties( ${MAYA_PLUGIN_NAME} PROPERTIES CLEAN_DIRECT_OUTPUT 1 )
if(WIN32)
set_target_properties( ${MAYA_PLUGIN_NAME} PROPERTIES SUFFIX ".mll" )
set_target_properties( ${MAYA_PLUGIN_NAME} PROPERTIES LINK_FLAGS "/export:initializePlugin /export:uninitializePlugin" )
endif(WIN32)

