#***************************************************************************
# Authors:     Oier Lauzirika Zarrabeitia (oierlauzi@bizkaia.eu)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# 02111-1307  USA
#
#  All comments concerning this program package may be sent to the
#  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************

cmake_minimum_required(VERSION 3.16)

# Define the project
project(
  xmipp-core
	VERSION 3.23.11
	LANGUAGES C CXX
)

include(GNUInstallDirs)

# Find dependencies
find_package(HDF5 1.10 REQUIRED)
find_package(TIFF REQUIRED)

# Register all source and header files
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/core/*.c)
file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/include/*.h ${PROJECT_SOURCE_DIR}/include/*.hpp)

# Create the shared library
add_library(${PROJECT_NAME} SHARED ${SOURCES})
set_target_properties(
	${PROJECT_NAME} PROPERTIES 
	SOVERSION ${PROJECT_VERSION}
	CXX_STANDARD 20
)
target_include_directories(
	${PROJECT_NAME} 
	PUBLIC 
	  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_precompile_headers(
	${PROJECT_NAME} 
	PRIVATE
	${HEADERS}
)
target_link_libraries(
	${PROJECT_NAME}
	PUBLIC
		HDF5::HDF5
		TIFF::TIFF
)

# Install library's binary files and headers
install(
	TARGETS ${PROJECT_NAME} 
	EXPORT ${PROJECT_NAME}-targets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(
	DIRECTORY ${PROJECT_SOURCE_DIR}/include/ 
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
