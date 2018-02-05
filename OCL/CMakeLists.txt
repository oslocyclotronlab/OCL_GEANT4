# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(OCL)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(OCL OCL.cc ${sources} ${headers})
target_link_libraries(OCL ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Scint. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(OCLSCRIPTS
  init.mac
  init_vis.mac
  vis.mac
  vis_VRML.mac
  run1.mac
  run100.mac
  run1779keV.mac
  runsims.sh
  PrintSimMacs.cpp
  manyruns.sh
  )

foreach(_script ${OCLSCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
#
# Checks whether a "data" directory exists at the place GEANT4 wants to save
# the root files. 
# If it does not exist, this will lead to errors when you run the simulation
#
IF(EXISTS "../data")
       message("-- Data directory exists (needed for the root outputs)")
else() 
       message(SEND_ERROR " \n WARNING \n \"../data\" directory DOES NOT exists (needed for the root outputs)")
endif()

##----------------------------------------------------------------------------
## For internal Geant4 use - but has no effect if you build this
## example standalone
##
#add_custom_target(ScintTest DEPENDS Scint)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS OCL DESTINATION bin)


