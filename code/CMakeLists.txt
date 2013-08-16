cmake_minimum_required( VERSION 2.6 )

project( transformation )

# Turn on BUILD_TESTS to enable the compilation of the test suite which will
# validate that the transformation class produces compatible results to the
# system OpenGL implementation. This requires GLUT in order to create an
# OpenGL context, without which the OpenGL functions will not function.
set( BUILD_TESTS TRUE )

# Turn on one of USE_EIGEN or USE_GMM to enable the use of the Eigen (resp. GMM++) 
# linear algebra libraries. These are used by the transformation::inverse() and
# transformation::glUnProject() methods, which will otherwise use a slower method
# based on Cramer's rule but which has no dependencies. If you enable GMM or
# Eigen you may need to set the user_Eigen3_INCLUDE_DIRS and user_GMM_INCLUDE_DIRS
# variables below should the find_package call not locate the libraries on your
# system
set( USE_EIGEN  FALSE )
set( USE_GMM    FALSE )

# set user-specific locations for Eigen and GMM if needed
# set( user_Eigen3_INCLUDE_DIRS path/to/eigen )
# set( user_GMM_INCLUDE_DIRS    path/to/gmm )

# conditional compilation based on USE_EIGEN
if( USE_EIGEN )

  # first try to find eigen using find_package()
  # otherwise fall back to using the user path
  find_package( Eigen3 )
  if( Eigen3_FOUND )
    include_directories( ${Eigen3_INCLUDE_DIRS} )
    add_definitions( -DTRANSFORMATION_USE_EIGEN )
  elseif( Eigen3_FOUND )
    include_directories( ${user_Eigen3_INCLUDE_DIRS} )
    add_definitions( -DTRANSFORMATION_USE_EIGEN )
  else( Eigen3_FOUND )
    message( FATAL_ERROR "Could not locate Eigen, set the variable user_Eigen3_INCLUDE_DIRS for your system in CMakeLists.txt!" )
  endif( Eigen3_FOUND )

endif( USE_EIGEN )

# conditional compilation based on USE_GMM
if( USE_GMM )
  # first try to find GMM using find_package()  
  # otherwise fall back to using the user path
  # first try to find eigen using find_package()
  # otherwise fall back to using the user path
  find_package( GMM )
  if( GMM_FOUND )
    include_directories( ${GMM_INCLUDE_DIRS} )
    add_definitions( -DTRANSFORMATION_USE_GMM )
  elseif( user_GMM_INCLUDE_DIRS )
    include_directories( ${user_GMM_INCLUDE_DIRS} )
    add_definitions( -DTRANSFORMATION_USE_GMM )
  else( GMM_FOUND )
    message( FATAL_ERROR "Could not locate GMM++, set the variable user_GMM_INCLUDE_DIRS for your system in CMakeLists.txt!" )
  endif( GMM_FOUND )

endif( USE_GMM )

# attempt to build the test executable
if( BUILD_TESTS )

  # try to find the OpenGL and GLUT packages
  find_package( OpenGL )
  find_package( GLUT )

  if( OPENGL_FOUND AND GLUT_FOUND )
    include_directories( ${OPENGL_INCLUDE_DIRS} )
    include_directories( ${GLUT_INCLUDE_DIRS} )
    add_definitions( -DTRANSFORMATION_USE_OPENGL )
  else( OPENGL_FOUND AND GLUT_FOUND )
    message("Could not find OpenGL and/or GLUT, disabling tests!" )
  endif( OPENGL_FOUND AND GLUT_FOUND )

endif( BUILD_TESTS )

# create the library target
add_library( transformation transformation.h transformation.cpp )
set_target_properties( transformation PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY         ${CMAKE_SOURCE_DIR}/../bin
  ARCHIVE_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/../bin 
  ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/../bin
)

#create the test executable target
if( BUILD_TESTS )
  add_executable( transformation_test test.cpp )
  add_dependencies( transformation_test transformation )
  target_link_libraries( transformation_test transformation ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

  set_target_properties( transformation_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY         ${CMAKE_SOURCE_DIR}/../bin
    RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/../bin 
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/../bin
  )
endif( BUILD_TESTS )



