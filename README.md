transformation

(c) James Gregson 2011-2013
==============

A basic 4x4 transformation class which duplicates much of the functionality of the OpenGL/GLU matrix functions with a very similar interface.  Written to make graphical debugging of geometric code easier, by allowing a uniform set of transformation code to be used. Moderately well tested to assure that the methods provided in the transformation class produce near identical results to the corresponding methods in the transformation class.

Licensed under the MIT license. See the LICENSE file for more information.

Build With CMake, from current directory run:

$> mkdir build

$> cd build

$> cmake -DCMAKE_BUILD_TYPE=Release ../code 

$> make

This should create a binary directory including the library file and a test executable (provided you haven't disabled this option). To build the Doxygen documentation run the following command from the top-level directory

$> doxygen Doxyfile

If you intend to use the transformation::inverse() or transformation::glUnProject() methods, you should probably enable either Eigen or GMM++ support. See the code/CMakeLists.txt file for more information on how to do this.