#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#ifdef TRANSFORMATION_USE_EIGEN
#ifdef Success
#undef Success
#endif
#include<Eigen/Dense>
#endif

#ifdef TRANSFORMATION_USE_GMM
#include<gmm/gmm.h>
#endif

/**
 @file transformation.h
 @brief Header file for the OpenGL compatible 4x4 homogeneous transformation library. This
 library allows duplication of the OpenGL transform pipeline via the transformation class
 and its static member equivalents to the OpenGL matrix functions. This allows software to
 use OpenGL style transformations without requiring OpenGL at compile time, or a valid
 OpenGL context at runtime.
 
 There are a few preprocessor flags which allow various options to be set:
 TRANSFORMATION_USE_OPENGL - include the OpenGL headers to allow the function
 transformation_test() to perform regression testing against OpenGL matrix functions
 
 One of the following should be defined:
 TRANSFORMATION_USE_EIGEN - include the Eigen linear algebra library to implement
 transformation::inverse(), which is used to invert transformations and in 
 transformation::gluUnProject()
 
 TRANSFORMATION_USE_GMM - include the gmm++ linear algebra library to implement
 transformation::inverse(), which is used to invert transformations and in 
 transformation::gluUnProject()
 */

/**
 @brief Unit test of the transformation class, to ensure that the library is functional
 and replicates the OpenGL functionality properly.
 
 For this test to function properly, the compiler flag TRANSFORMATION_USE_OPENGL must be
 defined, otherwise the OpenGL headers will not be included.  This function is the only
 reason why the OpenGL headers need to be included at all.
 */
bool transformation_test();

/**
 @brief A transformation class allowing general 4x4 linear transformations.  Includes a number
 of static member functions which duplicate (much of the) functionality of OpenGL matrices.
 Functions to convert between OpenGL and the transformation class are provided as well.
 
 This allows code to written using this class to be used with OpenGL applications and in
 general code, since there is no requirement for a valid OpenGL context, unlike hijacking
 the OpenGL matrix stack.
 
 Uses the open-source libraries Eigen or gmm++ to handle the linear algebra functions
 */
class transformation {
private:
	double	m_a[4][4]; /** matrix element, first index is row, second index is column */
public:
    /** default constructor */
	transformation();
    
    /** 
     @brief converts this transformation to OpenGL matrix
     @param[out] val column-major output (i.e. an OpenGL matrix)
     */
    void to_opengl( double *val );
    
    /**
     @brief initalizes this transformation from an OpenGL matrix 
     @param[in] val column-major input from which to set this transformation
     */
    void from_opengl( double *val );
    
#ifdef TRANSFORMATION_USE_EIGEN
    /**
     @brief stores the current transformation in the Eigen::Matrix4d input
     @param[in] mat matrix to store this transformation in
     */
    void to_eigen( Eigen::Matrix4d &mat );
    
    /** 
     @brief sets the current transformation from the input Eigen::Matrix4d input
     @param[in] mat matrix from which to load the transformation
     */
    void from_eigen( Eigen::Matrix4d &mat );
#endif    
    
#ifdef TRANSFORMATION_USE_GMM
    /**
     @brief stores the current transformation in the gmm::dense_matrix<double> input
     @param[in] mat matrix to store this transformation in
     */
    void to_gmm( gmm::dense_matrix<double> &mat );
    
    /** 
     @brief sets the current transformation from the input Eigen::Matrix4d input
     @param[in] mat matrix from which to load the transformation
     */
    void from_gmm( gmm::dense_matrix<double> &mat );
#endif  
    
    /**
     @brief returns the right, up and forward vectors for this transformation
     and the translation
     @param[out] r local "right" vector
     @param[out] u local "up" vector
     @param[out] f local "forward" vector
     @param[out] p positions/translation
     */
    void get_right_up_forward_position( double *r, double *u, double *f, double *p );
   
    /**
     @brief construct a transformation by setting the right, up, forward and 
     translation vectors, with the fourth row [0, 0, 0, 1]
     @param[in] r right vector
     @param[in] u up vector
     @param[in] f forward vector
     @param[in] t translation vector
    */
    void set_right_up_forward_position( double *r, double *u, double *f, double *p );
    
    /**
     @brief performs a Gram-Schmidt orthogonalization on the rotation matrix
     portion of the transformation matrix. This completely ignores (and so
     destroys) any scaling, shearing or perspective projection that has been
     applied, but can be used to help reduce the effect of repeated roundoff
     errors accumulating on repeated operations
    */
    void reorthogonalize_rotation();
    
    /**
     @brief returns a mutable reference to a matrix entry
     @param[in] r row of entry
     @param[in] c column of entry
     */
	double &operator()( int r, int c );
    
    /**
     @brief transforms a point by the current transformation
     and renormalizes the homogeneous coordinate to 1
     @param[in] pin pointer to input point [x,y,z,w]^T
     @param[out] pout pointer to output point [x,y,z,1]^T
     */
    void transform_homogeneous( double *pin, double *pout );
    
    /**
     @brief transforms a point by this transformation
     @param[in] pin pointer to input point [x,y,z]^T, homogeneous coordinate assumed to be 1
     @param[out] pout pointer to output point [x,y,z]^T, homogeneous coordinate set to be 1
     */
	void transform_point( double *pin, double *pout );
    
    /** 
     @brief transforms a vector by this transformation
     @param[in] vin pointer to input vector [x,y,z], homogeneous coordinate assumed to be 0
     @param[out] vout pointer to output vector [x,y,z], homogeneous coordinate set to be 0     
     */
	void transform_vector( double *vin, double *vout );
    
    /** 
     @brief returns the transpose of the current transformation. Does not change this transform.
     */
    transformation transpose();
    
    /**
     @brief returns this transformation post-multiplied by the input transformation. Does not change this transform.
     @param[in] in transformation to multiply the current transform by
     */
	transformation operator*( transformation in );
    
    /**
     @brief returns the inverse of this transform. Does not change this transform.
     */
	transformation inverse();
    
    /**
     @brief returns a transformation with all elements zeroed
     */
    static transformation zero();
    
    /**
     @brief static member function to return an identity matrix. Effectively glLoadIdentity
     */
	static transformation identity();
    
    /**
     @brief static member function to return a scaling matrix
     @param[in] x scale factor for x-axis
     @param[in] y scale factor for y-axis
     @param[in] z scale factor for z-axis
     */
    static transformation scale( double x, double y, double z );
    
    /**
     @brief static member function to return a translation matrix
     @param[in] x offset for x-axis
     @param[in] y offset for y-axis
     @param[in] z offset for z-axis
     */
    static transformation translate( double x, double y, double z );
    
    /**
     @brief static member function to return a rotation about the x-axis
     @param[in] theta rotation angle in radians
     */
	static transformation rotate_x( double theta );
    
    /** 
     @brief static member function to return a rotation about the y-axis
     @param[in] theta rotation angle in radians
     */
	static transformation rotate_y( double theta );
    
    /**
     @brief static member function to return a rotation about the z-axis
     @param[in] theta rotation angle in radians
     */
	static transformation rotate_z( double theta );
    
    /**
     @brief static member function to return a rotation about an arbitrary axis.
     @param[in] theta angle to rotate by, in RADIANS
     @param[in] x axis vector x component
     @param[in] y axis vector y component
     @param[in] z axis vector z component
     */
    static transformation rotate_axis( double theta, double x, double y, double z );
    
    /**
     @brief computes the closest rigid transformation between two corresponding
     point sets via singular-value decomposition. Requires Eigen.
     @param[in] npts number of points in pointsets
     @param[in] p0 packed array of points, i.e. [x0, y0, z0, x1, y1, z1, .... ]
     @param[in] p1 packed array of points
     @return rigid transformation which most closesly maps points in p0 to points in p1
    */
    static transformation best_rigid_transformation( int npts, double *p0, double *p1 );
    
    /**
     @brief static member function to return an identity tranformation. Convenience
     definition that simply wraps transformation::identity()
     */
    static transformation glIdentity();
    
    /**
     @brief static member function to return a scaling matrix. Convenience definition
     that simply wraps transformation::scale()
     
     @param[in] x scale factor for x-axis
     @param[in] y scale factor for y-axis
     @param[in] z scale factor for z-axis
     */
	static transformation glScale( double x, double y, double z );
    
    /**
     @brief static member function to return a translation. Convenience definition
     that simply wraps transformation::translate()
     
     @param[in] x translation along x-axis
     @param[in] y translation along y-axis
     @param[in] z translation along z-axis
     */
	static transformation glTranslate( double x, double y, double z );
    
    /**
     @brief static member function to return a rotation about an arbitrary axis. Convenience
     definition that wraps transformation::rotate_axis(), but handles conversion from degrees
     to radians.
     
     @param[in] theta angle to rotate by, in DEGREES
     @param[in] x axis vector x component
     @param[in] y axis vector y component
     @param[in] z axis vector z component
     */
    static transformation glRotate( double theta, double x, double y, double z );
    
    /**
     @brief static member function to return a orthographic projeciton matrix
     identical to that produced by OpenGL's glOrtho() function
     
     @param[in] left left clipping plane position
     @param[in] right right clipping plane position
     @param[in] bottom bottom clipping plane position
     @param[in] top top clipping plane position
     @param[in] near near clipping plane position
     @param[in] far far clipping plane position
     */
    static transformation glOrtho( double left, double right, double bottom, double top, double nearVal, double farVal );
    
    /**
     @brief static member fuction to return a perspective projection matrix 
     identical to that produced by OpenGL's glFrustum() function
     @param[in] left left clipping plane coordinate
     @param[in] right right clipping plane coordinate
     @param[in] bottom bottom clipping plane coordinate
     @param[in] top top clipping plane coordinate
     @param[in] nearVal near clipping plane distance, must be positive
     @param[in] farVal far clipping plane distance, must be positive
     */
    static transformation glFrustum( double left, double right, double bottom, double top, double nearVal, double farVal );
    
    /**
     @brief implementation of the glProject() function
     @param[in] objX x-coordinate of point on object
     @param[in] objY y-coordinate of point on object
     @param[in] objZ z-coordinate of point on object
     @param[in] model OpenGL formatted 'modelview' matrix
     @param[in] proj  OpenGL formatted projection matrix
     @param[in] view  OpenGL viewport [left, right, width, height]
     @param[out] winX output window x-coordinate
     @param[out] winY output window y-coordinate
     @param[out] winZ output window z-coordinate
     */    
    static void gluProject( double objX, double objY, double objZ, double *model, double *proj, int *view, double *winx, double *winY, double *winZ );
    
    /**
     @brief implementation of the gluUnProject() function
     @param[in] winX output window x-coordinate
     @param[in] winY output window y-coordinate
     @param[in] winZ output window z-coordinate
     @param[in] model OpenGL formatted 'modelview' matrix
     @param[in] proj  OpenGL formatted projection matrix
     @param[in] view  OpenGL viewport [left, right, width, height]
     @param[out] objX x-coordinate of point on object
     @param[out] objY y-coordinate of point on object
     @param[out] objZ z-coordinate of point on object
     */   
    static void gluUnProject( double winX, double winY, double winZ, double *model, double *proj, int *view, double *objX, double *objY, double *objZ );
    
    /**
     @brief constructs a transformation corresponding to an OpenGL modelview matrix
     which places the camera at [eyeX, eyeY, eyeZ]^T looking at [centerX, centerY, centerZ]^T
     with up vector [upX, upY, upZ]^T
     @param[in] eyeX eye position x-coordinate
     @param[in] eyeY eye position y-coordinate
     @param[in] eyeZ eye position z-coordinate
     @param[in] centerX viewed position x-coordinate
     @param[in] centerY viewed position y-coordinate
     @param[in] centerZ viewed position z-coordinate
     @param[in] upX up-vector x-coordinate
     @param[in] upY up-vector y-coordinate
     @param[in] upZ up-vector z-coordinate     
     */
    static transformation gluLookAt( double eyeX, double eyeY, double eyeZ, double centerX, double centerY, double centerZ, double upX, double upY, double upZ );
    
    /**
     @brief constructs the transformation corresponding to an OpenGL projection matrix
     generated by gluPerspective.
     @param[in] fovy field of view angle (degrees) in the y-direction
     @param[in] aspect aspect ratio, which determines the field of view angle in the x-direction
     @param[in] zNear near clip plane distance
     @param[in] zFar far clip plane distance
     */
    static transformation gluPerspective( double fovy, double aspect, double zNear, double zFar );
};

#endif
