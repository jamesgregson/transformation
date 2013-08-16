
#include"transformation.h"

// Do not include the OpenGL and GLU headers unless the flag 
// TRANSFORMATION_USE_OPENGL is defined.  This is ONLY used 
// by the library to run the regression test against the system
// OpenGL library
#ifdef TRANSFORMATION_USE_OPENGL
#ifdef __APPLE__
#include<OpenGL/GL.h>
#include<OpenGL/GLU.h>
#else
#include<GL/gl.h>
#include<GL/glu.h>
#endif
#endif

#ifdef TRANSFORMATION_USE_EIGEN
// includes required if the Eigen linear algebra library is
// to be used
#ifdef Success
#undef Success
#endif
#include<Eigen/Dense>
#elif TRANSFORMATION_USE_GMM
// includes required if the gmm++ linear algebra library is to be used
#include<gmm/gmm.h>
#else
// print a warning message during compilation about the lack of
// transformation::inverse() and the associated method 
// transformation::gluUnProject(.)
#warning No linear algebra library specified, library will have very slow transformation::inverse() and transformation::glUnProject()
#warning Define either TRANSFORMATION_USE_GMM (to use gmm++) or TRANSFORMATION_USE_EIGEN (to use Eigen) to specify a linear algebra package
#endif

#include<cmath>
#include<iostream>

#ifdef TRANSFORMATION_USE_OPENGL
void transformation_glErr(){
    int err = glGetError();
    if( err == GL_NO_ERROR ) return;
    std::cout << "OpenGL Error: ";
    switch( err ){
        case GL_INVALID_ENUM:
            std::cout << "GL_INVALID_ENUM";
            break;
        case GL_INVALID_VALUE:
            std::cout << "GL_INVALID_VALUE";
            break;
        case GL_INVALID_OPERATION:
            std::cout << "GL_INVALID_OPERATION";
            break;
        case GL_STACK_OVERFLOW:
            std::cout << "GL_STACK_OVERFLOW";
            break;
        case GL_STACK_UNDERFLOW:
            std::cout << "GL_STACK_UNDERFLOW";
            break;
        case GL_OUT_OF_MEMORY:
            std::cout << "GL_OUT_OF_MEMORY";
            break;
        case GL_TABLE_TOO_LARGE:
            std::cout << "GL_TABLE_TOO_LARGE";
            break;
        default:
            std::cout << "unknown error";
            break;
    }
    std::cout << std::endl;
    err=1;
}
#endif

transformation::transformation(){
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            m_a[i][j] = 0.0;
        }
    }
}

void transformation::to_opengl( double *val ){
    val[0]=m_a[0][0]; val[4]=m_a[0][1]; val[ 8]=m_a[0][2]; val[12]=m_a[0][3];
    val[1]=m_a[1][0]; val[5]=m_a[1][1]; val[ 9]=m_a[1][2]; val[13]=m_a[1][3];
    val[2]=m_a[2][0]; val[6]=m_a[2][1]; val[10]=m_a[2][2]; val[14]=m_a[2][3];
    val[3]=m_a[3][0]; val[7]=m_a[3][1]; val[11]=m_a[3][2]; val[15]=m_a[3][3];
}

void transformation::from_opengl( double *val ){
    m_a[0][0]=val[0]; m_a[0][1]=val[4]; m_a[0][2]=val[ 8]; m_a[0][3]=val[12];
    m_a[1][0]=val[1]; m_a[1][1]=val[5]; m_a[1][2]=val[ 9]; m_a[1][3]=val[13];
    m_a[2][0]=val[2]; m_a[2][1]=val[6]; m_a[2][2]=val[10]; m_a[2][3]=val[14];
    m_a[3][0]=val[3]; m_a[3][1]=val[7]; m_a[3][2]=val[11]; m_a[3][3]=val[15];
}

#ifdef TRANSFORMATION_USE_EIGEN
void transformation::to_eigen( Eigen::Matrix4d &mat ){
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            mat(i,j) = m_a[i][j];
        }
    }
}

void transformation::from_eigen( Eigen::Matrix4d &mat ){
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            m_a[i][j] = mat(i,j);
        }
    }
}
#endif    

#ifdef TRANSFORMATION_USE_GMM
void transformation::to_gmm( gmm::dense_matrix<double> &mat ){
    mat.resize(4,4);
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            mat(i,j) = m_a[i][j];
        }
    }
}

void transformation::from_gmm( gmm::dense_matrix<double> &mat ){
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            m_a[i][j] = mat(i,j);
        }
    }
}
#endif 

void transformation::get_right_up_forward_position( double *r, double *u, double *f, double *p ){
    r[0] = m_a[0][0]; u[0] = m_a[0][1]; f[0] = m_a[0][2]; p[0] = m_a[0][3];
    r[1] = m_a[1][0]; u[1] = m_a[1][1]; f[1] = m_a[1][2]; p[1] = m_a[1][3];
    r[2] = m_a[2][0]; u[2] = m_a[2][1]; f[2] = m_a[2][2]; p[2] = m_a[2][3];
}

double &transformation::operator()( int r, int c ){
    return m_a[r][c];
}

void transformation::transform_homogeneous( double *pin, double *pout ){
    double w = m_a[3][0]*pin[0] + m_a[3][1]*pin[1] + m_a[3][2]*pin[2] + m_a[3][3]*pin[3];
    pout[0] = (m_a[0][0]*pin[0] + m_a[0][1]*pin[1] + m_a[0][2]*pin[2] + m_a[0][3]*pin[3])/w;
    pout[1] = (m_a[1][0]*pin[0] + m_a[1][1]*pin[1] + m_a[1][2]*pin[2] + m_a[1][3]*pin[3])/w;
    pout[2] = (m_a[2][0]*pin[0] + m_a[2][1]*pin[1] + m_a[2][2]*pin[2] + m_a[2][3]*pin[3])/w;
    pout[3] = 1.0;
}

void transformation::transform_point( double *pin, double *pout ){
    double w = m_a[3][0]*pin[0] + m_a[3][1]*pin[1] + m_a[3][2]*pin[2] + m_a[3][3];
    pout[0] = (m_a[0][0]*pin[0] + m_a[0][1]*pin[1] + m_a[0][2]*pin[2] + m_a[0][3])/w;
    pout[1] = (m_a[1][0]*pin[0] + m_a[1][1]*pin[1] + m_a[1][2]*pin[2] + m_a[1][3])/w;
    pout[2] = (m_a[2][0]*pin[0] + m_a[2][1]*pin[1] + m_a[2][2]*pin[2] + m_a[2][3])/w;
}

void transformation::transform_vector( double *vin, double *vout ){
    vout[0] = (m_a[0][0]*vin[0] + m_a[0][1]*vin[1] + m_a[0][2]*vin[2]);
    vout[1] = (m_a[1][0]*vin[0] + m_a[1][1]*vin[1] + m_a[1][2]*vin[2]);
    vout[2] = (m_a[2][0]*vin[0] + m_a[2][1]*vin[1] + m_a[2][2]*vin[2]);		
}

transformation transformation::transpose(){
    transformation tmp;
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            tmp(i,j) = (*this)(i,j);
        }
    }
    return tmp;
}

transformation transformation::operator*( transformation in ){
    transformation tmp;
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            tmp(i,j) = 0.0;
            for( int k=0; k<4; k++ ){
                tmp(i,j) += (*this)(i,k)*in(k,j);	
            }
        }
    }
    return tmp;
}

transformation transformation::inverse(){
    transformation tmp;
#ifdef TRANSFORMATION_USE_EIGEN
    // codepath for using eigen
    Eigen::Matrix4d A, iA;
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            A(i,j) = m_a[i][j];
        }
    }
    iA = A.inverse();
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            tmp(i,j) = iA(i,j);
        }
    }
#elif TRANSFORMATION_USE_GMM
    // codepath for using gmm++
    gmm::dense_matrix<double> A(4,4);
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            A(i,j) = m_a[i][j];
        }
    }
    gmm::lu_inverse(A);
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            tmp(i,j) = A(i,j);
        }
    }
#else
    // first copy to a flat layout
    double m[16];
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            m[i+j*4] = m_a[i][j];
        }
    }
    // now invert, using code directly cribbed from Rodolphe Vaillant's homepage: http://www.irit.fr/~Rodolphe.Vaillant/?e=7
    double inv[16], det;
    inv[ 0] =  m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
    inv[ 4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
    inv[ 8] =  m[4] * m[ 9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[ 9];
    inv[12] = -m[4] * m[ 9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[ 9];
    inv[ 1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
    inv[ 5] =  m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
    inv[ 9] = -m[0] * m[ 9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[ 9];
    inv[13] =  m[0] * m[ 9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[ 9];
    inv[ 2] =  m[1] * m[ 6] * m[15] - m[1] * m[ 7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[ 7] - m[13] * m[3] * m[ 6];
    inv[ 6] = -m[0] * m[ 6] * m[15] + m[0] * m[ 7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[ 7] + m[12] * m[3] * m[ 6];
    inv[10] =  m[0] * m[ 5] * m[15] - m[0] * m[ 7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[ 7] - m[12] * m[3] * m[ 5];
    inv[14] = -m[0] * m[ 5] * m[14] + m[0] * m[ 6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[ 6] + m[12] * m[2] * m[ 5];
    inv[ 3] = -m[1] * m[ 6] * m[11] + m[1] * m[ 7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[ 9] * m[2] * m[ 7] + m[ 9] * m[3] * m[ 6];
    inv[ 7] =  m[0] * m[ 6] * m[11] - m[0] * m[ 7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[ 8] * m[2] * m[ 7] - m[ 8] * m[3] * m[ 6];
    inv[11] = -m[0] * m[ 5] * m[11] + m[0] * m[ 7] * m[ 9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[ 9] - m[ 8] * m[1] * m[ 7] + m[ 8] * m[3] * m[ 5];
    inv[15] =  m[0] * m[ 5] * m[10] - m[0] * m[ 6] * m[ 9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[ 9] + m[ 8] * m[1] * m[ 6] - m[ 8] * m[2] * m[ 5];
    
    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
    
    // error checking??? Nope.
    //if(det == 0)
    //    return false;
    
    det = 1.0 / det;

    // copy from flat layout to returned inverse
    for( int i=0; i<4; i++ ){
        for( int j=0; j<4; j++ ){
            tmp(i,j) = inv[i+j*4]*det;
        }
    }
    
#endif
    return tmp;
}

transformation transformation::zero(){
    transformation tmp;
    tmp(0,0)=tmp(0,1)=tmp(0,2)=tmp(0,3)=0.0;
    tmp(1,0)=tmp(1,1)=tmp(1,2)=tmp(1,3)=0.0;
    tmp(2,0)=tmp(2,1)=tmp(2,2)=tmp(2,3)=0.0;
    tmp(3,0)=tmp(3,1)=tmp(3,2)=tmp(3,3)=0.0;
    return tmp;
}

transformation transformation::identity(){
    transformation tmp;
    tmp(0,0)=tmp(1,1)=tmp(2,2)=tmp(3,3)=1.0;
    return tmp;
}

transformation transformation::rotate_x( double theta ){
    double s = sin( theta );
    double c = cos( theta );
    transformation tmp;
    tmp(0,0)=tmp(1,1)=tmp(2,2)=tmp(3,3)=1.0;
    tmp(1,1) =  c; tmp(1,2) = -s;
    tmp(2,1) =  s; tmp(2,2) =  c;
    return tmp;
}

transformation transformation::rotate_y( double theta ){
    double s = sin( theta );
    double c = cos( theta );
    transformation tmp;
    tmp(0,0)=tmp(1,1)=tmp(2,2)=tmp(3,3)=1.0;
    tmp(0,0) =  c; tmp(0,2) = s;
    tmp(2,0) = -s; tmp(2,2) = c;
    return tmp;
}

transformation transformation::rotate_z( double theta ){
    double s = sin( theta );
    double c = cos( theta );
    transformation tmp;
    tmp(0,0)=tmp(1,1)=tmp(2,2)=tmp(3,3)=1.0;
    tmp(0,0) = c; tmp(0,1) = -s;
    tmp(1,0) = s; tmp(1,1) = c;
    return tmp;
}

transformation transformation::scale( double x, double y, double z ){
    transformation tmp;
    tmp(0,0) = x;
    tmp(1,1) = y;
    tmp(2,2) = z;
    tmp(3,3) = 1.0;
    return tmp;
}

transformation transformation::translate( double x, double y, double z ){
    transformation tmp;
    tmp(0,0)=tmp(1,1)=tmp(2,2)=tmp(3,3)=1.0;
    tmp(0,3) = x;
    tmp(1,3) = y;
    tmp(2,3) = z;
    return tmp;
}

transformation transformation::rotate_axis( double theta, double x, double y, double z ){
    double qx, qy, qz, qw;
    double L = sqrt( x*x+y*y+z*z );
    double c=cos(theta/2.0), s=sin(theta/2.0);
    qx = s*x/L; 
    qy = s*y/L;
    qz = s*z/L;
    qw = c;
    transformation tmp;
    tmp(0,0) = 1.0-2.0*(qy*qy+qz*qz); tmp(1,0) = 2.0*(qx*qy+qw*qz);     tmp(2,0) = 2.0*(qx*qz-qw*qy);     tmp(3,0) = 0.0;
    tmp(0,1) = 2.0*(qx*qy-qw*qz);     tmp(1,1) = 1.0-2.0*(qx*qx+qz*qz); tmp(2,1) = 2.0*(qy*qz+qw*qx);     tmp(3,1) = 0.0;
    tmp(0,2) = 2.0*(qx*qz+qw*qy);     tmp(1,2) = 2.0*(qy*qz-qw*qx);     tmp(2,2) = 1.0-2.0*(qx*qx+qy*qy); tmp(3,2) = 0.0;
    tmp(0,3) = 0.0;                   tmp(1,3) = 0.0;                   tmp(2,3) = 0.0;                   tmp(3,3) = 1.0;
    return tmp;
}

transformation transformation::glIdentity(){
    return transformation::identity();
}

transformation transformation::glScale( double x, double y, double z ){
    return transformation::scale( x, y, z );
}

transformation transformation::glTranslate( double x, double y, double z ){
    return transformation::translate( x, y, z );
}

transformation transformation::glRotate( double theta, double x, double y, double z ){
    return transformation::rotate_axis( theta*M_PI/180.0, x, y, z );
}

transformation transformation::glOrtho( double left, double right, double bottom, double top, double nearVal, double farVal ){
    transformation tmp;
    double tx, ty, tz;
    tx = -(right+left)/(right-left);
    ty = -(top+bottom)/(top-bottom);
    tz = -(farVal+nearVal)/(farVal-nearVal);
    tmp(0,0) = 2.0/(right-left);
    tmp(1,1) = 2.0/(top-bottom);
    tmp(2,2) =-2.0/(farVal-nearVal);
    tmp(3,3) = 1.0;
    tmp(0,3) = tx;
    tmp(1,3) = ty;
    tmp(2,3) = tz;
    return tmp;
}


transformation transformation::glFrustum( double left, double right, double bottom, double top, double nearVal, double farVal ){
    double A, B, C, D;
    transformation tmp = transformation::identity();
    tmp(0,0) = 2.0*nearVal/(right-left);
    tmp(1,1) = 2.0*nearVal/(top-bottom);
    tmp(3,2) = -1.0;
    tmp(3,3) = 0.0;

    A =  (right+left)/(right-left);
    B =  (top+bottom)/(top-bottom);
    C = -(farVal+nearVal)/(farVal-nearVal);
    D = -(2*farVal*nearVal)/(farVal-nearVal);

    tmp(0,2) = A;
    tmp(1,2) = B;
    tmp(2,2) = C;
    tmp(2,3) = D;
    return tmp;
}

void transformation::gluProject( double objX, double objY, double objZ, double *model, double *proj, int *view, double *winX, double *winY, double *winZ ){
    transformation M, P, tmp;
    double pout[4], pin[] = { objX, objY, objZ, 1.0 };
    M.from_opengl( model );
    P.from_opengl( proj );
    tmp = P*M;
    tmp.transform_homogeneous( pin, pout );
    *winX = view[0] + 0.5*double(view[2])*( pout[0]+1.0 );
    *winY = view[1] + 0.5*double(view[3])*( pout[1]+1.0 );
    *winZ = 0.5*(pout[2]+1.0);
}

void transformation::gluUnProject( double winX, double winY, double winZ, double *model, double *proj, int *view, double *objX, double *objY, double *objZ ){
    transformation M, P, tmp;
    double pout[4], pin[] = { 2.0*(winX-double(view[0]))/double(view[2])-1.0, 2.0*(winY-double(view[1]))/double(view[3])-1.0, 2.0*winZ-1.0, 1.0 };
    M.from_opengl( model );
    P.from_opengl( proj );
    tmp = (P*M).inverse();
    tmp.transform_homogeneous( pin, pout );
    *objX = pout[0];
    *objY = pout[1];
    *objZ = pout[2];
}

transformation transformation::gluLookAt( double eyeX, double eyeY, double eyeZ, double centerX, double centerY, double centerZ, double upX, double upY, double upZ ){
    transformation tmp;
    double u[] = { upX, upY, upZ };
    double f[] = { centerX-eyeX, centerY-eyeY, centerZ-eyeZ };
    double L = sqrt( f[0]*f[0] + f[1]*f[1] + f[2]*f[2] );
    f[0]/=L; f[1]/=L; f[2]/=L;
    double r[] = { f[1]*u[2]-f[2]*u[1], f[2]*u[0]-f[0]*u[2], f[0]*u[1]-f[1]*u[0] };
    L = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
    r[0] /= L; r[1] /= L; r[2] /= L;
    u[0] = r[1]*f[2]-r[2]*f[1];
    u[1] = r[2]*f[0]-r[0]*f[2];
    u[2] = r[0]*f[1]-r[1]*f[0];
    tmp(0,0) = r[0]; tmp(0,1) = r[1]; tmp(0,2) = r[2];
    tmp(1,0) = u[0]; tmp(1,1) = u[1]; tmp(1,2) = u[2];
    tmp(2,0) =-f[0]; tmp(2,1) =-f[1]; tmp(2,2) =-f[2];
    tmp(3,3) = 1.0;
    return tmp * transformation::glTranslate( -eyeX, -eyeY, -eyeZ );
}

transformation transformation::gluPerspective( double fovy, double aspect, double zNear, double zFar ){
    transformation tmp;
    double f = 1.0/tan(fovy*M_PI/360.0);
    tmp(0,0) = f/aspect;
    tmp(1,1) = f;
    tmp(2,2) = (zFar+zNear)/(zNear-zFar); tmp(2,3) = 2.0*zFar*zNear/(zNear-zFar);
    tmp(3,2) = -1.0;
    return tmp;
}


// transformation class test suite
#ifndef TRANSFORMATION_USE_OPENGL
// TRANSFORMATION_USE_OPENGL not defined, print an error message and return failed
bool transformation_test(){
    std::cout << "Test suite cannot be run since TRANSFORMATION_USE_OPENGL is not defined." << std::endl;
    std::cout << "  returning FAILED" << std::endl;
    return false;
}
#else
bool transformation_test(){
    double eps    = 1e-6;
    int n_samples = 1000;
    
    std::cout << "running transformation_test()..." << std::endl;
    
    // test that a transformation loaded from an OpenGL matrix can
    // be re-output as an OpenGL matrix without modification 
    std::cout << "\ttesting transformation::to_opengl() & transformation::from_opengl()...";
    {
        double a[16], b[16];
        transformation tmp;
        for( int i=0; i<16; i++ ){
            a[i] = drand48();
        }
        tmp.from_opengl( a );
        tmp.to_opengl( b );
        for( int i=0; i<16; i++ ){
            if( fabs( a[i]-b[i] ) > eps ){
                std::cout << "failed." << std::endl;
            }
        }
    }
    std::cout << "passed." << std::endl;
    
    // test that transformation::glTranslate() produces the same
    // matrices as the OpenGL glTranslate()
    std::cout << "\ttesting transformation::glTranslate()...";
    {
        GLint mode;
        glGetIntegerv( GL_MATRIX_MODE, &mode ); transformation_glErr();
        glMatrixMode( GL_MODELVIEW ); transformation_glErr();
        glPushMatrix(); transformation_glErr();
        for( int i=0; i<n_samples; i++ ){
            transformation tmp;
            double x, y, z, ogl[16], mat[16];
            glLoadIdentity(); transformation_glErr();
            x = drand48()-0.5;
            y = drand48()-0.5;
            z = drand48()-0.5;
            glTranslated( x, y, z ); transformation_glErr();
            glGetDoublev( GL_MODELVIEW_MATRIX, ogl ); transformation_glErr();
            tmp = transformation::glTranslate( x, y, z );
            tmp.to_opengl( mat );
            for( int j=0; j<16; j++ ){
                if( fabs(mat[j]-ogl[j]) > eps){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
        glPopMatrix(); transformation_glErr();
        glMatrixMode( mode ); transformation_glErr();
    }
    std::cout << "passed." << std::endl;
    
    // test that transformation::glScale() produces the same
    // matrices as the OpenGL glScale()
    std::cout << "\ttesting transformation::glScale()...";
    {
        int mode;
        glGetIntegerv( GL_MATRIX_MODE, &mode ); transformation_glErr();
        glMatrixMode( GL_MODELVIEW ); transformation_glErr();
        glPushMatrix(); transformation_glErr();
        for( int i=0; i<n_samples; i++ ){
            transformation tmp;
            double x, y, z, ogl[16], mat[16];
            glLoadIdentity(); transformation_glErr();
            x = drand48()-0.5;
            y = drand48()-0.5;
            z = drand48()-0.5;
            glScaled( x, y, z ); transformation_glErr();
            glGetDoublev( GL_MODELVIEW_MATRIX, ogl ); transformation_glErr();
            tmp = transformation::glScale( x, y, z );
            tmp.to_opengl( mat );
            for( int j=0; j<16; j++ ){
                if( fabs(mat[j]-ogl[j]) > eps){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
        glPopMatrix(); transformation_glErr();
        glMatrixMode( mode ); transformation_glErr();
    }
    std::cout << "passed." << std::endl;
    
    // test that transformation::glRotate() produces the same
    // matrices as the OpenGL glRotate()
    std::cout << "\ttesting transformation::glRotate()...";
    {
        int mode;
        glGetIntegerv( GL_MATRIX_MODE, &mode ); transformation_glErr();
        glMatrixMode( GL_MODELVIEW ); transformation_glErr();
        glPushMatrix(); transformation_glErr();
        for( int i=0; i<n_samples; i++ ){
            transformation tmp;
            double theta, x, y, z, L, ogl[16], mat[16];
            glLoadIdentity(); transformation_glErr();
            theta = drand48()*360.0;
            x = drand48()-0.5;
            y = drand48()-0.5;
            z = drand48()-0.5;
            L = sqrt(x*x+y*y+z*z);
            x /= L; y /= L; z /= L;
            glRotated( theta, x, y, z ); transformation_glErr();
            glGetDoublev( GL_MODELVIEW_MATRIX, ogl ); transformation_glErr();
            tmp = transformation::glRotate( theta, x, y, z );
            tmp.to_opengl( mat );
            for( int j=0; j<16; j++ ){
                if( fabs(mat[j]-ogl[j]) > eps){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
        glPopMatrix(); transformation_glErr();
        glMatrixMode( mode ); transformation_glErr();
    }
    std::cout << "passed." << std::endl;
    
    // test that transformation::glOrtho() produces the same
    // matrices as the OpenGL glOrtho
    std::cout << "\ttesting transformation::glOrtho()...";
    {
        int mode;
        glGetIntegerv( GL_MATRIX_MODE, &mode ); transformation_glErr();
        glMatrixMode( GL_MODELVIEW ); transformation_glErr();
        glPushMatrix(); transformation_glErr();
        for( int i=0; i<n_samples; i++ ){
            transformation tmp;
            double ogl[16], mat[16], L, R, B, T, N, F;
            
            L = drand48();
            R = L + 1e-2 + drand48();
            B = drand48();
            T = B + 1e-2 + drand48();
            N = 1e-2 + drand48();
            F = N + 1e-2 + drand48();
            
            glLoadIdentity(); transformation_glErr();
            glOrtho( L, R, B, T, N, F ); transformation_glErr();
            glGetDoublev( GL_MODELVIEW_MATRIX, ogl ); transformation_glErr();
            tmp = transformation::glOrtho( L, R, B, T, N, F );
            tmp.to_opengl( mat );
            for( int j=0; j<16; j++ ){
                if( fabs(mat[j]-ogl[j]) > eps*100.0){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
        glPopMatrix(); transformation_glErr();
        glMatrixMode( mode ); transformation_glErr();       
    }
    std::cout << "passed." << std::endl;
    
    // test that transformation::glFrustum() produces the same
    // matrices as the OpenGL glFrustum
    std::cout << "\ttesting transformation::glFrustum()...";
    {
        int mode;
        glGetIntegerv( GL_MATRIX_MODE, &mode ); transformation_glErr();
        glMatrixMode( GL_PROJECTION ); transformation_glErr();
        glPushMatrix(); transformation_glErr();
        for( int i=0; i<n_samples; i++ ){
            transformation tmp;
            double ogl[16], mat[16], L, R, B, T, N, F;
            
            L = drand48();
            R = L + 1e-2 + drand48();
            B = drand48();
            T = B + 1e-2 + drand48();
            N = 1e-2 + drand48();
            F = N + 1e-2 + drand48();
            
            glLoadIdentity(); transformation_glErr();
            glFrustum( L, R, B, T, N, F ); transformation_glErr();
            glGetDoublev( GL_PROJECTION_MATRIX, ogl ); transformation_glErr();
            tmp = transformation::glFrustum( L, R, B, T, N, F );
            tmp.to_opengl( mat );
            for( int j=0; j<16; j++ ){
                if( fabs(mat[j]-ogl[j]) > eps*10.0){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
        glPopMatrix(); transformation_glErr();
        glMatrixMode( mode ); transformation_glErr();       
    }
    std::cout << "passed" << std::endl;
    
    // test that transformation::gluProject() produces the same
    // matrices as the OpenGL glProject
    std::cout << "\ttesting transformation::glProject()...";
    {
        for( int i=0; i<n_samples; i++ ){
            transformation M, P;
            int view[] = { drand48()*1000, drand48()*1000, drand48()*1000, drand48()*1000 };
            double model[16], proj[16], vtx[3], ogl[3];
            double pnt[] = { drand48(), drand48(), drand48() };
            
            double L, R, B, T, N, F;
            L = drand48();
            R = L + 1e-2 + drand48();
            B = drand48();
            T = B + 1e-2 + drand48();
            N = 1e-2 + drand48();
            F = N + 1e-2 + drand48();
            
            M = transformation::glRotate( drand48()*360.0, drand48()-0.5, drand48()-0.5, drand48()-0.5 );
            P = transformation::glFrustum( L, R, B, T, N, F );
            
            M.to_opengl( model );
            P.to_opengl( proj );
            
            gluProject( pnt[0], pnt[1], pnt[2], model, proj, view, &ogl[0], &ogl[1], &ogl[2] ); transformation_glErr();
            transformation::gluProject( pnt[0], pnt[1], pnt[2], model, proj, view, &vtx[0], &vtx[1], &vtx[2] );
            for( int j=0; j<3; j++ ){
                if( fabs( vtx[j]-ogl[j] ) > eps ){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
    }
    std::cout << "passed" << std::endl;
    
//#ifdef TRANSFORMATION_USE_GMM | TRANSFORMATION_USE_EIGEN
    // test that transformation::gluUnProject() produces the same
    // matrices as the OpenGL gluUnProject()
    std::cout << "\ttesting transformation::glUnProject()...";
    {
        for( int i=0; i<n_samples; i++ ){
            transformation M, P;
            int view[] = { drand48()*1000, drand48()*1000, drand48()*1000, drand48()*1000 };
            double model[16], proj[16], vtx[3], ogl[3];
            double pnt[] = { drand48()*1000, drand48()*1000, drand48()*1000 };
            
            double L, R, B, T, N, F;
            L = drand48();
            R = L + 1e-2 + drand48();
            B = drand48();
            T = B + 1e-2 + drand48();
            N = 1e-2 + drand48();
            F = N + 1e-2 + drand48();
            
            M = transformation::glRotate( drand48()*360.0, drand48()-0.5, drand48()-0.5, drand48()-0.5 );
            P = transformation::glFrustum( L, R, B, T, N, F );
            
            M.to_opengl( model );
            P.to_opengl( proj );
            
            gluUnProject( pnt[0], pnt[1], pnt[2], model, proj, view, &ogl[0], &ogl[1], &ogl[2] ); transformation_glErr();
            transformation::gluUnProject( pnt[0], pnt[1], pnt[2], model, proj, view, &vtx[0], &vtx[1], &vtx[2] );
            for( int j=0; j<3; j++ ){
                if( fabs( vtx[j]-ogl[j] ) > eps ){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
    }
    std::cout << "passed" << std::endl;
    
    // test that transformation::gluLookAt() produces the same
    // matrices as the OpenGL gluLookAt()
    std::cout << "\ttesting transformation::gluLookAt()...";
    {
        int mode;
        glGetIntegerv( GL_MATRIX_MODE, &mode ); transformation_glErr();
        glMatrixMode( GL_MODELVIEW ); transformation_glErr();
        glPushMatrix(); transformation_glErr();
        for( int i=0; i<n_samples; i++ ){
            double eye[] = { drand48()-0.5, drand48()-0.5, drand48()-0.5 };
            double cen[] = { drand48()-0.5, drand48()-0.5, drand48()-0.5 };
            double up[] = { drand48()-0.5, drand48()-0.5, drand48()-0.5 };
            double ogl[16], mat[16];
            
            glLoadIdentity(); transformation_glErr();
            gluLookAt( eye[0], eye[1], eye[2], cen[0], cen[1], cen[2], up[0], up[1], up[2] ); transformation_glErr();
            glGetDoublev( GL_MODELVIEW_MATRIX, ogl ); transformation_glErr();
            transformation tmp = transformation::gluLookAt( eye[0], eye[1], eye[2], cen[0], cen[1], cen[2], up[0], up[1], up[2] );
            tmp.to_opengl(mat);
            for( int j=0; j<16; j++ ){
                if( fabs( mat[j]-ogl[j] ) > eps ){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
        glPopMatrix(); transformation_glErr();
        glMatrixMode( mode ); transformation_glErr();
    }
    std::cout << "passed" << std::endl;
    
    // test that transformation::gluPerspective() produces the same
    // matrices as the OpenGL gluPerspective()
    std::cout << "\ttesting transformation::gluPerspective()...";
    {
        int mode;
        glGetIntegerv( GL_MATRIX_MODE, &mode ); transformation_glErr();
        glMatrixMode( GL_MODELVIEW ); transformation_glErr();
        glPushMatrix(); transformation_glErr();
        for( int i=0; i<n_samples; i++ ){
            double fovy = 5.0 + drand48()*85.0;
            double aspect = 0.1 + drand48()*1.8;
            double zNear = 1e-2 + drand48();
            double zFar = zNear + 1e-2 + drand48();
            double ogl[16], mat[16];
            
            glLoadIdentity(); transformation_glErr();
            gluPerspective(fovy, aspect, zNear, zFar ); transformation_glErr();
            glGetDoublev( GL_MODELVIEW_MATRIX, ogl ); transformation_glErr();
            transformation tmp = transformation::gluPerspective(fovy, aspect, zNear, zFar );
            tmp.to_opengl(mat);
            for( int j=0; j<16; j++ ){
                if( fabs( mat[j]-ogl[j] ) > eps*100.0 ){
                    std::cout << "failed." << std::endl;
                    return false;
                }
            }
        }
        glPopMatrix(); transformation_glErr();
        glMatrixMode( mode ); transformation_glErr();
    }
    std::cout << "passed" << std::endl;
    
    std::cout << "test suite passed!" << std::endl;    
    return true;
}
#endif

