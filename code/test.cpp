#ifdef __APPLE__
#include<GLUT/glut.h>
#else
#include<GL/glut.h>
#endif

#include"transformation.h"

int main( int argc, char **argv ){
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE );
    glutCreateWindow( "window" );
    glutInitWindowSize( 128, 128 );
    transformation_test();
    
    return 0;
}