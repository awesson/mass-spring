// TinkerToy.cpp : Defines the entry point for the console application.
//

#include "Particle.h"
#include "imageio.h"
#include "System.h"
#include "integrator.h"

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <GLUT/glut.h>

/* macros */


/* global variables */

static int N;
static float dt, d;
static int dsim;
static int dump_frames;
static int frame_number;

// static Particle *pList;
static std::vector<Particle*> pVector;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int mouse_release[3];
static int mouse_shiftclick[3];
static int omx, omy, mx, my;
static int hmx, hmy;
// flag for whether the mouse has already been registered as being down
static bool clicked;

static std::vector<SpringForce*> forceVector;
static std::vector<RodConstraint*> rodConstVector;
static std::vector<CircularWireConstraint*> wireConstVector;
static Integrator* integrator;
static System* sys = NULL;

/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/

static void free_data ( void )
{
    pVector.clear();
    forceVector.clear();
    wireConstVector.clear();
    rodConstVector.clear();
    delete sys;
    delete integrator;
}

static void clear_data ( void )
{
        sys->get_state( pVector );
	int ii, size = pVector.size();

	for(ii=0; ii<size; ii++){
		pVector[ii]->reset();
	}
}

static void init_system( void )
{
        clicked = false;
        const double dist = 0.1;
        const Vec2f center(-0.5, 0.5);
        const Vec2f x_offset(dist, 0.0);
        const Vec2f y_offset(0.0, -dist);

        // Create an array of 100 particles connected by warp, weft and shear springs.
        // Then connect these to two particles constrained to a circular wire

        // add particles
        for(int i = 0; i < 10; ++i){
            for(int k = 0; k < 10; ++k){
                pVector.push_back(new Particle(center + k*x_offset + i*y_offset, 1.f, i*10 + k));
            }
        }

        // the two anchor particles
        pVector.push_back(new Particle(center - 2*x_offset - y_offset, 3.f, 100));
        pVector.push_back(new Particle(center + 11*x_offset - y_offset, 3.f, 101));

        // the anchor particles constraints and spring forces
        wireConstVector.push_back(new CircularWireConstraint(pVector[100], center - 3*x_offset - y_offset, dist));
        wireConstVector.push_back(new CircularWireConstraint(pVector[101], center + 12*x_offset - y_offset, dist));
        forceVector.push_back(new SpringForce(pVector[100], pVector[0], 2*dist, 10.f, 1.f));
        forceVector.push_back(new SpringForce(pVector[101], pVector[9], 2*dist, 10.f, 1.f));

        // add a rod constraint or spring force between alternating particles in the first row
        for(int i = 0; i < 9; ++i){
            if(i%2 == 0)
                rodConstVector.push_back(new RodConstraint(pVector[i], pVector[i+1], dist));
            else
                forceVector.push_back(new SpringForce(pVector[i], pVector[i+1], dist, 4.f, 1.f));
        }

        for(int i = 0; i < 9; ++i){
            for(int k = 0; k < 10; ++k){
                if(k != 9){
                    // weft springs
                    forceVector.push_back(new SpringForce(pVector[10*i + k + 10], pVector[10*i + k + 11], dist, 4.f, 1.f));
                    // shear springs
                    forceVector.push_back(new SpringForce(pVector[10*i + k + 10], pVector[10*i + k + 1], sqrt(2)*dist, 4.f, 1.0f));
                }
                if(k != 0){
                    // shear springs
                    forceVector.push_back(new SpringForce(pVector[10*i + k + 10], pVector[10*i + k - 1], sqrt(2)*dist, 4.f, 1.f));
                }
                // warp springs
                forceVector.push_back(new SpringForce(pVector[10*i + k + 10], pVector[10*i + k], dist, 4.f, 1.f));

            }
        }
	
        sys = new System(pVector, forceVector, wireConstVector, rodConstVector);
}

/*
----------------------------------------------------------------------
OpenGL specific drawing routines
----------------------------------------------------------------------
*/

static void pre_display ( void )
{
	glViewport ( 0, 0, win_x, win_y );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();
	gluOrtho2D ( -1.0, 1.0, -1.0, 1.0 );
	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
	// Write frames if necessary.
	if (dump_frames) {
                const int FRAME_INTERVAL = 24;
		if ((frame_number % FRAME_INTERVAL) == 0) {
			const unsigned int w = glutGet(GLUT_WINDOW_WIDTH);
			const unsigned int h = glutGet(GLUT_WINDOW_HEIGHT);
                        unsigned char * buffer = (unsigned char *) malloc(w * h * 4 * sizeof(unsigned char));
			if (!buffer)
				exit(-1);
			// glRasterPos2i(0, 0);
			glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
			char filename[13];
			sprintf(filename, "img%.5i.png", frame_number / FRAME_INTERVAL);
			printf("Dumped %s.\n", filename);
			saveImageRGBA(filename, buffer, w, h);
		}
	}
	frame_number++;
	
	glutSwapBuffers ();
}

static void draw_particles ( void )
{
        sys->get_state( pVector );
	int size = pVector.size();

	for(int ii=0; ii< size; ii++)
	{
		pVector[ii]->draw();
	}
}

static void draw_forces ( void )
{
    sys->get_forces( forceVector );
    int size = forceVector.size();

    for(int ii=0; ii< size; ii++)
    {
        forceVector[ii]->draw();
    }
}

static void draw_constraints ( void )
{
    sys->get_rodConst( rodConstVector );
    sys->get_wireConst( wireConstVector );
    int rodSize = rodConstVector.size();
    int wireSize = wireConstVector.size();

    for(int ii=0; ii< rodSize; ii++)
    {
        rodConstVector[ii]->draw();
    }
    for(int ii=0; ii< wireSize; ii++)
    {
        wireConstVector[ii]->draw();
    }
}

/*
----------------------------------------------------------------------
relates mouse movements to tinker toy construction
----------------------------------------------------------------------
*/

static void get_from_UI ()
{
	int i, j;
	// int size, flag;
	int hi, hj;
	// float x, y;

	if ( !mouse_down[0] && !mouse_down[2] && !mouse_release[0] 
		&& !mouse_shiftclick[0] && !mouse_shiftclick[2] ) return;

	i = (int)((       mx /(float)win_x)*N);
	j = (int)(((win_y-my)/(float)win_y)*N);

	if ( i<1 || i>N || j<1 || j>N ) return;

	if ( mouse_down[0] ) {
            if(!clicked){
                // add a particle at the position of the mouse.
                // then add two spring particles to the top of the "cloth"
                std::vector<Particle*> p;
                p.resize(sys->size());
                sys->get_state( p );
                // have to convert the mouse location from pixel coordinates to window coordinates
                p.push_back(new Particle(Vec2f(2.f*mx/(double)win_x - 1.f,
                                               1.f - 2.f*my/(double)win_y), 1.f, sys->size()));
                sys->set_state( p );
                sys->add_springForce(new SpringForce(p[4], p[sys->size() - 1], 0.2f, 0.1f, 0.05f));
                sys->add_springForce(new SpringForce(p[5], p[sys->size() - 1], 0.2f, 0.1f, 0.05f));
                clicked = true;
            } else{
                // reposition the particle to the new location of the mouse
                std::vector<Particle*> p;
                p.resize(sys->size());
                sys->get_state( p );
                p[sys->size() - 1]->Position = Vec2f(2.f*mx/(double)win_x - 1.f,
                                                     1.f - 2.f*my/(double)win_y);
                sys->set_state( p );
            }
	}

	if ( mouse_down[2] ) {

	}

	hi = (int)((       hmx /(float)win_x)*N);
	hj = (int)(((win_y-hmy)/(float)win_y)*N);

        if( !mouse_down[0] && mouse_release[0] ) {
            if(clicked){
                // remove the particle and two springs
                clicked = false;
                std::vector<Particle*> p;
                p.resize(sys->size());
                sys->get_state( p );
                p.pop_back();
                sys->set_state( p );
                sys->pop_springForce();
                sys->pop_springForce();
            }
	}

	omx = mx;
	omy = my;
}

static void remap_GUI()
{
	int ii, size = pVector.size();
	for(ii=0; ii<size; ii++)
	{
		pVector[ii]->Position[0] = pVector[ii]->ConstructPos[0];
		pVector[ii]->Position[1] = pVector[ii]->ConstructPos[1];
	}
}

/*
----------------------------------------------------------------------
GLUT callback routines
----------------------------------------------------------------------
*/

static void key_func ( unsigned char key, int x, int y )
{
	switch ( key )
	{

        case ' ':
                dsim = !dsim;

	case 'c':
	case 'C':
		clear_data ();
		break;

	case 'd':
	case 'D':
		dump_frames = !dump_frames;
                frame_number = 0;
		break;

	case 'q':
	case 'Q':
		free_data ();
		exit ( 0 );
		break;
	}
}

static void mouse_func ( int button, int state, int x, int y )
{
	omx = mx = x;
	omx = my = y;

	if(!mouse_down[0]){hmx=x; hmy=y;}
	if(mouse_down[button]) mouse_release[button] = state == GLUT_UP;
	if(mouse_down[button]) mouse_shiftclick[button] = glutGetModifiers()==GLUT_ACTIVE_SHIFT;
	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y )
{
	mx = x;
	my = y;
}

static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}

static void idle_func ( void )
{
    if ( dsim ){
        integrator->integrate( *sys, dt );
        get_from_UI();
    }
    else        {remap_GUI();}

    glutSetWindow ( win_id );
    glutPostRedisplay ();
}

static void display_func ( void )
{
	pre_display ();

	draw_forces();
	draw_constraints();
	draw_particles();

	post_display ();
}


/*
----------------------------------------------------------------------
open_glut_window --- open a glut compatible window and set callbacks
----------------------------------------------------------------------
*/

static void open_glut_window ( void )
{
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( win_x, win_y );
	win_id = glutCreateWindow ( "Tinkertoys!" );

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

	pre_display ();

	glutKeyboardFunc ( key_func );
	glutMouseFunc ( mouse_func );
	glutMotionFunc ( motion_func );
	glutReshapeFunc ( reshape_func );
	glutIdleFunc ( idle_func );
	glutDisplayFunc ( display_func );
}


/*
----------------------------------------------------------------------
main --- main routine
----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{
	glutInit ( &argc, argv );

	if( argc < 2 ){
                printf("Usage: ./%s integrator=(1-euler, 2-RK2, 3-sympleticEuler, 4-RK4) [N dt d]", argv[0]);
		exit(0);
	}
	
        // decide which integrator to use
        switch ( argv[1][0] )
        {
        case '1':
            integrator = new EulerIntegrator();
            break;

        case '2':
            integrator = new RK2Integrator();
            break;

        case '3':
            integrator = new SymplecticEulerIntegrator();
            break;

        case '4':
            integrator = new RK4Integrator();
            break;

        default:
            integrator = new RK4Integrator();
        }
	
	if ( argc == 2 ) {
		N = 64;
                dt = 0.01f;
		d = 5.f;
		fprintf ( stderr, "Using defaults : N=%d dt=%g d=%g\n",
			N, dt, d );
	} else {
		N = atoi(argv[2]);
		dt = atof(argv[3]);
		d = atof(argv[4]);
	}

	printf ( "\n\nHow to use this application:\n\n" );
	printf ( "\t Toggle construction/simulation display with the spacebar key\n" );
	printf ( "\t Dump frames by pressing the 'd' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	dsim = 0;
	dump_frames = 0;
	frame_number = 0;
	
	init_system();
	
        win_x = 720;
        win_y = 720;
	open_glut_window ();

	glutMainLoop ();

	exit ( 0 );
}

