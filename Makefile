# $Id: gfx-config.in 343 2008-09-13 18:34:59Z garland $

CXX = g++
CXXFLAGS = -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H 
OBJS = Solver.o Particle.o TinkerToy.o RodConstraint.o SpringForce.o CircularWireConstraint.o imageio.o linearSolver.o System.o integrator.o

project1: $(OBJS)
	$(CXX) -o $@ $^ -lpng -framework GLUT -framework OpenGL
clean:
	rm $(OBJS) project1
