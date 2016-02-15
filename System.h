#pragma once

#include <gfx/vec2.h>
#include <vector>
#include <stdlib.h>
#include "CircularWireConstraint.h"
#include "RodConstraint.h"
#include "SpringForce.h"

#define G 0.003f
#define EPSILON 1.0e-30
#define Ks 100.0f
#define Kd 100.0f

class System
{
public:

        System(std::vector<Particle*> pVector, std::vector<SpringForce*> forceVector,
                std::vector<CircularWireConstraint*> wireConstVector,
                std::vector<RodConstraint*> rodConstVector);
        ~System(void);

        void deriv_eval(std::vector<Particle*> &);
        void get_state(std::vector<Particle*> &);
        void get_forces(std::vector<SpringForce*> &);
        void get_rodConst(std::vector<RodConstraint*> &);
        void get_wireConst(std::vector<CircularWireConstraint*> &);
        void set_state(std::vector<Particle*>);
        // allow for adding spring forces after initialization
        void add_springForce(SpringForce*);
        // allow for adding rod constraints after initialization
        void add_rodConst(RodConstraint*);
        // allow for adding wire constraints after initialization
        void add_wireConst(CircularWireConstraint*);
        void pop_springForce();
        void pop_rodConst();
        void pop_wireConst();
        int size();

private:

        std::vector<Particle*> pVector;
        std::vector<SpringForce*> forceVector;
        std::vector<CircularWireConstraint*> wireConstVector;
        std::vector<RodConstraint*> rodConstVector;
	
};
