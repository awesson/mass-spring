#include "System.h"
#include "linearSolver.h"


System::System(std::vector<Particle*> i_pVector, std::vector<SpringForce*> i_forceVector,
std::vector<CircularWireConstraint*> i_wireConstVector, std::vector<RodConstraint*> i_rodConstVector) :
	pVector(i_pVector),
	forceVector(i_forceVector),
	wireConstVector(i_wireConstVector), 
	rodConstVector(i_rodConstVector)
{
}

System::~System(void)
{
}

void System::deriv_eval(std::vector<Particle*>& o_pVector){
        int size = pVector.size();
        int num_f = forceVector.size();
	
        // reset forces to just gravity
        for(int i = 0; i < size; ++i){
            pVector[i]->forces = Vec2f(0.0, -G);
        }
	
        // add spring forces
        for(int i = 0; i < num_f; ++i){
                forceVector[i]->add_force();
	}

        int num_wC = wireConstVector.size();
        int num_const = num_wC + rodConstVector.size();
        // initialize matrix
        implicitMatrix* JWJ_t = new implicitMatrixImpl(wireConstVector, rodConstVector);

        double* lambda = (double*) malloc(num_const*sizeof(double));
        double* b = (double*) malloc(num_const*sizeof(double));

        // Calculate b
        for(int i = 0; i < num_const; ++i){
            lambda[i] = 0;

            // -(Jdot)*(qdot)
            if(i < num_wC){
                b[i] = -wireConstVector[i]->get_Jdot() * pVector[ wireConstVector[i]->get_id() ]->Velocity;
            } else{
                b[i] = -rodConstVector[i - num_wC]->get_Jdot() *
                       (pVector[ rodConstVector[i - num_wC]->get_id1() ]->Velocity -
                        pVector[ rodConstVector[i - num_wC]->get_id2() ]->Velocity);
            }

            // -JWQ
            if(i < num_wC){
                b[i] -= wireConstVector[i]->get_J() * pVector[ wireConstVector[i]->get_id() ]->forces;
            } else{
                b[i] -= rodConstVector[i - num_wC]->get_J() *
                       (pVector[ rodConstVector[i - num_wC]->get_id1() ]->forces -
                        pVector[ rodConstVector[i - num_wC]->get_id2() ]->forces);
            }

            // -ks*C
            if(i < num_wC){
                b[i] -= Ks * wireConstVector[i]->get_C();
            } else{
                b[i] -= Ks * rodConstVector[i - num_wC]->get_C();
            }

            // -kd*Cdot
            if(i < num_wC){
                b[i] -= Kd * wireConstVector[i]->get_Cdot();
            } else{
                b[i] -= Kd * rodConstVector[i - num_wC]->get_Cdot();
            }

        }

        int steps = MAX_STEPS;
        double err = ConjGrad(num_const, JWJ_t, lambda, b, EPSILON, &steps);
        if(err > EPSILON){
            printf("probably too many constraints to satisfy!!/n");
            exit(0);
        }

        // calculate J_t*lambda and add those constraint forces
        for(int i = 0; i < num_const; ++i){
            if(i < num_wC){
                pVector[ wireConstVector[i]->get_id() ]->forces += wireConstVector[i]->get_J() * lambda[i];
            } else{
                pVector[ rodConstVector[i - num_wC]->get_id1() ]->forces += rodConstVector[i - num_wC]->get_J() * lambda[i];
                pVector[ rodConstVector[i - num_wC]->get_id2() ]->forces -= rodConstVector[i - num_wC]->get_J() * lambda[i];
            }
        }

        // set the derivative of position to the velocity and
        // the derivative of the velocity to the total force divided by the mass
        for(int i = 0; i < size; ++i){
            pVector[i]->deriv_position = pVector[i]->Velocity;
            pVector[i]->deriv_velocity = pVector[i]->forces / pVector[i]->mass;
        }

        free(lambda);
        free(b);
        delete JWJ_t;

        // return particle vector with the derivatives
        o_pVector = pVector;
}

void System::get_state(std::vector<Particle*>& o_pVector){
        o_pVector = pVector;
}

void System::set_state(std::vector<Particle*> i_pVector){
        pVector = i_pVector;
}

void System::get_forces(std::vector<SpringForce*>& o_forceVector){
        o_forceVector = forceVector;
}

void System::get_rodConst(std::vector<RodConstraint*>& o_rodConstVector){
        o_rodConstVector = rodConstVector;
}

void System::get_wireConst(std::vector<CircularWireConstraint*>& o_wireConstVector){
        o_wireConstVector = wireConstVector;
}

void System::add_springForce(SpringForce* f){
    forceVector.push_back(f);
}

void System::pop_springForce(){
    forceVector.pop_back();
}

void System::add_rodConst(RodConstraint* rod){
    rodConstVector.push_back(rod);
}

void System::pop_rodConst(){
    rodConstVector.pop_back();
}

void System::add_wireConst(CircularWireConstraint* wire){
    wireConstVector.push_back(wire);
}

void System::pop_wireConst(){
    wireConstVector.pop_back();
}

int System::size(){
        return pVector.size();
}

