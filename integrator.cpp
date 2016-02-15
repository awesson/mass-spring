#include "integrator.h"
#include "System.h"

/**
 * Uses the basic Euler integration method, x' = x + dx/dt * dt.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void EulerIntegrator::integrate( System sys, float dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize( size );
    deriv_state.resize( size );

    // get the current state
    sys.get_state( state );

    // compute the current derivative
    sys.deriv_eval( deriv_state );

    // update the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += deriv_state[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j] * dt;
        }
    }

    // set the updated state
    sys.set_state( state );
}

/**
 * Uses the midpoint integration method.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void RK2Integrator::integrate( System sys, float dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize( size );
    mid_deriv.resize( size );
    deriv_state.resize( size );

    // get the current state
    sys.get_state( state );

    // compute the current derivative
    sys.deriv_eval( deriv_state );

    // get midpoint state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += deriv_state[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get derivative at the midpoint
    sys.set_state( state );
    sys.deriv_eval( mid_deriv );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] -= deriv_state[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] -= deriv_state[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get full point state using midpoint derivative
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += mid_deriv[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] += mid_deriv[i]->deriv_velocity[j] * dt;
        }
    }

    // set the state to (pos + pos' * dt, t + dt)
    sys.set_state( state );
}

/**
 * Uses the 4th order Runge-Kutta integration method.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void RK4Integrator::integrate( System sys, float dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize( size );
    k1.resize( size );
    k2.resize( size );
    k3.resize( size );
    k4.resize( size );

    // get the current state
    sys.get_state( state );
    sys.get_state( k1 );
    sys.get_state( k2 );
    sys.get_state( k3 );
    sys.get_state( k4 );

    // compute the current derivative
    sys.deriv_eval( k1 );

    // get midpoint state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += k1[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] += k1[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get derivative at the midpoint
    sys.set_state( state );
    sys.deriv_eval( k2 );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] -= k1[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] -= k1[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get second midpoint state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += k2[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] += k2[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get derivative at the new midpoint
    sys.set_state( state );
    sys.deriv_eval( k3 );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] -= k2[i]->deriv_position[j] * dt/2.f;
            state[i]->Velocity[j] -= k2[i]->deriv_velocity[j] * dt/2.f;
        }
    }

    // get final state using new midpoint derivative
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += k3[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] += k3[i]->deriv_velocity[j] * dt;
        }
    }

    // get derivative at the guess for the final state
    sys.set_state( state );
    sys.deriv_eval( k4 );

    // reset the state
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] -= k3[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] -= k3[i]->deriv_velocity[j] * dt;
        }
    }

    // get final state using all 4 derivatives
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += k1[i]->deriv_position[j] * dt/6.f + k2[i]->deriv_position[j] * dt/3.f
                                   + k3[i]->deriv_position[j] * dt/3.f + k4[i]->deriv_position[j] * dt/6.f;
            state[i]->Velocity[j] += k1[i]->deriv_velocity[j] * dt/6.f + k2[i]->deriv_velocity[j] * dt/3.f
                                   + k3[i]->deriv_velocity[j] * dt/3.f + k4[i]->deriv_velocity[j] * dt/6.f;
        }
    }

    // set the state to (pos + pos' * dt, t + dt)
    sys.set_state( state );
}

/**
 * Uses a symplectic euler integration method. First the postion is
 * calculated explicitly and the the velocity is calculated implicitly.
 * @param sys The system to integrate
 * @param dt The time step to integrate over
 */
void SymplecticEulerIntegrator::integrate( System sys, float dt ) const
{
    int size = sys.size();

    if (size == 0)
        return;
    state.resize(size);
    deriv_state.resize(size);

    // get the current state (pos, t)
    sys.get_state( state );

    // compute the current derivative (pos')
    sys.deriv_eval( deriv_state );

    // update the x component explicitly.
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            //state[i]->Position[j] += deriv_state[i]->deriv_position[j] * dt;
            state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j] * dt;
        }
    }

    // update the y component implicitly with the x component at t + dt
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < 2; ++j){
            state[i]->Position[j] += state[i]->Velocity[j] * dt;
            //state[i]->Velocity[j] += deriv_state[i]->deriv_velocity[j] * dt;
        }
    }

    // set the state to (pos + pos' * dt, t + dt)
    sys.set_state( state );
}

