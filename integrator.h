#pragma once

#include <vector>
#include "System.h"

/**
 * Interface for integrators that can step the simulation of a system.
 */
class Integrator
{
public:
    Integrator() { }
    virtual ~Integrator() { }

    /**
     * Step the simulation of the given system by the given timestep.
     * @param sys The system to integrate. It should be integrated starting
     *   from the system's current step.
     * @param dt The length of the time step to integrate.
     */
    virtual void integrate( System sys, float dt ) const = 0;

    // used for storing state vectors locally
    // without allocating memory every time.
    typedef std::vector<Particle*> StateList;
};

/**
 * Uses the basic Euler integration method, x' = x + dx/dt * dt.
 */
class EulerIntegrator : public Integrator
{
public:
    EulerIntegrator() { }
    virtual ~EulerIntegrator() { }
    virtual void integrate( System sys, float dt ) const;
private:
    mutable StateList state;
    mutable StateList deriv_state;
};

/**
 * Uses the midpoint integration method.
 */
class RK2Integrator : public Integrator
{
public:
    RK2Integrator() { }
    virtual ~RK2Integrator() { }
        virtual void integrate( System sys, float dt ) const;
private:
        mutable StateList state;
        mutable StateList mid_deriv;
        mutable StateList deriv_state;
};

/**
 * Uses the 4th order Runge-Kutta integration method.
 */
class RK4Integrator : public Integrator
{
public:
    RK4Integrator() { }
    virtual ~RK4Integrator() { }
    virtual void integrate( System sys, float dt ) const;
private:
    mutable StateList state;
    // the derivatives at each guess
    mutable StateList k1;
    mutable StateList k2;
    mutable StateList k3;
    mutable StateList k4;
};

/**
 * Uses a sympletic euler integration method, calculating
 * the position explicitly and then the velocity implicitly
 */
class SymplecticEulerIntegrator : public Integrator
{
public:
    SymplecticEulerIntegrator() { }
    virtual ~SymplecticEulerIntegrator() { }
    virtual void integrate( System sys, float dt ) const;
private:
	mutable StateList state;
	mutable StateList deriv_state;
};

