#pragma once

#include "Particle.h"
#define INF 1e5

class SpringForce {
 public:
  SpringForce(Particle* i_p1, Particle* i_p2, double i_dist, double i_ks, double i_kd);

  void add_force();
  void draw();

 private:

  Particle * const p1;   // particle 1
  Particle * const p2;   // particle 2
  double const dist;     // rest length
  double const ks, kd; // spring strength constants
};
