#include "linearSolver.h"

implicitMatrixImpl::implicitMatrixImpl(std::vector<CircularWireConstraint*> i_wireConstVector,
        std::vector<RodConstraint*> i_rodConstVector) : wireConstVector(i_wireConstVector),
        rodConstVector(i_rodConstVector){ }

implicitMatrixImpl::~implicitMatrixImpl(){ }

/**
 * Assumes that constraints are numbered so that wire constaints are lower numbers
 * and that they are numbered in the order of the respective vectors (wireConstVector and rodConstVector.
 **/
void implicitMatrixImpl::matVecMult(double x[], double r[]){
    int num_wC = wireConstVector.size();
    int num_rC = rodConstVector.size();
    int num_const = num_wC + num_rC;

    for(int i = 0; i < num_const; ++i){
        r[i] = 0;

        // adds constraint information for a wire constraint
        if(i < num_wC){
            // diagnal term
            r[i] += wireConstVector[i]->get_J() * wireConstVector[i]->get_J() * x[i]
                    / wireConstVector[i]->get_mass();

            // adds constraint information for particles with both a rod and wire constraint
            for(int k = 0; k < rodConstVector.size(); ++k){
                if(wireConstVector[i]->get_id() == rodConstVector[k]->get_id1()){
                    r[i] += wireConstVector[i]->get_J()*rodConstVector[k]->get_J() * x[k]
                            / wireConstVector[i]->get_mass();
                    // by symetry, can add this now
                    r[k] += wireConstVector[i]->get_J()*rodConstVector[k]->get_J() * x[i]
                            / wireConstVector[i]->get_mass();
                } else{
                    if(wireConstVector[i]->get_id() == rodConstVector[k]->get_id2()){
                        r[i] -= wireConstVector[i]->get_J()*rodConstVector[k]->get_J() * x[k]
                                / wireConstVector[i]->get_mass();
                        // by symetry, can add this now
                        r[k] -= wireConstVector[i]->get_J()*rodConstVector[k]->get_J() * x[i]
                                / wireConstVector[i]->get_mass();
                    }
                }
            }

        // adds constraint information for a rod constraint.
        //(rod-wire constraints will have already been added)
        }else{
            // diagonal term
            r[i] += rodConstVector[i - num_wC]->get_J() * rodConstVector[i - num_wC]->get_J() * x[i]
                    * (1.f / rodConstVector[i - num_wC]->get_mass1() + 1.f / rodConstVector[i - num_wC]->get_mass2());

            // add information about particles with two rod constraints
            for(int k = 0; k < rodConstVector.size(); ++k){
                if(k != i - num_wC){
                    if(rodConstVector[i - num_wC]->get_id1() == rodConstVector[k]->get_id1()){
                        r[i] += rodConstVector[i - num_wC]->get_J() * rodConstVector[k]->get_J() * x[i]
                                / rodConstVector[i - num_wC]->get_mass1();
                    } else{
                        if(rodConstVector[i - num_wC]->get_id2() == rodConstVector[k]->get_id2()){
                            r[i] += rodConstVector[i - num_wC]->get_J() * rodConstVector[k]->get_J() * x[i]
                                    / rodConstVector[i - num_wC]->get_mass2();
                        } else{
                            if(rodConstVector[i - num_wC]->get_id1() == rodConstVector[k]->get_id2()){
                                r[i] -= rodConstVector[i - num_wC]->get_J() * rodConstVector[k]->get_J() * x[i]
                                        / rodConstVector[i - num_wC]->get_mass1();
                            } else{
                                if(rodConstVector[i - num_wC]->get_id2() == rodConstVector[k]->get_id1()){
                                    r[i] -= rodConstVector[i - num_wC]->get_J() * rodConstVector[k]->get_J() * x[i]
                                            / rodConstVector[i - num_wC]->get_mass2();
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// vector helper functions

void vecAddEqual(int n, double r[], double v[])
{
  for (int i = 0; i < n; i++)
    r[i] = r[i] + v[i];
}

void vecDiffEqual(int n, double r[], double v[])
{
  for (int i = 0; i < n; i++)
    r[i] = r[i] - v[i];
}

void vecAssign(int n, double v1[], double v2[])
{
  for (int i = 0; i < n; i++)
    v1[i] = v2[i];
}

void vecTimesScalar(int n, double v[], double s)
{
  for (int i = 0; i < n; i++)
    v[i] *= s;
}

double vecDot(int n, double v1[], double v2[])
{
  double dot = 0;
  for (int i = 0; i < n; i++)
    dot += v1[i] * v2[i];
  return dot;
}

double vecSqrLen(int n, double v[])
{
  return vecDot(n, v, v);
}

double ConjGrad(int n, implicitMatrix *A, double x[], double b[], 
		double epsilon,	// how low should we go?
		int    *steps)
{
  int		i, iMax;
  double	alpha, beta, rSqrLen, rSqrLenOld, u;

  double *r = (double *) malloc(sizeof(double) * n);
  double *d = (double *) malloc(sizeof(double) * n);
  double *t = (double *) malloc(sizeof(double) * n);
  double *temp = (double *) malloc(sizeof(double) * n);

  vecAssign(n, x, b);

  vecAssign(n, r, b);
  A->matVecMult(x, temp);
  vecDiffEqual(n, r, temp);

  rSqrLen = vecSqrLen(n, r);

  vecAssign(n, d, r);

  i = 0;
  if (*steps)
    iMax = *steps;
  else
    iMax = MAX_STEPS;
		
  if (rSqrLen > epsilon)
    while (i < iMax) {	
      i++;
      A->matVecMult(d, t);
      u = vecDot(n, d, t);
      
      if (u == 0) {
	printf("(SolveConjGrad) d'Ad = 0\n");
	break;
      }
      
      // How far should we go?
      alpha = rSqrLen / u;
      
      // Take a step along direction d
      vecAssign(n, temp, d);
      vecTimesScalar(n, temp, alpha);
      vecAddEqual(n, x, temp);
      
      if (i & 0x3F) {
	vecAssign(n, temp, t);
	vecTimesScalar(n, temp, alpha);
	vecDiffEqual(n, r, temp);
      } else {
	// For stability, correct r every 64th iteration
	vecAssign(n, r, b);
	A->matVecMult(x, temp);
	vecDiffEqual(n, r, temp);
      }
      
      rSqrLenOld = rSqrLen;
      rSqrLen = vecSqrLen(n, r);
      
      // Converged! Let's get out of here
      if (rSqrLen <= epsilon)
	break;			    
      
      // Change direction: d = r + beta * d
      beta = rSqrLen/rSqrLenOld;
      vecTimesScalar(n, d, beta);
      vecAddEqual(n, d, r);
    }
  
  // free memory

  free(r);
  free(d);
  free(t);
  free(temp);
		
  *steps = i;
  return(rSqrLen);
}


