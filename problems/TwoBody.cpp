#include "TwoBody.h"

#include <math.h>
#include <stdexcept>
#include <iostream>
using namespace std;

#define CUBE(x) ((x)*(x)*(x))
#define SQUARE(x) ((x)*(x))
#define SQR(x) SQUARE(x)

#define Ith(v, i) (NV_Ith_S((v), (i)-1))
#define IJth(v, i, j) (DENSE_ELEM(v, (i)-1, (j)-1))

#define CAPTURE 101

int TwoBody::rhsFunc( 
		realtype t, N_Vector z, N_Vector z_prime, void *myself)
{
	TwoBody* me = (TwoBody*) myself;

	static realtype rho_x, rho_y, rho;
	static realtype x_1, y_1, x_2, y_2;

	x_1 = NV_Ith_S(z, 0);
	y_1 = NV_Ith_S(z, 2);
	x_2 = NV_Ith_S(z, 1);
	y_2 = NV_Ith_S(z, 3);

	rho_x = x_2 - x_1;
	rho_y = y_2 - y_1;

	rho = sqrt(rho_x*rho_x + rho_y*rho_y);
	
	if (rho < me->_capture_radius) {
		me->setError(CAPTURE, "Bodies captured by each other.");
		return -1;
	}
	
	static realtype rho3;
	rho3 = rho*rho*rho; 

	static realtype m1, m2;
	m1 = me->_mass_1;
	m2 = me->_mass_2;
	
	Ith(z_prime,1) = Ith(z,5);
	Ith(z_prime,2) = Ith(z,6);	
	Ith(z_prime,3) = Ith(z,7);	
	Ith(z_prime,4) = Ith(z,8);	
	Ith(z_prime,5) = - m2*(x_1 - x_2)/rho3;	
	Ith(z_prime,6) =   m1*(x_1 - x_2)/rho3;	
	Ith(z_prime,7) = - m2*(y_1 - y_2)/rho3;	
	Ith(z_prime,8) =   m1*(y_1 - y_2)/rho3;	
	
	return 0; 
}

int TwoBody::denseJac( long int N, DenseMat J, realtype t, 
	N_Vector z, N_Vector fy, void *jac_data, N_Vector tmp1, 
	N_Vector tmp2, N_Vector tmp3 )
{
	TwoBody * me = (TwoBody *) jac_data;

	static realtype x_1, y_1, x_2, y_2, rho_x, rho_y, rho;
		
	x_1 = NV_Ith_S(z, 0);
	y_1 = NV_Ith_S(z, 1);
	x_2 = NV_Ith_S(z, 2);
	y_2 = NV_Ith_S(z, 3);

	rho_x = x_2 - x_1;
	rho_y = y_2 - y_1;

	rho = sqrt(rho_x*rho_x + rho_y*rho_y);
	
	static realtype rho3;
	rho3 = rho*rho*rho; 
	static realtype rho5;
	rho5 = rho3*rho*rho;

	static realtype m1, m2;
	m1 = me->_mass_1;
	m2 = me->_mass_2;

	IJth(J,1,1) = 0;
	IJth(J,1,2) = 0;
	IJth(J,1,3) = 0;
	IJth(J,1,4) = 0;
	IJth(J,1,5) = 1;
	IJth(J,1,6) = 0;
	IJth(J,1,7) = 0;
	IJth(J,1,8) = 0;

	IJth(J,2,1) = 0;
	IJth(J,2,2) = 0;
	IJth(J,2,3) = 0;
	IJth(J,2,4) = 0;
	IJth(J,2,5) = 0;
	IJth(J,2,6) = 1;
	IJth(J,2,7) = 0;
	IJth(J,2,8) = 0;

	IJth(J,3,1) = 0;
	IJth(J,3,2) = 0;
	IJth(J,3,3) = 0;
	IJth(J,3,4) = 0;
	IJth(J,3,5) = 0;
	IJth(J,3,6) = 0;
	IJth(J,3,7) = 1;
	IJth(J,3,8) = 0;

	IJth(J,4,1) = 0;
	IJth(J,4,2) = 0;
	IJth(J,4,3) = 0;
	IJth(J,4,4) = 0;
	IJth(J,4,5) = 0;
	IJth(J,4,6) = 0;
	IJth(J,4,7) = 0;
	IJth(J,4,8) = 1;

	IJth(J,5,1) = -m2/rho3 + 3*m2*SQR(Ith(z,1) - Ith(z,2))/rho5;
	IJth(J,5,2) = m2/rho3 - 3*m2*SQR(Ith(z,1) - Ith(z,2))/rho5;
	IJth(J,5,3) = 3*m2*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,5,4) = -3*m2*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,5,5) = 0;
	IJth(J,5,6) = 0;
	IJth(J,5,7) = 0;
	IJth(J,5,8) = 0;

	IJth(J,6,1) = m1/rho3 - 3*m1*SQR(Ith(z,1) - Ith(z,2))/rho5;
	IJth(J,6,2) = -m1/rho3 + 3*m1*SQR(Ith(z,1) - Ith(z,2))/rho5;
	IJth(J,6,3) = -3*m1*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,6,4) = 3*m1*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,6,5) = 0;
	IJth(J,6,6) = 0;
	IJth(J,6,7) = 0;
	IJth(J,6,8) = 0;

	IJth(J,7,1) = 3*m2*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,7,2) = -3*m2*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,7,3) = -m2/rho3 + 3*m2*SQR(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,7,4) = m2/rho3 - 3*m2*SQR(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,7,5) = 0;
	IJth(J,7,6) = 0;
	IJth(J,7,7) = 0;
	IJth(J,7,8) = 0;

	IJth(J,8,1) = -3*m1*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,8,2) = 3*m1*(Ith(z,1) - Ith(z,2))*(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,8,3) = m1/rho3 - 3*m1*SQR(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,8,4) = -m1/rho3 + 3*m1*SQR(Ith(z,3) - Ith(z,4))/rho5;
	IJth(J,8,5) = 0;
	IJth(J,8,6) = 0;
	IJth(J,8,7) = 0;
	IJth(J,8,8) = 0;

	return 0; 
}

bool TwoBody::handleOption(int arg)
{
	switch (arg)
	{
		case 'R': // initial separation (\rho_0) 
			setInitialSeparation(atof(optarg));
			return true;

		case 'V': // launch speed in inertial reference frame
			setDeltaV(atof(optarg));
			return true;

		case 'c': // capture radius
			setCaptureRadius(atof(optarg));
			return true;

		case 'f': // mass fraction
			setMassFraction(atof(optarg));
			return true;
	}
}

void TwoBody::setInitialSeparation(realtype r)
{
	_initial_separation = r;
	updateInitialVector();
}

void TwoBody::setDeltaV(realtype s)
{
	_delta_v = s;
	updateInitialVector();
}

/**
 * Update the initial velocity vector to reflect the initial conditions
 * specified by _delta_v and _initial_separation.
 */
void TwoBody::updateInitialVector() 
{
	vector<realtype> v(getVectorSize());
	
	// v = (x_1, x_2, y_1, y_2, x_1', x_2', y_1', y_2')
	v.at(0) = 0;
	v.at(1) = _initial_separation;
	v.at(2) = 0;
	v.at(3) = 0;

	v.at(4) = 0;
	v.at(5) = 0;
	v.at(6) = - (_delta_v * _mass_2);
	v.at(7) =   (_delta_v * _mass_1);

	setInitialVector(v);
}

void TwoBody::setCaptureRadius(realtype r)
{
	_capture_radius = r;
}

void TwoBody::setMassFraction(realtype f)
{
	if (f < 0 || f > 1)
		throw invalid_argument("TwoBody::setMassFraction(): mass fraction must be between 0 and 1.");
	
	_mass_1 = f;
	_mass_2 = 1 - f;
}

void TwoBody::recordDataInit()
{
	cout << "t\tx_1\ty_1\tx_2\ty_2" << endl;
}

void TwoBody::recordData(realtype t, N_Vector z)
{
	cout << t << "\t" 
		 << NV_Ith_S(z, 0) << "\t" 
		 << NV_Ith_S(z, 2) << "\t"
		 << NV_Ith_S(z, 1) << "\t"
		 << NV_Ith_S(z, 3) << endl;
}
	
