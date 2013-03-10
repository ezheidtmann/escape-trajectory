#include "SingleBody.h"

#include <math.h>
#include <stdexcept>
#include <iostream>
using namespace std;

#define CUBE(x) ((x)*(x)*(x))
#define SQUARE(x) ((x)*(x))

#define CAPTURE 101

int SingleBody::rhsFunc( 
		realtype t, N_Vector z, N_Vector z_prime, void *myself)
{
	SingleBody* me = (SingleBody*) myself;

	// x' = v
	// v' = x'' = x^-2
	
	static realtype x;
	x = NV_Ith_S(z, 0);

	if (x < me->_capture_radius) {
		me->setError(CAPTURE, "Projectile captured.");
		return -1;
	}
	
	NV_Ith_S(z_prime, 0) = NV_Ith_S(z, 1);
	NV_Ith_S(z_prime, 1) = RCONST(-1.0) / (SQUARE(x));

	return 0; 
}

void SingleBody::setLaunchRadius(realtype r)
{
	_launch_radius = r;
	updateInitialVector();
}

int SingleBody::denseJac( long int N, DenseMat J, realtype t, 
	N_Vector z, N_Vector fy, void *jac_data, N_Vector tmp1, 
	N_Vector tmp2, N_Vector tmp3 )
{
	return 0; 
}

bool SingleBody::handleOption(int arg)
{
	switch (arg)
	{
		case 'R': // launch radius
			setLaunchRadius(atof(optarg));
			return true;

		case 'V': // launch speed in inertial reference frame
			setLaunchSpeed(atof(optarg));
			return true;

		case 'c': // capture radius
			setCaptureRadius(atof(optarg));
			return true;
	}
}

void SingleBody::setLaunchSpeed(realtype s)
{
	_launch_speed = s;
	updateInitialVector();
}

void SingleBody::updateInitialVector() 
{
	vector<realtype> v(getVectorSize());
	v.at(0) = _launch_radius;
	v.at(1) = _launch_speed;
	setInitialVector(v);
}

void SingleBody::setCaptureRadius(realtype r)
{
	_capture_radius = r;
}

void SingleBody::recordData(realtype t, N_Vector z)
{
	cout << t << "\t" 
		 << NV_Ith_S(z, 0) << "\t" 
		 << NV_Ith_S(z, 1) << endl;
}
	
