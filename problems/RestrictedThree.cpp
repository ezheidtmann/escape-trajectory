#include "RestrictedThree.h"

#include <math.h>
#include <stdexcept>
#include <iostream>
using namespace std;

#define CUBE(x) ((x)*(x)*(x))
#define SQUARE(x) ((x)*(x))

#define SIGN(x) (((x) > 0) - ((x) < 0))

#define CAPTURE_EARTH 101
#define CAPTURE_SUN   102

int RestrictedThree::rhsFunc( 
		realtype t, N_Vector z, N_Vector z_prime, void *myself)
{
	RestrictedThree* me = (RestrictedThree*) myself;

	// All distances are measured from the center of mass
	// of the two-body system. Unless otherwise specified,
	// these are distances of the particle from the other
	// objects.
	
	static realtype x;
	static realtype y;
	
	static realtype x_component;
	static realtype y_component;
	static realtype dist_x_earth;
	static realtype dist_y_earth;
	static realtype dist_earth;
	static realtype dist_x_sun;
	static realtype dist_y_sun;
	static realtype dist_sun;
	static realtype earth_multiplier;
	static realtype sun_multiplier;

	x_component = cos(t);
	y_component = sin(t);
	
	x = NV_Ith_S(z, 0);
	y = NV_Ith_S(z, 1);
	
	dist_x_earth = 
			x - me->_dist_com_earth*x_component;
	dist_y_earth = 
			y - me->_dist_com_earth*y_component;
	
	dist_x_sun = 
			x + me->_dist_com_sun*x_component;
	dist_y_sun =	
			y + me->_dist_com_sun*y_component;

	dist_earth = sqrt(SQUARE(dist_x_earth) + SQUARE(dist_y_earth));
	dist_sun = sqrt(SQUARE(dist_x_sun) + SQUARE(dist_y_sun));

	// If we're within a capture radius, it's time to abort.
	if (dist_earth < me->_capture_earth)
	{
		me->setError(CAPTURE_EARTH, "Projectile captured by EARTH");
		return -1; // unrecoverable error
	}
	else if (dist_sun < me->_capture_sun)
	{
		me->setError(CAPTURE_SUN, "Projectile captured by SUN");
		return -1;
	}

	// Record delta-epsilon / epsilon if requested
	if (me->_shouldRecordPseudoEnergy) 
	{
		// TODO: It might be more efficient to only do this calculation
		// when it's needed by recordData(). Investigate this if the 
		// simulation isn't fast enough.
		static realtype kineticEnergy;
		static realtype gravEnergy;

		// Eg ~= Me / Re + Ms / Rs
		gravEnergy = me->_mass_earth / dist_earth +
					 me->_mass_sun   / dist_sun;
		// Ek ~= 0.5 * v^2 + omega*(Vx*y - Vy*x)
		kineticEnergy = 0.5 * 
				(
					SQUARE(NV_Ith_S(z, 2)) + SQUARE(NV_Ith_S(z, 3))
				) 
				+
				( // omega == 1
				 	NV_Ith_S(z, 2) * y // Vx * y
					-
					NV_Ith_S(z, 3) * x // Vy * x
				);

		me->_pseudoEnergy = kineticEnergy - gravEnergy + 1.5;
		me->_pseudoEnergyIsDefined = true;

		if (t == me->getInitialTime())
			me->_initialPseudoEnergy = me->_pseudoEnergy;
	}

	earth_multiplier = me->_mass_earth / CUBE( dist_earth );
	sun_multiplier = me->_mass_sun / CUBE( dist_sun );

	// z_prime = (x', y', x'', y'')
	NV_Ith_S(z_prime, 2) = - earth_multiplier * dist_x_earth 
						   - sun_multiplier * dist_x_sun; 
	NV_Ith_S(z_prime, 3) = - earth_multiplier * dist_y_earth 
						   - sun_multiplier * dist_y_sun;
	NV_Ith_S(z_prime, 0) = NV_Ith_S(z, 2);
	NV_Ith_S(z_prime, 1) = NV_Ith_S(z, 3);

	return 0; // FIXME: is it true that we always succeed?
}

void RestrictedThree::setMassFraction(realtype earth_mass_fraction)
{
	_dist_com_earth = 1 - earth_mass_fraction;
	_dist_com_sun   = earth_mass_fraction;
	_mass_earth     = earth_mass_fraction;
	_mass_sun       = 1 - earth_mass_fraction;
}

void RestrictedThree::setLaunchRadius(realtype r)
{
	_launch_radius = r;
	updateInitialVectorPolar();
}

int RestrictedThree::denseJac( long int N, DenseMat J, realtype t, 
	N_Vector z, N_Vector fy, void *jac_data, N_Vector tmp1, 
	N_Vector tmp2, N_Vector tmp3 )
{

	RestrictedThree *me = (RestrictedThree*) jac_data;
		
	// Begin copy-and-paste from rhsFunc()
	

	static realtype x_component;
	static realtype y_component;
	static realtype dist_x_earth;
	static realtype dist_y_earth;
	static realtype dist_x_sun;
	static realtype dist_y_sun;
	static realtype earth_multiplier;
	static realtype sun_multiplier;
	
	x_component = cos(t);
	y_component = sin(t);
	
	dist_x_earth = 
			NV_Ith_S(z, 1) - me->_dist_com_earth*x_component;
	dist_y_earth = 
			NV_Ith_S(z, 2) - me->_dist_com_earth*y_component;
	
	dist_x_sun = 
			NV_Ith_S(z, 1) + me->_dist_com_sun*x_component;
	dist_y_sun =	
			NV_Ith_S(z, 2) + me->_dist_com_sun*y_component;

	// End copy-and-paste
	
	static realtype dist_earth;
	static realtype dist_sun;

	dist_earth = sqrt( 
					SQUARE(dist_x_earth) + SQUARE(dist_y_earth) 
				);
	dist_sun = sqrt( 
					SQUARE(dist_x_sun) + SQUARE(dist_y_sun) 
				);
	
	earth_multiplier = me->_mass_earth / 
			CUBE( dist_earth );
	sun_multiplier = me->_mass_sun / 
			CUBE( dist_sun );


#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) // IJth numbers rows,cols 1..NEQ 

	IJth(J,1,1) = 0;	IJth(J,1,2) = 0;	IJth(J,1,3) = 1;	IJth(J,1,4) = 0;
	IJth(J,2,1) = 0;	IJth(J,2,2) = 0;	IJth(J,2,3) = 0;	IJth(J,2,4) = 1; 
	IJth(J,3,1) = (3 * me->_mass_earth * SQUARE(dist_x_earth)) /
			(CUBE(dist_earth) * SQUARE(dist_earth)) - sun_multiplier * dist_x_sun;
	IJth(J,3,2) = (3 * me->_mass_earth * dist_x_earth * dist_y_earth) / 
			(CUBE(dist_earth) * SQUARE(dist_earth)) - sun_multiplier * dist_x_sun;
	IJth(J,3,3) = 0;
	IJth(J,3,4) = 0;
	IJth(J,4,1) = (3 * me->_mass_earth * dist_x_earth * dist_y_earth) /
			(CUBE(dist_earth) * SQUARE(dist_earth)) - sun_multiplier * dist_y_sun;
	IJth(J,4,2) = (3 * me->_mass_earth * SQUARE(dist_y_earth)) / 
			(CUBE(dist_earth) * SQUARE(dist_earth)) - sun_multiplier * dist_y_sun;
	IJth(J,4,3) = 0;
	IJth(J,4,4) = 0;

	return 0; // FIXME: Is it true that we always succeed?
}

bool RestrictedThree::handleOption(int arg)
{
	switch (arg)
	{
		case 'f': // earth mass fraction
			setMassFraction(atof(optarg));
			return true;

		case 'R': // launch radius
			setLaunchRadius(atof(optarg));
			return true;

		case 'V': // launch speed in inertial reference frame
			setLaunchSpeed(atof(optarg));
			return true;

		case 'q': // launch angle (q = theta in Mathematica)
			setLaunchAngle(atof(optarg));
			return true;

		case 'X': // launch speed in x-direction
			setLaunchSpeedX(atof(optarg));
			return true;

		case 'Y': // launch speed in y-direction
			setLaunchSpeedY(atof(optarg));
			return true;

		case 'M': // velocity parameter magnification
			setLaunchSpeedMagnification(atof(optarg));
			return true;

		case 'S': // minimum velocity in component-wise spec
			setLaunchSpeedMinimum(atof(optarg));
			return true;

		case 'c': // capture radius
			setCaptureRadius(atof(optarg));
			return true;

		case 'E': // record pseudoenergy
			enablePseudoEnergy();
			return true;

		case 'D': // record total distance
			enableTotalDistance();
			return true;

		case 'W': // record winding number
			enableWindingNumber();
			return true;

		case 'P': // speed is given as p = v^2 - V_E^2
			_speed_parameter = SPEED_PARAMETER;
			return true;

		case 'G': // speed is given as pseudo-enerGy
			_speed_parameter = SPEED_ENERGY;
			return true;
	}
}

void RestrictedThree::setLaunchSpeed(realtype s)
{
	_launch_speed = s;
	updateInitialVectorPolar();
}
	
void RestrictedThree::setLaunchAngle(realtype theta)
{
	_launch_angle = theta;
	updateInitialVectorPolar();
}

void RestrictedThree::setLaunchSpeedMagnification(realtype m)
{
	_v_map_slope = m;
	updateInitialVectorComp();
}

void RestrictedThree::setLaunchSpeedMinimum(realtype vs)
{
	_min_launch_speed = vs;
	updateInitialVectorComp();
}

void RestrictedThree::setLaunchSpeedX(realtype vx)
{
	_launch_speed_x = vx;
	updateInitialVectorComp();
}

void RestrictedThree::setLaunchSpeedY(realtype vy)
{
	_launch_speed_y = vy;
	updateInitialVectorComp();
}

void RestrictedThree::updateInitialVectorPolar() 
{
	// This code implicity assumes the initial time is t = 0.
	// At this time, the earth is at (_dist_com_earth, 0)
	// and moving at (0, _dist_com_earth).
	//
	// Details on the velocity:
	//
	// The earth travels around a circle of radius _dist_com_earth
	// in time 2*pi (because omega = 1). Such a circle has circumference
	// 2*pi * _dist_com_earth, so (ignoring units)
	//   v = d/t = (2*pi*r)/(2*pi) = r = _dist_com_earth
	//
	
	// If speed is given as a "parameter", we treat it as the experimentally 
	// determined "parameter" p = v^2 -v_E^2. The instance variable
	// _speed_parameter is set by the -P command-line option. If it is
	// not set, the speed is treated as the value v/v_E.
	realtype speed;
	if (_speed_parameter)
	{
		speed = sqrt(_launch_speed + 2*_mass_earth/_launch_radius);
		if (speed != speed) // speed == NaN, which means sqrt(-1)
			speed = 0;
	}
	else
		speed = _launch_speed * (sqrt(2*_mass_earth/_launch_radius));
	
	vector<realtype> v(getVectorSize());
	realtype theta = M_PI * _launch_angle / RCONST(180.0);
	v.at(0) = _launch_radius * cos(theta) + _dist_com_earth; 
	v.at(1) = _launch_radius * sin(theta); 
	v.at(2) = speed          * cos(theta); 
	v.at(3) = speed          * sin(theta) + _dist_com_earth; 
	setInitialVector(v);
}

/**
 * Update the initial vector to match launch radius and velocity
 * in X and Y directions (the "component"-wise representation of
 * initial velocity). If speed is given as a parameter, then the
 * following mapping function is applied:
 *
 * Given a parameter speed in components \f$P_x\f$ and \f$P_y\f$,
 * the actual velocity is found by components:
 *
 * \f{eqnarray*}
 *  V_x &=& \frac{V_x}{m} + V_s \cos \Theta\\
 *  v    &=& \frac{V_x}{m} + V_s \cos \left( \arccos \left( 
 *      	\frac{P_x}{\sqrt{P_x^2 + P_y^2}} \right) \right)\\
 *   v   &=& \frac{V_x}{m} 
 *      	+ V_s \left( \frac{P_x}{ \sqrt{P_x^2 + P_y^2} } \right)\\
 *  V_y &=& \frac{V_y}{m}
 *  		+ V_s \left( \frac{P_y}{ \sqrt{P_x^2 + P_y^2} } \right) 
 * \f}
 *
 * Internally the variables are represented like this:
 *
 *  - \c _launch_speed_x: \f$P_x\f$
 *  - \c _launch_speed_y: \f$P_y\f$
 *  - \c _v_map_slope: \f$m\f$
 *  - \c vx: \f$V_x\f$
 *  - \c vy: \f$V_y\f$
 *
 * If the speed is given as pseudo-energy (if the \c -G option is present),
 * then the parameters \c _launch_speed_x and \c _launch_speed_y are used
 * to calculate a launch angle and then the magnitude is interpreted as 
 * pseudo-energy. The launch speed is then found from the formula for
 * pseudo-energy.
 *
 * @see setLaunchSpeedX()
 * @see setLaunchSpeedY()
 * @see setLaunchSpeedMagnification()
 * @see setLaunchSpeedMinimum()
 */
void RestrictedThree::updateInitialVectorComp()
{
	realtype vx, vy, unmappedx, unmappedy, unmapped_mag;

	realtype mag, cosq, sinq;
	mag = sqrt(SQUARE(_launch_speed_x) + SQUARE(_launch_speed_y));

	cosq = _launch_speed_x / mag;
	sinq = _launch_speed_y / mag;

	unmappedx = (_launch_speed_x) * _v_map_slope + _min_launch_speed * cosq;
	unmappedy = (_launch_speed_y) * _v_map_slope + _min_launch_speed * sinq;
	unmapped_mag = mag * _v_map_slope + _min_launch_speed;

    vector<realtype> v(getVectorSize());
    v.at(0) = _launch_radius * cosq + _dist_com_earth;
    v.at(1) = _launch_radius * sinq;

	if (_speed_parameter == SPEED_PARAMETER)
	{
		// Calculate the square of escape speed from earth alone
		realtype ve2;
		ve2 = 2*_mass_earth/_launch_radius;

		vx = sqrt( unmappedx*cosq + ve2*cosq*cosq ) * SIGN(cosq);
		vy = sqrt( unmappedy*sinq + ve2*sinq*sinq ) * SIGN(sinq);
	}
	else if (_speed_parameter == SPEED_ENERGY)
	{
		realtype gpe, x, y;
		
		x = v.at(0);
		y = v.at(1);
		
		gpe = _mass_earth / _launch_radius +
				_mass_sun / 
					sqrt(
						SQUARE(_dist_com_sun + x) +
						SQUARE(y)
					);
			
		realtype a, b, c; // a,b,c in ax^2 + bx + c = 0

		a = 0.5;
		b = _mass_sun * sinq + y * cosq - x * sinq;
		c = _mass_sun*_mass_sun / 2 - x * _mass_sun - 
				gpe - unmapped_mag + 1.5;
		
		realtype v0;

		v0 = (-b + sqrt(b*b - 4*a*c)) / (2*a); // quadratic formula

		vx = v0 * cosq;
		vy = v0 * sinq;
	}
	else
	{
		vx = vx * (sqrt(2*_mass_earth/_launch_radius));
		vy = vy * (sqrt(2*_mass_earth/_launch_radius));
	}

    v.at(2) = vx;
    v.at(3) = vy + _dist_com_earth;

    setInitialVector(v);

}

void RestrictedThree::setCaptureRadius(realtype r)
{
	_capture_sun = _capture_earth = r;
}

void RestrictedThree::recordDataInit()
{
	cout << "#t\tx\ty";

	if (_shouldRecordPseudoEnergy)
		cout << "\te";
	if (_shouldRecordDistance)
		cout << "\td";
	if (_shouldRecordWinding)
		cout << "\twE\twS\twC";

	cout << endl;
}

void RestrictedThree::recordData(realtype t, N_Vector z)
{
	cout << t << "\t" 
		 << NV_Ith_S(z, 0) << "\t" 
		 << NV_Ith_S(z, 1);

#define ABS(x) (x) > 0 ? (x) : -(x)
	if (_shouldRecordPseudoEnergy)
	{
		if (_pseudoEnergyIsDefined)
			cout << "\t" << ABS(_pseudoEnergy - _initialPseudoEnergy) 
				/ _initialPseudoEnergy;
		else
			cout << "\t" << 0;
	}

	if (_shouldRecordDistance)
	{
		cout << "\t" << _totalDistance;
	}
	
	if (_shouldRecordWinding)
	{
		cout << "\t" << _earthWinding.getWindingNumber()
			 << "\t" << _sunWinding.getWindingNumber()
			 << "\t" << _comWinding.getWindingNumber();
	}

	cout << endl;
}

void RestrictedThree::postStep(realtype t, N_Vector z)
{
	// This looks like a reimplementation of rhsFunc(),
	// and it pretty much is. But copy-paste is easier (and faster)
	// than modularizing this computation.

	static realtype x, y;
	static realtype x_component;
	static realtype y_component;
	static realtype dist_x_earth;
	static realtype dist_y_earth;
	static realtype dist_x_sun;
	static realtype dist_y_sun;

	static bool first_time = true;

	x = NV_Ith_S(z, 0);
	y = NV_Ith_S(z, 1);

	if (_shouldRecordWinding)
	{
		x_component = cos(t);
		y_component = sin(t);
	
		dist_x_earth = 
			x - _dist_com_earth * x_component;
		dist_y_earth =
			y - _dist_com_earth * y_component;
	
		dist_x_sun = 
			x + _dist_com_sun * x_component;
		dist_y_sun =
			y + _dist_com_sun * y_component;

		// Now for the actual recording.
	
		_earthWinding.addStep(dist_x_earth, dist_y_earth);
		_sunWinding.addStep(dist_x_sun, dist_y_sun);
		_comWinding.addStep(x, y);
	}

	if (_shouldRecordDistance) 
	{
		if (!first_time)
		{
			_totalDistance += sqrt(
				SQUARE( x - _lastx ) +
				SQUARE( y - _lasty )
			);
		}
		
		_lastx = x;
		_lasty = y;
		first_time = false;
	}
}

