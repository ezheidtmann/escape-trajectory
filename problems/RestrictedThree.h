#include "AbstractProblem.h"
#include "Winding.h"

#include <math.h>

/**
 * Implements the restricted three-body problem.
 *
 * In this problem, there are two bodies of comparable mass orbiting
 * their common center of mass. Due to carefully chosen initial conditions,
 * we can assume that they orbit the center of mass in concentric
 * circles. We call these two bodies Earth and Sun, but the mass ratio
 * is not normally as extreme as in the actual solar system.
 * 
 * We then launch a projectile from Earth within the plane of the massive
 * bodies and calculate its trajectory. The equations of motion of the 
 * projectile are:
 *
 * \f{eqnarray*}
 * 	\ddot{x} & = & -\frac{GM_E}{d_E^3}\left(x-L_E\cos\omega t\right)-
 * 		\frac{GM_S}{d_S^3}\left(x+L_S\cos\omega t\right)\\
 * 	\ddot{y} & = & -\frac{GM_E}{d_E^3}\left(y-L_E\sin\omega t\right)-
 * 		\frac{GM_S}{d_S^3}\left(y+L_S\sin\omega t\right)
 * \f}
 *
 * where \f$M_E\f$ is the mass of the earth, \f$M_S\f$ is the mass of
 * the sun, \f$L_E\f$ and \f$L_S\f$ are the distances from the center of
 * mass to the earth and sun respectively, and \f$d_E\f$ and \f$d_S\f$ 
 * are the distances from the projectile to the earth and sun
 * respectively. 
 */
class RestrictedThree : public AbstractProblem {
  public:
	/**
	 * The CVODE callback function implementing the ODEs that represent
	 * this problem. See rhsFuncAbstract() for details on the behavior
	 * and expectation of this function.
	 *
	 * @param t The current time in the integration
	 * @param z The current "position" vector, having values
	 * 			\f$(x, y, x^\prime, y^\prime)\f$ or
	 * 			\f$(x, y, V_x, V_y)\f$.
	 * @param z_prime The current "velocity" vector, having values
	 * 			\f$(x^\prime, y^\prime, 
	 * 				x^{\prime\prime}, y^{\prime\prime})\f$ or
	 * 			\f$(x^\prime, y^\prime, V_x^\prime, V_y^\prime)\f$.
	 * @param myself A pointer to a RestrictedThree object, which is
	 * 			treated as if it were the "this" pointer.
	 */
	static int rhsFunc( realtype t, N_Vector z, N_Vector z_prime, void *myself);

	/**
	 * The dense Jacobian function. I don't pretend to understand
	 * its implementation.
	 */
	static int denseJac( long int N, DenseMat J, realtype t, 
		N_Vector y, N_Vector fy, void *jac_data, N_Vector tmp1, 
		N_Vector tmp2, N_Vector tmp3 );
	
	/**
	 * Set Earth's mass fraction. This determines the mass and
	 * orbital radius for both Earth and Sun.
	 *
	 * The only other constant in this system is the angular 
	 * velocity of the two primaries. This is constant by Kepler's
	 * third law: it depends only on the total mass and maximum
	 * separation of the two bodies.
	 *
	 * @param earth_mass_fraction The fraction of the total mass
	 * 			given to earth. Default: 0.3
	 */
	void setMassFraction(realtype earth_mass_fraction);
	
	/**
	 * Set launch radius. The infinitesimal projectile is 
	 * launched not from the center of the earth (where the 
	 * gravitatonal force is infinite), but from some small
	 * but finite distance away from the center of the Earth.
	 * This distance is the launch radius.
	 *
	 * @param r Launch radius, as a fraction of the maximum
	 * 			separation between Earth and Sun.
	 */
	void setLaunchRadius(realtype r);
	
	/**
	 * This problem consists of 4 ODEs, so this function always
	 * returns 4.
	 */
	unsigned int getVectorSize() { return 4; }

	/**
	 * Additional command-line options specific to this problem.
	 *
	 * Unless otherwise noted, every option takes as an argument
	 * a single floating-point value.
	 *
	 *  - \c -f: earth mass fraction
	 *  - \c -R: launch radius
	 *  - \c -V: launch speed
	 *  - \c -q: launch angle (q = theta)
	 *  - \c -c: capture radius
	 *  - \c -X: launch speed in the X direction (\c -P is ignored).
	 *  	This or -Y must appear after -R if specified.
	 *  - \c -Y: launch speed in the Y direction (\c -P is ignored).
	 *  	This or -X must appear after -R if specified.
	 *  - \c -M: magnification factor for velocity parameter
	 *  - \c -S: minimum ("starting") velocity for component-wise
	 *  	velocity parameter
	 *  - \c -P: treat launch speed (\c -V) as the "parameter"
	 *   \f$p = v^2 - v_E^2\f$, rather than as the relative speed
	 *   \f$v_R = v/v_E\f$. This option does not take an
	 *   argument.
	 *  - \c -G: treate launch speed (\c -X or \c -Y) as pseudo-energy.
	 *  This argument does not take an argument.
	 *  - \c -E: enable recording pseudo-energy. This option
	 *   does not take an argument.
	 *  - \c -D: enable recording total distance covered
	 *  - \c -W: enable recording "winding number", or total net
	 *   radians (measured from the center of mass) covered. Large 
	 *   winding numbers correlate with loops; small numbers mean
	 *   straight shots or switchbacks.
	 *
	 * @see handleOption()
	 * @return An option string representing these options, in the format
	 * required by \c getopt(3).
	 */
	string getOptionString() { return "f:R:V:q:c:X:Y:M:S:PEDWG"; }

	/**
	 * Handle command-line options. 
	 *
	 * @see getOptionString()
	 */
	bool handleOption( int arg );

	/**
	 * Set launch speed. If speed is to be interpreted absolutely (which
	 * is the case if the \c -A option was passed to the program), this is
	 * the initial speed \f$v\f$ relative to the earth.
	 *
	 * If speed is interpreted relatively (the default), this parameter 
	 * is actually the value \f$p = v^2 - v_E^2\f$, where \f$v_E\f$ is
	 * the escape velocity from earth alone, given by
	 * 
 	 * \f{eqnarray*}
	 *   \frac{1}{2}mv_0^2 & = & \frac{GMm}{R_L}\\
	 *   v_0^2 & = & \frac{2GM}{R_L}\\
	 *   v_0 & = & \sqrt{\frac{2M}{R_L}} = v_E
	 * \f}
	 *
	 * Here, \f$R_L\f$ is the launch radius, \f$M\f$ is the mass of the
	 * earth, and \f$G = 1\f$ in our chosen system of units.
	 *
	 * @note This method of specifying launch parameters has been
	 * 		superseded by the component-wise method. See 
	 * 		setLaunchSpeedX(), setLaunchSpeedY(), and
	 * 		updateInitialVectorComp().
	 */
	void setLaunchSpeed(realtype s);

	/**
	 * Set launch angle. The angle is measured counter-clockwise from
	 * the positive x-axis, as is standard in mathematics.
	 *
	 * @note This method of specifying launch parameters has been
	 * 		superseded by the component-wise method. See 
	 * 		setLaunchSpeedX(), setLaunchSpeedY(), and
	 * 		updateInitialVectorComp().
	 *
	 * @param angle Launch angle in degrees, CCW from horizontal.
	 */
	void setLaunchAngle(realtype angle);

	/**
	 * Set launch speed in the X direction
	 *
	 * @param vx X-component of the velocity
	 */
	void setLaunchSpeedX(realtype vx);

	/**
	 * Set launch speed in the Y direction
	 *
	 * @param vy Y-component of the velocity
	 */
	void setLaunchSpeedY(realtype vy);

	/**
	 * Set launch speed magnification. This has an effect on 
	 * both polar and component-wise launch speed specification,
	 * but is primarily of use when speed is specified by 
	 * components. Usually it is safe to leave this at its default
	 * value of 1.
	 *
	 * @param m Slope of mapping function
	 */
	void setLaunchSpeedMagnification(realtype m);

	/**
	 * Set minimum launch speed. This has an effect only when launch
	 * speed is specified by components.
	 *
	 * @param vs The minimum launch speed
	 */
	void setLaunchSpeedMinimum(realtype vs);
	
	/**
	 * Write out the header line.
	 *
	 * @see The meanings of the column labels are described in the 
	 * 		documentation for recordData().
	 */
	void recordDataInit();
	
	/**
	 * Record data to standard output.
	 *
	 * Each line consists of 3 to 8 columns separated by tab characters.
	 * The first line begins with a "#" character and indicates the
	 * meaning of each of the columns:
	 *
	 *  - \c t: time, where \f$t = 2\pi\f$ is one complete revolution
	 *  	of the two-body system (one "year")
	 *  - \c x: x-coordinate, relative to the center of mass
	 *  - \c y: y-coordinate, relative to the center of mass
	 *  - \c e (\f$\Delta\epsilon / \epsilon_0\f$): relative change in the
	 *  	pseudo-energy of the system
	 *  - \c d: total distance covered by the projectile over the
	 *  	integration period
	 *  - \c wE: winding number, with center of rotation at the earth
	 *  - \c wS: winding number, with center of rotation at the sun
	 *  - \c wC: winding number, with center of rotation at the center
	 *  	of mass.
	 *  	
	 * @param t Time value at this data point
	 * @param z Value of the data vector at this point
	 * @see enablePseudoEnergy()
	 * @see enableWindingNumber()
	 * @see enableTotalDistance()
	 * @see recordDataInit()
	 */
	void recordData(realtype t, N_Vector z);
	
	/**
	 * Perform the winding and distance calcuations (if requested)
	 * after each CVode() call. 
	 *
	 * @param t Time value at this data point
	 * @param z Value of the data vector at this point
	 */
	void postStep(realtype t, N_Vector z);
	
	/**
	 * Set capture radius. If the projectile comes within this distance
	 * of either the Earth or Sun, it is considered to have collided 
	 * and the simulation comes to a halt. Specifically, rhsFunc() sets
	 * an error message and returns an unrecoverable error to CVODE.
	 */
	void setCaptureRadius(realtype r);

	/**
	 * Enable recording pseudo-energy. This is the Jacobi constant
	 * of the motion, which must be constant for the trajectory to
	 * be valid. Enabling this option causes the relative change 
	 * versus the initial value to be recorded at each data point, 
	 * which increases computation time by approximately 20%.
	 *
	 * @todo Make this calculation more efficient, possibly by
	 * 		only performing it every so often.
	 * @see recordData()
	 * @see rhsFunc()
	 */
	void enablePseudoEnergy() { _shouldRecordPseudoEnergy = true; }
	
	/**
	 * Enable recording winding number.
	 *
	 * @see Winding
	 */
	void enableWindingNumber() { _shouldRecordWinding = true; }

	/**
	 * Enable recording total distance covered
	 */
	void enableTotalDistance() { _shouldRecordDistance = true; }
	
	/**
	 * Default constructor that picks some reasonable values for the 
	 * constants. 
	 */
	RestrictedThree():
		AbstractProblem(),
		_launch_radius(0.1),
		_capture_earth(0.005),
		_capture_sun(0.005),
		_shouldRecordPseudoEnergy(false),
		_pseudoEnergyIsDefined(false),
		_speed_parameter(SPEED_RELATIVE),
		_shouldRecordDistance(false),
		_shouldRecordWinding(false),
		_totalDistance(RCONST(0.0))
	{ 
		setMassFraction(0.3);
		setLaunchAngle(0.0);
		setLaunchSpeed(1.0);
		setTolerance(RCONST(1e-12), RCONST(1e-14));
		setInitialTime(0.0);
		setFinalTime(32.0);
		setTimestep(0.0001);
	}

  private:
	/**
	 * Update the initial vector to match launch radius,
	 * speed, and angle (the "polar" representation of initial velocity)
	 *
	 * @note This method of specifying launch parameters has been
	 * 		superseded by the component-wise method. See 
	 * 		setLaunchSpeedX(), setLaunchSpeedY(), and
	 * 		updateInitialVectorComp().
	 */
	void updateInitialVectorPolar();

	void updateInitialVectorComp();
	
	// "The Jacobi integral constant of the motion"
	// If the boolean is set, the realtype variable is maintained
	// by rhsFunc() and output by recordData()
	realtype _pseudoEnergy, _initialPseudoEnergy;
	bool _shouldRecordPseudoEnergy;
	bool _pseudoEnergyIsDefined;

	// Winding number (rotation about center of mass)
	bool _shouldRecordWinding;
	Winding _earthWinding, _sunWinding, _comWinding;

	// Total distance (length of line traced by projectile)
	bool _shouldRecordDistance;
	realtype _totalDistance;

	// Last x and y coords, needed for winding and distance
	realtype _lastx, _lasty;
	
	//// CONSTANTS ////
	// These values are not to be modified for the duration
	// of the integration. They are not declared const because
	// they may be modified in the lifetime of the object.
	
	// Radial distance of earth and sun (aka the primaries) from the
	// center of mass
	realtype _dist_com_earth, _dist_com_sun;
	
	// Mass of the primaries
	realtype _mass_earth, _mass_sun;

	// Radial distance from Earth where the projectile starts
	realtype _launch_radius;

	// Launch speed and angle (angle in degrees)
	realtype _launch_speed, _launch_angle;

	// How should the launch speed be interpreted?
	enum SpeedParam {
			SPEED_RELATIVE, // speed relative to earth, as a fraction
							// of escape speed from earth alone.
			SPEED_PARAMETER, // p = v^2 - v_E^2
			SPEED_ENERGY    // \epsilon (pseudo-energy)
	};
	SpeedParam _speed_parameter;

	// Launch speed in x and y dimension
	realtype _launch_speed_x, _launch_speed_y;

	// Slope of component-wise v-mapping function
	// @see updateInitialVectorComp()
	realtype _v_map_slope;

	// Minimum launch speed. Component-wise speed parameter values of 0
	// are assumed to be this value.
	realtype _min_launch_speed;

	// Radius of earth and sun
	realtype _capture_sun, _capture_earth;

};

