#include "AbstractProblem.h"

#include <math.h>

/**
 * Implements a two-body problem where one of the bodies is of
 * infinitesimal mass.
 */
class SingleBody : public AbstractProblem {
  public:
	/**
	 * The CVODE callback function implementing the ODEs that represent
	 * this problem. See rhsFuncAbstract() for details on the behavior
	 * and expectation of this function.
	 *
	 * @param t The current time in the integration
	 * @param z The current "position" vector, having values
	 * 			\f$(x, x^\prime)\f$ or
	 * 			\f$(x, V_x)\f$.
	 * @param z_prime The current "velocity" vector, having values
	 * 			\f$(x^\prime, x^{\prime\prime})\f$ or
	 * 			\f$(x^\prime, V_x^\prime)\f$.
	 * @param myself A pointer to a SingleBody object, which is
	 * 			treated as if it were the "this" pointer.
	 */
	static int rhsFunc( realtype t, N_Vector z, N_Vector z_prime, void *myself);

	/**
	 * The dense Jacobian function. This is not actually implemented
	 * in this problem.
	 *
	 * @return Always returns 0
	 */
	static int denseJac( long int N, DenseMat J, realtype t, 
		N_Vector y, N_Vector fy, void *jac_data, N_Vector tmp1, 
		N_Vector tmp2, N_Vector tmp3 );
	
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
	 * This problem consists of 2 ODEs, so this function always
	 * returns 2.
	 */
	unsigned int getVectorSize() { return 2; }

	/**
	 * Additional command-line options specific to this problem.
	 *
	 * Unless otherwise noted, every option takes as an argument
	 * a single floating-point value.
	 *
	 *  - \c -R: launch radius
	 *  - \c -V: launch speed
	 *  - \c -c: capture radius
	 *
	 * @return An option string representing these options, in the format
	 * required by \c getopt(3).
	 */
	string getOptionString() { return "R:V:c:"; }

	/**
	 * Handle command-line options. 
	 *
	 * @see getOptionString()
	 */
	bool handleOption( int arg );

	/**
	 * Set launch speed. This speed is relative to the massive body.
	 */
	void setLaunchSpeed(realtype s);

	/**
	 * Record data to stdout. The data produced is suitable for parametric
	 * plotting with gnuplot.
	 *
	 * Each line consists of 3 columns separated by tab characters.
	 *
	 *  - t: time, where \f$t = 2\pi\f$ is one complete revolution
	 *  	of the two-body system (one "year")
	 *  - x: distance of the projectile from the massive body.
	 *  - Vx: speed of the projectile
	 *  	
	 * @todo Reduce the number of times this is called --
	 * 		removing the body of this function decreases runtime
	 * 		by almost a factor of 7.
	 * @param t Time value at this data point
	 * @param z Value of the data vector at this point
	 * @see enablePseudoEnergy()
	 */
    void recordData(realtype t, N_Vector z);

	/**
	 * Set capture radius. If the projectile comes within this distance
	 * of either the massive body, it is considered to have collided 
	 * and the simulation comes to a halt. Specifically, rhsFunc() sets
	 * an error message and returns an unrecoverable error to CVODE.
	 */
	void setCaptureRadius(realtype r);

	/**
	 * Default constructor that picks some reasonable values for the 
	 * constants. 
	 */
	SingleBody():
		AbstractProblem(),
		_launch_radius(0.1),
		_capture_radius(0.005)
	{ 
		setLaunchSpeed(1.0);
		setTolerance(RCONST(1e-12), RCONST(1e-14));
		setInitialTime(0.0);
		setFinalTime(32.0);
		setTimestep(0.0001);
	}

  private:
	/**
	 * Update the initial vector to match launch radius,
	 * speed, and angle.
	 */
	void updateInitialVector();
	
	//// CONSTANTS ////
	// These values are not to be modified for the duration
	// of the integration. They are not declared const because
	// they may be modified in the lifetime of the object.
	
	// Radial distance from the massive body where the projectile 
	// starts
	realtype _launch_radius;

	// Launch speed 
	realtype _launch_speed;

	// Radius of massive body 
	realtype _capture_radius;

};

