#include "AbstractProblem.h"

#include <math.h>

/**
 * Implements the general two-body problem in a curved 2-dimensional universe.
 *
 * This problem object was developed by Martha I. Roseberry, with programming
 * assistance from Evan Heidtmann.
 *
 * @todo Add proper documentation.
 *
 * @author Martha Roseberry
 * @author Evan Heidtmann
 */
class CurvedTwoBody: public AbstractProblem {
  public:
	/**
	 * The CVODE callback function implementing the ODEs that represent
	 * this problem. See rhsFuncAbstract() for details on the behavior
	 * and expectation of this function.
	 *
	 * @param t The current time in the integration
	 * @param z The current "position" vector
	 * @param z_prime The current "velocity" vector
	 * @param myself A pointer to a TwoBodyCurved object, which is
	 * 			treated as if it were the "this" pointer.
	 */
	static int rhsFunc( realtype t, N_Vector z, N_Vector z_prime, void *myself);

	/**
	 * The dense Jacobian function. 
	 */
	static int denseJac( long int N, DenseMat J, realtype t, 
		N_Vector y, N_Vector fy, void *jac_data, N_Vector tmp1, 
		N_Vector tmp2, N_Vector tmp3 );
	
	/**
	 * This problem consists of 8 ODEs, so this function always
	 * returns 8.
	 */
	unsigned int getVectorSize() { return 8; }

	/**
	 * Additional command-line options specific to this problem.
	 *
	 * Unless otherwise noted, every option takes as an argument
	 * a single floating-point value.
	 *
	 *  - \c -P: initial separation \f$\rho_0\f$
	 *  - \c -R: universe radius
	 *  - \c -V: difference in the initial velocities or particles one and two
	 *  - \c -c: capture radius
	 *  - \c -f: mass fraction
	 *
	 * @return An option string representing these options, in the format
	 * required by \c getopt(3).
	 */
	string getOptionString() { return "P:R:V:c:f:"; }

	/**
	 * Handle command-line options. 
	 *
	 * @see getOptionString()
	 */
	bool handleOption( int arg );

	/**
	 * Set \f$\Delta v\f$
	 */
	void setDeltaV(realtype s);

	/**
	 * Write the first line of output.
	 */
	void recordDataInit();
	
	/**
	 * Record data to stdout. 
	 *
	 * Each line consists of five columns separated by tab characters.
	 *
	 *  - t: time
	 *  - j1: \f$\phi\f$ for particle one
	 *  - q1: \f$\theta\f$ for particle one
	 *  - j2: \f$\phi\f$ for particle two
	 *  - q2: \f$\theta\f$ for particle two
	 *  	
	 * @param t Time value at this data point
	 * @param z Value of the data vector at this point
	 */
    void recordData(realtype t, N_Vector z);

	/**
	 * Set capture radius. If the bodies come within this distance
	 * of each other, they are considered to have collided and the 
	 * simulation comes to a halt. Specifically, rhsFunc() sets
	 * an error message and returns an unrecoverable error to CVODE.
	 */
	void setCaptureRadius(realtype r);

	/**
	 * Set initial separation. This is the value \f$\rho_0\f$.
	 */
	void setInitialSeparation(realtype r);
	
	/**
	 * Set sphere radius. This is the radius of the universe.
	 */
	void setSphereRadius(realtype r);

	/**
	 * Set fraction of the total mass given to the first body.
	 *
	 * @param f Mass fraction, must be between 0 and 1.
	 * @throw invalid_argument If f > 1 or f < 0
	 */
	void setMassFraction(realtype f);

	/**
	 * Default constructor that picks some reasonable values for the 
	 * constants. 
	 */
	CurvedTwoBody():
		AbstractProblem(),
		_capture_radius(0.005),
		_mass_1(0.3),
		_mass_2(0.7),
		_initial_separation(1.0),
		_delta_v(1.0),
		_sphere_radius(1.0)
	{ 
		updateInitialVector();
		setTolerance(RCONST(1e-12), RCONST(1e-14));
		setInitialTime(0.0);
		setFinalTime(32.0);
		setTimestep(0.01);
	}

  private:
	/**
	 * Update the initial vector to match initial separation,
	 * velocity difference, and mass fraction.
	 */
	void updateInitialVector();
	
	//// CONSTANTS ////
	// These values are not to be modified for the duration
	// of the integration. They are not declared const because
	// they may be modified in the lifetime of the object.
	
	// Initial separation between the two masses
	realtype _initial_separation;

	// The difference between the inital velocities
	realtype _delta_v;

	// Capture radius
	realtype _capture_radius;

	// Mass of first and second body, respectively
	realtype _mass_1, _mass_2;
	
	//Radius of the 2 sphere
	realtype _sphere_radius;

};

