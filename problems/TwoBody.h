#include "AbstractProblem.h"

#include <math.h>

/**
 * Implements the general two-body problem. In this problem, two bodies
 * of mass \f$m_1\f$ and \f$m_2\f$ exist in a flat 2-dimensional universe.
 * They are initially separated by a distance \f$\rho_0\f$ and are moving
 * at velocities \f$v_1\f$ and \f$v_2\f$, directed perpendicular to the 
 * line connecting the two bodies and in opposite directions. The magnitude
 * of the initial velocities are related by conservation of linear momentum:
 * 
 * \f{eqnarray*}
 * 	m_1 v_1 &=& m_2 v_2
 * \f}
 *
 * This particular implementation is parameterized by the initial
 * separation \f$\rho_0\f$ and the value \f$\Delta v = v_1 - v_2\f$
 */
class TwoBody : public AbstractProblem {
  public:
	/**
	 * The CVODE callback function implementing the ODEs that represent
	 * this problem. See rhsFuncAbstract() for details on the behavior
	 * and expectation of this function.
	 *
	 * @param t The current time in the integration
	 * @param z The current "position" vector, having values
	 * 			\f$(x_1, x_2, y_1, y_2, x_1^\prime, x_2^\prime,
	 * 				y_1^\prime, y_2^\prime)\f$.
	 * @param z_prime The current "velocity" vector, having values
	 * 			\f$(x_1^\prime, x_2^\prime, y_1^\prime, y_2^\prime,
	 * 				x_1^{\prime\prime}, x_2^{\prime\prime},	
	 * 				y_1^{\prime\prime}, y_2^{\prime\prime})\f$.
	 * @param myself A pointer to a TwoBody object, which is
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
	 *  - \c -R: \f$\rho_0\f$: The initial separation of the bodies
	 *  - \c -V: \f$\Delta v = v_1 - v_2\f$
	 *  - \c -c: capture radius
	 *
	 * @see setInitialSeparation()
	 * @see setCaptureRadius()
	 * @return An option string representing these options, in the format
	 * required by \c getopt(3).
	 */
	string getOptionString() { return "R:V:c:f:"; }

	/**
	 * Handle command-line options. 
	 *
	 * @see getOptionString()
	 */
	bool handleOption( int arg );

	/**
	 * Set launch speed. This speed is relative to the massive body.
	 */
	void setDeltaV(realtype s);

	/**
	 * Write the first line of output.
	 */
	void recordDataInit();
	
	/**
	 * Record data to stdout. The data produced is suitable for parametric
	 * plotting with gnuplot.
	 *
	 * Each line consists of 3 columns separated by tab characters.
	 *
	 *  - t: time
	 *  - x_1: X-coordinate of the first body
	 *  - y_1: Y-coordinate of the first body
	 *  - x_2: X-coordinate of the second body
	 *  - y_2: Y-coordinate of the second body
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
	TwoBody():
		AbstractProblem(),
		_capture_radius(0.005),
		_mass_1(0.3),
		_mass_2(0.7),
		_initial_separation(1.0),
		_delta_v(1.0)
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

	/// The difference between the inital velocities:
	/// \f$\Delta v = v_1 - v_2\f$
	realtype _delta_v;

	// Capture radius
	realtype _capture_radius;

	// Mass of first and second body, respectively
	realtype _mass_1, _mass_2;

};

