#ifndef LOGTIME_H
#define LOGTIME_H

#include <cvode.h>
#include <math.h>

/**
 * Class implementing logarithmic time counting. This is simply an 
 * implementation of this recurrence relation:
 * 
 * \f[
 *   t_n = r t_{n-1}
 * \f]
 *
 * The factor \f$r\f$ is set by setFactor() and indirectly by setFactorByN(),
 * while \f$t_0\f$ is set by setFirstTime().
 */
class LogTime {
  public:
	/**
	 * Construct a default LogTime object. By default, \f$r=2\f$ and
	 * \f$t_0=1\f$.
	 */
	LogTime():
		_factor(2.0),
		_next_time(1.0)
	{}

	// mutators
	
	/**
	 * Set the first time value \f$t_0\f$.
	 */
	void setFirstTime( realtype t ) { _next_time = t; }

	/**
	 * Set the factor \f$r\f$ by which time increases each time.
	 */
	void setFactor( realtype r ) { _factor = r; }

	/**
	 * Set \f$r\f$ indirectly by the parameter \f$n\f$ in the relation
	 * \f[
	 *   r = 2^{1 / n}
	 * \f]
	 */
	void setFactorByN( realtype n ) { _factor = pow(2.0, 1.0 / n); }

	// methods for use by integrator
	
	/**
	 * Find out if the given time value is at or beyond the next timestep.
	 */
	bool atNextTimestep( realtype t ) { return t >= _next_time; }

	/**
	 * Increment the internal time counter. This should be called whenever
	 * the timed event occurs.
	 */
	void increment() { _next_time *= _factor; }

  private:
	realtype _factor;
	realtype _next_time;

};

#endif

