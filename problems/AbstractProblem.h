#ifndef ABSTRACTPROBLEM_H
#define ABSTRACTPROBLEM_H

#include <unistd.h>
#include <stdlib.h>

#include <cvode.h>
#include <nvector_serial.h>
#include <sundials_dense.h>

#include <vector>
#include <string>
using namespace std;

/**
 * Abstract parent class for a CVODE problem. 
 *
 * All CVODE problems should be implemented as a child class of this
 * class. Although all methods are virtual and thus can be overridden,
 * those functions which are more commonly reimplemented in child classes
 * are grouped at the top of this file.
 *
 * The most important methods are rhsFuncAbstract(), denseJacAbstract(), 
 * and recordData() and friends. Child classes will typically also add
 * some methods to set parameters and options.
 */
class AbstractProblem {
  public:
	
	/// @name Methods which child classes should override
	//@{
	
	/**
	 * The right-hand side of the vector ODE, a callback from CVODE.
	 * Classes implementing this function should name it \em rhsFunc().
	 *
	 * This function matches the prototype of \c CVRhsFn:
	 * 
	 * <tt>typedef int (*CVRhsFn) (realtype t, N_Vector y, N_Vector ydot, 
	 *   						 void *f_data)</tt>
	 *
	 * This is the function \f$f(t,z)\f$ in the expression 
	 * \f$z^\prime = f(t,z)\f$ as required by CVODE. Both \f$z^\prime\f$
	 * and \f$z\f$ are vectors.
	 *
	 * CVODE requires that this be a static function, but the last parameter
	 * allows this to act like a regular member function. The "myself" 
	 * pointer is treated as the "this" pointer is within a regular member 
	 * function.
	 *
	 * @param t The current value of the independent variable, time in this 
	 * 			case.
	 * @param z The input vector, as generated by initial conditions or
	 * 			by CVODE
	 * @param z_prime The output vector, to be filled by this function.
	 * @param myself A pointer to an AbstractProblem object where the function
	 * 			should get its configuration values. 
	 */
	static int rhsFuncAbstract( realtype t, N_Vector z, N_Vector z_prime, void *myself );

	/**
	 * The dense Jacobian function. Child classes should name this function
	 * \em denseJac().
	 */
	static int denseJacAbstract( long int N, DenseMat J, realtype t, 
		N_Vector y, N_Vector fy, void *jac_data, N_Vector tmp1, 
		N_Vector tmp2, N_Vector tmp3 );

	/**
	 * Get the size of the vector used by this problem.
	 */
	virtual unsigned int getVectorSize() = 0;

	/**
	 * Get additional options supported by this problem. The option
	 * string should be in the format required by getopt(3).
	 */
	virtual string getOptionString() { return ""; }

	/**
	 * Handle options specified by getOptionString().
	 *
	 * @param arg The option argument, as returned by getopt_long().
	 * @return True if the option was successfully handled,
	 * 			false otherwise.
	 */
	virtual bool handleOption(int arg) { return false; }

	/**
	 * Initialize data recording. Concrete classes might use this
	 * method to output a header line or initalize a buffer.
	 *
	 * @see recordData()
	 */
	virtual void recordDataInit() {};
	
	/**
	 * Record data point.
	 *
	 * @param t Time value at this data point
	 * @param z Value of the data vector at this point
	 */
	virtual void recordData(realtype t, N_Vector z);

	/**
	 * Finish up recording data. Concrete classes might use this
	 * method to flush a buffer.
	 *
	 * @see recordData()
	 */
	virtual void recordDataFinish() {}

	/**
	 * Called before each CVode() call. Implementation is optional.
	 */
	virtual void preStep(realtype t, N_Vector z) {}

	/**
	 * Called after each CVode() call. Implementation is optional.
	 */
	virtual void postStep(realtype t, N_Vector z) {}
	
	//@}
	
	/**
	 * Set tolerance for the integration.
	 *
	 * @param relativeTolerance The scalar relative tolerance
	 * @param absoluteTolerances The absolute tolerances for all the 
	 * 			vector elements. The size of this vector needs to
	 * 			be equal to the value returned by getVectorSize().
	 */
	virtual void setTolerance( realtype relativeTolerance, 
			const vector<realtype> &absoluteTolerances );

	/**
	 * Set tolerance for the integration. A convenience function.
	 *
	 * @param relativeTolerance The scalar relative tolerance
	 * @param absoluteTolerance The scalar absolute tolerance. 
	 * 			This value is simply copied into a vector of
	 * 			the proper length.
	 */
	virtual void setTolerance( realtype relativeTolerance,
			realtype absoluteTolerance );

	/**
	 * Set relative tolerance
	 */
	virtual void setRelativeTolerance( realtype relativeTolerance );

	/**
	 * Set absolute tolerance. This function sets all elements of
	 * the absolute tolerances vector to the same value.
	 *
	 * @param absoluteTolerance The absolute tolerance. Smaller values 
	 * 		mean better precision. Default: 1.0e-4
	 */
	virtual void setAbsoluteTolerance( realtype absoluteTolerance );

	/**
	 * Set absolute tolerances. This allows each element of the
	 * vector to have a different absolute tolerance.
	 */
	virtual void setAbsoluteTolerance( 
			const vector<realtype> &absoluteTolerances );
	
	/**
	 * Get relative tolerance.
	 */
	virtual realtype getRelativeTolerance() { return _relTol; }

	/**
	 * Get absolute tolerances. Modifications made to the returned
	 * value will change absolute tolerances for the integration,
	 * if the modifications are made before Integrator::initIntegration()
	 * is called.
	 */
	virtual N_Vector getAbsoluteTolerances() { return _absTol; }
	
	/**
	 * Set the initial value for t in this problem.
	 */
	virtual void setInitialTime(realtype t) { _t_0 = t; }
	
	/**
	 * Get the initial value for t in this problem.
	 */
	virtual realtype getInitialTime() { return _t_0; }

	/**
	 * Set the final value for t in this problem
	 */
	virtual void setFinalTime(realtype t) { _t_f = t; }
	
	/**
	 * Get the final value for t in this problem
	 */
	virtual realtype getFinalTime() { return _t_f; }

	/**
	 * Get vector of initial conditions in this problem.
	 */
	virtual N_Vector getInitialVector() { return _initialVector; }

	/**
	 * Set vector of initial conditions in this problem.
	 */
	virtual void setInitialVector(const vector<realtype> &v);

	/**
	 * Set vector of initial conditions in this problem. The values
	 * are copied into the internal vector, so the parameter may be
	 * safely deallocated after calling.
	 */
	virtual void setInitialVector(N_Vector v);
	
	/**
	 * Set timestep. Reducing this value only produces more recorded
	 * data points; the precision of the integration should be independent
	 * of this parameter. However, it has been observed that small values
	 * produce different trajectories than large values.
	 *
	 * @param step Difference between successive time values requested
	 * 		from CVODE.
	 */
	virtual void setTimestep(realtype step) { _timestep = step; }

	/**
	 * Get timestep
	 */
	virtual realtype getTimestep() { return _timestep; }

	/**
	 * Set error message. This should be called from rhsFunc just before
	 * it returns an unrecoverable error.
	 *
	 * @param code Numerical error code. Generally, this becomes the return
	 * 		value of the program itself.
	 * @param msg String describing the error
	 */
	virtual void setError(int code, const string &msg) 
		{ _error_code = code; _error_msg = msg; }

	/**
	 * Get error message.
	 */
	virtual string getErrorMessage() { return _error_msg; }

	/**
	 * Get error code. This should be zero if nothing went 
	 * wrong, and nonzero if something did go wrong
	 *
	 * @return The error code, suitable for the program's
	 * 		return value.
	 */
	virtual int getErrorCode() { return _error_code; }
	
	
	/**
	 * Default constructor.
	 *
	 * Subclass constructors need to run the line of code which is
	 * commented out below, or something like it. It can't be done
	 * here because it relies on the pure virtual getVectorSize()
	 * function.
	 */
	AbstractProblem() :
		_relTol(RCONST(0.0)), _absTol(NULL),
		_initialVector(NULL), _timestep(RCONST(1.0)), 
		_t_0(RCONST(0.0)), _t_f(RCONST(1000.0)),
		_error_code(0)
	{
		// Subclasses should run this to prevent
		// null tolerances being passed to CVODE.
		//
		//setTolerance(RCONST(1e-4), RCONST(1e-3));
	}

	/**
	 * Destructor. Try to free some memory.
	 */
	~AbstractProblem() 
	{
		if (_initialVector != NULL)
		{
			N_VDestroy(_initialVector);
			_initialVector = NULL;
		}
		
		if (_absTol != NULL)
		{
			N_VDestroy(_absTol);
			_absTol = NULL;
		}
	}
	
  private:
	realtype _relTol;
	N_Vector _absTol;

	N_Vector _initialVector;

	realtype _timestep;
	realtype _t_0, _t_f;

	string _error_msg;
	int _error_code;
	
	//// CONSTANTS //// ( to be added by child classes )
	// These values are not to be modified for the duration
	// of the integration. They are not declared const because
	// they may be modified in the lifetime of the object.
};

#endif // ABSTRACT_PROBLEM_H

