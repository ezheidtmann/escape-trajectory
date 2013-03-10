#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "LogTime.h"

#include <cvode.h>
#include <nvector_serial.h>
#include <cvode/cvode_dense.h>

#include <stdexcept>
#include <iostream>
using namespace std;

/**
 * A CVODE integrator, parameterized by problem type. The template
 * parameter should be a class implementing most of the functions
 * specified in the AbstractProblem class. In particular, the rhsFunc()
 * and denseJac() functions are required for CVODE to do its job.
 *
 * The modifiers may be called any number of times before calling 
 * initIntegration(), but certain parameters both in the integrator
 * and in the problem itself cannot be changed after initIntegration()
 * is run.
 *
 * The following example shows how this class is typically used:
 *
\code
Integrator<PROBLEM_TYPE>* integrator;
integrator = new Integrator<PROBLEM_TYPE>(p);

p.setFinalTime(23.723);
p.setTimestep(0.01);
integrator->setPointCount(350);
integrator->setCvodeMaxInternalSteps(1e8);
// ... other calls to change settings or parameters

integrator->initIntegration();
while (!integrator->done()) // done() returns false iff t > 23.723
{
    // move forward 0.01 time units
    // this function will also produce data output
	// (approximately 350 points as set above)
	integrator->runStep();
}
integrator->cleanup(); // free some memory.
\endcode
				
 *
 * @see The main() function shows how one might use this template class.
 */
template<typename P> 
class Integrator {
  public:
	
	/**
	 * Construct an integrator that integrates over the 
	 * given problem.
	 */
	Integrator(P &problem);

	/**
	 * Destructor
	 */
	~Integrator();
	
	/**
	 * Get a reference to the problem this integrator integrates.
	 */
	P& getProblem() { return _problem; }
	
	/**
	 * Set the integration method and iterator.
	 *
	 * @param integrationMethod A valid CVODE integration method: one
	 * 			of CV_ADAMS or CV_BDF.
	 * @param iterator A valid CVODE iterator: one of CV_NEWTON or
	 * 			CV_FUNCTIONAL.
	 */
	void setParameters(int integrationMethod, int iterator);

	/**
	 * Set the maximum number of steps that CVODE is allowed to take.
	 *
	 * Bigger values for this parameter tend to reduce "too much work"
	 * errors and have no negative effect on performance.
	 *
	 * @param mxsteps Maximum number of internal steps per timestep. 
	 * 		If not specified or negative, the old value is used.
	 */
	void setCvodeMaxInternalSteps(long int mxsteps = -1);

	/**
	 * Set the integration method.
	 *
	 * @param integrationMethod A valid CVODE integration method: one
	 *          of CV_ADAMS or CV_BDF.
	 */
	void setIntegrationMethod(int integrationMethod);

	/**
	 * Set the iterator.
	 *
	 * @param iterator A valid CVODE iterator: one of CV_NEWTON or
	 * 			CV_FUNCTIONAL.
	 */
	void setIterator(int iterator);

	/**
	 * Set approximate number of data points to record for standard
	 * linear time data recording.
	 *
	 * The first and last data points are always recorded. The 
	 * integrator then tries to record at least \arg count points, but
	 * may not be able to do it exactly. In particular, if \arg count
	 * is greater than the number of times runStep() is run, then 
	 * only that many points are recorded.
	 * than the value set here.
	 *
	 * Because IO is the most expensive part of most problems, this 
	 * can really reduce runtime when trajectories need to be known
	 * only roughly.
	 *
	 * @param count The desired number of data points. Default: 1000
	 */
	void setPointCount(unsigned int count) { _numPoints = count; }

	/**
	 * Enable logarithmic time for data recording. The parameters
	 * are set by accessing the logtime member variable.
	 */
	void enableLogTime() { _shouldUseLogTime = true; }

	LogTime logtime;
	
	/**
	 * Set up CVODE and get ready to integrate.
	 *
	 * @throws logic_error When a CVODE function fails
	 */
	void initIntegration();

	/**
	 * Run one step of the simulation. The size of the step
	 * is defined by the particular problem type P in use.
	 *
	 * @throws logic_error If the CVode() call fails
	 */
	void runStep();

	/**
	 * Get the current state vector. It is a pointer, but due to
	 * the structure of CVODE, any modifications to this vector 
	 * will not affect the integrator.
	 *
	 * @return An N_Vector object (which is actually a pointer)
	 * 		representing the current state of the system
	 * 		being integrated.
	 */
	N_Vector getCurrentVector() { return _v; }

	/**
	 * Get the current time of integration.
	 *
	 * @return The current value for \f$t\f$ in the integration
	 */
	realtype getCurrentTime() { return _currentTime; }

	/**
	 * Check if the integration is finished
	 *
	 * @return True if it's done, false if there's still
	 * 			more work to do.
	 */
	bool done() 
	{
		return _currentTime >= _problem.getFinalTime();
	}

	/**
	 * Do some cleanup stuff. This includes making sure the last data point
	 * was recorded and trying to free some memory.
	 *
	 * @todo Make sure this does a good job of freeing memory. Currently, 
	 * we expect most programs to quit soon after finishing the integration,
	 * so we don't worry too much about memory leaks.
	 */
	int cleanup();

	/**
	 * Construct a copy of the integrator. This may only be called
	 * before initIntegration(); an exception will be thrown otherwise.
	 */
	void makeCopy(Integrator<P>& otherInt);

	/**
	 * An error handler function for CVODE. The last parameter is 
	 * expected to be a pointer to an Integrator<P> object.
	 */
	static void cvodeErrorHandler 
		(int error_code, const char *module, const char *function,
		 char *msg, void *eh_data);
	
  private:
	P &_problem;

	N_Vector _v;
	
	int _integrationMethod;
	int _iterator;

	long int _cvode_mxsteps;

	void * _cvode_mem;

	realtype _currentTime;
	
	unsigned int _numPoints;
	/**
	 * The change in time between recorded data points
	 */
	realtype _recordingDeltaT;
	realtype _lastRecordedTime;

	bool _shouldUseLogTime;
};

template<typename P>
Integrator<P>::Integrator (P &problem) :
		_problem(problem),
		_currentTime(problem.getInitialTime()),
		_integrationMethod(CV_BDF),
		_iterator(CV_NEWTON),
		_cvode_mxsteps(1000000),
		_numPoints(1000),
		_shouldUseLogTime(false),
		_v(NULL),
		_cvode_mem(NULL)
{}
	
template<typename P>
void Integrator<P>::setIntegrationMethod( int integrationMethod )
{
	if (integrationMethod != CV_ADAMS && integrationMethod != CV_BDF)
		throw invalid_argument("Integrator<P>::setParameters(): Invalid integration method " + integrationMethod);
	_integrationMethod = integrationMethod;
}

template<typename P>
void Integrator<P>::setIterator( int iterator )
{
	if (iterator != CV_NEWTON && iterator != CV_FUNCTIONAL)
		throw invalid_argument("Integrator<P>::setParameters(): Invalid iterator " + iterator);
	_iterator = iterator;
}


template<typename P>
void Integrator<P>::setParameters(int integrationMethod, int iterator) 
{
	setIntegrationMethod(integrationMethod);
	setIterator(iterator);
}

template<typename P>	
void Integrator<P>::setCvodeMaxInternalSteps(long int mxsteps)
{
	if (mxsteps > 0)
		_cvode_mxsteps = mxsteps;

	if (_cvode_mem != NULL)
		CVodeSetMaxNumSteps(_cvode_mem, _cvode_mxsteps);
}

template<typename P>
void Integrator<P>::initIntegration() 
{
	int retval;
	int vector_size;
	N_Vector initial_vector;

	_cvode_mem = CVodeCreate(_integrationMethod, _iterator);

#ifdef DEBUG
	cout << "Ran CVodeCreate(" << _integrationMethod << ", " << _iterator << ") = " << _cvode_mem << endl;
#endif
	if (_cvode_mem == NULL) {
		throw logic_error("Integrator<P>::initIntegration(): CVodeCreate() failed.");
	}
	
	retval = CVodeSetFdata(_cvode_mem, (void*) &_problem);

	if (retval != CV_SUCCESS) {
		throw logic_error("Integrator<P>::initIntegration(): CVodeSetFdata() failed with status "+ retval);
	}

	// Copy the problem's initial conditions vector into our own
	// vector so we can use the same vector in CVodeMalloc() and
	// CVode()
	vector_size = _problem.getVectorSize();
	initial_vector = _problem.getInitialVector();
	if (_v == NULL)
		_v = N_VNew_Serial(vector_size);
	for (int i = 0; i < vector_size; i++)
		NV_Ith_S(_v, i) = NV_Ith_S(initial_vector, i);

#ifdef DEBUG
	cout << "Initial conditions vector:" << endl;
	N_VPrint_Serial(_problem.getInitialVector());
#endif

	// Prototype: CVodeMalloc(handle, function, t_0, initialVector,
	// 						  relTol, absTolVector)
	retval = CVodeMalloc(_cvode_mem, P::rhsFunc, _problem.getInitialTime(),
					_v, CV_SV,
					_problem.getRelativeTolerance(), 
 					_problem.getAbsoluteTolerances());

#ifdef DEBUG
	cout << "Ran CVodeMalloc(...) = " << retval << endl;
#endif

	if (retval != CV_SUCCESS) {
		throw logic_error("Integrator<P>::initIntegration(): CVodeMalloc() failed with status "+ retval);
	}

	retval = CVodeSetErrHandlerFn(
		_cvode_mem, Integrator<P>::cvodeErrorHandler, (void*) this);

	if (retval != CV_SUCCESS) {
		throw logic_error("Integrator<P>::initIntegration(): CVodeSetErrHandlerFn() failed with status "+ retval);
	}

	// If we're using the Newton iterator, we also need to use the
	// Jacobian. 
	if (_iterator == CV_NEWTON) {
		retval = CVDense(_cvode_mem, _problem.getVectorSize());

		if (retval != CVDENSE_SUCCESS) {
			throw logic_error("Integrator<P>::initIntegration(): CVDense() failed with status "+ retval);
		}

		retval = CVDenseSetJacFn(_cvode_mem, P::denseJac, (void*) &_problem);

		if (retval != CVDENSE_SUCCESS) {
			throw logic_error("Integrator<P>::initIntegration(): CVDenseSetJacFn() failed with status "+ retval);
		}
	}

	if (_numPoints != 0)
	{
	// Write the first line of output
		_problem.recordDataInit();

	// Record the first data point, because runStep() only does 
	// this after running a step.
		_problem.recordData(_problem.getInitialTime(), 
						_problem.getInitialVector());
		_lastRecordedTime = _problem.getInitialTime();
	
	// Set how often data points are recorded
		_recordingDeltaT = 
			(_problem.getFinalTime() - _problem.getInitialTime()) 
				/ (realtype) _numPoints;
	}
	else
	{
		_lastRecordedTime = _problem.getInitialTime();
		_recordingDeltaT = 2 * (_problem.getFinalTime() 
				- _problem.getInitialTime());
	}
	
	// Set the number of internal steps -- high values good!
	setCvodeMaxInternalSteps();
}

template<typename P>
void Integrator<P>::runStep()
{
	static int retval; // static for speed

#ifdef DEBUG2
	cout << "Running CVode(" << _cvode_mem << ", " << _currentTime + _problem.getTimestep()
			<< ", " << _v << ", " << &_currentTime << ", " << CV_NORMAL << ")" << endl;
#endif

	_problem.preStep(_currentTime, _v);
	
	// Advance one timestep, where timestep is defined by the
	// problem definition.
	static realtype newTime;
	newTime = _currentTime + _problem.getTimestep();
	newTime = newTime > _problem.getFinalTime() ? 
			_problem.getFinalTime() : newTime;
	retval = CVode(_cvode_mem, newTime, _v, &_currentTime, CV_NORMAL);

	_problem.postStep(_currentTime, _v);
	
	// If it's past time, record a data point
	if (_shouldUseLogTime)
	{
		if (logtime.atNextTimestep(_currentTime))
		{
			logtime.increment();
			_problem.recordData(_currentTime, _v);
		}
	}
	else
	{
		if (_currentTime >= _lastRecordedTime + _recordingDeltaT) {
			_problem.recordData(_currentTime, _v);
			_lastRecordedTime = _currentTime;
		}
	}

	// If the right-hand function failed, 
	// record the dying data point. The caller can get the error
	// message and error code from the problem.
	if (retval == CV_RHSFUNC_FAIL) {
		if (_currentTime != _lastRecordedTime && _numPoints != 0)
			_problem.recordData(_currentTime, _v);
	}
	if (retval != CV_SUCCESS) {
		throw logic_error(string("Integrator<P>::runStep(): CVode() failed"));
	}
}

template<typename P>
int Integrator<P>::cleanup()
{
	if (_currentTime != _lastRecordedTime && _numPoints != 0) 
	{
		_problem.recordData(_currentTime, _v);
		_lastRecordedTime = _currentTime;
	}

	if (_cvode_mem != NULL)
	{
		CVodeFree(&_cvode_mem);
		_cvode_mem = NULL;
	}
}

template<typename P>
Integrator<P>::~Integrator()
{
	cleanup();
	if (_v != NULL)
	{
		N_VDestroy(_v);
		_v = NULL;
	}
}

template<typename P>
void Integrator<P>::cvodeErrorHandler 
	(int error_code, const char *module, const char *function,
	 char *msg, void *eh_data)
{
	// We known the rhsFunc may fail, and already take care 
	// of it.
	if (error_code == CV_RHSFUNC_FAIL)
		return;

	// Other errors, on the other hand, are not supposed to
	// happen.
	
	cerr << "CVODE error code " << error_code << " in:" << endl;
	cerr << "\tModule: " << module << endl;
	cerr << "\tFunction: " << function << endl;
	cerr << "\tMessage: " << msg << endl;
}

#endif

