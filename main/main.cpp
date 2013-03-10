#include "Integrator.h"
#include "RestrictedThree.h"
#include "SingleBody.h"
#include "TwoBody.h"
#include "CurvedTwoBody.h"

#include <getopt.h>
#include <math.h>

#include <iostream>
#include <stdexcept>
using namespace std;

#ifndef PROBLEM_TYPE
/**
 * The name of the class implementing the problem to be integrated.
 * This can be changed at compile time to create a different executable
 * without modifying the source directly.
 */
#define PROBLEM_TYPE RestrictedThree
#endif

/** @file main.cpp
 *
 * This file includes the main() function, which parses options and
 * starts the integration.
 */

/**
 * \mainpage
 *
 * \section intro_sec Introduction
 *
 * The CVODE Trajectory Generator is a command-line program based on
 * the CVODE library which
 * integrates a set of differential equations to calculate a trajectory 
 * in phase space. The initial conditions for the system are given on
 * the command line and the differential equations are provided by a 
 * "problem" object. The supported problem types are listed at the end 
 * of this section.
 *
 * The Generator is capable of integrating the differential equations
 * for only one set of initial conditions per execution. The initial 
 * conditions are specified on the command-line, either by the 
 * \c --initial-vector option or by options specific to a problem. The
 * format of output data may vary with the problem being used, but all
 * current problems output columns of data (separated by tabs)
 * representing the state of 
 * the system at many intermediate points during the integration.
 *
 * \section example_sec The Escape problem: an example
 *
 * The Generator was originally written for the Escape problem. In this
 * problem, we simulate the path taken by a small particle of infinitesimal
 * mass under the Newtonian gravitational influence of two massive
 * bodies which orbit each other in concentric circles. The problem is
 * a restricted version of the classical three-body problem; it is 
 * implemented in the RestrictedThree class.
 *
 * In the current build system, \c escape and \c rthree are identical 
 * executables which integrate the RestrictedThree problem. With the 
 * default parameters and initial conditions, \c escape produces output
 * like this:
 * 
\verbatim
 $ ./escape
#t      x       y
0       0.8     0
0.0321  0.867814        0.0224674
0.0641  0.923033        0.0448505
0.0961  0.97147 0.0672085
0.1281  1.01532 0.0895327
0.1602  1.05587 0.111886
0.1923  1.09369 0.134191
0.2244  1.12926 0.156445
0.2565  1.16294 0.178643
...
\endverbatim
 *
 * The first column is the system time, where \f$t = 0\f$ is the start
 * of integration, and \f$t = 2\pi\f$ marks the time when the 2-body system 
 * completes its first revolution. The next two columns give the familiar
 * Cartesian coordinates of the projectile's position in the plane. When
 * additional options are specified, additional columns may appear:
 *
\verbatim
 $ ./escape -E -G -X -0.24 -Y -0.584 
#t      x       y       e
0       0.661989        -0.092494       0
0.0321  0.638866        -0.1253 8.50784e-11
0.0641  0.619952        -0.145911       1.07883e-10
0.0961  0.603042        -0.159546       1.13102e-10
0.1281  0.587256        -0.16827        1.22007e-10
0.1602  0.572063        -0.173198       1.23826e-10
0.1923  0.557199        -0.174981       1.24331e-10
0.2244  0.542441        -0.174081       1.24602e-10
0.2565  0.527617        -0.170828       1.2514e-10
...
\endverbatim
 * 
 * Here the initial conditions are specified by the \c -X and \c -Y
 * options. The \c -G option specifies how the \c -X and \c -Y options 
 * are to be interpreted, and the \c -E option adds an additional column 
 * giving the deviation from initial pseudoenergy. (all these options are 
 * documented in RestrictedThree::getOptionString())
 *
 * The frequency of output is controlled with the \c -N option and the
 * upper bound of integration is specified with \c -e. For example:
\verbatim
 $ ./escape -e 44 -N 4
#t      x       y
0       0.8     0
11.0001 2.28949 3.84313
22.0002 -0.592399       0.147827
33.0002 -4.77448        0.0824081
44      -3.82037        -1.39882
 $
\endverbatim
 * 
 * The \c escape executable and others like it are infrequently used
 * alone like this. Usually they are used in the so-called "dispatch" 
 * scripts developed for the Escape problem, which can be used to run
 * \c escape millions of times with different
 * intial conditions and distribute the computational task across
 * many computers using Apple's Xgrid software.
 * 
 * \section param_sec Command-line parameters
 *
 * The generator supports many command-line parameters. The following
 * parameters are valid for all problems:
 *
 * - \c -p: \n
 *   Select problem to simulate. This is currently not supported; for now
 *   each problem produces its own executable.
 *   
 * - \c -I, \c --integration-method: \n
 *   Select integrator. This should be  either "bdf" or "adams".
 *   - CVODE Manual: <a href="http://www.llnl.gov/casc/sundials/documentation/cv_guide/node6.html#SECTION00651000000000000000">5.5.1 CVODE initialization and deallocation functions</a> - CVodeCreate()
 *   
 * - \c -i, \c --iterator:\n
 *   Select iterator. This should be either "newton" or "functional".
 *   - CVODE Manual: <a href="http://www.llnl.gov/casc/sundials/documentation/cv_guide/node6.html#SECTION00651000000000000000">5.5.1 CVODE initialization and deallocation functions</a> - CVodeCreate()
 *   
 * - \c -v, \c --initial-vector: \n
 *   A vector of initial conditions, separated by commas. The
 *   list must be of the length required by the problem in use.
 *   
 * - \c -r, \c --relative-tolerance: \n
 *   Relative tolerance. Values around \c 1e-8 or \c 1e-12 have
 *   been found to work well. 
 *   - CVODE Manual: <a href="http://www.llnl.gov/casc/sundials/documentation/cv_guide/node6.html#SECTION00652000000000000000">5.5.2 Advice on choice and use of tolerances</a>
 *   
 * - \c -T:\n
 *   Absolute tolerances, as a comma-separated list. The list
 *   must be of the length required by the problem in use.
 *   - CVODE Manual: <a href="http://www.llnl.gov/casc/sundials/documentation/cv_guide/node6.html#SECTION00652000000000000000">5.5.2 Advice on choice and use of tolerances</a>
 *   
 * - \c -t, \c --absolute-tolerance: \n
 *   Absolute tolerance, as a scalar. The list of absolute tolerances
 *   is filled with copies of this value.
 *   - CVODE Manual: <a href="http://www.llnl.gov/casc/sundials/documentation/cv_guide/node6.html#SECTION00652000000000000000">5.5.2 Advice on choice and use of tolerances</a>
 *   
 * - \c -b, \c --begin: \n
 *   Beginning time. The unit of time is defined by the
 *   problem in use. It is generally safe to leave this at its default;
 *   in fact, the generator may have trouble if this parameter is used.
 *   Default: 0
 *   
 * - \c -e, \c --end: \n
 *   Ending time.
 *
 * - \c -s, \c --timestep: \n
 *   Internal timestep for CVODE integration. Smaller values allow the
 *   integrator to go farther without dying with a "too much work" error.
 *   But make this too small, and you gain no precision but lose 
 *   performance. 
 *  
 * - \c -N, \c --number-points, \c --data-length: \n
 *   Number of points to record. Default: 1000
 *
 * - \c -L, \c --logtime-start: \n
 *   Enable logarithmic time for data recording and set \f$t_0\f$.
 *
 * - \c -n, \c --logtime-n: \n
 *   Enable logarithmic time for data recording and set \f$n\f$ so
 *   \f$r=2^{1/n}\f$
 *
 * - \c -m, \c --max-internal-steps: \n
 *   Maximum number of internal CVODE steps to take per timestep. Bigger
 *   values for this parameter tend to eliminate "too much work" errors
 *   and seem to have no negative effect on performance. Default: 1000000
 *
 * Other options are specified by the particular problem in use.
 *
 * - RestrictedThree::getOptionString()
 * - SingleBody::getOptionString()
 * - TwoBody::getOptionString()
 * - CurvedTwoBody::getOptionString()
 */

vector<realtype> split_realtype(const string &sep,string text);

/**
 * Parse command-line options, initialize problem and integrator, and
 * finally run the integration.
 */
int main (int argc, char** argv) 
{
	PROBLEM_TYPE p;

	try {
		Integrator<PROBLEM_TYPE>* integrator;
		integrator = new Integrator<PROBLEM_TYPE>(p);

		// Let's do some option parsing!
		int arg, opt_index;
		string optstring = "p:I:i:v:r:T:t:b:e:s:N:m:L:n:";
		optstring += p.getOptionString();
		while (1) 
		{
			static struct option long_options[] = {
				{"integration-method", 1, 0, 'I'},
				{"iterator",           1, 0, 'i'},
				{"initial-vector",     1, 0, 'v'},
				{"relative-tolerance", 1, 0, 'r'},
				{"absolute-tolerance", 1, 0, 't'},
				{"begin",              1, 0, 'b'},
				{"end",                1, 0, 'e'},
				{"timestep",           1, 0, 's'},
				{"number-points",      1, 0, 'N'},
				{"data-length",        1, 0, 'N'},
				{"max-internal-steps", 1, 0, 'm'},
				{"logtime-start",      1, 0, 'L'},
				{"logtime-n",          1, 0, 'n'}
			};

			arg = getopt_long(argc, argv, optstring.c_str(), long_options, &opt_index);
			if (arg == -1)
				break;

			switch (arg) {
				case 'p': // problem
					cout << "-p: Sorry, choosing problem at runtime isn't supported. (yet)" << endl;
					break;

				case 'I': // Integrator
					if (string(optarg) == string("bdf"))
						integrator->setIntegrationMethod(CV_BDF);
					else if (string(optarg) == string("adams"))
						integrator->setIntegrationMethod(CV_ADAMS);
					else
						cerr << "-I: Valid integration methods are \"bdf\" and \"adams\"" << endl;
					break;

				case 'i': // iterator
					if (string(optarg) == string("newton"))
						integrator->setIterator(CV_NEWTON);
					else if (string(optarg) == string("functional"))
						integrator->setIterator(CV_FUNCTIONAL);
					else
						cerr << "-i: Valid iterators are \"newton\" and \"functional\"" << endl;
					break;

				case 'v': // initial conditions vector
					p.setInitialVector(split_realtype(",", optarg));
					break;

				case 'r': // relative tolerance (float scalar)
					p.setRelativeTolerance(atof(optarg));
					break;

				case 'T': // absolute tolerance (comma-separated vector)
					p.setAbsoluteTolerance(split_realtype(",", optarg));
					break;

				case 't': // absolute tolerance (float scalar)
					p.setAbsoluteTolerance(atof(optarg));
					break;

				case 'b': // beginning time
					p.setInitialTime(atof(optarg));
					break;

				case 'e': // ending time
					p.setFinalTime(atof(optarg));
					break;

				case 's': // timestep
					p.setTimestep(atof(optarg));
					break;

				case 'N':
					integrator->setPointCount(atoll(optarg));
					break;

				case 'm':
					integrator->setCvodeMaxInternalSteps(atoll(optarg));
					break;

				case 'L':
					integrator->enableLogTime();
					integrator->logtime.setFirstTime(atof(optarg));
					break;

				case 'n':
					integrator->enableLogTime();
					integrator->logtime.setFactorByN(atof(optarg));
					break;
				
				case ':':
				case '?':
					cerr << "Error: Something messed up in the options." << endl;
					return 1;
					
				default:
					if (!p.handleOption(arg))
					{
						cerr << "Error: Don't know what to do with option: -" 
							 << arg << endl;
					}
			}
					
		}

		// End option parsing.
		
#ifdef DEBUG
		cerr << "Address of initial vector: " << p.getInitialVector() << endl;
#endif
		integrator->initIntegration();
		while (!integrator->done()) 
		{
			integrator->runStep();
		}
		integrator->cleanup();
	}
	catch (exception &e) {
		// Only complain if the problem isn't aware of it.
		if (p.getErrorCode() == 0)
		{
			cerr << "Caught exception!" << endl;
			cerr << e.what() << endl;
		}
	}

	// The problem's error code is 0 if nothing bad happened, and nonzero
	// if something did go awry. So let's use it as the program's return
	// value!
	
	return p.getErrorCode();
}

/**
 * Split a string into a vector of \c realtype values. 
 *
 * In the following example, the vector \a v is filled with the four
 * values 1.0, 2.7, 3.2, and 4.03:
 * 
\code
vector<realtype> v;
v = split_realtype(",", "1.0,2.7,3.2,4.03");
\endcode
 *
 * The implementation does not try to be particularly efficient. It
 * does, however, do a pretty good job of being correct.
 *
 * @param sep Delimeter on which to split. \a sep is stripped
 * 		from the output.
 * @param text Text to split.
 */
vector<realtype> split_realtype(const string &sep,string text)
{
    vector<realtype> numbers;
    string::size_type end;
    do
    {
        end = text.find(sep);
        if (end == string::npos)
            end = text.length() + 1;

        numbers.push_back(atof(text.substr(0,end).c_str()));
        text.replace(0,end+sep.length(),"");

    } while (text.length());
    return numbers;
}

