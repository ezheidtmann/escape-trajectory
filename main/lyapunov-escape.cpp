#include "Integrator.h"
#include "RestrictedThree.h"

#include <getopt.h>
#include <math.h>

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
using namespace std;

#ifndef PROBLEM_TYPE
#define PROBLEM_TYPE RestrictedThree
#endif

/**
 * \mainpage
 *
 * \section intro_sec Introduction
 *
 * TODO
 *
 * \section param_sec Command-line parameters
 *
 * The executable supports many command-line parameters. The following
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
 * - \c -T: \n
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
 *   problem in use. It is generally safe to leave this at its default.
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

void setup_integrator(Integrator<PROBLEM_TYPE>* integrator, int argc, char** argv)
{

	PROBLEM_TYPE &p = integrator->getProblem();
	// Let's do some option parsing!
	int arg, opt_index;
	string optstring = "p:I:i:v:r:T:t:b:e:s:N:m:L:n:";
	optstring += p.getOptionString();

	optind = 0; // Reset internal getopt variable
				// so option parsing works more than once.
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
			{"logtime-n",          1, 0, 'n'},
			{0, 0, 0, 0}
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
				exit(1);
				
			default:
				if (!p.handleOption(arg))
				{
					cerr << "Error: Don't know what to do with option: -" 
						 << arg << endl;
				}
		} // end switch
				
	} // end while
}

realtype exponential_regression(vector<realtype>& v, realtype delta_t,
				realtype* correlation_square = NULL)
{
	// Exponential regression is hard. Let's do linear instead!
	
	realtype sxx = 0, 
			 sx = 0,
			 sy = 0, 
			 sxy = 0,
			 syy = 0;
	realtype N = (realtype) v.size();

	#ifdef DEBUG
	cout << "Regression got vector: ";
	copy(v.begin(), v.end(), ostream_iterator<realtype>(cout, " "));
	cout << endl;
	#endif
	
	sx = (N) * (N+1) / 2.0 * delta_t;
	sxx = (N) * (N+1) * (2*N+1) / 6.0 * delta_t;
	
	int n = 1; // This must be initially 1 to match the closed-form
			   // expressions above.
	
	for (vector<realtype>::iterator i = v.begin(); 
					i != v.end(); i++, n++)
	{
		*i = log(*i); // natural log

		sy += *i;
		sxy += n * delta_t * *i;
		syy += *i * *i;
	}
	
	realtype lambda;
	lambda = - (-N * sxy + sx * sy) / (-sx*sx + N*sxx);
	
	if (correlation_square != NULL)
		*correlation_square = (N * sxy - sx * sy) * (N * sxy - sx * sy)/
			( (N * sxx - sx*sx) * (N * syy - sy*sy) );
	
	#ifdef DEBUG
	cout << "Regression found lambda: " << lambda << endl;
	cout << "Regression found R^2 = " << *correlation_square << endl;
	#endif
	
	return lambda;
}

realtype compute_average(vector<realtype>& v)
{
	return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

realtype compute_stddev(vector<realtype>& v, realtype avg)
{
	realtype square_avg = 0.0;

	for (vector<realtype>::iterator i = v.begin();
					i != v.end(); i++)
	{
		square_avg += *i * *i;
	}
	square_avg /= v.size();

	return sqrt( square_avg - avg*avg );
}

int main (int argc, char** argv) 
{
	PROBLEM_TYPE *orig_p, *mod_p;
	Integrator<PROBLEM_TYPE> *orig_i, *mod_i;

	realtype initial_deflection, hypercube_radius;
	realtype norm, lambda, rsquared;

	int i, steps_before_readjust, total_deflections;
	steps_before_readjust = 50;

	N_Vector ov, mv, diffv;

	unsigned int vector_size;

	vector<realtype> diff_vec, lambda_vec, rsquared_vec;

	// TODO: replace with getopt parameter
	hypercube_radius = 1e-9;
	total_deflections = 50;

	try {
		
		for (int delta_num = 0; delta_num < total_deflections; delta_num++)
		{
			lambda_vec.clear();
			rsquared_vec.clear();
			orig_p = new PROBLEM_TYPE;
			mod_p  = new PROBLEM_TYPE;
			orig_i = new Integrator<PROBLEM_TYPE>(*orig_p);
			mod_i  = new Integrator<PROBLEM_TYPE>(*mod_p);
			
			// setup both identically
			setup_integrator(orig_i, argc, argv);
			setup_integrator(mod_i, argc, argv);

			orig_i->setPointCount(0);
			mod_i->setPointCount(0);

			vector_size = orig_p->getVectorSize();
	
			orig_i->initIntegration();
			int steps = 0;
			ov = orig_i->getCurrentVector();
			mv = N_VNew_Serial(vector_size);
			// add deflection: within a hypercube with side length 2*hypercube_radius
			diffv = N_VNew_Serial(vector_size);
			for (int i = 0; i < vector_size; i++)
			{
				NV_Ith_S(diffv, i) = ((realtype) rand()) / ((realtype) RAND_MAX) * hypercube_radius;
			}
			initial_deflection = sqrt(N_VDotProd(diffv, diffv));
			N_VLinearSum(1, diffv, 1, ov, mv);
			mod_p->setInitialVector(mv);
			mod_i->initIntegration();
			N_VDestroy(mv);
			mv = mod_i->getCurrentVector();
#ifdef DEBUG2
			cout << "initial difference: ";
			N_VPrint_Serial(diffv);
			cout << "initial original: ";
			N_VPrint_Serial(ov);
			cout << "initial modded: "; 
			N_VPrint_Serial(mv);
#endif
			diff_vec.clear();
			while (!orig_i->done()) 
			{
				if (steps >= steps_before_readjust)
				{
					steps = -1;
					// compute diffv = mv - ov
					//N_VLinearSum(1, mv, -1, ov, diffv); // (already computed below)
					// compute norm = ||diffv||
					//norm = sqrt(N_VDotProd(diffv, diffv)); // (already computed below)
					// change length
					N_VScale(initial_deflection/norm, diffv, diffv);
					// calculate new state vector for modified intergrator
					N_VLinearSum(1, diffv, 1, ov, mv);
					// and restart the intergration from that point
					mod_i->cleanup();
					mod_p->setInitialVector(mv);
					mod_p->setInitialTime(orig_i->getCurrentTime());
					mod_i->initIntegration();
					lambda = exponential_regression(diff_vec, 
									mod_p->getTimestep(), &rsquared);
					lambda_vec.push_back(lambda);
					rsquared_vec.push_back(rsquared);
					diff_vec.clear();
				}
					
				orig_i->runStep();
				mod_i->runStep();
				++steps;
				N_VLinearSum(1.0, mv, -1.0, ov, diffv);
				norm = sqrt(N_VDotProd(diffv, diffv));
				diff_vec.push_back(norm);

#ifdef DEBUG2
				cout << "original vector: " << endl;
				N_VPrint_Serial(ov);
				cout << "modified vector: " << endl;
				N_VPrint_Serial(mv);
				cout << "Diff vector: " << endl;
				N_VPrint_Serial(diffv);
#endif
			}

			orig_i->cleanup();
			mod_i->cleanup();
			N_VDestroy(diffv);

			realtype avg = compute_average(lambda_vec);
			realtype stddev = compute_stddev(lambda_vec, avg);
			realtype avg_r2 = compute_average(rsquared_vec);
			realtype stddev_r2 = compute_stddev(rsquared_vec, avg_r2);

			cout << avg << "\t" 
				 << stddev << "\t" 
				 << avg_r2 << "\t"
				 << stddev_r2 << endl;
		}
	}
	catch (exception &e) {
		// Only complain if the problem isn't aware of it.
		if (mod_p->getErrorCode() == 0 && orig_p->getErrorCode() == 0)
		{
			cerr << "Caught exception!" << endl;
			cerr << e.what() << endl;
		}
	}

	// The problem's error code is 0 if nothing bad happened, and nonzero
	// if something did go awry. So let's use it as the program's return
	// value!
	
	return mod_p->getErrorCode();
}

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

