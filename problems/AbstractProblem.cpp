#include "AbstractProblem.h"

#include <iostream>
#include <stdexcept>
using namespace std;

void AbstractProblem::setTolerance( realtype relativeTolerance,
	const vector<realtype> &absoluteTolerances )
{
	setRelativeTolerance(relativeTolerance);
	setAbsoluteTolerance(absoluteTolerances);
}

void AbstractProblem::setTolerance( realtype relativeTolerance,
	realtype absoluteTolerance )
{
	setRelativeTolerance(relativeTolerance);
	setAbsoluteTolerance(absoluteTolerance);
}

void AbstractProblem::setRelativeTolerance( realtype relativeTolerance )
{
	_relTol = relativeTolerance;
}

void AbstractProblem::setAbsoluteTolerance( const vector<realtype> &absoluteTolerances )
{
	if (absoluteTolerances.size() != getVectorSize()) {
		throw invalid_argument("AbstractProblem::setTolerance(): absoluteTolerances vector is of improper length.");
//		throw invalid_argument("AbstractProblem::setTolerance(): absoluteTolerances vector is of length "+ absoluteTolerances.size() +", should be "+ getVectorSize());
	}
	
	if (_absTol != NULL) {
		N_VDestroy_Serial(_absTol);
		_absTol = NULL;
	}
	
	_absTol = N_VNew_Serial(getVectorSize());
	realtype * data = NV_DATA_S(_absTol);
	for (int i = 0; i < getVectorSize(); i++) {
		data[i] = absoluteTolerances.at(i);
	}
	
}

void AbstractProblem::setAbsoluteTolerance( realtype absoluteTolerance )
{
	setAbsoluteTolerance(vector<realtype>(getVectorSize(), absoluteTolerance));
}

void AbstractProblem::setInitialVector( const vector<realtype> &v ) 
{
	if (v.size() != getVectorSize()) {
		throw invalid_argument("AbstractProblem::setInitialVector(): vector is of improper length.");
		//throw invalid_argument("AbstractProblem::setInitialVector(): vector is of length "+ v.size() +", should be "+ getVectorSize());
	}

	if (_initialVector != NULL) {
		N_VDestroy_Serial(_initialVector);
		_initialVector = NULL;
	}

	_initialVector = N_VNew_Serial(getVectorSize());
	realtype * data = NV_DATA_S(_initialVector);
	for (int i = 0; i < getVectorSize(); i++) {
		data[i] = v.at(i);
	}
}
	
void AbstractProblem::setInitialVector(N_Vector v)
{
	if (v == NULL)
		throw invalid_argument("AbstractProblem::setInitialVector(): vector is NULL!");
	
	if (NV_LENGTH_S(v) != getVectorSize()) {
		throw invalid_argument("AbstractProblem::setInitialVector(): vector is of improper length.");
	}

	if (_initialVector != NULL) {
		N_VDestroy_Serial(_initialVector);
		_initialVector = NULL;
	}

	_initialVector = N_VNew_Serial(getVectorSize());
	for (int i = 0; i < getVectorSize(); i++) {
		NV_Ith_S(_initialVector, i) = NV_Ith_S(v, i);
	}
}

void AbstractProblem::recordData(realtype t, N_Vector z) 
{
	cout << "Time: " << t << endl;
	N_VPrint_Serial(z);
}
