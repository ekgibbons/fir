#include <complex>
#include <stdexcept>
#include <iostream>
#include <vector> 
#include <algorithm>

#include <math.h>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "numpy_wrappers.hpp"

extern "C" {
#include "remez_scipy.h"
#include "remez_jano.h"
#include "four1.h"
}

#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

#if _MSC_VER
using boost::uint8_t;
#endif

namespace py = boost::python;
namespace np = boost::python::numpy;

int CheckPad(int num)
{

    int result;
    if (num % 2)
    {
	num -= (num %2);
    }

    double remainder = log10((double)num)/log10(2.0);
    remainder = ceil(remainder);
    result = (int)pow(2,remainder);
    
    return result;
}


np::ndarray RemezScipy(const int &numTaps,
		       const np::ndarray &bandsNp,
		       const np::ndarray &desiredNp,
		       const np::ndarray &weightNp,
		       const int &type,
		       const int &maxiter,
		       const int &gridDensity)
{

    int numBands = (int)bandsNp.shape(0);
    int numDesired = (int)desiredNp.shape(0);
    int numWeight = (int)weightNp.shape(0);

    
    double *h = static_cast<double *>(malloc(numTaps*sizeof(double)));

    for (int i = 0; i < numTaps; i++)
    {
    	h[i] = 0;
    }

    double *bands = static_cast<double *>(malloc(numBands*sizeof(double)));
    double *desired = static_cast<double *>(malloc(numDesired*sizeof(double)));
    double *weight = static_cast<double *>(malloc(numWeight*sizeof(double)));
    
    Numpy2Array(bandsNp, bands, numBands);
    Numpy2Array(desiredNp, desired, numDesired);
    Numpy2Array(weightNp, weight, numWeight);

    numBands /= 2; /* need to divide by two, the bands array is just the corners */
    
    /* remez function call */
    remez_scipy(h, numTaps, numBands, bands, desired, weight, type,
    		maxiter, gridDensity);

    np::ndarray result = Array2Numpy(h,numTaps);
    
    free(h);
    free(bands);
    free(desired);
    free(weight);       
    
    return result;
}


np::ndarray fft(const np::ndarray &inputNp,
		int nPoints=0)
{
    int numOrig = (int)inputNp.shape(0);

    if (numOrig > nPoints)
    {
	nPoints = numOrig;
    }
    nPoints = CheckPad(nPoints);

    double *transformArray = static_cast<double *>(calloc(2*(size_t)nPoints,
							  sizeof(double)));
	
    std::string dtype = py::extract<std::string>(py::str(inputNp.get_dtype()));
    
    if (dtype.compare("complex128") == 0)
    {
	CNumpy2CArray(inputNp, transformArray, numOrig);
    }
    else if (dtype.compare("float64") == 0)
    {
	Numpy2CArray(inputNp, transformArray, numOrig);
    }
    
    four1(transformArray-1,nPoints,1);

    np::ndarray result = CArray2CNumpy(transformArray,nPoints);

    free(transformArray);

    return result;    
}

np::ndarray ifft(const np::ndarray &inputNp,
		 int nPoints=0)
{
    int numOrig = (int)inputNp.shape(0);

    if (numOrig > nPoints)
    {
	nPoints = numOrig;
    }
    nPoints = CheckPad(nPoints);

    double *transformArray = static_cast<double *>(calloc(2*(size_t)nPoints,
							  sizeof(double)));

    std::string dtype = py::extract<std::string>(py::str(inputNp.get_dtype()));
    
    if (dtype.compare("complex128") == 0)
    {
	CNumpy2CArray(inputNp, transformArray, numOrig);
    }
    else if (dtype.compare("float64") == 0)
    {
	Numpy2CArray(inputNp, transformArray, numOrig);
    }

    four1(transformArray-1,nPoints,-1);

    for (int ii = 0; ii < nPoints; ii++)
    {
	transformArray[ii] /= nPoints;
    }
    
    np::ndarray result = CArray2CNumpy(transformArray,nPoints);

    free(transformArray);

    return result;    
}

   

BOOST_PYTHON_MODULE(sigtools)
{
    Py_Initialize();
    np::initialize();
    py::def("RemezScipy",RemezScipy);
    py::def("fft",fft);
    py::def("ifft",ifft);

}
