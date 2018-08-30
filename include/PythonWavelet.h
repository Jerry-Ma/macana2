#ifndef PYTHONWAVELET_H
#define PYTHONWAVELET_H
#include <Python.h>
#include <string>
#include <cstdlib>

#include "nr3.h"

class PythonWavelet{
	private:
//		size_t nSamples;
//		PyObject *waveModule;
		PyObject *waveFunc;
		PyObject *splineFunc;
//		PyObject *waveFunc2d;
		std::string waveName;
		double cutStd;
	public:
		PythonWavelet ();
		void filter_data(double *data, size_t nDetectors, size_t nSamples);
//		void filter_data2D (double *data, size_t nx, size_t ny);
		void setWaveName (std::string wavelet);
		std::string getWaveName ();
		void setCutStd (double cutStd);
		double getCutStd();
		~PythonWavelet();

		MatDoub fit2dspline (MatDoub signal, VecDoub ra, VecDoub dec, MatDoub weight, double coverage,  size_t cellSize);

};

#endif
