#ifndef DPI_H
#define DPI_H

// ***INCLUDES***
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

using namespace std;
using namespace Eigen;

// ***DEFINES***
#define FS 360.0
#define FC 8.0
#define PI 3.14159265358979323846

// ***FUNCTIONS***
vector<int> dpi_based_qrs_detector(VectorXf signal,float fs,float wnd, float p);
MatrixXf readRecording(const char*);
VectorXf hpf(VectorXf signal, float fc, float fs);
VectorXf getHalfWaveOfSignal(VectorXf signal);
MatrixXf getDenominators(int wnd, float p);
VectorXf convolve(VectorXf u, VectorXf v, string mode);
VectorXf smooth(VectorXf signal);
VectorXf derivative(VectorXf signal);
VectorXi findIndices(VectorXi logicVector);
tuple<VectorXi,VectorXi> zeroCrossing(VectorXf der, float threshold);
int swing(VectorXf dpi, VectorXi indPos, VectorXi indNeg);
VectorXf index(VectorXf vec, VectorXi ind);
int improveComplex(int indexStart, int nPoints, VectorXf signal);
#endif
