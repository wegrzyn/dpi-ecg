#ifndef DPI_H
#define DPI_H

// ***INCLUDES***
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

using namespace std;
using namespace Eigen;

// ***DEFINES***
#define FS 250.0
#define FC 8.0
#define PI 3.14159265358979323846

// ***FUNCTIONS***
MatrixXf readRecording(const char*);
VectorXf hpf(VectorXf signal, float fc, float fs);
VectorXf getHalfWaveOfSignal(VectorXf signal);
MatrixXf getDenominators(int wnd, float p);
#endif
