#ifndef DPI_H
#define DPI_H

// ***INCLUDES***
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// ***DEFINES***
#define FS 250.0
#define FC 8.0


// ***FUNCTIONS***
MatrixXd readRecording(const char*);
MatrixXd hpf(MatrixXd signal, double fc, double fs);

#endif
