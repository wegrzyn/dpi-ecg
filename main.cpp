#include "dpi.h"

int main(){

  //Import signal form text file to MatrixXf,
  //where columns are: time[probes], lead1, lead2
  MatrixXf signal;
  signal = readRecording("100.txt");

  // Get only first lead to ECG analysis
  VectorXf lead1(VectorXf::Map(signal.col(1).data(), signal.rows()));

  // Highpass filtering with cut-off frequency equal to 8Hz
  lead1 = hpf(lead1,FC,FS);
  lead1(0) = -100;
  cout << lead1.head(10) << endl;
  // Half-wave signal
  lead1 = getHalfWaveOfSignal(lead1);

  cout << lead1.head(10) << endl;

  getDenominators(1800, 5.0);
}
