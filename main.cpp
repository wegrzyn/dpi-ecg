#include "dpi.h"

int main(){
  MatrixXd signal;
  signal = readRecording("100.txt");
  signal = hpf(signal,FC,FS);

  cout << signal.block(0,0,10,3) << endl;
}
