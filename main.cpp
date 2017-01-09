#include "dpi.h"

int main(){
//    Import signal form text file to MatrixXf,
//    where columns are: time[probes], lead1, lead2
    MatrixXf signal;
    VectorXi annotation;
    signal = readRecording("100.txt");
    // Get only first lead to ECG analysis
    VectorXf lead1(VectorXf::Map(signal.col(1).data(), signal.rows()));
//    cout << "lead1:\n" << lead1.head(10) << endl;
    vector<int> qrs;
    qrs = dpi_based_qrs_detector(lead1,FS, 1800.0, 5.0);

    for( size_t i = 0; i < qrs.size(); i++ )
        printf("%d, ", qrs[i]);
//        Get vector of annotations from file
    annotation = readAnnotation("ann.txt");
//    cout << annotation << endl;

}
