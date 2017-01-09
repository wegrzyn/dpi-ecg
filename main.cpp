#include "dpi.h"

int main(){
//    Import signal form text file to MatrixXf,
//    where columns are: time[probes], lead1, lead2
    VectorXf time,lead1,lead2;
    VectorXi ann,qrs;

    tie(time,lead1,lead2) = readRecording("100.txt");
//    cout << "lead1:\n" << lead1.head(10) << endl;
    qrs = dpi_based_qrs_detector(lead1,FS, 1800.0, 5.0);
//    cout << qrs.transpose() << endl;
//        Get vector of annotations from file
    ann = readAnnotation("ann100.txt");
//    cout << ann.transpose() << endl;

    float sensitivity,precision,accuracy;
    int wnd = int(FS*0.150);
    tie(accuracy,sensitivity,precision) = validateDetector(ann, qrs, wnd);
}
