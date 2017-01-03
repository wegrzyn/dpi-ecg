#include "dpi.h"


MatrixXd readRecording(const char* pathToEcgFile){
    string line;
    unsigned int d1,d2;
    vector<double> vsignal,vlead1,vlead2;
    std::ifstream myfile(pathToEcgFile);
    if (myfile.is_open()){
            while (getline (myfile,line)){
            line.erase(remove(line.begin(), line.end(),' '), line.end());
            d1 = line.find('\t');
            d2 = line.find('\t',d1+1);
            vsignal.push_back(atof(line.substr(0,d1).c_str()));
            vlead1.push_back(atof(line.substr(d1,d2).c_str()));
            vlead2.push_back(atof(line.substr(d2,line.length()).c_str()));
        }
    myfile.close();
    }
    else cout << "Unable to open file";
    vsignal.insert(vsignal.end(), vlead1.begin(), vlead1.end());
    vsignal.insert(vsignal.end(), vlead2.begin(), vlead2.end());
    MatrixXd signal(vlead1.size(),3);
    signal = MatrixXd::Map(vsignal.data(),vlead1.size(),3);
//    cout << signal.block(0,0,10,3) << endl;
    return signal;
};

MatrixXd hpf(MatrixXd signal, double fc, double fs){
return signal;
}
