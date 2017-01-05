#include "dpi.h"

MatrixXf readRecording(const char* pathToEcgFile){
    string line;
    unsigned int d1,d2;
    vector<float> vsignal,vlead1,vlead2;
    // Importing text file to Eigen Array:
//    std::ifstream myfile(pathToEcgFile);
//    unsigned int numberOfLines=0;
//
////    ArrayXd signal;
//    while (getline (myfile,line)){
//        numberOfLines++;
//    }
//    ArrayXd signal(numberOfLines,1);
//    myfile.clear();
//    myfile.seekg(0, ios::beg);
//    for(unsigned int i=0;getline(myfile,line); i++){
//        line.erase(remove(line.begin(), line.end(),' '), line.end());
//        d1 = line.find('\t');
//        d2 = line.find('\t',d1+1);
//        signal.row(i) << atof(line.substr(d1,d2).c_str());
//    }
//
//    myfile.close();
// Importing text file to Eigen Matrix
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
    MatrixXf signal(vlead1.size(),3);
    signal = MatrixXf::Map(vsignal.data(),vlead1.size(),3);
    cout << signal.block(0,0,10,3) << endl;
    return signal;
};

VectorXf hpf(VectorXf signal, float fc, float fs){

    FFT<float> fft;
    size_t nRows = signal.size();
    float freqPoint = fs/nRows;
    int nPointsToFc = floor(fc/freqPoint);

    VectorXf raised_cosine(nPointsToFc);
    raised_cosine.setLinSpaced(nPointsToFc,1/float(nPointsToFc),1);
    raised_cosine = 0.5-0.5*(raised_cosine.array()*PI).cos();
    VectorXf freqDomainFilter(nRows);
    freqDomainFilter.setOnes(nRows);
    freqDomainFilter.head(nPointsToFc) = raised_cosine;
    freqDomainFilter.tail(nPointsToFc) = raised_cosine.reverse();

    VectorXcf spectrum(nRows);
    fft.fwd( spectrum, signal);
    spectrum = spectrum.cwiseProduct(freqDomainFilter);
    fft.inv( signal,spectrum);
    return signal;
}

VectorXf getHalfWaveOfSignal(VectorXf signal){
    Array<bool, Dynamic,1> result = signal.array()>0.0;
    signal = signal.array() * result.cast<float>();
    return signal;
}

MatrixXf getDenominators(int wnd, float p){
    VectorXf vec(wnd);
    vec.setLinSpaced(wnd+1,1.0,float(wnd)+1.0);
    vec = 1.0/vec.array();
    vec = vec.array().pow(1.0/p);
    MatrixXf dpi_denom(wnd,wnd);
    dpi_denom = vec * VectorXf::Ones(wnd).transpose();
    dpi_denom = dpi_denom.triangularView<Lower>();

    cout << "\nDPI Denominators:" << endl;
    cout << dpi_denom.block(0,0,10,10) << endl;
    return dpi_denom;
}


