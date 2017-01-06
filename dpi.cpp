#include "dpi.h"

vector<int> dpi_based_qrs_detector(VectorXf signal,float fs, float wnd, float p){

    int n285 = floor(0.285*fs);         //The period of the highest possible heart rate (210 BPM)
    int nWnd = floor(wnd/1000.0*fs);    //Computation window
    int nWnd285 = n285 + nWnd;       //Computation window with additional 285 ms
    int currQrs = 1,prevQrs = 0;
    int indPrevQrs = 4;
    vector<int> qrs;

    VectorXf hhecg(nWnd),dpi(nWnd),der(nWnd);
    MatrixXf dpi_denom(nWnd,nWnd);

    // Highpass filtering in the frequency domain (fc = 8 Hz)
    signal = hpf(signal,FC,fs);

    //Prepare triangle matrix of DPI denominators:
    dpi_denom = getDenominators(nWnd,p);



    while (indPrevQrs < signal.size()-nWnd285){
        // Half wave of filtered signal
        hhecg = getHalfWaveOfSignal(signal.segment(indPrevQrs,nWnd));
        // Dynamic Plosion Index in window
        dpi = (dpi_denom * hhecg).cwiseInverse();
        // Smoothing DPI
        dpi = smooth(dpi);
        // Derivative computation
        der = derivative(dpi);



        indPrevQrs +=200;
        qrs.push_back(indPrevQrs);
        prevQrs++;
        currQrs++;
    }
    return qrs;
}


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

    VectorXf raisedCosine(nPointsToFc);
    raisedCosine.setLinSpaced(nPointsToFc,1/float(nPointsToFc),1);
    raisedCosine = 0.5-0.5*(raisedCosine.array()*PI).cos();
    VectorXf freqDomainFilter(nRows);
    freqDomainFilter.setOnes(nRows);
    freqDomainFilter.head(nPointsToFc) = raisedCosine;
    freqDomainFilter.tail(nPointsToFc) = raisedCosine.reverse();

    VectorXcf spectrum(nRows);
    fft.fwd(spectrum, signal);
    spectrum = spectrum.cwiseProduct(freqDomainFilter);
    fft.inv(signal,spectrum);
    return signal;
}

VectorXf getHalfWaveOfSignal(VectorXf signal){
    Array<bool, Dynamic,1> result = signal.array()>0.0;
    signal = signal.array().abs() * result.cast<float>();
    return signal;
}

MatrixXf getDenominators(int wnd, float p){
    VectorXf vec(wnd);
    vec.setLinSpaced(wnd,1.0,float(wnd));
    vec = 1.0/vec.array();
    vec = vec.array().pow(1.0/p);
    MatrixXf dpi_denom(wnd,wnd);
    dpi_denom = vec * RowVectorXf::Ones(wnd);
    dpi_denom = dpi_denom.triangularView<Lower>();
//    cout << "\nDPI Denominators:" << endl;
//    cout << dpi_denom.block(0,0,10,10) << endl;
    return dpi_denom;
}

VectorXf convolve(VectorXf vec1, VectorXf vec2, string mode){
    size_t n1 = vec1.size(), n2 = vec2.size();
    VectorXf vec1Wide = VectorXf::Zero(n1+2*n2-2);
    VectorXf conv = VectorXf::Zero(n1+n2-1);
    vec1Wide.head(n2-1).fill(vec1(0));
    vec1Wide.tail(n2-1).fill(vec1(n1-1));
    vec1Wide.segment(n2-1,n1) = vec1;
    for (size_t i = 0; i < (n1+n2-1); i++){
        conv(i) = (vec1Wide.segment(i,n2).cwiseProduct(vec2.reverse())).sum();
    };
    if (mode.compare("same")==0){
        return conv.segment(floor(n2/2),n1);
    }
    else if (mode.compare("valid")==0){
        return conv.segment(n2-1,n1-n2+1);
    }
    else{  // mode: 'full'
        return conv;
    }
}

VectorXf smooth(VectorXf signal){
    VectorXf s(5);
    s << 0.2,0.2,0.2,0.2,0.2;
    return convolve(signal, s, "same");
}

VectorXf derivative(VectorXf signal){
    VectorXf d(5);
    d << -1,0,0,0,1;
    return convolve(signal,d,"same");
}


