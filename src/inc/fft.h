#ifndef FFT_H
#define FFT_H
#include <iostream>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>

/* Note:
 * Please, refer to links below to see how it should works:
 * https://leonidov.su/fft-lection-notes/ (Russian)
 * https://www.ee.iitm.ac.in/~csr/teaching/pg_dsp/lecnotes/fft.pdf (English)
 * */

using std::vector;
using std::complex;

#ifndef PI
#define PI 3.141592653589793L
#endif

/* Declares type of output vector in read operations */
enum class OutputVector {
    Complex,
    Absolute
};

class FourierTransform {
public:
    /* Sets input vector for forward Fourier transform calculation */
    void setForwardInputVector(vector<complex<double>>&);
    /* Sets input vector for inverse Fourier transform calculation */
    void setInverseInputVector(vector<complex<double>>&);
    /* Calculates Fourier transform in forward direction */
    void forwardCalculation();
    void forwardCalculation(vector<complex<double>>&);
    /* Calculates Fourier transform in inverse direction */
    void inverseCalculation();
    void inverseCalculation(vector<complex<double>>&);
    /* Reads output vector */
    vector<complex<double>> readForwardOutputVector(OutputVector);
    vector<complex<double>> readInverseOutputVector(OutputVector);
private:
    union Flag {
        unsigned char all;
        struct Bit {
            unsigned char forwardInputVectorDefined :1;
            unsigned char inverseInputVectorDefined :1;
            unsigned char forwardOutputVectorDefined :1;
            unsigned char inverseOutputVectorDefined :1;
            unsigned char RESERVED :4;
        } bit;
    } flag;
    vector<complex<double>> forwardInputVector;
    vector<complex<double>> inverseInputVector;
    vector<complex<double>> forwardOutputVector;
    vector<complex<double>> inverseOutputVector;
    void checkSize(vector<complex<double>>&);
    unsigned binaryInverse(unsigned, unsigned);
    void forward();
    void inverse();
};

#endif /* FFT_H */
