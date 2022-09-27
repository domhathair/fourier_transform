#include "fft.h"

void FourierTransform::setForwardInputVector(vector<complex<double>> &input) {
    forwardInputVector.assign(input.begin(), input.end());
    flag.bit.forwardInputVectorDefined = true;
    flag.bit.forwardOutputVectorDefined = false;
}

void FourierTransform::setInverseInputVector(vector<complex<double>> &input) {
    inverseInputVector.assign(input.begin(), input.end());
    flag.bit.inverseInputVectorDefined = true;
    flag.bit.inverseOutputVectorDefined = false;
}

void FourierTransform::forwardCalculation() {
    if (!flag.bit.forwardInputVectorDefined)
        forward();
    else
        std::cerr << "Error: forward input vector is not defined.\n";
}

void FourierTransform::forwardCalculation(vector<complex<double>> &input) {
    setForwardInputVector(input);
    forward();
    flag.bit.forwardOutputVectorDefined = true;
}

void FourierTransform::inverseCalculation() {
    if (flag.bit.inverseInputVectorDefined) {
        inverse();
    } else
        std::cerr << "Error: inverse input vector is not defined.\n";
}

void FourierTransform::inverseCalculation(vector<complex<double>> &input) {
    setInverseInputVector(input);
    inverse();
    flag.bit.inverseOutputVectorDefined = true;
}

vector<complex<double>> FourierTransform::readForwardOutputVector(OutputVector type) {
    if (flag.bit.forwardOutputVectorDefined)
        switch (type) {
        case OutputVector::Complex:
            return forwardOutputVector;
        case OutputVector::Absolute:
            vector<complex<double>> realForwardOutputVector(forwardOutputVector.size());
            for (unsigned cnt = 0; cnt < forwardOutputVector.size(); cnt++)
                realForwardOutputVector[cnt] = abs(forwardOutputVector[cnt]) / forwardOutputVector.size();
            return realForwardOutputVector;
        }
    std::cerr << "Error: first calculate the forward Fourier transform.\n";
    return vector<complex<double>>(0.0, 0.0);
}

vector<complex<double>> FourierTransform::readInverseOutputVector(OutputVector type) {
    if (flag.bit.inverseOutputVectorDefined)
        switch (type) {
        case OutputVector::Complex:
            return inverseOutputVector;
        case OutputVector::Absolute:
            vector<complex<double>> realInverseOutputVector(inverseOutputVector.size());
            for (unsigned cnt = 0; cnt < inverseOutputVector.size(); cnt++)
                realInverseOutputVector[cnt] = abs(inverseOutputVector[cnt]) / inverseOutputVector.size();
            return realInverseOutputVector;
        }
    std::cerr << "Error: first calculate the inverse Fourier transform.\n";
    return vector<complex<double>>(0.0, 0.0);
}

void FourierTransform::checkSize(vector<complex<double>> &input) {
    /* Make sure that vector size is power of two */
    if (log2(input.size()) != unsigned(log2(input.size()))) {
        unsigned size = pow(2, unsigned(log2(input.size())) + 1);
        input.resize(size);
    }
}

unsigned FourierTransform::binaryInverse(unsigned input, unsigned size) {
    unsigned backward = 0;
    for (unsigned shift = 0; shift < size; shift++)
        backward |= ((input >> shift) & 1) << (size - 1 - shift);
    return backward;
}

void FourierTransform::forward() {
    checkSize(forwardInputVector);
    forwardOutputVector.resize(forwardInputVector.size());
    unsigned sizeBit = trunc(log2(forwardOutputVector.size() - 1) / log2(2)) + 1;

    for (unsigned cnt = 0; cnt < forwardOutputVector.size(); cnt++)
        forwardOutputVector[cnt] = forwardInputVector[binaryInverse(cnt, sizeBit)];

    unsigned numberOfTreeLevels = log2(forwardOutputVector.size()) / log2(2);
    for (unsigned fcnt = 0; fcnt < numberOfTreeLevels; fcnt++) {
        unsigned powOfTwo = pow(2, fcnt + 1);
        complex<double> Wm(cos(2 * PI / powOfTwo), -sin(2 * PI / powOfTwo));
        for (unsigned scnt = 0; scnt < forwardOutputVector.size(); scnt += powOfTwo) {
            complex<double> W = 1;
            for (unsigned tcnt = 0; tcnt < powOfTwo / 2; tcnt++) {
                complex<double> odd = W * forwardOutputVector[scnt + tcnt + powOfTwo / 2];
                complex<double> even = forwardOutputVector[scnt + tcnt];
                forwardOutputVector[scnt + tcnt] = even + odd;
                forwardOutputVector[scnt + tcnt + powOfTwo / 2] = even - odd;
                W *= Wm;
            }
        }
    }
}

void FourierTransform::inverse() {
    checkSize(inverseInputVector);
    inverseOutputVector.resize(inverseInputVector.size());
    unsigned sizeBit = trunc(log2(inverseOutputVector.size() - 1) / log2(2)) + 1;

    for (unsigned cnt = 0; cnt < inverseOutputVector.size(); cnt++)
        inverseOutputVector[cnt] = inverseInputVector[binaryInverse(cnt, sizeBit)];

    unsigned numberOfTreeLevels = log2(inverseOutputVector.size()) / log2(2);
    for (unsigned fcnt = 1; fcnt <= numberOfTreeLevels; fcnt++) {
        unsigned powOfTwo = pow(2, fcnt);
        complex<double> Wm(cos(2 * PI / powOfTwo), sin(2 * PI / powOfTwo));
        for (unsigned scnt = 0; scnt < inverseOutputVector.size(); scnt += powOfTwo) {
            complex<double> W = 1;
            for (unsigned tcnt = 0; tcnt < powOfTwo / 2; tcnt++) {
                complex<double> odd = W * inverseOutputVector[scnt + tcnt + powOfTwo / 2];
                complex<double> even = inverseOutputVector[scnt + tcnt];
                inverseOutputVector[scnt + tcnt] = even + odd;
                inverseOutputVector[scnt + tcnt + powOfTwo / 2] = even - odd;
                W *= Wm;
            }
        }
    }
}
