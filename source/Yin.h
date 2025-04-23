//
// Created by Milan Sastry on 3/15/25.
//

#ifndef YIN_H
#define YIN_H

#include <vector>
//int channel = 0; //use the left channel

class Yin {
    public:
        Yin(float rate, int size, double thresh);
        float getPitch(float* buffer);

    private:
        float sampleRate;
        int bufferSize;
        double threshold;
        std::vector<float> yinBuffer;

        void difference(float* buffer);
        void normalizeDifference();
        int absoluteThreshold();
        float parabolicInterpolation(int tauEstimate);
};

#endif //YIN_H
