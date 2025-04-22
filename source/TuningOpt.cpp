//
// Created by Milan Sastry on 4/21/25.
//
#include "../third_party/LEAF/leaf/leaf.h"
#include "Yin.h"
#include "sfx.h"
#include "Birl.h"
#include "Tune.cpp"
#include <vector>
#include "TuneA.cpp"

int main(){
    LEAF leaf;
    int sampleRate = 44100;
    LEAF_init(&leaf, sampleRate, birl::medium_memory, 500000, []() {return (float)rand() / RAND_MAX; });
    tMempool_init(&birl::smallPool, birl::small_memory, 80328, &leaf);
    tMempool_init(&birl::largePool, birl::large_memory, 33554432, &leaf);
    birl::initGlobalSFXObjects(leaf);
    birl::SFXPhysicalModelPMAlloc(leaf);

    for (int i =0; i < NUM_OF_TONEHOLES; i++)
    {
        printf ("desired %d: %f\n",i,birl::desiredFrequencies[i]);
    }
    float* freqs = getFreqs();
    printDiffs (birl::desiredFrequencies, freqs);
    spsaGradientDescent (birl::desiredFrequencies, 5e-4, 0.01, 150, 1e-13,0.8);
    // multiSampleSPSA  (birl::desiredFrequencies, 5e-4, 0.01,0.9,150, 1e-6,2,0.80);
    // tuneMethodA (birl::desiredFrequencies, 31);
    freqs = getFreqs();
    printDiffs (birl::desiredFrequencies, freqs);
    }