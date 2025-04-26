#include "../third_party/LEAF/leaf/leaf.h"
#include "Yin.h"
#include "sfx.h"
#include "Birl.h"
#include "Tune.cpp"
#include "TuneA.cpp"
#include "Plot.cpp"

static constexpr int MIN_FUNDAMENTAL = 100;
static constexpr int NUM_STEPS = 400;
static constexpr int NUM_HOLES = 10;

// Allocate a 2D array [numSteps][numHoles]
float** allocateDiffs(int numSteps, int numHoles) {
    float** diffs = new float*[numSteps];
    for (int i = 0; i < numSteps; ++i) {
        diffs[i] = new float[numHoles];
    }
    return diffs;
}

// Unified evaluation for initial loss or GD-adjusted loss
float** evalLosses(
    int startFund,
    int stepSize,
    int numSteps,
    bool useGradientDescent)
{
    float** diffs = allocateDiffs(numSteps, NUM_HOLES);

    for (int idx = 0; idx < numSteps; ++idx) {
        float fundamental = startFund + idx * stepSize;
        birl::SFXPhysicalModelTune(fundamental);

        if (useGradientDescent) {
            // parameters: desired freqs, gain, perturbation, iterations, tolerance, momentum
            spsaGradientDescent(birl::desiredFrequencies, 5e-4, 0.01, 200, 1e-13, 0.8);
        }

        float* freqs = getFreqs();
        for (int j = 0; j < NUM_HOLES; ++j) {
            float desired = birl::desiredFrequencies[j];
            diffs[idx][j] = freqs[j] - desired;
        }
    }
    return diffs;
}

// Plotting helper that builds axis, calls eval, and exports CSV
void plotLosses(
    int startFund,
    int stepSize,
    int numSteps,
    bool useGD,
    const std::string& filename)
{
    // Evaluate losses
    float** diffs = evalLosses(startFund, stepSize, numSteps, useGD);

    // Build frequency axis
    float* frequencies = new float[numSteps];
    for (int i = 0; i < numSteps; ++i) {
        frequencies[i] = startFund + i * stepSize;
    }

    // Export CSV
    if (exportFrequencyTuningToCSV(diffs, filename, frequencies, numSteps, NUM_HOLES)) {
        std::cout << "Data successfully exported to " << filename << "\n";
    }

    // Cleanup
    for (int i = 0; i < numSteps; ++i) {
        delete[] diffs[i];
    }
    delete[] diffs;
    delete[] frequencies;
}

int main() {
    LEAF leaf;
    int sampleRate = 44100;
    LEAF_init(&leaf, sampleRate, birl::medium_memory, 500000, birl::myRandom);
    tMempool_init(&birl::smallPool, birl::small_memory, 80328, &leaf);
    tMempool_init(&birl::largePool, birl::large_memory, 33554432, &leaf);
    birl::initGlobalSFXObjects(leaf);
    birl::SFXPhysicalModelPMAlloc(leaf);

    // Plot initial losses (step size = 1)
    // plotLosses(MIN_FUNDAMENTAL, 1, NUM_STEPS, false, "initial_loss_data_centz.csv");

    // Plot GD-adjusted losses (step size = 10)
    // plotLosses(MIN_FUNDAMENTAL, 4, 100, true, "gd_loss_data_cents.csv");

    birl::SFXPhysicalModelTune(394);
    float* freqs = getFreqs();
    for (int i = 0; i < NUM_HOLES; ++i) {
        printf("hole %d: %f\n", i, LEAF_midiToFrequency(freqs[i]));
    }
    printDiffs(birl::desiredFrequencies, freqs);
    // multiSampleSPSA(birl::desiredFrequencies, 1e-3, 0.01, 200, 1e-13, 5, 0.7);
    // spsaGradientDescent(birl::desiredFrequencies, 5e-3, 0.005, 300, 1e-13, 0.8);
    // Print desired vs actual diffs
    for (int i = 0; i < NUM_HOLES; ++i) {
        printf("desired %d: %f\n", i, birl::desiredFrequencies[i]);
    }
    freqs = getFreqs();
    printDiffs(birl::desiredFrequencies, freqs);

    return 0;
}
