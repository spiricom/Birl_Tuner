#pragma once
#include <iostream>
#include <fstream>
#include <string>

// Export CSV with dynamic number of steps and holes
bool exportFrequencyTuningToCSV(
    float** freqData,
    const std::string& filename,
    const float* frequencies,
    int numSteps,
    int numHoles)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    // Write headers
    file << "Frequency";
    for (int hole = 0; hole < numHoles; ++hole) {
        file << ",Hole_" << hole;
    }
    file << "\n";

    // Write rows
    for (int i = 0; i < numSteps; ++i) {
        file << frequencies[i];
        for (int j = 0; j < numHoles; ++j) {
            file << "," << freqData[i][j];
        }
        file << "\n";
    }

    file.close();
    return true;
}
