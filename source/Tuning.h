#ifndef TUNING_H_
#define TUNING_H_

#include "Birl.h"
#include "sfx.h"
// ONEHOLE. TWO NOTES.
const int NUM_NOTES = 10;

// CONSTS.
const double MIN_D1 = 0.5;
const double MIN_DH = 0.03;
const double DH_FIRST_GUESS = 0.03;

const double BORE_DIAMETER = 1.89;
const double TONEHOLE_DIAMETER = 0.4765;
const double TONEHOLE_HEIGHT = 0.34;
const double MINIMUM_CHIMNEY = 0.34;

const double FICTITIOUS_DUCT = (MINIMUM_CHIMNEY + TONEHOLE_DIAMETER) * pow ((BORE_DIAMETER / TONEHOLE_DIAMETER), 2.0) - 0.45 * BORE_DIAMETER;


// = {2.519840, 2.244920, 2.000000, 1.887750, 1.681790, 1.498310, 1.334830, 1.259920, 1.122460, 1.000000, 0.9439};
//const double tuning[] = {10.0/4.0, 18.0/8.0, 2.0/1.0, 15.0/8.0, 5.0/3.0, 3.0/2.0, 4.0/3.0, 5.0/4.0, 9.0/8.0, 1.0, 15.0/16.0};
const double midiTuning[] = {2.519842,2.244924,2.0,1.887749,1.681793,1.498307,1.33484, 1.259921,1.122462,1.0};
static inline double convertTocm(double samps) {
    return (samps * (C_cm / (sRate_*OVERSAMPLE)));
}

static double convertToSamples(double cm) {
    return (cm * ((sRate_*OVERSAMPLE) / C_cm));
}


static double calcg(int thNum) {
    // CHANGED THIS FROM NUM_NOTES - 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (thNum < 0 || thNum > NUM_NOTES) {
        printf("calcg: thNum out of bounds: %d\n", thNum);
        return 0.0;
    }
    return (midiTuning[thNum] / midiTuning[thNum+1]) - 1;
}

static double calcd1(int LC, double LS) {
    return (LS - (double) LC)/(0.3 * convertToSamples(2.54));
    /* return (LS - (double) LC)/0.3; */
}

// In samples.
static double calcLC(double LS, double LC) {
    double d1 = calcd1(LC, LS);
    double cutAmount = 0.3f * d1;
    return LS - cutAmount;

}

// In samples.
static double calcLS(double Fc) {
    return (sRate_*OVERSAMPLE)/(4.0 * Fc);
}

static double calclL(double d1, int thNum, double LS) {
    double dH = TONEHOLE_DIAMETER;
    double g = calcg(thNum);
    double LSh = (1.0/midiTuning[thNum]) * LS;
    double LBh = TONEHOLE_HEIGHT + dH * ((d1*d1)/(dH*dH)) - 0.45*d1;
    double z = 0.5 * g * sqrt(1 + 4*(LBh/(g*LS))) - 0.5*g;
    return (LSh - (z*LSh));
}

// In samples iff LSh is in samples.
static double calcLBh(int thNum, double LSh, double lLint) {
    double g = calcg(thNum);
    double lL = (double) lLint;

    double gLSh = g * LSh;
    double nmrtr1 = (LSh + 0.5*gLSh - lL) * (LSh + 0.5*gLSh - lL);
    nmrtr1 /= gLSh;
    return nmrtr1 - (gLSh / 4.0);
}

// static double calcdH(int thNum, double d1, double LS, double lL) {
//     double LSh = (1.0/midiTuning[thNum]) * LS;
//     double LBh = calcLBh(thNum, LSh, lL);
//     return (d1*d1) / (LBh + 0.45 * d1);
// }

static double checkTuning(double d1, double dH, double LSh, double lL, double g) {
    double LBh = dH * ((d1*d1)/(dH*dH)) - 0.45*d1;
    double z = 0.5 * g * sqrt(1 + 4*(LBh/(g*LSh))) - 0.5*g;
    return (sRate_*OVERSAMPLE)/(4 * (lL + (z*LSh)));
}

static double calc_acoustic_length(double fundamental_frequency) {
    return C_cm / (4.0 * fundamental_frequency); // returns length in centimeters
}
static double calc_g(int tonehole_number)
{
    return midiTuning[tonehole_number]/midiTuning[tonehole_number-1] - 1.0;
    // return pow (2.0, (1.0/12.0)) - 1;

}

static float calc_cut_length(double fundamental_frequency, double bore_diameter)
{
    return calc_acoustic_length(fundamental_frequency) - 0.3 * bore_diameter;
}

static double acoustic_length_at_tonehole (int tonehole_number, double desired_frequency)
{
    return C_m / (4 * desired_frequency) * 100.0;
}

/* Tests different values of tonehole diameter in order to find the cut length.*/

static double cut_length_at_tonehole (int tonehole_number, double desired_frequency, double g)
{
    /* Now let's cut a little of the tube. We'll call this delta(l_H). */
    double cut_guess = -1.0; // start guessing here

    /* Based on Equation (3), we'll estimate a z value. */
    double z_guess = cut_guess / acoustic_length_at_tonehole (tonehole_number, desired_frequency); //cm

    /* We want to find a cut_length that results in what z should be. */
    double z_should_be = 0.5 * g * sqrt(1 + (4 * FICTITIOUS_DUCT) / (g * acoustic_length_at_tonehole (tonehole_number, desired_frequency))) - 0.5 * g;
    /* Now we need to loop different values of cut_guess until z_guess = z based upon Equation (4). */
    while (z_should_be > z_guess)
    {
        cut_guess += 0.00001;

        z_guess = cut_guess / acoustic_length_at_tonehole (tonehole_number, desired_frequency);
    }
    /* Now with the new value of z, use Equation (3) to find the proper cut length. */
    printf("z_guess = %f, z_should_be = %f\n", z_guess, z_should_be);
    double cut_length_for_frequency = acoustic_length_at_tonehole (tonehole_number, desired_frequency) * (1 - z_guess);
    return cut_length_for_frequency;
}

#endif


