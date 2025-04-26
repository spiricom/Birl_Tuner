#include "sfx.h"
#include "leaf.h"
#include "Tuning.h"
namespace birl
{

char small_memory[SMALL_MEM_SIZE];
char medium_memory[MED_MEM_SIZE];
char large_memory[LARGE_MEM_SIZE];


tMempool smallPool;
tMempool largePool;

float defaultControlKnobValues[ControlNil][NUM_SYNTH_KNOB_VALUES];
float controlKnobValues[ControlNil][NUM_SYNTH_KNOB_VALUES];
//uint8_t knobActive[NUM_ADC_CHANNELS];
float prevDisplayValues[NUM_SYNTH_KNOB_VALUES];


#define EXP_BUFFER_SIZE 128
float expBuffer[EXP_BUFFER_SIZE];
float expBufferSizeMinusOne = EXP_BUFFER_SIZE - 1;

#define DECAY_EXP_BUFFER_SIZE 512
float decayExpBuffer[DECAY_EXP_BUFFER_SIZE];
float decayExpBufferSizeMinusOne = DECAY_EXP_BUFFER_SIZE - 1;

float fingers[NUM_OF_TONEHOLES];
float maxToneholeArg[NUM_OF_TONEHOLES];
bool buttons[NUM_OF_BUTTONS];

uint8_t dummyMemory[128];


float randomCall()
{
    return 0.0f;
}
float myRandom() {
    return (float) rand() /RAND_MAX;
}

void initGlobalSFXObjects(LEAF &leaf)
{
    LEAF_init(&leaf, 44100.0f, medium_memory, MED_MEM_SIZE, myRandom);
    tMempool_init (&smallPool, small_memory, SMALL_MEM_SIZE, &leaf);
    tMempool_init (&largePool, large_memory, LARGE_MEM_SIZE, &leaf);

    LEAF_generate_exp(expBuffer, 1000.0f, -1.0f, 0.0f, -0.0008f, EXP_BUFFER_SIZE);
    // exponential buffer rising from 0 to 1
    LEAF_generate_exp(decayExpBuffer, 0.001f, 0.0f, 1.0f, -0.0008f, DECAY_EXP_BUFFER_SIZE);
    // exponential decay buffer falling from 1 to 0
    for (int i = 0; i < NUM_OF_TONEHOLES; i++)
    {
        maxToneholeArg[i] = 1.0f;
    }
    defaultControlKnobValues[PhysicalModelPM][0] = 0.2f;        // gain
    defaultControlKnobValues[PhysicalModelPM][1] = 100.0f;      // fundamental
    defaultControlKnobValues[PhysicalModelPM][2] = 9.0f;        // num_fingers
    defaultControlKnobValues[PhysicalModelPM][3] = 0.995f;      // dcblocker1
    defaultControlKnobValues[PhysicalModelPM][4] = 0.995f;      // dcblocker2
    defaultControlKnobValues[PhysicalModelPM][5] = 0.169301f;   // biquad_coeff1
    defaultControlKnobValues[PhysicalModelPM][6] = 0.338601f;   // biquad_coeff2
    defaultControlKnobValues[PhysicalModelPM][7] = 0.169301f;   // biquad_coeff3
    defaultControlKnobValues[PhysicalModelPM][8] = -0.482013f;  // biquad_coeff4
    defaultControlKnobValues[PhysicalModelPM][9] = 0.186622f;   // biquad_coeff5
    defaultControlKnobValues[PhysicalModelPM][10] = 2000.0f;    // pf1_cutoff
    defaultControlKnobValues[PhysicalModelPM][11] = 0.5f;       // pf1_q
    defaultControlKnobValues[PhysicalModelPM][12] = 1000.0f;    // pf2_cutoff
    defaultControlKnobValues[PhysicalModelPM][13] = 1.0f;       // pf2_q
    defaultControlKnobValues[PhysicalModelPM][14] = 5000.0f;    // lp1_cutoff
    defaultControlKnobValues[PhysicalModelPM][15] = 0.5f;       // lp1_q
    defaultControlKnobValues[PhysicalModelPM][16] = 5000.0f;    // lp2_cutoff
    defaultControlKnobValues[PhysicalModelPM][17] = 0.5f;       // lp2_q
    defaultControlKnobValues[PhysicalModelPM][18] = 0.0f;       // shaper drive
    defaultControlKnobValues[PhysicalModelPM][19] = 0.0f;       // shaper mix
    defaultControlKnobValues[PhysicalModelPM][20] = 0.2f;       // noise gain
    defaultControlKnobValues[PhysicalModelPM][21] = 16000.0f;   // noise_bp_cutoff
    defaultControlKnobValues[PhysicalModelPM][22] = 1.0f;       // noise_bp_q

    for (int c = 0; c < ControlNil; c++)
    {
        for (int v = 0; v < sizeof(controlKnobValues[c]); v++) {
            controlKnobValues[c][v] = defaultControlKnobValues[c][v];
        }
    }
}

/* physical model internal poly */
Tube tubes[MAX_TONEHOLES + 1];
tPoleZero toneHoles[MAX_TONEHOLES];
//DCFilter *dcblocker1;
//DCFilter *dcblocker2;
tHighpass dcblocker1;
tHighpass dcblocker2;
tBiQuad biquad;
tSVF pf1;
tSVF pf2;
tSVF lp1;
tSVF lp2;
tSVF noiseBP;

double pam_[MAX_TONEHOLES];
double pbp_[MAX_TONEHOLES];
double pthp_[MAX_TONEHOLES];
double scatter_[MAX_TONEHOLES];
double thCoeff_[MAX_TONEHOLES];

double rb_;
double outputGain;
double noiseGain;
double breathPressure;
double prevBreathPressure;
float breathArray[2];

float fundamental;
int tubeIndex;
int numToneholes;

float count;
float min;
float max;
float mDrive;
float shaperMix;
    float desiredFrequencies[NUM_OF_TONEHOLES+1];

void SFXPhysicalModelPMAlloc(LEAF &leaf)
{
    leaf.clearOnAllocation = 1;
    for (int i = 0; i < MAX_TONEHOLES; i++) {
        tPoleZero_initToPool(&toneHoles[i], &smallPool);
    }

    outputGain = defaultControlKnobValues[PhysicalModelPM][0];
    noiseGain = defaultControlKnobValues[PhysicalModelPM][20];

    fundamental = defaultControlKnobValues[PhysicalModelPM][1];
    tubeIndex = 0;
    numToneholes = defaultControlKnobValues[PhysicalModelPM][2];

    breathPressure = 0.0;
    prevBreathPressure = 0.0;
    count = 0;
    min = 0.0;
    max = 0.0;
    mDrive = 0.0;
    shaperMix = 0.0;
    sRate_ = leaf.sampleRate;

    for (int i = 0; i < MAX_TONEHOLES + 1; i++) {
//        tubes[i] = initTube(3); // IDK???
        tLinearDelay_init(&tubes[i].upper, 100, 512, &leaf);
        tLinearDelay_init(&tubes[i].lower, 100, 512, &leaf);
    }
//    dcblocker1 = initDCFilter(defaultControlKnobValues[PhysicalModelPM][3]);
//    dcblocker2 = initDCFilter(defaultControlKnobValues[PhysicalModelPM][4]);
    tHighpass_initToPool(&dcblocker1, 13.0, &smallPool);
    tHighpass_initToPool(&dcblocker2, 13.0, &smallPool);
    // choice of 13 Hz for cutoff due to translation from DC Blocker to LEAF.

    tBiQuad_initToPool(&biquad, &smallPool);
    tBiQuad_setCoefficients(biquad, defaultControlKnobValues[PhysicalModelPM][5], defaultControlKnobValues[PhysicalModelPM][6], defaultControlKnobValues[PhysicalModelPM][7], defaultControlKnobValues[PhysicalModelPM][8], defaultControlKnobValues[PhysicalModelPM][9]);

    tSVF_initToPool(&pf1,     SVFTypePeak,      defaultControlKnobValues[PhysicalModelPM][10], defaultControlKnobValues[PhysicalModelPM][11], &smallPool);
    tSVF_initToPool(&pf2,     SVFTypePeak,      defaultControlKnobValues[PhysicalModelPM][12], defaultControlKnobValues[PhysicalModelPM][13], &smallPool);
    tSVF_initToPool(&lp1,     SVFTypeLowpass,   defaultControlKnobValues[PhysicalModelPM][14], defaultControlKnobValues[PhysicalModelPM][15], &smallPool);
    tSVF_initToPool(&lp2,     SVFTypeLowpass,   defaultControlKnobValues[PhysicalModelPM][16], defaultControlKnobValues[PhysicalModelPM][17], &smallPool);
    tSVF_initToPool(&noiseBP, SVFTypeBandpass,  defaultControlKnobValues[PhysicalModelPM][21], defaultControlKnobValues[PhysicalModelPM][22], &smallPool);

    //SFXPhysicalModelTune(fundamental);
}

void SFXPhysicalModelSetToneholeRadius(int index, float radius) {
    if (radius < MIN_TONEHOLE_RADIUS || radius > MAX_TONEHOLE_RADIUS) {
        // juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon,
        //     "Radius out of range",
        //     "Error: Radius is out of range",
        //     "OK");
        return;
    }
    rth_[index] = radius;
    scatter_[index] = -pow(radius,2) / ( pow(radius,2) + 2*pow(rb_,2) );

    // Calculate toneHole coefficients.
    double te = radius;    // effective length of the open hole
    thCoeff_[index] = (te*2*(sRate_*OVERSAMPLE) - C_m) / (te*2*(sRate_*OVERSAMPLE) + C_m);
}

float SFXPhysicalModelGetToneholeRadius(int  index) {
    return rth_[index];
}
void SFXPhysicalModelSetTonehole(int index, float newValue) {
    double new_coeff;
    newValue = 1.0 - newValue;
    if (newValue <= 0.0)
        new_coeff = 0.999999995;
    else if (newValue >= 1.0)
        new_coeff = thCoeff_[index];
    else
        new_coeff = newValue * (thCoeff_[index] - 0.999999995) + 0.999999995;
    tPoleZero_setA1(toneHoles[index], -new_coeff);
    tPoleZero_setB0(toneHoles[index], new_coeff);
}
void SFXPhysicalModelSetBreathPressure(float input) {
    breathPressure = input;
    //printf("%9.9f \n", breathPressure);
}


void SFXPhysicalModelSetTubeLength(int index, double newLength)
{
    tubeLengths_[index] = newLength;
    tLinearDelay_setDelay(tubes[index].upper, tubeLengths_[index]);
    tLinearDelay_setDelay(tubes[index].lower, tubeLengths_[index]);

}
double SFXPhysicalModelGetTubeLength(int index)
{
    return tubeLengths_[index];
}
void SFXPhysicalModelCalcTHCoeffs() {
        // Calculate initial tone hole three-port scattering coefficients
        for (int i = 0; i < MAX_TONEHOLES; i++) {
            scatter_[i] = -pow(rth_[i],2) / ( pow(rth_[i],2) + 2*pow(rb_,2) );

            // Calculate toneHole coefficients and set for initially open.
            thCoeff_[i] = (rth_[i]*2*(sRate_*OVERSAMPLE) - C_m) / (rth_[i]*2*(sRate_*OVERSAMPLE) + C_m);


            // Initialize fingers.
            tPoleZero_setA1(toneHoles[i], -thCoeff_[i]);
            tPoleZero_setB0(toneHoles[i], thCoeff_[i]);
            tPoleZero_setB1(toneHoles[i], -1.0);
        }
}
void SFXPhysicalModelTune(float fundamental) {
    double lowpassDelay = 0;
    double effectiveLength = calcLS(fundamental) - lowpassDelay;
    float filterDelay = 0;
    float scale = 1.0f;



        //double tempy = calclL(BORE_DIAMETER, i, effectiveLength);
        double effective_length_tube = calc_acoustic_length (fundamental);
        printf("effective length tube: %f\n", effective_length_tube);
        double cut_length_tube = calc_cut_length (fundamental, BORE_DIAMETER);
        printf("cut length tube: %f\n\n", cut_length_tube);
        double length_of_tube_so_far = 0.0;

        for (int i = 0; i < MAX_TONEHOLES + 1 ; i++)
        {
            double desired_frequency = fundamental * midiTuning[i];
            printf("i = %i, tuning: %f, desired frequency: %f\n", i, midiTuning[i], desired_frequency);
            double g = calc_g(i);
            tubeLengths_[i] = cut_length_at_tonehole (i, desired_frequency, g) - length_of_tube_so_far; //cm
            printf("calc_g at %i = %f\n", i, g);
            length_of_tube_so_far += tubeLengths_[i];
            printf("length of tube so far at index %i: %f, after adding a segment of %f\n", i, length_of_tube_so_far, tubeLengths_[i]);
            tLinearDelay_setDelay(tubes[i].upper, tubeLengths_[i]);
            tLinearDelay_setDelay(tubes[i].lower, convertToSamples(tubeLengths_[i] * 0.01));
            printf("convertToSamples(tube_segments) at %i = %f\n\n", i, convertToSamples(tubeLengths_[i] * 0.01));


            // tLinearDelay_setDelay(tubes[i].upper, tubeLengths_[i]);
            // tLinearDelay_setDelay(tubes[i].lower, tubeLengths_[i]);

            //prevlL += tubeLengths_[i] + filterDelay;
            //previousCut = correction;
            printf("th %d: lL = %f\n", i, tubeLengths_[i]);
            //printf("prev1L %f\n", prevlL);
        }



    double lL = tubeLengths_[0];

    for (int i = 0; i < NUM_OF_TONEHOLES; i++) {
        originalRth_[i] = TONEHOLE_DIAMETER / 200.0f;//convertTocm(calcdH(i, BORE_DIAMETER, effectiveLength, lL))/200.0; //dividing by 200 because we converted to cm, and rth should be in meters, but also wants radius not diameter - so divide by 100 then divide by 2
        rth_[i] = originalRth_[i];
        lL += tubeLengths_[i+1];
    }

    rb_ = BORE_DIAMETER / 200.0f;
    printf("rb: %f\n", rb_);

    lL = 0.0;
    for (int i = 0; i < NUM_OF_TONEHOLES; i++) {
        lL += tubeLengths_[i];
        double LSh = (1.0/midiTuning[i]) * effectiveLength;
        float freq = checkTuning(BORE_DIAMETER, convertToSamples(rth_[i]*200.0), LSh, lL, calcg(i));
        desiredFrequencies[i] = LEAF_frequencyToMidi (midiTuning[i]*fundamental);
        printf("th %d rth: %f m, output freq when open: %f\n", i, rth_[i], LEAF_midiToFrequency(desiredFrequencies[i]));
    }
    desiredFrequencies[9] = LEAF_frequencyToMidi (fundamental);
    // Calculate the tonehole coefficients.
    SFXPhysicalModelCalcTHCoeffs();

}


float SFXPhysicalModelInterpolateLinear(float a, float b, float alpha) {
    return (alpha * a) + ((1.0-alpha) * b);
}





float SFXPhysicalModelPMTick() {
    double sample = 0.0f;
    double bellReflected;
    mDrive = birl::controlKnobValues[0][18];
    shaperMix = birl::controlKnobValues[0][19];

    double breath = breathPressure;
    double noise = (double) rand() / (double) RAND_MAX;

    int numHoles = 9;

    noise = noiseGain * tSVF_tick(noiseBP, noise);
    breath += breath * noise;

    double pressureDiff = tLinearDelay_tickOut(tubes[0].lower) - breath;
    double reedLookup = pressureDiff * reedTable (pressureDiff);

    breath = LEAF_clip(-1.0, breath + reedLookup, 1.0);
    breath = SFXPhysicalModelInterpolateLinear(shaper(breath, mDrive), breath, shaperMix);
    //breath = tSVF_tick(pf1, breath);
    //breath = tSVF_tick(lp1, breath);
    breath = tHighpass_tick(dcblocker1, breath);

//    breath = inputSVFPeak(pf_, breath);
//    breath = inputSVFLP(lp_, breath);
//    breath = inputDCFilter(dcBlocker_, breath);

    // tone-hole scatter
    for (int i = 0; i < numHoles; i++)
    {
        double pap = tLinearDelay_tickOut(tubes[i].upper);
        double pbm = tLinearDelay_tickOut(tubes[i+1].lower);
        double pthm = toneHoles[i]->lastOut;

        double scatter = scatter_[i] * (pap + pbm - (2.0*pthm));
        pbp_[i] = pap + scatter;
        pam_[i] = pbm + scatter;
        pthp_[i] = pap + scatter + pbm - pthm;

        if (i == 0)
        {
            sample = pap + pam_[i];
        }
    }

    /* bell filters */
    double bell = tLinearDelay_tickOut(tubes[numHoles].upper);

    // Reflection = Inversion + gain reduction + lowpass filtering.
    //bell = tSVF_tick(pf2, bell);
   //bell = tSVF_tick(lp2, bell);
    bell = tHighpass_tick(dcblocker2, bell);
    //    bell = inputDCFilter(dcBlocker2, bell);

    bellReflected = bell * -0.9999995;

    // tone-hole update
    for (int i = 0; i < numHoles; i++)
    {
        tPoleZero_tick (toneHoles[i], pthp_[i]);
        tLinearDelay_tickIn(tubes[i].lower, pam_[i]);
        tLinearDelay_tickIn(tubes[i+1].upper, pbp_[i]);
        float fingerScaled = LEAF_clip(0.0f, fingers[i] * 2.0f, 1.0f);
        SFXPhysicalModelSetTonehole(i, fingerScaled);
    }

    tLinearDelay_tickIn(tubes[0].upper, breath);
    tLinearDelay_tickIn(tubes[numHoles].lower, bellReflected);

    sample = tanhf(sample);
    return sample;
}

void SFXPhysicalModelPMFree(void) {
    for (int i = 0; i < MAX_TONEHOLES; i++) {
        tPoleZero_free(&toneHoles[i]);
    }
    for (int i = 0; i < MAX_TONEHOLES - 1; i++)
        freeTube(&tubes[i]);
    tHighpass_free(&dcblocker1);
    tHighpass_free(&dcblocker2);
    tBiQuad_free(&biquad);
    tSVF_free(&pf1);
    tSVF_free(&pf2);
    tSVF_free(&lp1);
    tSVF_free(&lp2);
    tSVF_free(&noiseBP);
}


}
