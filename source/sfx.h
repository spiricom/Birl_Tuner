#ifndef SFX_H
#define SFX_H


//#include "PluginEditor.h"
#include "ui.h"
#include "Filters.h"
#include "Tube.h"
#include "Tuning.h"
#include "../third_party/LEAF/leaf/leaf.h"


#define SMALL_MEM_SIZE 80328
#define MED_MEM_SIZE 519000
#define LARGE_MEM_SIZE 33554432 // size of SDRAM IC




namespace birl {

extern char small_memory[SMALL_MEM_SIZE];
extern char medium_memory[MED_MEM_SIZE];
extern char large_memory[LARGE_MEM_SIZE];


extern tMempool smallPool;
extern tMempool largePool;

/* physical model */
extern Tube tubes[MAX_TONEHOLES + 1];
extern tPoleZero toneHoles[MAX_TONEHOLES];
//extern DCFilter *dcblocker1;
//extern DCFilter *dcblocker2;
extern tHighpass dcblocker1;
extern tHighpass dcblocker2;
extern tBiQuad biquad;
extern tSVF pf1;
extern tSVF pf2;
extern tSVF lp1;
extern tSVF lp2;
extern tSVF noiseBP;


extern float controlKnobValues[ControlNil][NUM_SYNTH_KNOB_VALUES];
extern uint8_t knobActive[NUM_SYNTH_KNOB_VALUES];
extern float fingers[NUM_OF_TONEHOLES];
extern float maxToneholeArg[NUM_OF_TONEHOLES];
extern bool buttons[NUM_OF_BUTTONS];

extern float breathArray[2];
extern float desiredFrequencies[NUM_OF_TONEHOLES+1];


//extern PlayMode samplerMode;
extern float sampleLength;

extern uint32_t freeze;

void initGlobalSFXObjects(LEAF &leaf);

/* physical model physical model */


void SFXPhysicalModelSetTubeLength(int index, double newLength);
double SFXPhysicalModelGetTubeLength(int index);
void SFXPhysicalModelPMAlloc(LEAF &leaf);
float SFXPhysicalModelPMTick();
float myRandom();
void SFXPhysicalModelPMFree(void);
void SFXPhysicalModelSetToneholeRadius(int index, float radius);
void SFXPhysicalModelSetTonehole(int index, float newValue);
void SFXPhysicalModelSetBreathPressure(float input);
void SFXPhysicalModelCalcTHCoeffs();
void SFXPhysicalModelTune(float fundamental);
void SFXPhysicalModelRetune(float fundamental);
float SFXPhysicalModelInterpolateLinear(float a, float b, float alpha);
float SFXPhysicalModelGetToneholeRadius(int index);


}

#endif
