#ifndef _SOURCECNTRL_H_
#define _SOURCECNTRL_H_

#include "lattikem.h"

enum SourceControlStates {VOLTAGE_CONTROL, CURRENT_CONTROL};
enum Waves {SQUARE_WAVE, SINE_WAVE, TRIANGLE_WAVE};

void traj_wave(Trajectory *traj, int stepsPerSwitch, int thisExternalCondition, int waveType, double offset, double modulation, int stepDelay);
void sourceControl_Trajectory(Trajectory *traj);
void sourceControl_RegisterSettings();
int sourceControl_getControlMode();

#endif
