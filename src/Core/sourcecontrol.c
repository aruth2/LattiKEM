#include "sourcecontrol.h"

int stepsPerSwitch;//switches from forward bias to reverse bias with this frequency.
int potentialStepDelay;
double offset;
double modulation;
int controlMode;
int waveMode;

void traj_wave(Trajectory *traj, int stepsPerSwitch, int thisExternalCondition, int waveType, double offset, double modulation, int stepDelay)
{ 
	printf("Initializing a wave type %d in slot %d with %d steps per switch a delay of %d an offset of %g and a modulation of %g\n",waveType,thisExternalCondition,stepsPerSwitch,stepDelay, offset,modulation);
    //printf("External Conditions are saved to %d\n",traj->externalConditions);
    //int iWave,iStep;
    int usedStepsPerSwitch;
    int numSteps = getNumSteps();
    int waveSteps = numSteps-stepDelay;
    if(stepsPerSwitch > waveSteps)
	usedStepsPerSwitch = waveSteps;
	else
	usedStepsPerSwitch=stepsPerSwitch;
    

	double value;
	int iStep,iStepWave,iWave;
	for(iStep=0;iStep<fmin(stepDelay,numSteps);iStep++)
	{
		
		*(traj->externalConditions+iStep * traj->numExternalConditions+thisExternalCondition) = offset;
		//printf("On step %d value is %g\n",iStep,*(traj->externalConditions+iStep * traj->numExternalConditions+thisExternalCondition));
	
	}
	for(iStepWave=0;iStepWave<stepsPerSwitch;iStepWave++)
	{
		switch (waveType)
		{
			case SQUARE_WAVE:
			if(iStepWave < stepsPerSwitch/2)
			value = modulation + offset;
			else
			value = offset;
			break;
			case TRIANGLE_WAVE:
			if(iStepWave < stepsPerSwitch/2)
			value = (1 - fabs(iStepWave - stepsPerSwitch/4.0)/(stepsPerSwitch/4.0)) * modulation + offset;
			else
			value = -1 * (1 - fabs(iStepWave - 3.0*stepsPerSwitch/4.0)/(stepsPerSwitch/4.0)) * modulation + offset;
			break;
			case SINE_WAVE:
			value = (sin(2*M_PI*iStepWave/stepsPerSwitch)) * modulation + offset;
			break;
		}
		for(iWave=0,iStep=iWave*stepsPerSwitch+iStepWave+stepDelay;iStep<numSteps;iWave++,iStep=iWave*stepsPerSwitch+iStepWave+stepDelay)
		{
			
			*(traj->externalConditions+(iStep) * traj->numExternalConditions+thisExternalCondition) = value;
			//if(offset != 0 || modulation != 0)
			//if(iStep == traj->numSteps-1)
			//printf("On step %d value is %g\n",iStep,value);
		}
	}

}

void sourceControl_Trajectory(Trajectory *traj)
{
	traj_wave(traj,stepsPerSwitch,0,waveMode,offset,modulation,potentialStepDelay);
}

void sourceControl_RegisterSettings()
{
	registerDouble(&offset,"offset",0);
	registerDouble(&modulation,"modulation",0);
	registerInt(&stepsPerSwitch,"stepsPerSwitch",100);
	registerInt(&potentialStepDelay,"potentialStepDelay",0);
	registerEnum(2,&controlMode,"controlMode",VOLTAGE_CONTROL,"voltage","current");
	registerEnum(3,&waveMode,"waveMode",SQUARE_WAVE,"square","sine","triangle");
	lattikem_registerSettings();
}

int sourceControl_getControlMode()
{
	return controlMode;
}
