#include "sourcecontrol.h"

int stepsPerSwitch;//switches from forward bias to reverse bias with this frequency.
int stepDelay;
double offset;
double modulation;
double dutyCycle;
double phaseShift;
int controlMode;
int waveMode;
int decaySteps;

void traj_wave(Trajectory *traj, int stepsPerSwitch, int thisExternalCondition, int waveType, double offset, 
double modulation, int stepDelay, double dutyCycle,double phaseShift, int stepDecay)
{  
	int stepShift = phaseShift/(2*M_PI)*stepsPerSwitch;
	printf("Initializing a wave type %d in slot %d with %d steps per switch a delay of %d a step shift of %d \n"
	"an offset of %g a modulation of %g, and %d decay steps\n",waveType,thisExternalCondition,stepsPerSwitch,stepDelay, stepShift, offset ,modulation, decaySteps);
    //printf("External Conditions are saved to %d\n",traj->externalConditions);
    //int iWave,iStep;
    //int usedStepsPerSwitch;
    int numSteps = getNumSteps();
    //int waveSteps = numSteps-stepDelay;
    //if(stepsPerSwitch > waveSteps)
	//usedStepsPerSwitch = waveSteps;
	//else
	//usedStepsPerSwitch=stepsPerSwitch;
    

	double value;
	int iStep,iStepWave,iWave;
	for(iStep=0;iStep<fmin(stepDelay,numSteps);iStep++)
	{
		
		*(traj->externalConditions+iStep * traj->numExternalConditions+thisExternalCondition) = offset;
		//printf("On step %d value is %g\n",iStep,*(traj->externalConditions+iStep * traj->numExternalConditions+thisExternalCondition));
	
	}
	int cycleStep;
	for(iStepWave=0;iStepWave<stepsPerSwitch;iStepWave++)
	{
		cycleStep = iStepWave+stepShift;
		switch (waveType)
		{
			case SQUARE_WAVE:
				if(cycleStep < stepsPerSwitch * dutyCycle)
			//if(iStepWave < stepsPerSwitch/2)
					value = modulation + offset;
				else
					value = offset;
			break;
			case TRIANGLE_WAVE:
				if(cycleStep < stepsPerSwitch/2)
					value = (1 - fabs(cycleStep - stepsPerSwitch/4.0)/(stepsPerSwitch/4.0)) * modulation + offset;
				else
					value = -1 * (1 - fabs(cycleStep - 3.0*stepsPerSwitch/4.0)/(stepsPerSwitch/4.0)) * modulation + offset;
			break;
			case SINE_WAVE:
				value = (sin(2*M_PI*cycleStep/stepsPerSwitch)) * modulation + offset;
			break;
			case SQUARE_DECAY:
				if(cycleStep < stepsPerSwitch * dutyCycle)
					value = modulation + offset;
				else
					if(decaySteps == 0)
						value = offset;
					else
						value = modulation* exp( - (cycleStep - stepsPerSwitch * dutyCycle)/decaySteps) + offset;
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
	traj_wave(traj,stepsPerSwitch,0,waveMode,offset,modulation,stepDelay,dutyCycle,phaseShift, decaySteps);
}

void sourceControl_RegisterSettings()
{
	registerDouble(&offset,"offset",0);
	registerDouble(&modulation,"modulation",0);
	registerDouble(&dutyCycle,"dutyCycle",0.5);
	registerDouble(&phaseShift,"phaseShift",0.0);
	registerInt(&stepsPerSwitch,"stepsPerSwitch",100);
	registerInt(&stepDelay,"stepDelay",0);
	registerInt(&decaySteps,"decaySteps",0);
	registerEnum(2,&controlMode,"controlMode",VOLTAGE_CONTROL,"voltage","current");
	registerEnum(4,&waveMode,"waveMode",SQUARE_WAVE,"square","sine","triangle","squareDecay");
	//lattikem_registerSettings();
}

int sourceControl_getControlMode()
{
	return controlMode;
}
