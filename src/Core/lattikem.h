#ifndef _LATTIKEM_H_
#define _LATTIKEM_H_

#include "pkmc.h"
#include "perovskite.h"
//#include "lightcube.h"
#include <unistd.h>


//How many steps in between saves of moves and energy

enum JOB_TYPE {JOB_AUTO, JOB_RESTART, JOB_CONTINUE, JOB_POSTPROCESS, JOB_MAIN};
enum AVERAGING_METHOD {STEP_AVERAGE, TIME_AVERAGE, NOT_COMBINED};
enum EQUIVOCATION_METHOD {STEP_SUM, TIME_ELAPSED};
 
void saveTrajectoryCrystals(Trajectory *traj);
//void saveTrajectoryLightCube(Trajectory *traj);
void performTrajectory(Trajectory *traj);
void simulateTrajectories();
void stepAverage(Trajectory *trajectories);
void timeAverage(Trajectory *trajectories);
void postProcess(Trajectory *trajectories);
void lattikem_registerSettings();
void lattikem_startupMessage();
char *getDir();
#endif
