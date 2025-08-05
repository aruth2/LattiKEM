#ifndef _SETTINGS_H_
#define _SETTINGS_H_

#include "supp.h"
#include <stdarg.h>

#define maxSettings 1000
#define stringSize 100
#define maxEnumSize 10
#define maxListSize 1000

void allocateSettings();
void registerInt(int *pointer, char *descriptor, int defaultValue);
void registerDouble(double *pointer, char *descriptor, double defaultValue);
void registerString(char *pointer, char *descriptor, char *defaultValue);
void registerEnum( int num, int *pointer, char *descriptor, int defaultValue,  ... );
int * registerIntList( char *descriptor, int *numEntriesPointer);
void loadSettings(FILE *settingsFile);
void saveSettings(FILE *settingsFile);

#endif
