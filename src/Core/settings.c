#include "settings.h"

char *commentFlags = "!@#$%^&*`~';:";

int **intSettings;
char **intDescriptors;
int *intDefaults;
int numIntSettings = 0;

double **doubleSettings;
char **doubleDescriptors;
double *doubleDefaults;
int numDoubleSettings = 0; 

char **stringSettings;
char **stringDescriptors;
char **stringDefaults;
int numStringSettings = 0;
 
int **enumSettings;
char **enumDescriptors;
int *enumSizes;
char ***enumSelectors;
int *enumDefaults;
int numEnumSettings = 0;

int  **intListSettings;
int **numIntListValues;
char **intListDescriptors;
int numIntListSettings = 0;


int saveDefaults;
void allocateSettings()
{
	int i,j;
	
	intSettings = (int **)calloc(maxSettings,sizeof(int *));
	intDescriptors = (char **)calloc(maxSettings,sizeof(char *));
	intDefaults = (int *)calloc(maxSettings,sizeof(int));
	
	doubleSettings = (double **)calloc(maxSettings,sizeof(double *));
	doubleDescriptors = (char **)calloc(maxSettings,sizeof(char *));
	doubleDefaults = (double *)calloc(maxSettings,sizeof(double));
	
	stringSettings = (char **)calloc(maxSettings,sizeof(char *));
	stringDescriptors = (char **)calloc(maxSettings,sizeof(char *));
	stringDefaults = (char **)calloc(maxSettings,sizeof(char *));
	
	enumSettings = (int **)calloc(maxSettings,sizeof(int *));
	enumSizes = (int *)calloc(maxSettings,sizeof(int));
	enumDefaults = (int *)calloc(maxSettings,sizeof(int));
	enumDescriptors = (char **)calloc(maxSettings,sizeof(char *));
	enumSelectors = (char ***)calloc(maxSettings,sizeof(char **));
	
	intListSettings = (int **)calloc(maxSettings,sizeof(int *));
	numIntListValues = (int **)calloc(maxSettings,sizeof(int *));
	intListDescriptors = (char **)calloc(maxSettings,sizeof(char *));
	
	for(i=0;i<maxSettings;i++)
	{
		*(intDescriptors+i) = (char *)calloc(stringSize,sizeof(char));
		*(doubleDescriptors+i) = (char *)calloc(stringSize,sizeof(char));
		*(stringDescriptors+i) = (char *)calloc(stringSize,sizeof(char));
		*(enumDescriptors+i) = (char *)calloc(stringSize,sizeof(char));
		*(intListDescriptors+i) = (char *)calloc(stringSize,sizeof(char));
		
		*(stringSettings+i) = (char *)calloc(stringSize,sizeof(char));

		*(stringDefaults+i) = (char *)calloc(stringSize,sizeof(char));
		*(intListSettings+i) = (int *)calloc(maxListSize,sizeof(int));
		*(enumSelectors+i) = (char **)calloc(maxEnumSize,sizeof(char *));
		for(j=0;j<maxEnumSize;j++)
		{
			*(*(enumSelectors+i)+j) = (char *)calloc(stringSize,sizeof(char));
		}
		
	}
	
	
	registerInt(&saveDefaults,"saveDefaults",0);
}

void registerInt(int *pointer, char *descriptor, int defaultValue)
{
	*(intSettings+numIntSettings) = pointer;
	strcpy(*(intDescriptors+numIntSettings),descriptor);
	*(intDefaults+numIntSettings) = defaultValue;
	*pointer = defaultValue;
	numIntSettings++;
}

void registerDouble(double *pointer, char *descriptor, double defaultValue)
{
	*(doubleSettings+numDoubleSettings) = pointer;
	strcpy(*(doubleDescriptors+numDoubleSettings),descriptor);
	*(doubleDefaults+numDoubleSettings) = defaultValue;
	*pointer = defaultValue;
	numDoubleSettings++;
}

void registerString(char *pointer, char *descriptor, char *defaultValue)
{
	*(stringSettings+numStringSettings) = pointer;
	strcpy(*(stringDescriptors+numStringSettings),descriptor);
	strcpy(*(stringDefaults+numStringSettings), defaultValue);
	strcpy(pointer, defaultValue);
	numStringSettings++;
}

int * registerIntList( char *descriptor, int *numEntriesPointer)
{
	int *pointer = (int *)calloc(maxListSize,sizeof(int));
	*(intListSettings + numIntListSettings) = pointer;
	*(numIntListValues + numIntListSettings) = numEntriesPointer;
	strcpy(*(intListDescriptors+numIntListSettings),descriptor);
	*numEntriesPointer = 0;
	numIntListSettings++;
	
	return pointer;
}

void registerEnum( int num, int *pointer, char *descriptor, int defaultValue,  ... )
{
    va_list arguments;                     	
    va_start ( arguments, defaultValue );
    
    *(enumSettings+numEnumSettings) = pointer;
	strcpy(*(enumDescriptors+numEnumSettings),descriptor);
	*(enumDefaults+numEnumSettings) = defaultValue;
	*pointer = defaultValue;
	*(enumSizes+numEnumSettings) = num;
	
    int i;
    for ( i = 0; i < num; i++ )        
    {
        strcpy(*(*(enumSelectors+numEnumSettings)+i), va_arg ( arguments, char * )); 
    }
    numEnumSettings++;
    
    va_end ( arguments ); 
}

void loadSettings(FILE *settingsFile)
{
	printf("\nLoading Settings \n\n");
	
	int i,j;
	for(i=0;i<numIntSettings;i++)
		readInt(settingsFile,*(intDescriptors+i),*(intSettings+i),commentFlags);
	for(i=0;i<numDoubleSettings;i++)
		readDouble(settingsFile,*(doubleDescriptors+i),*(doubleSettings+i),commentFlags);
		
		
	for(i=0;i<numStringSettings;i++)
		readString(settingsFile,*(stringDescriptors+i),*(stringSettings+i),commentFlags);
	
	for(i=0;i<numIntListSettings;i++)
		readIntList(settingsFile,*(intListDescriptors+i),*(intListSettings+i),*(numIntListValues+i),commentFlags);
		
	char stringBuffer[1000];
	int enumFound;
	for(i=0;i<numEnumSettings;i++)
	{
		strcpy(stringBuffer,*(*(enumSelectors+i)+**(enumSettings+i)));
		readString(settingsFile,*(enumDescriptors+i),stringBuffer,commentFlags);
		enumFound = 0;
		for(j = 0; j<*(enumSizes+i); j++)
		{
			if(!strcmp(stringBuffer,*(*(enumSelectors+i)+j)))
			{
				//printf("Enum %s is %s\n",*(enumDescriptors+i),*(*(enumSelectors+i)+j));
				**(enumSettings+i) = j;
				enumFound = 1;
				break;
			}
		}
		if(enumFound == 0)
			printf("%s does not match any allowed value for %s. Using default value of %s\n",stringBuffer,*(enumDescriptors+i),*(*(enumSelectors+i)+**(enumSettings+i)));
	}
}

void saveSettings(FILE *settingsFile)
{
	//printf("There are %d ints, %d doubles, %d strings, and %d enums to be saved\n",numIntSettings,numDoubleSettings,numStringSettings,numEnumSettings);
	int i;
	for(i=0;i<numIntSettings;i++)
		if(**(intSettings + i) != *(intDefaults+i) || saveDefaults)
			saveInt(settingsFile,*(intDescriptors+i),**(intSettings+i));

	for(i=0;i<numDoubleSettings;i++)
		if(**(doubleSettings + i) != *(doubleDefaults+i) || saveDefaults)
			saveDouble(settingsFile,*(doubleDescriptors+i),**(doubleSettings+i));
					
	for(i=0;i<numStringSettings;i++)
		if(strcmp(*(stringSettings+i),*(stringDefaults+i)) || saveDefaults)
			saveString(settingsFile,*(stringDescriptors+i),*(stringSettings+i));
	
	for(i=0;i<numEnumSettings;i++)
		if(**(enumSettings+i) != *(enumDefaults+i) || saveDefaults)
			saveString(settingsFile,*(enumDescriptors+i),*(*(enumSelectors+i)+**(enumSettings+i)));
	
	for(i=0;i<numIntListSettings;i++)
		saveIntList(settingsFile,*(intListDescriptors+i),*(intListSettings+i), **(numIntListValues+i));	
}
