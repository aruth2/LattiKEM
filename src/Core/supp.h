#ifndef _SUPP_H_
#define _SUPP_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define numsigma 3
#define numgausspoints 100

#define numthetapoints 100
#define thetawidth 2.0

double generaterandom();
double listmin(double *list, int num);
double listmax(double *list, int num);
char *strcat2(char *str1, char *str2);
void mkdir2(char *name);
int strip(FILE *infile, char *search, char *output,char *delimeter, int fill, int dorewind);
int intcontains(int *list, int search, int numitems);
int binarysearch(int *list, int search, int numitems);
int binarysearchDouble(double *list, double search, double epsilon, int numitems);
void listinsert(int *list, int numitems,int location, int value);
void listinsertdouble(double *list, int numitems,int location, double value);
void printlist(int *list, int numitems);
int factorial(int i);
int choose(int a, int b);
void combination(int* c,int n,int p, int x);
void readString(FILE *infile, char *descriptor, char *value, char *commentFlags);
void readInt(FILE *infile, char *descriptor, int *value, char *commentFlags);
void readDouble(FILE *infile, char *descriptor, double *value, char *commentFlags);
void convolve(double *xout, double *yout, double *xina, double *yina, double *xinb, double *yinb, int numout, int numa, int numb, double llimit, double ulimit, int unweighted);
void gaussianconvolution(double *xout, double *yout, double *xin, double *yin, int numout, int numin, double llimit, double ulimit, double sigma, int weighted);
void stepconvolution(double *xout, double *yout, double *xin, double *yin, int numout, int numin, double llimit, double ulimit, int reverse, int weighted);
int countLines(FILE *infile);
int countColumns(FILE *infile);
void weightToProbabilityRange(double *weight, double *probRange, int numItems);
int chooseItem(double *probRange, int numItems);
char *formatStep(int step, int maxSteps);
void listFiles(char *dir, char **filenames, int *numfiles);
void freeFileNames(char **filenames, int numfiles);
double dotProduct(double *v1, double *v2, int nDim);
double magnitude(double *v, int nDim);
void derivativeArray(double *x, double *y, double *dydx, int numPoints);
double avgWithNAN(double value1, double value2);
void padListLeft(double *list, int *numitems);
void padListRight(double *list, int *numitems);
void zeroPoint(double *series, int numPoints, int zeroStep);
double seriesAt(double uEvaluate, double *u, double *v, int numItems);
double interpolate(double uEvaluate, double u1, double u2, double v1, double v2);
void linearInterpolationWeights(double *list, double search, int numItems, int *offset, double *weight1, double *weight2);
void saveString(FILE *outfile, char *descriptor, char *value);
void saveInt(FILE *outfile, char *descriptor, int value);
void saveDouble(FILE *outfile, char *descriptor, double value);
double *unit(int num);
void normalize(double *vector, int len, double value);
void  SetBit( uint32_t A[],  int k );
int TestBit( uint32_t A[],  int k );
int roundup(int number, int multiplicity);

#endif
