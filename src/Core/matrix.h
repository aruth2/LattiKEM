#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* matrix.h - Created by Anthony Ruth on June 17, 2017
 * 
 * Contains functions for basic matrix operations and formating matrices into strings
 * 
 * */
void matrixmultiplication(double *output, double *leftmatrix, double *rightmatrix, int m, int n, int o);
double *transpose(double *matrix, int a, int b);
double singlecofactor(double *matrix, int size, int a, int b);
double *cofactor (double *matrix, int size);
double determinant(double *matrix, int size);
double *inversematrix(double *matrix, int size);
void makematrix(char *output,double *matrix, int numberOfRows, int numberOfColumns, char *delimiter);
void printmatrix(double *matrix, int numberOfRows, int numberOfColumns, char *delimiter);
#endif
