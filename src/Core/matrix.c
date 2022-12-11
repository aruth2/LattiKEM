#include "matrix.h"

/* matrix.c - Created by Anthony Ruth on June 17, 2017
 * 
 * Contains functions for basic matrix operations and formating matrices into strings
 * 
 * */


void matrixmultiplication(double *output, double *leftmatrix, double *rightmatrix, int m, int n, int o)
{
	/**********************************************/
	/* multiplies a mxn matrix by a nxo matrix    */
	/* it is safe for output and rightmatrix to be*/
	/* the same or output and leftmatrix to be the*/
	/* same.
	/*
	/**********************************************/
	
	int im,in,io;
	double *values = malloc(m*o*sizeof(double));
	
	for(im=0;im<m;im++)
	for(io=0;io<o;io++)
	{
	for(in=0,*(values+im*o+io)=0;in<n;in++)
	*(values+im*o+io) += *(leftmatrix + im*n+in) * *(rightmatrix+in*o+io);
	}
	for(im=0;im<m;im++)
	for(io=0;io<o;io++)
	*(output+im*o+io) = *(values+im*o+io);
	free(values);
	
}

double *transpose(double *matrix, int m, int n)
{
    /*
     * Transposed and mxn matrix into its nxm adjoint
     * */
     
	double *transposedmatrix = malloc(m * n * sizeof(double));
	int im,in;
	for(im = 0;im<m;im++)
	for(in = 0;in<n;in++)
	*(transposedmatrix + in*m + im) = *(matrix +im * n + in);
	return transposedmatrix;
}

double singlecofactor(double *matrix, int size, int a, int b)
{
	//printf("Finding single cofactor %d x %d of matrix\n",a,b);
	//printmatrix(matrix,size,size);
			int k,l, kprime, lprime; //The indices of the submatrix
			double *submatrix = malloc((size-1)*(size-1)*sizeof(double));
			double value;
			
			for(k = 0;k<size-1;k++)
			{
				if(k>=a)
				kprime = k +1;
				else
				kprime = k;
			for(l = 0;l<size-1;l++)
			{
				if(l>=b)
				lprime = l+1;
				else
				lprime = l;
				*(submatrix+k*(size-1) +l) = *(matrix+kprime*(size) + lprime);
				}
			}
			value =determinant(submatrix,size-1);
			free(submatrix); 
			return value;
}

double *cofactor (double *matrix, int size)
{
	//Constructs a cofactor matrix from a matrix
	
	
		double *outputmatrix = malloc(size*size*sizeof(double));
		int i,j; //The indices of the matrix to use
		for(i = 0;i<size;i++)
		for(j = 0;j<size;j++)
		{
			*(outputmatrix + i*size+j) = pow(-1,i+j) * singlecofactor(matrix,size,i,j);
		}
	//printf("Cofactor Matrix\n");
	//printmatrix(outputmatrix,size,size);
	return outputmatrix;
}

double determinant(double *matrix, int size)
{
	//printf("Finding determinant of matrix\n");
	//printmatrix(matrix,size,size);
	if(size == 1)
	return *matrix;
	
	else
	{
		int i;
		double value = 0;
		for(i = 0;i<size;i++)
		value += pow(-1,i) * *(matrix+i) * singlecofactor(matrix,size,0,i);
		return  value;
	}

}

double *inversematrix(double *matrix, int size)
{
	//printf("Inverting matrix \n");
	//printmatrix(matrix,size,size);
	double *cofact = cofactor(matrix,size);
	double *outputmatrix = transpose(cofact,size,size);
	free(cofact);
	int i,j;
	double det = determinant(matrix,size);
	for(i = 0;i<size;i++)
	for(j = 0;j<size;j++)
	*(outputmatrix +i*size+j) /= det;
	return outputmatrix;
}



void makematrix(char *output,double *matrix, int numberOfRows, int numberOfColumns, char *delimiter)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	char buffer[100000];
	sprintf(output,"");
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
			sprintf(buffer,"%s",output);
			if ((i+1) % numberOfColumns == 0)
			{ 
				sprintf(output,"%s%.16g",buffer,*(matrix+i));  
				sprintf(buffer,"%s",output);
				sprintf(output,"%s\n",buffer);
			}
			else
			sprintf(output,"%s%.16g%s",buffer,*(matrix+i),delimiter);
			
	}
}

void printmatrix(double *matrix, int numberOfRows, int numberOfColumns, char *delimiter)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
			if ((i+1) % numberOfColumns == 0)
			{ 
				printf("%.16g",*(matrix+i));  
				
				printf("\n");
			}
			else
			printf("%.16g%s",*(matrix+i),delimiter);
			
	}
}
