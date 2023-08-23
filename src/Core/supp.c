#include "supp.h"

/* supp.h - Created by Anthony Ruth on June 17, 2017
 * Support library for kinetic monte carlo
 * 
 * 
 * */
 
double generaterandom()
{	
/******************************************/
/* Generates a random double from 0       */
/* to 1. The number of different numbers  */
/* allowed is equal to RAND_MAX           */
/******************************************/
	double random;
	random = ((double)rand())/((double)RAND_MAX);	
	return(random);
}

void padListLeft(double *list, int *numitems)
{
	listinsertdouble(list,*numitems,0,2* *(list) - *(list+1));
	(*numitems)++;
}

void padListRight(double *list, int *numitems)
{
	listinsertdouble(list,*numitems,*numitems,2.0 * *(list+*numitems-1) - *(list + *numitems - 2));
	(*numitems)++;
}

//Not in use
double listmin(double *list, int num)
{
	double min = *list;
	int i;
	for(i = 1;i<num;i++)
	if(*(list+i) < min)
	min = *(list+i);
	return min;
}

//Not in use
double listmax(double *list, int num)
{
	double max = *list;
	int i;
	for(i = 1;i<num;i++)
	if(*(list+i) > max)
	max = *(list+i);
	return max;
}
char *strcat2(char *str1, char *str2)
{
	//like strcat, but it puts the output into a new string.
	int l1 = strlen(str1);
	int l2 = strlen(str2);
	char *output = malloc((l1+l2+1)*sizeof(char));//1 extra char for null terminator
	strcpy(output,str1);
	strcat(output,str2);
	return output;
}

void mkdir2(char *name)
{
	char buffer[1000];
	sprintf(buffer,"mkdir %s",name);
	printf("%s\n",buffer);
	system(buffer);
}

int strip(FILE *infile, char *search, char *output,char *delimeter, int fill, int dorewind)
{
	/**************************************************************/
	/* finds an line that fits the mold                         */
	/* search(white space)delimeter(white space)output(white space)
	 * if fill = 1, fill in output. Else return after search is found
	 * */
	 
	 //printf("Searching for %s\n",search);
	 char line[1000];
	 char formatstring[1000];
     int error = 1;
	 while ( 1 == fscanf( infile, " %999[^\n]%*[\n]", line ) )
	 {
		 if(strstr(line,search))
		 {
             //printf("found line \n%s\n",line);
			 if(fill)
			 {
				 sprintf(formatstring,"%%*[^%s]%%*[%s]%%[^%s]",delimeter,delimeter,delimeter);
			 //printf("scan string \n%s\n format string \n%s\n",strstr(line,search),formatstring);
			 sscanf(strstr (line, search),formatstring,output);
			 //printf("%s\n",output);
			}
            error = 0;//clear error when found
			break; 
		 }
	 }
	 if(dorewind)
     rewind(infile);
           
     return error;
}

int intcontains(int *list, int search, int numitems)
{
	if(numitems == 0)
	return 0;
	
	int i;
	for(i = 0;i<numitems;i++)
	if(*(list+i) == search)
	return 1;
	
	return 0;
}

int binarysearch(int *list, int search, int numitems)
{
	//returns n if the list does not contain the search where n would be the array index of the search after an insertionsort
	//returns -n-1 if the list contains the search where n is the location in the list
	//Assumes the list is sorted in ascending order

//	printf("Searching for element %d in list:\n",search);
	//printlist(list,numitems);
	
	int min = 0;
	int max = numitems-1;
	int middle = (max+min)/2;
	int contains = 0;
	
	for(;min<=max;middle = (max+min)/2)
	{
		if(search < *(list+middle))
            {
            max = middle-1;
            //printf("Search is less than %d setting max to %d\n",*(list+middle),max);
            }
        else if (search > *(list+middle))
            {
            min = middle+1;
            //printf("Search is greater than %d setting max to %d\n",*(list+middle),min);
            }
        else
            {
            contains = 1;
            break;
            }
	}
	if(contains)
        return -middle-1;
	else
        return min;
}

int binarysearchDouble(double *list, double search, double epsilon, int numItems)
{
	//returns n if the list does not contain the search where n would be the array index of the search after an insertionsort
	//returns -n-1 if the list contains the search where n is the location in the list
	//Assumes the list is sorted in ascending order
//	printf("Searching for element %d in list:\n",search);
	//printlist(list,numItems);
	
	int min = 0;
	int max = numItems-1;
	int middle = (max+min)/2;
	int contains = 0;
	
	for(;min<=max;middle = (max+min)/2)
	{
		if(search < *(list+middle)-epsilon)
            {
            max = middle-1;
            //printf("Search is less than %d setting max to %d\n",*(list+middle),max);
            }
        else if (search > *(list+middle)+epsilon)
            {
            min = middle+1;
            //printf("Search is greater than %d setting max to %d\n",*(list+middle),min);
            }
        else
            {
            contains = 1;
            break;
            }
	}
	if(contains)
        return -middle-1;
	else
        return min;
}

void listinsert(int *list, int numItems,int location, int value)
{
	//Inserts value at location by bubbling up all data from [location:numItems-1]
	//printf("List before adding %d to spot %d\n",value,location);
	//printlist(list,numItems);
	int i;
	for(i=numItems;i>location;i--)
        *(list+i) = *(list+i-1);
	*(list+location) = value;
	//printf("List after adding\n");
	//printlist(list,numItems);
	
}

void listinsertdouble(double *list, int numItems,int location, double value)
{
	//Inserts value at location by bubbling up all data from [location:numItems-1]
	//printf("List before adding %d to spot %d\n",value,location);
	//printlist(list,numItems);
	int i;
	for(i=numItems;i>location;i--)
        *(list+i) = *(list+i-1);
	*(list+location) = value;
	//printf("List after adding\n");
	//printlist(list,numItems);
	
}

void printlist(int *list, int numItems)
{
	int i;
	for(i=0;i<numItems;i++)
        printf("%d %d\n",i,*(list+i));
}

int factorial(int i)
{
if(i == 0)
return 1;
else
return i * factorial(i-1);
}

int choose(int a, int b)
{
    if (b == 0) 
    return 1;
    return (a * choose(a - 1, b - 1)) / b;
}


void combination(int* c,int n,int p, int x){
 /** [combination c n p x]
 * get the [x]th lexicographically ordered set of [p] elements in [n]
 * output is in [c], and should be sizeof(int)*[p] */
    //printf("Finding the %dth set of %dC%d\n",x,n,p);
    if(p == 0)
    return;

    x++; //The way this is written it is 1 based instead of 0 based    
    if(p==1)
    {
    c[0] = x;
	return;
	}
    
    int i,r,k = 0;
    for(i=0;i<p-1;i++){
        c[i] = (i != 0) ? c[i-1] : 0;
        do {
            c[i]++;
            r = choose(n-c[i],p-(i+1));
            k = k + r;
        } while(k < x);
        k = k - r;
    }
    c[p-1] = c[p-2] + x - k;
}

void clipComments(char *string, char *commentFlags)
{
	if(commentFlags == NULL)
		return;
		
	int iChar;
	char *match;
	for (iChar=0;iChar<strlen(commentFlags);iChar++)
	{
		if ((match = strchr(string,*(commentFlags+iChar))) != NULL)
			*(match) = '\0';
	}
}
/*
void clipComments(char *string, char *commentFlags)
{
	if(commentFlags == NULL)
		return;
		
	int iChar, iMatch;
	for (iChar=0;iChar<strlen(commentFlags);iChar++)
	{
		if ((iMatch = strchr(string,*(commentFlags+iChar))) != NULL)
			*(string+iMatch) = '\0';
	}
}*/

void readInt(FILE *infile, char *descriptor, int *value, char *commentFlags)
{
    char buffer[1000];
    if(strip(infile,descriptor,buffer,"= ",1,1))
    printf("Nothing found for %s using default value of %d\n",descriptor,*value);
    else
    {
	clipComments(buffer,commentFlags);
    *value = atoi(buffer);
	printf("%s is %d\n",descriptor,*value);
	}
}

void readString(FILE *infile, char *descriptor, char *value, char *commentFlags)
{
    char buffer[1000];
    strcpy(buffer,"");
    if(strip(infile,descriptor,buffer,"= ",1,1))
    {
    printf("Nothing found for %s, using default value of %s\n",descriptor,value);
	}
    else
    {
	clipComments(buffer,commentFlags);
    printf("%s is %s\n",descriptor,buffer);
    strcpy(value,buffer);
	}
}

void readDouble(FILE *infile, char *descriptor, double *value, char *commentFlags)
{
    char buffer[1000];
    if(strip(infile,descriptor,buffer,"= ",1,1))
    printf("Nothing found for %s, using default value of %g\n",descriptor,*value);
    else
    {
	clipComments(buffer,commentFlags);
    *value = atof(buffer);
    printf("%s is %g\n",descriptor,*value);
	}
}

void saveString(FILE *outfile, char *descriptor, char *value)
{
    fprintf(outfile,"%s = %s\n",descriptor,value);
}

void saveInt(FILE *outfile, char *descriptor, int value)
{
    fprintf(outfile,"%s = %d\n",descriptor,value);
}

void saveDouble(FILE *outfile, char *descriptor, double value)
{
    fprintf(outfile,"%s = %g\n",descriptor,value);
}

void convolve(double *xout, double *yout, double *xina, double *yina, double *xinb, double *yinb, int numout, int numa, int numb, double llimit, double ulimit, int weighted)
{
    /* Convolves F(x) = int(G(a)H(x-a)da)
     * This is solved numerically by linearly interpolating H for each value x in F.
     * F(x) = Sum(G(a_i) * (H(b_j)*(b_{j+1}-x)+ H(b_{j+1})*(x-b_j))/(b_{j+1}-b_j))
     * Where here b_j means the nearest value b_j in H less than (x-a_i)
     * and b_{j+1} means the nearest value b_j+1 in H greater than (x-a_i)
     * This function assumes that xinb is sorted in ascending order
     * */
    //printf("Dimensions of the convolution are %d x %d x %d\n",numout,numa,numb);
    //printf("Performing a convolution with output limits %g to %g, input limit 1 %g to %g and input limit 2 %g to %g\n",llimit,ulimit,*xina,*(xina+numa-1),*xinb,*(xinb+numb-1)); 
    int ina, inb, out;
    double ai, bj, bjplus, x;
    for(out=0;out<numout;out++)
    {
        *(yout + out) = 0;
        x = *(xout + out) = llimit + (ulimit-llimit) * out /(numout-1);
        for(ina=0;ina<numa;ina++)
        {
            ai = *(xina+ina);
            inb = binarysearchDouble(xinb,x-ai,0,numb);
                if (inb == numb || inb == 0) //There is no value b_j which will suffice so go on to the next value of a_i
                continue;
            //This means the number exactly matchs the search criteria. The linear interpolation should still work
            if(inb < 0)
            inb = -inb -1;    
            bj = *(xinb+inb);
            bjplus = *(xinb+inb+1);
            if(weighted)
            *(yout+out) += *(yina+ina) * (*(yinb+inb) * (bjplus-(x-ai)) + *(yinb+inb+1) * ((x-ai)-bj)) / (bjplus-bj);            
            else
            *(yout+out) += (*(yinb+inb) * (bjplus-(x-ai)) + *(yinb+inb+1) * ((x-ai)-bj)) / (bjplus-bj);

        }
    }
}


void gaussianconvolution(double *xout, double *yout, double *xin, double *yin, int numout, int numin, double llimit, double ulimit, double sigma, int weighted)
{
    double gausspoints[numgausspoints], gaussian[numgausspoints];
    int i;
    for(i=0;i<numgausspoints;i++)
    {
        gausspoints[i] = (sigma*numsigma)*(2*(double)i/(numgausspoints-1) - 1.0);
        gaussian[i] =  1/sqrt(M_PI_2)/sigma * exp(-pow(gausspoints[i]/sigma,2)/2);
    }
    convolve(xout, yout, xin, yin, gausspoints, gaussian, numout, numin, numgausspoints, llimit, ulimit,weighted);
}

void stepconvolution(double *xout, double *yout, double *xin, double *yin, int numout, int numin, double llimit, double ulimit, int reverse, int weighted)
{
    //This calculates the convolution of a function yi(xi) with Theta(x-xi)
    //Using the reverse option this changes to Theta(xi-x)
    double thetapoints[numthetapoints], theta[numthetapoints];
    int i;
    for(i=0;i<numthetapoints;i++)
    {
        thetapoints[i] = thetawidth*(2*(double)i/(numthetapoints-1.0) - 1.0);
        if(reverse)
        theta[i] =  thetapoints[i]<0;
        else
        theta[i] =  thetapoints[i]>0;
    }
    convolve(xout, yout, xin, yin, thetapoints, theta, numout, numin, numthetapoints, llimit, ulimit, weighted);
}

int countLines(FILE *infile)
{
	char *line=NULL;
	size_t len=0;
	int numLines=0;
	while(getline(&line,&len,infile) != -1)
		numLines++;
	rewind(infile);
	return numLines;
}

int countColumns(FILE *infile)
{
	char *line=NULL;
	size_t len=0;
	int numColumns=1;
	//char *offset=0;
	getline(&line,&len,infile);
	while((line = strstr(line," ")) != NULL)
	{
		numColumns++;
		line = line + 1;
	}
	rewind(infile);
	return numColumns;
}

void weightToProbabilityRange(double *weight, double *probRange, int numItems)
{
	int i;
	double sum=0;
	for(i=0;i<numItems;i++)
		sum+=*(weight+i);
	
	*probRange = *weight/sum;
	for(i=1;i<numItems;i++)
		*(probRange+i) = *(probRange+i-1) + *(weight+i)/sum;

}

int chooseItem(double *probRange, int numItems)
{
	double rand = generaterandom();
	int choice;
	for(choice=0;choice<numItems;choice++)
	{
		//printf("rand %g rate %g\n",rand,hoprates[chosenhop]);
		if(rand < *(probRange+choice))
		break;
	}
	return choice;
}

//This function has caused a memory leak before. Be careful about calling it too many times without freeing the string
char *formatStep(int step, int maxSteps)
{
	int numformat = log10(maxSteps)+1;
	char formatString[100];
	sprintf(formatString,"%%0%dd",numformat);
	//printf("format string: %s\n",formatString);
	
	char *stepString = calloc((numformat+1),sizeof(char));
	sprintf(stepString,formatString,step);
	return stepString;
}

void listFiles(char *dir, char **filenames, int *numfiles)
{
	char command[1000];
	sprintf(command,"ls %s",dir);
	FILE *p = popen(command,"r");
	//char **filenames = malloc(maxfiles*sizeof(char *));
	
	char *line=NULL;
	size_t len=0;
	*numfiles = 0;
	while(getline(&line,&len,p) != -1)
	{
	filenames[*numfiles] = calloc(strlen(line),sizeof(char));
	strncpy(filenames[*numfiles],line,strlen(filenames[*numfiles]));//Remove the new line char
	printf("%s\n",filenames[*numfiles]);
	(*numfiles)++;
	}
	pclose(p);
}

void freeFileNames(char **filenames, int numfiles)
{
	int i;
	for(i=0;i<numfiles;i++)
	free(filenames[i]);
}

double dotProduct(double *v1, double *v2, int nDim)
{
	double sum;
	int i;
	for(i=0,sum=0;i<nDim;i++)
	sum+=*(v1+i) * *(v2+i);
	
	return sum;
}

double magnitude(double *vector, int nDim)
{
	double sum;
	int i;
	//printf("Vector is:");
	for(i=0,sum=0;i<nDim;i++)
	{
	//printf("%g\t",*(vector+i));
	sum+=*(vector+i) * *(vector+i);
	}
	//printf("\n");
	return sqrt(sum);
}

void normalize(double *vector, int len, double value)
{
	//printf("Normalizing to %g\n",value);
	double sum=0;
	int i;
	for(i=0;i<len;i++)
		sum += *(vector+i);
	//printf("Factor is to %g\n",value/sum);
	for(i=0;i<len;i++)
		if(sum == 0)
			*(vector + i) = 0;
		else
			*(vector + i) *= value/sum;
}

void derivativeArray(double *x, double *y, double *dydx, int numPoints)
{
	/* Given an array of y values and an array of x values, this calculates an array of derivatives dydx corresponding to each x value.
	 * At the first and last endpoint a forward difference and a backwards difference are used. 
	 * In the middle a central difference is used.
	 * */
	int i;
	for(i=0;i<numPoints;i++)
	{
			if(i==0)
			*(dydx) = -(*(y+1)-*(y))/(*(x+1)-*(x));	
			else if(i ==numPoints-1)
			*(dydx+numPoints-1) = -(*(y+numPoints-1)-*(y+numPoints-2))/(*(x+numPoints-1)-*(x+numPoints-2));	
			else
			*(dydx+i) = -(*(y+i+1)-*(y+i-1))/(*(x+i+1)-*(x+i-1));	
	}
}

double avgWithNAN(double value1, double value2)
{
	/* This allows averaging of two values one or both of which may be nan which in this case indicates no data rather than zero.
	 * 
	 * */
	 double sum=0;
	 int num=0;
	 if(!isnan(value1))
	 {
		 sum+= value1;
		 num++;
	 }
	 if(!isnan(value2))
	 {
		 sum+= value2;
		 num++;
	 }
	 if(num==0)
	 return NAN;
	 else
	 return sum/num;
}

void zeroPoint(double *series, int numPoints, int zeroStep)
{
	int iPoint;
	double seriesZero = *(series+zeroStep);
	for(iPoint=0;iPoint<numPoints;iPoint++)
	*(series+iPoint) -= seriesZero;
}

double seriesAt(double uEvaluate, double *u, double *v, int numItems)
{
	//This evaluates the expected value of v(ue) when given vi(ui).
	//ue is found by binary search of u. If it lies within the range
	//of ui, a linear interpolation is performed between the two surrounding
	//indicies, vi and vi+1. If it is not found, then a NAN is returned by
	//this function.
	int index = binarysearchDouble(u,uEvaluate,0,numItems);
	if((index == 0) || (index == numItems))//The value being search for lies outside of the range of the array u.
	return NAN;
	
	return interpolate(uEvaluate,*(u+index-1),*(u+index),*(v+index-1),*(v+index));
}

double interpolate(double uEvaluate, double u1, double u2, double v1, double v2)
{
	double value;
	if (!((uEvaluate >= u1 && uEvaluate <= u2) || (uEvaluate <= u1 && uEvaluate >= u2)))
	printf("This is extrapolation instead of interpolation. Limits %g and %g evaluating at %g\n",u1,u2,uEvaluate);

	//perform some error check here. This is extrapolation not interpolation.
	if (u2 == u1)//This is a divide by zero;, no need for interpolation anyway
	return (v1+v2)/2.0;
	value = v1 + (uEvaluate-u1)*(v2-v1)/(u2-u1);
	return value;
}

void linearInterpolationWeights(double *list, double search, int numItems, int *offset, double *weight1, double *weight2)
{
	/* This searches for an item in an ascending order list. If the search is between two items, the offset will be set to the index of the lowered-number item.
	 * The weights are proportional to how close the search is to each item and the weights sum to 1.
	 * If the search is before the first item, the offset is set to zero, weight1 = 1, weight2 = 0
	 * if the search is after the last item, the offset is numItems-2, weight1 = 0, weight2 = 1
	 * 
	 * */
	 int n = binarysearchDouble(list, search, 0, numItems);
	 if(n < 0)
	 n = -n - 1;

	 
	 if(n == 0)
	 {
		 search = *(list);
	 }
	 else
	 if(n > numItems-2)
	 {
		 n = numItems-2;
		 search = *(list+numItems-1);
	 }
	 else
	 n--;
	 
	 *(offset) = n;
	 
	 double difference = *(list+n+1) - *(list+n);
	 if(difference == 0)
	 {
	 *(weight1) = 0.5;
	 *(weight2) = 0.5;
	}
	else
	{
	 *(weight1) = (*(list+n+1)-search)/difference;
	 *(weight2) = (search-*(list+n))/difference;
	}
 }

double *unit(int num)
{
	//Returns an array of all 1s
	int i;
	double *array = malloc(num*sizeof(double));
	for(i=0;i<num;i++)
	*(array+i)=1;
	return array;
}

void  SetBit( uint32_t A[],  int k )
{
    A[k/32] |= 1 << (k%32);  // Set the bit at the k-th position in A[i]
}

int TestBit( uint32_t A[],  int k )
{
    return ( (A[k/32] & (1 << (k%32) )) != 0 );     
}
int roundup(int number, int multiplicity)
{
	if(number % multiplicity == 0)
	return number;
	else
	return (number/multiplicity+1)*multiplicity;
}
