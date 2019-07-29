/******************************************/
/***        Author: Pei-Wei Tsai        ***/
/***   E-mail: pwtsai@bit.kuas.edu.tw   ***/
/***         Version: 2.0.0.0           ***/
/***   Released on: November 1, 2009    ***/
/******************************************/

/* The "Benchmark.h" and the "Benchmark.cpp" composes a library contains i_FuncNum test functions. */
#include "Benchmark.h"			// Include its header file.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "fstream"
using namespace std;

double **p;
int NumDimofP;
int NumDataofP;
double* T;
int numNeuron;
int epoch;

double MF(double point , double center , double spread){
  double b = 1;
  double r3;
//  cout<<"point:"<<point<<endl;
 // cout<<"center:"<<center<<endl;
 // cout<<"spread:"<<spread<<endl;
  //cout<<"result"<<-1 * ( pow ( (double)(point - center) / ( spread +  0.1) , 2 * b )) <<endl;
  //r3 = exp(-1 * ( pow ( (double)(point - center) / ( spread  ) , 2 * b )) );
    r3 = 1/(1+ pow ((double)(point - center) / ( spread  ) , 2)) ;
 // cout<<"r3:"<<r3<<endl;
  return r3;
}

void InitializeF(){

  double numData = 100;
  
  p = new double* [4];
  for(int i = 0 ; i < 4 ; i++)
    p[i] = new double[(int)numData];

    
  
    string yasin = "yasin";
    ifstream fin;
    fin.open(yasin.c_str());
  
  NumDimofP = 4;
  NumDataofP = numData;
  
  T = new double[100]; 
  for(int i = 0 ; i < numData ; i++ ){
    fin>>T[i];
	}

  for(int i = 4 ; i < numData ; i++ ){
     p[0][i] = T[i-18];
     p[1][i] = T[i-12];
     p[2][i] = T[i-6];
     p[3][i] = T[i];
  }  
  
  numNeuron = 2;
  epoch = 5;
}

double min(double* numbers , int indices){
  double min = numbers[0];
  for(int i = 0 ; i < indices ; i++)
    if(min > numbers[i]){
      min = numbers[i];
    }
  return min;    
}

double max(double* numbers , int indices){
  double max = numbers[0];
  for(int i = 0 ; i < indices ; i++)
    if(max < numbers[i]){
      max = numbers[i];
    }
  return max;    
}

double sum(double* numbers , int indices){
  double result = 0 ;
  //cout<<"numbers[0]"<<numbers[0]<<endl;
  for(int i = 0 ; i < indices ; i++)
    result += numbers[i];
  
  return result;
}

double Benchmark(int i_fn,int i_dim,double d_pos[])
{

	int i, j;
	double d_tmp[4];	// variables for temporary storage.
	double d_result = 0;	// cotains the result to return to the main program.
  
	/* initialize the parameters - Start */
	for(i=0;i<4;i++)
		d_tmp[i]=0.0;
	/* initialize the parameters - End   */
	if(i_fn==1)	// Anfis Function
	{
	  double z = 0.0001;
	  int numDim = NumDimofP;
	  int numData = NumDataofP;
	  //
#include <stdint.h>
	  uint32_t NumRule = pow(numNeuron , numDim);
	  //
	  
	  double **a;
	  a = new double *[NumRule];
	  for(int i = 0 ; i < NumRule ; i++)
	    a[i] = new double[numDim + 1];
	  
	  unsigned int biggestInteger = -1;
	  srand(time(0));
	    
	  // defination of a;
	  for(int i = 0 ; i < NumRule ; i++)
	    for(int j = 0 ; j < numDim + 1 ; j++)
	      a[i][j] = 2 * ((double)(rand())/biggestInteger) - 1;
	  //
	  double *lowc = new double[numDim];
	  double *highc = new double[numDim];
	  double *NC = new double[numDim];
	  
	  for(int i = 0 ; i < numDim ; i++){
	    lowc[i] = min(p[i] , numDim);
	      highc[i] = max(p[i] , numDim);
	      NC[i] = (highc[i] - lowc[i]) / ( numNeuron + 1 );
	  }   
	  
	  // defination of c
	  double **c;
	  c = new double *[numNeuron];
	  for(int i = 0 ; i < numNeuron ; i++)
	    c[i] = new double[numDim];
	  
	  for(int i = 0 ; i < numNeuron ; i++)
	    for(int j = 0 ; j < numDim ; j++)
	      c[i][j] = lowc[j] + NC[j] * ( i + 1 );
	  
	   double *cen = new double[NumRule];
	  double *variance = new double[NumRule];
	  for(int i = 0 ; i < NumRule ; i++){
	    cen[i] = 0;
	    variance[i] = 1;
	  }
	  
	  //
	  // defination of s
	  double *s;
	  s = new double[numDim];
	  for(int i = 0 ; i < numDim ; i++)
	    s[i] = (highc[i] - lowc[i]) / pow(numNeuron , 2);
	  //
	  int max = numNeuron;
	  
	  // defination of Coff
	  double **Coff;
	  int index = pow(max,numDim);
	  Coff = new double *[index];
	  
	  for(int i = 0 ; i < index ; i++)
	    Coff[i] = new double[numDim];
	  for(int i = 0 ; i < numDim ; i++)
	    Coff[0][i] = 1;
	  int indexI = 0 ;
	  int indexK = 1;
	  int indexJ = numDim - 1;

	  // Calculating Coff
	  while(indexJ >= 0){
	    while( Coff[indexI][indexJ] < max ){
	      for(int counter = 0 ; counter < numDim ; counter++){
		Coff[indexK][counter] = Coff[indexI][counter];
		if(counter == indexJ)
		  Coff[indexK][counter] = Coff[indexI][indexJ] + 1;
	      }
	      indexK++;
	      indexI++;
	    }
	    indexI = 0;
	    indexJ--;
	  }
	  // end Calculating Coff
	  double error;
	  double e;
	  
	  // defination f1
	  double **f1;
	  f1 = new double*[numNeuron];
	  for(int i = 0 ; i < numNeuron ; i++)
	    f1[i] = new double[numDim];
	  
	  // defination f2
	  double* f2;
	  f2 = new double[(int)pow(numNeuron,numDim)];
	  
	  // defination f3;
	  double* f3;
	  f3 = new double[NumRule];
	  
	  // defination f4;
	  double* f4;
	  f4 = new double[NumRule];
	  
	  double sum1 = 0 ;
	  
	  // defination f5
	  double f5;
	  double RoundE;
   
	  double **RoundA;
	  RoundA = new double*[NumRule];
	  for(int i = 0 ; i < NumRule ; i++)
	    RoundA[i] = new double[numDim+1];

	  
	  
	  for(int ep = 0 ; ep < epoch ; ep++){
	    error = 0;
	    for(int j = 0 ; j < numData ; j++){
	      for(int i  = 0 ; i < numNeuron ; i++){
		for(int k = 0 ; k < numDim ; k++){
		  f1[i][k] = MF(p[k][j] , c[i][k] , s[k]);
		}
	      }
	      
	      for(int i = 0 ; i < pow(numNeuron,numDim) ; i++){
		f2[i] = 1;
		for(int k = 0 ; k < numDim ; k++){
		  f2[i] = ((double)f1[(int)Coff[i][k] - 1][k]) * f2[i];
		}
	      }
	      
	      for(int i = 0 ; i < NumRule ; i++){
		//cout<<"man injam:" << sum(f2,pow(numNeuron,numDim));
		f3[i] = f2[i] / sum(f2,pow(numNeuron,numDim));
	      }
	      
	      
	      double yi;
	      for(int i = 0 ; i < NumRule ; i++){
		if( i % 2 == 0){
		  yi = f3[i] *( cen[i] - variance[i] * sqrt((double)1/f3[i] - 1)); 

		}
		else{
		  yi = f3[i] *( cen[i] + variance[i] * sqrt((double)1/f3[i] - 1)); 

		}
		f4[i] = f3[i] * yi;

	      }
	      
	      int k;
	      f5 = sum(f4,NumRule);
	      e = 0.5 * pow( (T[j] - f5 ), 2);
	      error = error + e;
	      
	      for(int i = 0 ; i < numNeuron ; i++){
		 for(int k = 0 ; k < numDim ; k++){
		    c[i][k] = highc[k] - (pow(numNeuron , 2) * d_pos[numDim * (i) + k]) + NC[k] * ( i + 1 );

		 }
	    }
	      double* dposS = new double[numDim];
	      int base = numNeuron * numDim;
	      for(int i = 0; i < numDim ;i++)
		dposS[i] = d_pos[i+base];
		
		for(int i = 0 ; i < numDim ; i++){
		  s[i] = dposS[i];
		}
		
	      double* dposCen = new double[NumRule];
	      base += numDim;
	      for(int i = 0; i < NumRule ;i++)
		dposCen[i] = d_pos[i+base];
		
		for(int i = 0 ; i < NumRule ; i++)
		  cen[i]= dposCen[i];
		
		 double* dposVar = new double[NumRule];
		base += NumRule;
		for(int i = 0; i < NumRule ;i++)
		    dposVar[i] = d_pos[i+base];
 		
		for(int i = 0 ; i < NumRule ; i++)
 		  variance[i] = dposVar[i];	     
	    }
	  }
	  d_result = pow ( ( (double)1 / numData ) * error , 0.5);
	
	}
		
	
	
	
	
	
	
	else if(i_fn==2)	// Sum of different power functions
	{
		for(j=0;j<i_dim;j++)
		{
			d_result += pow(abs(d_pos[j]),j+1);
		}
	}
	else if(i_fn==3)	// Ackley functions
	{
		for(j=0;j<i_dim;j++)
		{
			d_tmp[0] += pow(d_pos[j],2);
			d_tmp[1] += cos(d_PII*d_pos[j]);
		}
		d_tmp[0] = d_tmp[0]/i_dim;
		d_tmp[0] = -1 * 0.2 * sqrt(d_tmp[0]);
		
		d_tmp[1] = d_tmp[1]/i_dim;

	
		d_result= -1 * 20 * exp(d_tmp[0]) - exp(d_tmp[1]) + 20 + exp(1);
	}
/*
	else if(i_fn==4)	// Test Function 4
	{
	}
	else if(i_fn==5)  // Test Function 5
	{
	}
	else // Test Function 6
	{
	}
*/
  return d_result;
}
