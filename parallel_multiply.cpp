#include<iostream>
#include<xmmintrin.h>
#include "multiply.h"
#include<omp.h>
#include<stdio.h>
using namespace std;


void matrix_multiplication_SSE(float* A, float* B, float* result, int n1, int m1, int m2, int lda, int ldb){
    int tid, nthreads, i, j, k, chunk;
    chunk=10;
    //cout<<"n1:"<<n1<<" "<<"m1:"<<m1<<" "<<"m2:"<<m2<<"\n";
    #pragma omp parallel shared(A,B,result,nthreads,chunk) private(tid,i,j,k) 
    {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
                nthreads = omp_get_num_threads();
		int maxNumThreads = omp_get_max_threads();
    		//printf("Maximum number of threads for this machine: %i\n", maxNumThreads);
                //printf("Starting matrix multiple with %d threads\n",nthreads);
         }       
	//printf("Thread %d starting matrix multiply...\n",tid);
 	#pragma omp for schedule (static, chunk) collapse(2)
  	for (i=0; i<n1; i++)    
   	{
    	//	printf("Thread=%d did row=%d\n",tid,i);
   		for(j=0; j<m2; j++){
		//#pragma omp parallel shared(i,j,A,B,result,n1,m1,m2,lda,ldb) private(k) num_threads(5)
		//{     
     		 for (k=0; k<m1; k++)
       			 *(result+i*lda+j) += A[i*ldb+k] * B[k*lda+j];
    		}
		}
	}   
}
