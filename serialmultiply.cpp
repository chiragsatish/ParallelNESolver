#include<iostream>
#include<stdio.h>
#include<xmmintrin.h>
#include "multiply.h"
using namespace std;

void matrix_multiplication_SSE(float* A, float* B,float * result,int n1,int m1,int m2,int lda,int ldb,float* Y, float *result_xty,int m3)
{
    //double st = omp_get_wtime();
    //cout<<"n1"<<n1<<" "<<"m1"<<m1<<" "<<"m2"<<m2<<"\n";
    int i,j,k;
    for(i = 0; i < n1; ++i){
        for(j = 0; j < m2; ++j){
            for(k = 0; k < m1; ++k)
            {
                *(result+i*lda+j) += A[i*ldb+k] * B[k*lda+j];
                if(j==0)
                        *(result_xty+i*lda+j) += A[i*ldb+k] * Y[k*lda+j];
            }
        }
    }
}
