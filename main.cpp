#include<iostream>
#include<chrono>
#include<xmmintrin.h>
#include<fstream>
#include "transpose.h"
#include "multiply.h"
#include "GaussElimination.h"
#include "GaussSiedel.h"
using namespace std;

#define ROUND_UP(x, s)  (((x)+((s)-1)) & -(s))
//set these values depending on matrix
const int n = 50000;
const int m = 5;
const int block_size=16;
int lda = ROUND_UP(m, 16); //to ensure that m, n are converted into multiples of block size
int ldb = ROUND_UP(n, 16);
//Variables for initial SSE Operations
float *A = (float*)_mm_malloc(sizeof(float)*lda*ldb, 64);//X matrix
float *B = (float*)_mm_malloc(sizeof(float)*lda*ldb, 64);//X transpose
float *Y = (float*)_mm_malloc(sizeof(float)*lda*ldb,64);
float *result_xxt = (float*)_mm_malloc(sizeof(float)*lda*ldb, 64);//A
float *result_xty = (float*)_mm_malloc(sizeof(float)*lda*ldb,64); //B
//Variables for finding X
int main(){
        int i,j;
        ifstream fp1("A.txt");
        ifstream fp2("B.txt");
        if(!fp1 || !fp2){
        	cout<<"Error opening the file\n";
        	return 255;
        }
	auto start = chrono::high_resolution_clock::now(); 
        for(i=0;i<n;i++){
                for(j=0;j<m;j++){
                        //*(A+i*lda+j)=rand()%10;
                        fp1>>*(A+i*lda+j);
                }
        }
        for(i=0;i<n;i++){
        		//*(Y+i*lda)=rand()%10;
        		fp2>>*(Y+i*lda);
        }
        transpose_block_SSE4x4(A,B,n,m,lda,ldb,block_size);
        //Print matrices
        /*for(i=0;i<n;i++){
                for(j=0;j<m;j++){
                        cout<<*(A+i*lda+j)<<" ";
                }
                cout<<"\n";
        }
        for(i=0;i<n;i++)
        	cout<<*(Y+i*lda)<<"\n";
        cout<<"----------------------------------\n";
        for(i=0;i<m;i++){
                for(j=0;j<n;j++){
                        cout<<*(B+i*ldb+j)<<" ";
                }
                cout<<"\n";
        }
        cout<<"-------------------------------------\n";*/
        matrix_multiplication_SSE(B, A,result_xxt,m,n,m,lda,ldb);
        matrix_multiplication_SSE(B,Y,result_xty,m,n,1,lda,ldb);
        /*for(i=0;i<m;i++){
                for(j=0;j<m;j++){
                        cout<<*(result_xxt+i*lda+j)<<" ";
                }
                cout<<"\n";
        }
        cout<<"----------------------------------------\n";
        for(i=0;i<m;i++){
        	cout<<*(result_xty+i*lda)<<"\n";
        }*/
        cout<<"GaussElimination\n";
	cout<<"------------------------------------------\n";
        GaussElimination(result_xxt,result_xty,m,lda,ldb);
	cout<<"------------------------------------------\n";
        int flag=0;
	cout<<"\nGaussSiedel\n";
	cout<<"------------------------------------------\n";
        flag = GaussSiedel(result_xxt,result_xty,m,lda,ldb);
	cout<<"------------------------------------------\n";
	auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start); 
	cout<<"Time taken:"<<duration.count()<<"ns"<<"\n";	
	FILE *fp = fopen("runtime","a");
	fprintf(fp,"%d\n",duration);
}
