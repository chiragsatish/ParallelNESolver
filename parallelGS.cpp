#include<iostream>
#include<iomanip>
#include<cmath>
#include "GaussSiedel.h"
#include<time.h>
#include<omp.h>
//#include<chrono>
using namespace std;

int GaussSiedel(float* xxt, float *xty, int n,int lda, int ldb)
{
    //cout<<"In parallel GS\n";
    cout.precision(4);
    cout.setf(ios::fixed);
    int i,j,k,flag=0,count=0;
    double Augmat[n][n+1];          
    double XVec[n];                
    double res,y;
    float maxit = 2*10000*10000;
    #pragma omp parallel shared(Augmat,n) private(i,j) if(n>100) 
    {
    #pragma omp for schedule(static) collapse(2) if(n>100)
    for(i=0;i<n;i++){
            for(j=0;j<n+1;j++){
                    if(j==n)
                        Augmat[i][j]=*(xty+i*lda);
                    else
                        Augmat[i][j]=*(xxt+i*lda+j);
            }
    }
    }
    /*cout<<"AUGMAT\n";
        for(i=0;i<n;i++){
                for(j=0;j<n+1;j++)
                        cout<<Augmat[i][j]<<" ";
                cout<<"\n";
        }*/
    cout<<"Enter the initial values of the variables:\n";
    for (i=0;i<n;i++)
	XVec[i]=0;
        //cin>>XVec[i];
    cout<<"Enter the accuracy upto which you want the solution:\n";
    //cin>>res;
    res=0.000001;
    clock_t t;
    //auto start = chrono::high_resolution_clock::now(); 
    do{                           
        for (i=0;i<n;i++){               
            y=XVec[i];
            XVec[i]=Augmat[i][n];
	    #pragma omp parallel shared(XVec,Augmat,i) private(j) if(n>100) 
	    {
	    #pragma omp for schedule(static) if(n>100)
            for (j=0;j<n;j++){
		#pragma omp critical
		{
                if (j!=i){
                XVec[i]=XVec[i]-Augmat[i][j]*XVec[j];
            }
	    }
	    }
	    }
            XVec[i]=XVec[i]/Augmat[i][i];
            if (abs(XVec[i]-y)<=res){           
                flag++;
            }
	}
        count++;   
    }while(flag<n && count<=maxit);//If the values of all the variables don't differ from their previous values with error more than res, stop the loop
    //auto stop = chrono::high_resolution_clock::now(); 
    //auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start); 
    //cout<<"Time taken:"<<duration.count()<<"ns"<<"\n";
    if(flag<n){
        cout<<"Did not converge\n";
    }
    cout<<"Iterations:"<<count<<"\n";
    for (i=0;i<n;i++)
        cout<<"x"<<i<<" = "<<XVec[i]<<endl; 
    return flag;      
}
