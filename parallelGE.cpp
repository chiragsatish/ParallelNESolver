#include<iostream>
#include<cmath>
#include<omp.h>
//#include<chrono>
#include "GaussElimination.h"
using namespace std;

void GaussElimination(float* xxt, float *xty, int n,int lda, int ldb){
	//cout<<"In parallel GE\n";
        double Augmat[n][n+1];
        double XVec[n];
        int i,j,k,maxRow;
	double maxelement;      
        //clock_t start_t, middle_t, end_t, first_t, second_t, total_t; 
        //auto start = chrono::high_resolution_clock::now();
        #pragma omp parallel shared(i,Augmat,n,maxelement,maxRow) private(j,k) if(n>100)
		{
		#pragma omp for schedule(static) collapse(2) //if(n>100)
        	for(j=0;j<n;j++){
                	for(k=0;k<n+1;k++){
                        	if(k==n)
                            		Augmat[j][k]=*(xty+j*lda);
                       		else
                            		Augmat[j][k]=*(xxt+j*lda+k);
               		 }
       		}
	}
	/*cout<<"AUGMAT\n";
	for(i=0;i<n;i++){
		for(j=0;j<n+1;j++)
			cout<<Augmat[i][j]<<" ";
		cout<<"\n";
	}*/
	//#pragma omp parallel shared(Augmat,n,maxelement,maxRow,i) private(j,k) num_threads(5)
	//{
        for (i=0; i<n; i++) {
        // Search for maximum in this column
           maxelement = abs(Augmat[i][i]);
           maxRow = i;
		#pragma omp parallel shared(Augmat,n,maxelement,maxrow,i) private(j,k)  if(n>100)
		{
		#pragma omp for schedule(static) if(n>100)
            	for (k=i+1; k<n; k++) {
			#pragma omp critical
			{
                		if (abs(Augmat[k][i]) > maxelement) {
                   			maxelement = abs(Augmat[k][i]);
                    			maxRow = k;
                		}
			}
       		}
        	// Swap maximum row with current row (column by column)
        	#pragma omp for schedule(static) if(n>100)
       		 for (k=i; k<n+1;k++) {
           		 double tmp = Augmat[maxRow][k];
           		 Augmat[maxRow][k] = Augmat[i][k];
           		 Augmat[i][k] = tmp;
       		 }
        // Make all rows below this one 0 in current column
        	 #pragma omp for schedule(static) collapse(2) if(n>100) 
       		 for (k=i+1; k<n; k++) {
           		 double c = -Augmat[k][i]/Augmat [i][i];
           		 for (j=i; j<n+1; j++) {
               			 if (i==j) {
                   			 Augmat[k][j] = 0;
               			 } 
				else {
                   			 Augmat[k][j] += c * Augmat[i][j];
               			 }
            		}
        	}
		}
	}    
    //auto backfillstart = chrono::high_resolution_clock::now();
    for (i=n-1; i>=0; i--) {
        XVec[i] = Augmat[i][n]/Augmat[i][i];
        for (k=i-1;k>=0; k--) {
            Augmat[k][n] -= Augmat[k][i] * XVec[i];
        }
    }
    //auto backfillend = chrono::high_resolution_clock::now();
    //auto end = chrono::high_resolution_clock::now();
    //auto duration = chrono::duration_cast<chrono::nanoseconds>(end-start);
    //auto backfillduration = chrono::duration_cast<chrono::nanoseconds>(backfillend-backfillstart);
    //cout<<"Time taken:"<<duration.count()<<"ns"<<"\n";
    //cout<<"Backfill duration:"<<backfillduration.count()<<"ns"<<"\n";
    for (i=0;i<n;i++)
        cout<<"x"<<i<<" = "<<XVec[i]<<endl;
}
