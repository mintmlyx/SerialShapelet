//
//  main.c
//  ComputeGap
//
//  Created by Mushroom on 4/12/13.
//  Copyright (c) 2013 Self. All rights reserved.
//

#include <stdio.h>
#include "serial_compute.h"
#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

int main(int argc, const char * argv[])
{

    // insert code here...
    int shapelet_len = 1000;
    double* shapelet = (double*) malloc(sizeof(double)*shapelet_len);
    
    
    //generate the shapelet
    printf("shapelet: ");
    for (int i=0; i<shapelet_len; i++) {
        
        printf(" , ");
        shapelet[i] = pow(((double)i)*0.83, 2.0)+((double)shapelet_len)*0.4;
        printf("%f", shapelet[i]);
        
    }
    
    int set_no = 100;
    double** dataset = (double**) malloc(sizeof(double*)*20);
    int* data_len = (int*) malloc(sizeof(int)*set_no);
    
    for (int i=0; i<set_no; i++) {
    
        //set the first one exactly the same as the second one
        if(i<5){
            
            data_len[i] = shapelet_len;
            
        }else{
            
            data_len[i] = rand()%10000 + 1000;
            
        }
        
    }
    
    //generate the shapelet for each set
    for (int i=0; i<set_no; i++) {
        
        dataset[i] = (double*) malloc(sizeof(double)*data_len[i]);
        
        
        printf("\nDataset[%d]: ", i);
        for (int j=0; j<data_len[i]; j++) {
            
            if (i<5) {
                
                dataset[i][j] = (double) (shapelet[j]+i);
                
            }else{
                
                dataset[i][j] = ((double)(rand()%1000))*0.89;//log(((double)(i+j))*0.45+pow(j, 2.0));
            }
            
            //generate the functions
            printf(" , ");
            printf("%f", dataset[i][j]);
            
        }
        
    }
    
    time_t start,end;
    start = clock();
    //test the computeGap
    double* threshold_dt = (double*) malloc(sizeof(double));
    double gap = computeGap(threshold_dt, 4, shapelet, shapelet_len, dataset, set_no, data_len);
    
    printf("\nGap: %f", gap);
    printf("\nDt : %f", (*threshold_dt));
    end = clock();
    
    long execution_time = (end-start)/CLOCKS_PER_SEC;
    printf("Execution time: %ld", execution_time);
    
    return 0;
}



