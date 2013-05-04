//
//  serial_compute.c
//  ComputeGap
//
//  Created by Mushroom on 4/12/13.
//  Copyright (c) 2013 Self. All rights reserved.
//

#include <stdio.h>
#include "serial_compute.h"

int compare (const void * a, const void * b)
{
    return CompareDoubles2(*(double*)a, *(double*)b);
}

int CompareDoubles2 (double A, double B)
{
    if ((nearbyint(A*1000000.0)-nearbyint(B*1000000.0))>0.0) {
        return 1;
    }else if ((nearbyint(A*1000000.0)-nearbyint(B*1000000.0))==0.0){
        return 0;
    }else{
        return -1;
    }
    
}


//return the gap and set the dt in the input parameter
//cluster_no is the target number of clusters of the total dataset
double computeGap(double* threshold_dt, int cluster_no, double* shapelet, int shapelet_len, double** dataset, int set_no,
                  int* data_len){
    
    double* distance_array = (double*) malloc(sizeof(double)*set_no);
    
    //compute the distance of the shapelet to each data sample
    //printf("Display distance array...\n");
    for (int i=0; i<set_no; i++) {
        
        distance_array[i] = computeDistance(shapelet, shapelet_len, dataset[i], data_len[i]);
        //printf("Distance of dataset[%i]: %f\n",i,distance_array[i]);
        
    }
    /*
     printf("\n Unsorted distance array\n");
     for (int i=0; i<set_no; i++) {
     printf("%f , ",distance_array[i]);
     }*/
    
    //sort the distance array in ascending order
    qsort(distance_array, set_no, sizeof(double), compare);
    //sorted distance array
    /*printf("\nSorted distance array\n");
     for (int i=0; i<set_no; i++) {
     printf("%f , ",distance_array[i]);
     }*/
    
    double maxGap = 0.0;
    double dt = 0.0;
    
    for (int i= 0; i<set_no-1; i++) {
        
        double tmp_dist = (distance_array[i]+distance_array[i+1])/2.0;
        
        //mark which cluster one sample is in
        //1 as D_a 0 as D_b
        int* data_cluster = (int*) malloc(sizeof(int)*set_no);
        
        //record the number of the points in each cluster
        int cluster_a = 0;
        
        //check which cluster it is in
        for (int j=0; j<set_no; j++) {
            
            if (CompareDoubles2(distance_array[j],tmp_dist)<0) {
                
                cluster_a ++;
                data_cluster[j] = 1;
                
            }else{
                
                data_cluster[j] = 0;
                
            }
            
        }
        
        //balance the set ratio
        //printf("cluster_a no: %d", cluster_a);
        
        if(cluster_a>0 && cluster_a<set_no){
            
            double ratio = (double)cluster_a/(double)(set_no-cluster_a);
            
            //printf("\nratio: %f 1/cluster_no: %f \n",ratio, 1.0/(double)cluster_no);
            
            if ( (CompareDoubles2((1.0/(double)cluster_no),ratio)<0) && (CompareDoubles2(ratio, 1.0 - (1.0/(double)cluster_no))<0) ){
                
                
                double cluster_a_mean = 0.0;
                double cluster_a_std = 0.0;
                double cluster_b_mean = 0.0;
                double cluster_b_std = 0.0;
                
                //calculate the distance mean of cluster_a and cluster_b
                for (int j=0; j<set_no; j++) {
                    
                    if (data_cluster[j]==1) {
                        
                        //in cluster a
                        cluster_a_mean += distance_array[j];
                        //cluster_a_std stores the sum of the power of 2
                        cluster_a_std = cluster_a_std +  pow(distance_array[j], 2.0);
                        
                    }else{
                        
                        //in cluster b
                        cluster_b_mean += distance_array[j];
                        cluster_b_std += pow(distance_array[j], 2.0);
                        
                    }
                }
                
                //cluster a
                cluster_a_mean = cluster_a_mean/(double)cluster_a;//normalize the sum to the mean
                if (CompareDoubles2((cluster_a_std/(double)cluster_a) - pow(cluster_a_mean, 2.0),0.0)==0) {
                    cluster_a_std = 0.000001;
                }
                else{
                    cluster_a_std = sqrt((cluster_a_std/(double)cluster_a) - pow(cluster_a_mean, 2.0));
                }
                
                //cluster b
                cluster_b_mean = cluster_b_mean / (double)(set_no - cluster_a);//normalize the sum to the mean
                if(CompareDoubles2(( cluster_b_std / (double)(set_no - cluster_a) ) - pow(cluster_b_mean, 2.0), 0.0)==0){
                    cluster_b_std = 0.000001;
                }else{
                    cluster_b_std = sqrt(( cluster_b_std / (double)(set_no - cluster_a) ) - pow(cluster_b_mean, 2.0));
                }
                
                double tmp_gap = cluster_b_mean - cluster_b_std - cluster_a_mean - cluster_a_std;
                //printf("\n%f, %f, %f, %f\n", cluster_b_mean , cluster_b_std,cluster_a_mean,cluster_a_std);
                //printf("tmp_gap at sample %d: %f\n", i ,tmp_gap);
                
                if ( CompareDoubles2(tmp_gap,maxGap)>0 ) {
                    
                    maxGap = tmp_gap;
                    dt = tmp_dist;
                    
                }
                
            }//if for ratio
            
        }//if cluster check
        
        free(data_cluster);
        
    }
    
    free(distance_array);
    
    //set the threshold dt
    *threshold_dt = dt;
    
    //return the maximum gap
    return maxGap;
    
    
}


//return the distance of each sample
double computeDistance(double* shapelet, int shapelet_len, double* data, int data_len){
    
    //z-normal
    double* shapelet_normalized;
    shapelet_normalized = zNormal(shapelet, shapelet_len);
    
    double min_dist = INFINITY;
    
    for (int j=0; j<(data_len-shapelet_len+1); j++) {
        
        //normalize the target comparing array
        double* comp_shapelet = (double*) malloc(sizeof(double)*shapelet_len);
        
        //create the subarray
        for (int k=0; k<shapelet_len; k++) {
            
            comp_shapelet[k] = data[j+k];
            
        }
        
        //znormal the subarray
        double* norm_shapelet = zNormal(comp_shapelet, shapelet_len);
        free(comp_shapelet);
        
        //calcuate the euclidean distance
        double dist = euclideanDistance(shapelet_normalized, norm_shapelet, shapelet_len);
        
        if(CompareDoubles2(dist,min_dist)<0){
            
            min_dist = dist;
            
        }
        
        free(norm_shapelet);
        
    }
    
    free(shapelet_normalized);
    
    return (min_dist/(double)shapelet_len);
}

double* zNormal(double* array, int len){
    
    //calculate the mean and std
    double tmp = 0.0;
    
    for (int i=0; i<len; i++) {
        
        tmp = tmp + array[i];
        
    }
    
    double mean = tmp/(double)len;
    double std;
    
    //reset tmp to zero
    tmp = 0.0;
    for (int i=0; i<len; i++) {
        
        tmp = tmp + pow((array[i]-mean), 2.0);
        
    }
    std = sqrt(tmp/(double)len);
    
    if (CompareDoubles2(std, 0.0) == 0) {
        std = 0.000001;//avoid zero
    }
    
    //calculate z value
    double* z_array = (double*) malloc(sizeof(double)*len);
    
    for (int i=0; i<len; i++) {
        
        z_array[i] = (array[i]-mean)/std;
        
    }
    
    return z_array;
}

double euclideanDistance(double* array1, double* array2, int len){
    
    double sum = 0.0;
    for (int i=0; i<len; i++) {
        
        sum = sum + pow((array1[i]-array2[i]), 2.0);
        
    }
    
    return sum;
    
}
