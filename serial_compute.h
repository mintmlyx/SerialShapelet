//
//  serial_compute.h
//  ComputeGap
//
//  Created by Mushroom on 4/12/13.
//  Copyright (c) 2013 Self. All rights reserved.
//

#ifndef ComputeGap_serial_compute_h
#define ComputeGap_serial_compute_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//function declaration
double* zNormal(double* array, int len);
double euclideanDistance(double* array1, double* array2, int len);
double computeDistance(double* shapelet, int shapelet_len, double* data, int data_len);
int compare (const void * a, const void * b);
double computeGap(double* threshold_dt, int cluster_no, double* shapelet, int shapelet_len, double** dataset, int set_no,int* data_len);

int CompareDoubles2 (double A, double B);

#endif
