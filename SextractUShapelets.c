#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "serial_compute.h"

#define UPPER 50
#define LOWER 45
#define STEP 5

void compute_mean_stddev(int *dataset_A, int dataset_Alen, double *dist, int newsize, double *mean, double *stddev)
{
    
	double sum = 0.0, avg = 0.0, dev = 0.0;
	int i;
    
	//for(i=0 ; i<dataset_Alen; i++) printf("%d ", dataset_A[i]);
    
	//for(i=0 ; i<newsize; i++) printf("%f ", dist[i]);
	//printf("\n");
    
	for(i = 0; i < dataset_Alen; i++) {
        
		sum += dist[dataset_A[i]];
		
	}
    
	printf("sum is %f\n", sum);
    
	avg = sum/dataset_Alen;
    
	for(i=0; i < dataset_Alen; i++) {
        
        dev += pow((dist[dataset_A[i]] - avg),  2.0);
        
	}
    
	*mean = avg;
	*stddev = sqrt(dev / dataset_Alen);
}

int clustered(int newcluster, int newsize) {
    
	//printf("Newsize is %d Newcluster is %d\n", newsize, newcluster);
	return ((newsize == newcluster) || (newcluster < 2));
}

int max_index(double *gap, int total)
{
	int max = 0;
	int i;
	for(i=1; i< total; i++) {
		//printf("gap of i is %f max is %f indx %d\n ", gap[i], gap[max], i);
		if(CompareDoubles2(gap[max],  gap[i]) < 0) {
			
			max = i;
		}
	}
    
	printf("\n\nmax is %d\n\n ", max);
	return max;
}
void extractU_Shapelets(double **pd_Dataset, int* ds_len, int n_sample, int sLen, int app_no, char app[][100], char file[][100], char*outputname, int start_ts_id)
{
    
	int sl;
	int i;
	int iter = 0;
	int cluster_id[n_sample];
	double *ts;
	int ts_len;
	int cnt;
	int newsize;
	int k,j,l;
	int index, index2;
	double mean, stddev, range;
    
	/*Empty discriminatory ushapelets*/
	double* ushapelet[n_sample];
	int ushapelet_len[n_sample];
    
	ts = pd_Dataset[start_ts_id];
	ts_len = ds_len[start_ts_id];
    
	printf("\nextractU_Shapelet\n");
    
	memset(cluster_id, -1, n_sample *sizeof(int));
    
    FILE *fp;
    
	fp = fopen(outputname, "a");
    
	if(!fp) {
		perror(outputname);
		return;
	}
    
    fprintf(fp, "\n\n\nNew Clustering Batch------------------------------\n");
    fprintf(fp,"************************\n");
    fprintf(fp, "Iteration %d\n", iter);
    fprintf(fp, "Application %s\n Dataset path %s\n", app[start_ts_id], file[start_ts_id]);
    
	while(1) {
		
		cnt = 0;
		newsize = 0;
        
		for(sl=sLen-LOWER; sl <= sLen+UPPER && sl <= ts_len; sl+=STEP) {
			cnt += (ts_len - sl+1);
            //			printf("sl = %d ts_len =%d cnt = %d\n", sl, ts_len, cnt);
		}
        
		double* p_subseq[cnt];
		int ps_len[cnt];
		double gap[cnt];
		double dt[cnt];
        
        
		for (k = 0, j= 0 ; k< n_sample; k++) {
            
			if(cluster_id[k] == -1)
				newsize++;
		}
        
		double* n_dataset[newsize];
		int n_datalen[newsize];
		int k,j;
        
		int old_id[n_sample];
		memset(old_id, 0, newsize*sizeof(int));
        
		/*create new unclustered data set*/
        
		for (k = 0, j= 0 ; k< n_sample; k++) {
            
			if(cluster_id[k] == -1) {
                
				// copy the data set and len and keep a mapping
                
				n_dataset[j] = pd_Dataset[k];
				n_datalen[j] = ds_len[k];
				old_id[j] = k;
				j++;
			}
		}
        
		
		double dist[newsize];
        
		/*index of the distance  within threshold*/
		int dataset_A[newsize];//--------Check: declared in both gobal and local
		int dataset_Alen;//--------Check: declared in both gobal and local
        
        
		memset(dataset_A, -1, newsize *sizeof(int));
        
		/*********VERIFY*******/
        
		if(ts_len < sLen) {
			
			printf(" The time series is too short to classify\n");
			
			break; // break?
		}
		int cluster_no = app_no;
        
		/*****VERIFY END****/
        
        
		cnt = 0;
        
		/*For all possible subsequences for a timeseries from the new dataset*/
		for(sl=sLen-LOWER; sl <= sLen+UPPER && sl <= ts_len; sl+=STEP) {
            
			//printf("sl is %d ts_len -sl +1 is %d \n\n", sl, ts_len - sl + 1);
            
			for(i=0; i< ts_len - sl +1; i++) {
                
				p_subseq[cnt] = ts + i;
				ps_len[cnt] = sl;
				/*Compute the gap and threshold for each of the subsequence*/
				gap[cnt] = computeGap(&dt[cnt], cluster_no, p_subseq[cnt], ps_len[cnt], n_dataset, newsize , n_datalen);
				//printf("i=%d sl=%d ps_len=%d cnt = %d\n", i, sl, ps_len[cnt], cnt);
				//printf("gap is %2.6f dt is %2.6f\n", gap[cnt], dt[cnt]);
				cnt++;
			}
		}
        
		/*Find the subsequence which gives the maximum gap for the dataset*/
		printf("max gap is %d\n", cnt);
		index = max_index(gap, cnt);
        
		/*Add the discriminatory subsequence to the ushapelet list*/
		printf("Discovered ushapelet gap is %2.6f dt is %2.6f index %d len %d\n", gap[index], dt[index], index, ps_len[index]);
        
        fprintf(fp, "Shapelet: ");
		for(k=0; k<ps_len[index]; k++)
			fprintf(fp,"%2.2f ", p_subseq[index][k]);
		printf("\n");
        
        
		ushapelet[iter] = p_subseq[index];
		ushapelet_len[iter] = ps_len[index];
        fprintf(fp, "Shapelet len: %d\n",ps_len[index]);
		
		dataset_Alen = 0;
		j=0;
		for(l=0; l<newsize; l++) {
            
			/*Compute the minimum distance of the shapelet from each of the dataset */
			dist[l]= computeDistance(p_subseq[index],  ps_len[index], n_dataset[l], n_datalen[l]);
            
			/*If the computed distance is less than threshold then add to Dataset A*/
			if (CompareDoubles2(dist[l], dt[index]) <= 0) {
				//printf("distance within threshold %2.2f %2.2f %d\n", dist[l], dt[index], l);
				dataset_A[j] = l;
				j++;
				dataset_Alen++;
			}
		}
		
        
		if(clustered(dataset_Alen, newsize)) break;
		else {
            
			mean = 0.0;
			stddev = 0.0;
			range = 0.0;
            
			/*Compute the mean standard deviation and range of the Dataset A*/
			compute_mean_stddev(dataset_A, dataset_Alen, dist, newsize, &mean, &stddev);
            
			range = mean + stddev;
            
			//printf("%2.2f is the range mean %2.2f stddev %2.2f \n", range, mean, stddev);
            
			/*Exclude all the dataset within the range by marking it as clustered*/
            
			for (k = 0, j= 0 ; k< newsize; k++) {
                
				if(CompareDoubles2(dist[k], range) <= 0) {
					//printf("Clustered dataset %d\n", old_id[k]);
                    fprintf(fp, "Appname: %s Filename: %s\n", app[old_id[k]], file[old_id[k]]);
					cluster_id[old_id[k]] = iter;
				}
			}
            
            /*Find the dataset far away from the ushapelet*/
			//printf("max distance is at %d\n", newsize);
			index2 = max_index(dist, newsize);
			ts = n_dataset[index2];
            
			//printf("Finding next data set is at  %d\n", index2);
            
            fprintf(fp,"************************\n");
            fprintf(fp, "Iteration %d\n", iter+1);
            fprintf(fp, "Application %s\n Dataset path %s\n", app[old_id[index2]], file[old_id[index2]]);
            
		}
        
		++iter;
		
	}
    
    fprintf(fp,"************************\n");
    fprintf(fp, "Remaining set\n");
	for( int z=0; z< n_sample; z++) {
		
        if(cluster_id[z]==-1)
            fprintf(fp, "Appname: %s Filename: %s\n", app[z], file[z]);
        
	}
    
	fclose(fp);
	printf("\n");
    
}

//sine wave generation
double sin(double x)
{
    double res=0, pow=x, fact=1;
    
    for(int i=0; i<5; ++i)
    {
        res+=pow/fact;
        pow*=x*x;
        fact*=(2*(i+1))*(2*(i+1)+1);
    }
    
    return res;
}

double* generate_sinewave(double phase, double amp, double freq, int len){
    
    double* result = (double*) malloc(sizeof(double)*len);
    
    for (int i=0; i<len; i++) {
        
        result[i] = amp*sin(phase+(double)i/freq);
        
    }
    
    return result;
    
}

double* generate_stepwave(int step_len, int len, double amp){
    
    double* result = (double*) malloc(sizeof(double)*len);
    
    for (int i=0; i<len; i++) {
        
        result[i] = (double)(((i/step_len)%2)*2 - 1)*amp;
    }
    
    return result;
    
}

int main(int argc, char *argv[])
{
    int set_no = 30;//even number better
 	double* pd_Dataset[set_no];
	int ds_len[set_no];
    char appname[set_no][100];
    char inputfile[set_no][100];
	
    FILE* serial_synthesis_file = fopen("serial_synthesis_input.txt", "w");
    
    //three sets of data
    //set 1 sine
    for (int i=0; i<set_no/2; i++) {
        
        ds_len[i] = rand()%5 + 100;
        double phase = ((double) (rand()%10))/3.0;
        double amp = ((double) (rand()%10))/3.0+1.0;
        double freq = (double) (rand()%10)+1.0;
        pd_Dataset[i] = generate_sinewave(phase, amp, freq, ds_len[i]);
        sprintf(appname[i], "sinewave");
        sprintf(inputfile[i], "sinewave_%d.txt", i);
        
    }
    
    //set 2 step
    for (int i=set_no/2; i<set_no; i++) {
        
        ds_len[i] = rand()%5 + 100;
        int step = rand()%5+1;
        double amp = ((double) (rand()%10))/3.0+1.0;
        pd_Dataset[i] = generate_stepwave(step, ds_len[i],amp);
        sprintf(appname[i], "stepwave");
        sprintf(inputfile[i], "stepwave_%d.txt", i);
        
    }
    
    printf("\nInput set generated......");
    fprintf(serial_synthesis_file, "Input set:\n");
    for (int i=0; i<set_no; i++) {
        
        for (int j=0; j<ds_len[i]; j++) {
            fprintf(serial_synthesis_file, "%f  ", pd_Dataset[i][j]);
        }
        
        fprintf(serial_synthesis_file, "\n\n");
        
    }
    
    fclose(serial_synthesis_file);
	
    
    char* outputfile = "serial_synthesis_output.txt";
    
    for (int i=0; i<30; i++) {
        extractU_Shapelets(pd_Dataset, ds_len, set_no, 50, 3, appname, inputfile,outputfile,i);
    }
    
    
	for (int i = 0; i< 30; i++)
		free(pd_Dataset[i]);
    
}
