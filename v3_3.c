#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "coo2csc.h"
#include "timediff.h"


extern const char *__progname;


int main(int argc, char const *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nnz;
    uint32_t *I, *J;
    uint32_t *c3;
    uint32_t isOneBased=0;
    struct timespec tStart,tEnd;
    uint32_t i=0,j=0,k=0,l=0;
    
    //--------------------------------------------------file reading-----------------------------------------------

    if(argc!=3){
        printf("\nusage: ./v3_3 ./*.mtx times_to_run\n");
        exit(1);
    }
    if(atoi(argv[2])<0 || atoi(argv[2])>100){
        printf("1 <= num_times <= 100\n");
        exit(1);
    }
    int num_times=atoi(argv[2]);

    if ((f = fopen(argv[1], "r")) == NULL) 
        exit(1);
    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    if(!mm_is_symmetric(matcode)){
        printf("matrix %s isnt symmetric\n",argv[1]);
        exit(1);
    }
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nnz)) !=0)
        exit(1);
    I = (uint32_t *) malloc(nnz * sizeof(uint32_t));
    J = (uint32_t *) malloc(nnz * sizeof(uint32_t));
    for (i=0; i<nnz; i++)
    {
        fscanf(f, "%u %u\n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    if (f !=stdin) fclose(f);

    //------------------------------------------coo2csc, memmory allocation and variable creation-----------------------------------------------

    c3      =(uint32_t *)malloc(M * sizeof(uint32_t));

    uint32_t * csc_row = (uint32_t *)malloc((M + 1) * sizeof(uint32_t));
    uint32_t * csc_col = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
    
    coo2csc(csc_col,csc_row,I,J,nnz,M,isOneBased);
    
    free(I);
    free(J);

    long long triangles;
    double minTime = 100;


    //--------------------------------------------------main algorithm-----------------------------------------------
    for(int times=0;times<num_times;times++){

        for(i=0;i<M;i++) c3[i]=0;
        
        clock_gettime(CLOCK_MONOTONIC, &tStart);
        
        for(i=0;i<M-2;i++){                                          
            for(j=csc_row[i];j<csc_row[i+1]-1;j++){
                for(l=csc_row[csc_col[j]];l<csc_row[csc_col[j]+1];l++){
                    for(k=j+1;k<csc_row[i+1];k++){                              //col[k] = k v2 
                        if(csc_col[l]==csc_col[k]){
                            c3[i]++; c3[csc_col[j]]++; c3[csc_col[k]]++;
                        }
                    }
                }
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &tEnd);
        struct timespec tResult = diff(tStart,tEnd);
        
        double timeR = (double)(tResult.tv_sec+(double)tResult.tv_nsec/1000000000);
        if(timeR<minTime)
            minTime=timeR;

        triangles=0;
        for(i = 0; i<M;i++)
            triangles+=c3[i];

        

    }

    printf("%s,%d,%.9lf,%lld\n",__progname,1,minTime, triangles/3);
    
    free(csc_col);
    free(csc_row);
    free(c3);

    return 0;
}
