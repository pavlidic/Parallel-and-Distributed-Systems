#include <stdio.h>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "mmio.h"
#include "coo2csc.h"
#include "timediff.h"
#include <string.h>

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


    //--------------------------------------------------file reading-----------------------------------------------

    __cilkrts_end_cilk();

    if(argc!=4){
        printf("\nusage: ./v3_3mergeCilk ./*.mtx num_threads times_to_run\n");
        exit(1);
    }
    if(atoi(argv[2])>__cilkrts_get_nworkers() || atoi(argv[2])<1){
        printf("num_threads should be between 1 and %d on this machine\n",__cilkrts_get_nworkers());
        exit(1);
    }
    __cilkrts_set_param("nworkers",argv[2]);
    int threads =__cilkrts_get_nworkers();
    __cilkrts_init();

    if(atoi(argv[3])<0 || atoi(argv[3])>100){
        printf("1 <= num_times <= 100\n");
        exit(1);
    }
    int num_times=atoi(argv[3]);
    

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
    for (uint32_t i=0; i<nnz; i++)
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
    
    uint32_t i,j;
    uint32_t n1,n2,t1,t2;

    long long triangles;
    double minTime = 100;
    
    //------------temporary threadprivate c3 array------------

    uint32_t *c3temp[threads];

    //--------------------------------------------------main algorithm-----------------------------------------------

    for(int times=0;times<num_times;times++){
        
        for(i=0;i<M;i++)
            c3[i]=0;

        
        for(i=0;i<threads;i++){
            c3temp[i]=(uint32_t*)malloc(M*sizeof(uint32_t));
            memcpy(c3temp[i],c3,M*sizeof(uint32_t)); 
        }   
        

        clock_gettime(CLOCK_MONOTONIC, &tStart);
        
        cilk_for (uint32_t i=0;i<M-2;i++){
            cilk_for (uint32_t j=csc_row[i];j<csc_row[i+1]-1;j++){ 
                int threadid=__cilkrts_get_worker_number();
                t1=csc_row[csc_col[j]];
                t2=j+1;
                n1=csc_row[csc_col[j]+1];
                n2=csc_row[i+1];

                while(t1<n1 && t2<n2){

                    if(csc_col[t1]==csc_col[t2]){
                        
                        c3temp[threadid][i]++; c3temp[threadid][csc_col[j]]++; c3temp[threadid][csc_col[t2]]++;
                        
                        t1++;
                        t2++;
                    }else if(csc_col[t1]<csc_col[t2]){
                        t1++;
                    }else{
                        t2++;
                    }

                }
            }
        }

        for(i=0;i<threads;i++){
            for (j=0;j<M;j++){
                c3[j]+=c3temp[i][j];
            }
            free(c3temp[i]);
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

    printf("%s,%d,%.9lf,%lld\n",__progname,threads,minTime, triangles/3);


    free(csc_col);
    free(csc_row);
    
    free(c3);

    return 0;
}
