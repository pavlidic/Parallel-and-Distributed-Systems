#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "coo2csc.h"
#include "timediff.h"
#include "M312.h"

extern const char *__progname;


int main(int argc, char const *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nnz;
    uint32_t *I, *J;
    uint32_t *c3,*values;
    uint32_t isOneBased=0;
    struct timespec tStart,tEnd;
    
    
    //--------------------------------------------------file reading-----------------------------------------------

    if(argc!=3){
        printf("\nusage: ./v4 ./*.mtx times_to_run \n");
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
    
    I = (uint32_t *) malloc(nnz *2 * sizeof(uint32_t));
    J = (uint32_t *) malloc(nnz *2 * sizeof(uint32_t));

    for (int i=0; i<nnz; i++)
    {
        fscanf(f, "%u %u\n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
        I[i+nnz]=J[i];
        J[i+nnz]=I[i];
    }
    if (f !=stdin) fclose(f);

    //------------------------------------------coo2csc, memmory allocation and variable creation-----------------------------------------------

    c3      =(uint32_t *)malloc(M   * sizeof(uint32_t));
    values  =(uint32_t *)malloc(nnz * sizeof(uint32_t));

    uint32_t * csc_row  = (uint32_t *)malloc((M + 1) * sizeof(uint32_t));   //entire A matrix 
    uint32_t * csc_col  = (uint32_t *)malloc(nnz *2  * sizeof(uint32_t));
    uint32_t * csc_rowA = (uint32_t *)malloc((M + 1) * sizeof(uint32_t));   //uper triangular part of A
    uint32_t * csc_colA = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
    
    coo2csc(csc_col, csc_row, I,J,nnz*2,M,isOneBased);
    coo2csc(csc_colA,csc_rowA,I,J,nnz,  M,isOneBased);
    
    free(I);
    free(J);
    
    for(int i=0;i<M;i++)
        mergeSort(csc_col,csc_row[i],csc_row[i+1]-1);       //sorting of every row


    uint32_t i,j,r;
    uint32_t sum;
    uint32_t n1,n2,t1,t2;

    long long triangles;
    double minTime = 100;

    uint32_t* e=(uint32_t*)malloc(M*sizeof(uint32_t));
    for(int i=0;i<M;i++){
        e[i]=1;
    }

    //--------------------------------------------------main algorithm-----------------------------------------------
    
    for(int times=0;times<num_times;times++){
        
        for(i=0;i<nnz;i++)  values[i]=0;
        
        for(i=0;i<M;i++)    c3[i]=0;
        
        clock_gettime(CLOCK_MONOTONIC, &tStart);
        
        for(i=0;i<M;i++){
            for(j=csc_rowA[i];j<csc_rowA[i+1];j++){
                r=csc_colA[j];
                sum=0;
                t1=csc_row[i];
                t2=csc_row[r];
                n1=csc_row[i+1];
                n2=csc_row[r+1];
                
                while(t1<n1 && t2<n2){
                    if(csc_col[t1]==csc_col[t2]){
                        sum++;
                        t1++;
                        t2++;
                    }else if(csc_col[t1]<csc_col[t2]){
                        t1++;
                    }else{
                        t2++;
                    }
                }

                values[j]=sum;
                
            }
        }


        //----------------matrix multiplication with e (half312) and /2-------------------------------
        half312(M,csc_colA,csc_rowA,e,c3,values);
        for(int i=0;i<M;i++)
            c3[i]/=2;

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
 


    free(csc_row);
    free(csc_col);
    free(csc_rowA);
    free(csc_colA);
    free(c3);
    free(values);
    free(e);

}