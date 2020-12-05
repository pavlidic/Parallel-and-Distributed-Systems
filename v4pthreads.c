#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "mmio.h"
#include "coo2csc.h"
#include "timediff.h"
#include "M312.h"




struct thread_data{
    uint32_t  thread_id, M, maxThreads;
    uint32_t* csc_col;
    uint32_t* csc_colA;
    uint32_t* csc_row;
    uint32_t* csc_rowA;
    uint32_t* values;
};

void *triangleCounting(void* t){
    struct thread_data *mydata;
    mydata = (struct thread_data*)t;
    uint32_t j,i,r;
    uint32_t  thread_id=mydata->thread_id;
    uint32_t sum=0;
    uint32_t t1,t2,n1,n2;
    uint32_t M=mydata->M;
    uint32_t maxThreads=mydata->maxThreads;

    //uint32_t counter=0;
    //long long counter2=0;

    for(i=0;i<M;i++)
        for(j=mydata->csc_rowA[i];j<mydata->csc_rowA[i+1];j++){
            //printf("%d\n",(j&i)%maxThreads);
            if((j)%maxThreads==thread_id){         //<------hash function---------
                //counter++;
                sum=0;
                r=mydata->csc_colA[j];
                t1=mydata->csc_row[i];
                t2=mydata->csc_row[r];
                n1=mydata->csc_row[i+1];
                n2=mydata->csc_row[r+1];

                while(t1<n1 && t2<n2){
                    if(mydata->csc_col[t1]==mydata->csc_col[t2]){
                        sum++;
                        t1++;
                        t2++;
                    }else if(mydata->csc_col[t1]<mydata->csc_col[t2]){
                        t1++;
                    }else{
                        t2++;
                    }
                }
                //counter2+=sum;
                mydata->values[j]=sum;

            }
        }
    //printf("thread: %d, did %d of nnz, sum= %lld\n",thread_id,counter,counter2);

    
    pthread_exit(NULL);
}


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

    if(argc!=4){
        printf("\nusage: ./v4pthreads ./*.mtx num_threads times_to_run\n");
        exit(1);
    }
    if(atoi(argv[2])<1){
        printf("num_threads should be bigger than 0\n");
        exit(1);
    }
    int maxThreads=atoi(argv[2]);

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
    
    coo2csc(csc_col,csc_row,I,J,nnz*2,M,isOneBased);
    coo2csc(csc_colA,csc_rowA,I,J,nnz,M,isOneBased);
    
    free(J);
    free(I);
    
    long long triangles;
    double minTime = 100;

    uint32_t* e=(uint32_t*)malloc(M*sizeof(uint32_t));
    for(int i=0;i<M;i++){
        e[i]=1;
    }
    
    for(int i=0;i<M;i++)
        mergeSort(csc_col,csc_row[i],csc_row[i+1]-1);       //sorting of every row
    
    //-----------------------------------------thread creation----------------------------------------------------------------------------

    pthread_t threads[maxThreads];
    struct thread_data thread_data_array[maxThreads];
    
    uint32_t i;
    int rc;
    void *status;
    
    //--------------------------------------------------main algorithm-----------------------------------------------

    for(int times=0; times<num_times;times++){
        
        for(i=0;i<M;i++)    c3[i]=0;
        
        for(i=0;i<nnz;i++)  values[i]=0;


        clock_gettime(CLOCK_MONOTONIC, &tStart);

        for(i=0;i<maxThreads;i++){
            
            //-------------------saving data for each thread-------------------
            
            thread_data_array[i].csc_colA=csc_colA;
            thread_data_array[i].csc_rowA=csc_rowA;
            thread_data_array[i].csc_col=csc_col;
            thread_data_array[i].csc_row=csc_row;
            thread_data_array[i].values=values;
            thread_data_array[i].M=M;
            thread_data_array[i].thread_id=i;
            thread_data_array[i].maxThreads=maxThreads;

            rc = pthread_create(&threads[i],NULL,triangleCounting,(void*)&thread_data_array[i]);
        

            if (rc) {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }

        }

        //---------------------waiting for all threads to finish----------------------
        
        for(i=0;i<maxThreads;i++){
            rc = pthread_join(threads[i], &status);

            if (rc) {
                printf("ERROR; return code from pthread_join() is %d\n", rc);
                exit(-1);
            }
        }

        //----------------matrix multiplication with e (half312) and /2-------------------------------
        half312(M,csc_colA,csc_rowA,e,c3,values);
        for(i=0;i<M;i++)
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

    printf("%s,%d,%.9lf,%lld\n",__progname,maxThreads,minTime, triangles/3);


    free(csc_row);
    free(csc_col);
    free(csc_rowA);
    free(csc_colA);
    free(c3);
    free(values);
    free(e);
    pthread_exit(NULL);

}