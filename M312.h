#ifndef M312_DOT_H
#define M312_DOT_H

#include <stdint.h>
#include <stdbool.h>

//matrix multiplication when we have the whole matrix in csc form
void whole312(uint32_t M,uint32_t* col,uint32_t* row,uint32_t* u,uint32_t* c,uint32_t* val){
    uint32_t i,j;
    for(i=0;i<M;i++)
        c[i]=0;
    for(i=0;i<M;i++)
        for(j=row[i];j<row[i+1];j++)
            c[i]+=u[col[j]]*val[j];
}

//matrix multiplication when we only have half the matrix in csc form
void half312(uint32_t M,uint32_t* col,uint32_t* row,uint32_t* u,uint32_t* c,uint32_t* val){
    uint32_t i,j;
    for(i=0;i<M;i++)
        c[i]=0;
    for(i=0;i<M;i++)
        for(j=row[i];j<row[i+1];j++){
            c[i]+=u[col[j]]*val[j];
            if(i!=col[j]){
                c[col[j]]+=u[i]*val[j];
            }
        }
}

uint32_t exists(uint32_t* col,uint32_t* row,uint32_t i,uint32_t j){
    uint32_t result=0;

    if(j<i){
        int temp = i;
        i=j;
        j=temp;
    }

    for(uint32_t k=row[i]; k<row[i+1]; k++)
        if(col[k]==j){
            result=1;
            break;
        }
    
    return result;
}

void hadamard(uint32_t* col, uint32_t* row, long long* values, uint32_t M, uint32_t nnz){
    uint32_t i,j,k,count=0,r;
    long long sum;

    for(i=0;i<nnz;i++)
        values[i]=0;
    for(i=0;i<M;i++){
        //printf("%d\n ",i);
        for(j=row[i];j<row[i+1];j++){
            r=col[j];
            sum=0;
            for(k=0;k<M;k++)
                if(exists(col,row,i,k) && exists(col,row,k,r))
                    sum++;
            values[count]=sum;
            count++;
        }
    }


}

void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
 
    /* create temp arrays */
    int L[n1], R[n2];
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
 
/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(int arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
 
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
 
        merge(arr, l, m, r);
    }
}




#endif
