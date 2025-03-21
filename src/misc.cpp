#include "misc.h"

size_t mod_fwrite (const void *array, unsigned long long size, unsigned long long count, FILE *stream){
    unsigned long long pos, tot_size, pos_ct;
    const unsigned long long block_size = 4*512*512*512;
    
    tot_size = size*count; // total size of buffer to be written
    
    //check if the buffer is smaller than our pre-defined block size
    if (tot_size <= block_size)
        return fwrite(array, size, count, stream);
    else{
        pos = 0;
        pos_ct=0;
        // buffer is large, so write out in increments of block_size
        while ( (pos+block_size) < tot_size ){
            fwrite((char *)array+pos, block_size, 1, stream);
            pos_ct++;
            pos = block_size*pos_ct;
        }
        return fwrite((char *)array+pos, tot_size-pos, 1, stream);
    }
}


size_t mod_fread(void * array, unsigned long long size, unsigned long long count, FILE * stream){
    unsigned long long pos, tot_size, pos_ct;
    const unsigned long long block_size = 512*512*512*4;
    
    tot_size = size*count; // total size of buffer to be written
    //check if the buffer is smaller than our pre-defined block size
    if (tot_size <= block_size)
        return fread(array, size, count, stream);
    else{
        pos = 0;
        pos_ct=0;
        // buffer is large, so write out in increments of block_size
        while ( (pos+block_size) < tot_size ){
            fread((char *)array + pos, block_size, 1, stream);
            pos_ct++;
            pos = block_size*pos_ct;
        }
        return fread((char *)array+pos, tot_size-pos, 1, stream);
    }
}

int compare_floats (const void *a, const void *b)
{
    const float *da = (const float *) a;
    const float *db = (const float *) b;
    
    return (*da > *db) - (*da < *db);
}

float erfcc(float x)
{
    double t,q,ans;
    
    q=fabs(x);
    t=1.0/(1.0+0.5*q);
    ans=t*exp(-q*q-1.2655122+t*(1.0000237+t*(0.374092+t*(0.0967842+
                                                         t*(-0.1862881+t*(0.2788681+t*(-1.13520398+t*(1.4885159+
                                                                                                      t*(-0.82215223+t*0.17087277)))))))));
    return x >= 0.0 ? ans : 2.0-ans;
}

double FMAX(double a, double b) {
    if (a > b)
        return a;
    else
        return b;
}

