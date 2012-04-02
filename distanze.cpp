#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include "strutture.h"


extern options opts;
extern double *mylog;

void distance::allocate(int n) {
    if (opts.da_calcolare & (SHAN | RID)) {
        common_factor = new int[n];
        reduced1 = new int[n];
        reduced2 = new int[n];
        product_reduced = new int[n];
    }
    product = new u_int64_t[n];
}

distance::distance(int n) {
    N = n;
    allocate(N);
}

distance::distance(const distance &d1) {
    N = d1.N;
    allocate(N);
}

distance::~distance() {
    if (N) {
        if (opts.da_calcolare & (SHAN | RID)) {
            delete[] common_factor;
            delete[] reduced1;
            delete[] reduced2;
            delete[] product_reduced;
        }
        delete[] product;
    }
    N = 0;
}

void distance::fill(const linear_partition& e1, const linear_partition& e2) {
    this->binary_partition(e1, e2);
}

void print_binary_partition(int*p, int N) {
    // |...|...||...|...
    for (int i = 0; i < N; i++)
        printf("%c", (p[i]) ? '|' : '.');
    printf("\n");
}

double entropy_binary_partition(const int *p, int n) {
    //p: binary partition
    //n: total length
    int i = 0;
    int mu;
    int begin;
    double H = 0;

    //the first position always starts an atom
    begin = 0;
    for (i = 1; i < n; i++) {
        //whenever we find 1 (a new atom)
        if (p[i]) {
            //the closed (old)atom's length is calculated
            mu = i - begin;
            //the new one is ready to go
            begin = i;
            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            H += (double) mu * mylog[mu];
        }
    }
    //we check the last one, in case it's left hanging
    mu = n - begin;
    H += mu * mylog[mu];

    //proper entropy normalization
    H = -H / n + mylog[n];
    return (H);
}

//function giving the distance between 2 partitions, by intersecting 
//and calculating the relevant entropy

void distance::binary_partition(const linear_partition &first, const linear_partition &second) {

    int N = first.N;

#ifdef RIDUZIONE
    //the common factor of the partition, AND'ing
    for (int i = 0; i < N; i++)
        common_factor[i] = first.binary[i] & second.binary[i];
    common_factor[0] = 1;

    //reducing both the partitions XOR'ing with the common factor
    for (int i = 0; i < N; i++)
        reduced1[i] = common_factor[i] ^ first.binary[i];
    reduced1[0] = 1;
    for (int i = 0; i < N; i++)
        reduced2[i] = common_factor[i] ^ second.binary[i];
    reduced2[0] = 1;

    //intersection of the reduced partitions, OR'ing
    for (int i = 0; i < N; i++)
        product_reduced[i] = reduced1[i] | reduced2[i];
    //  alternatively, one could XOR the unreduced partitions
    //  intersection_reduced[i] =first->binary[i] ^ second->binary[i];
    product_reduced[0] = 1;
    
    double  hr1 = entropy_binary_partition(reduced1, N),
            hr2 = entropy_binary_partition(reduced2, N),
            hr12 = entropy_binary_partition(product_reduced, N);
    
    this->dist_s_r = 2 * hr12 - hr1 - hr2;
    
#else

    //the partition product/intersection, OR'ing
    for (int i = 0; i < N; i++)
        product[i] = first.binary[i] | second.binary[i];
    product[0] = 1;

    //calculating ALL the entropies (3 per non-reduced Rohlin dist, 3 per reduced)
    double
	    h1 = first.entropia_shannon,
            h2 = second.entropia_shannon,
            h12 = entropy_binary_partition(product, N);


    //packing the output tuple with the distances
    this->dist_s = 2 * h12 - h1 - h2;
    
#endif
}


void distance::fill(const general_partition& e1, const general_partition& e2) {

    if (opts.da_calcolare & GENERAL_RID) {
        partizione_comune.linear_intersection(e1, e2);
        //partizione_comune.trivial(e1.N);
        ridotto1.reduce(e1, partizione_comune);
        ridotto2.reduce(e2, partizione_comune);

     
        linear_product_sorted(ridotto1,ridotto2);
        dist_fuzzy_r = dist_fuzzy;
        dist_fuzzy_r_t = dist_fuzzy_t;
    }

    if (opts.da_calcolare & GENERAL) {
        switch (opts.alg) {
            case PMATRIX:
            case SORTED:
            default:
                linear_product_sorted(e1, e2);
                break;
        }
    }
}

inline int compare (const void * a, const void * b){
  return ( *(u_int64_t*)a - *(u_int64_t*)b );
}


void distance::linear_product_sorted(const general_partition &p1,const  general_partition &p2) {
    int i;
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;

    for (i = 0; i < N; i++){
        u_int64_t temp1=p1.labels[i];
        u_int64_t temp2=p2.labels[i];
    
        product[i] = (temp1<<32) | temp2;
    }
    
    
    qsort(product,N,sizeof(u_int64_t),compare);

    //the first position always starts an atom
    begin = 0;
    label_count=1;
    u_int64_t old_val=product[0];
    for (i = 1; i < N; i++) {
        //whenever we find 1 (a new atom)
        if (product[i]!=old_val) {
            //a new atom starts
            label_count++;
            //the closed (old)atom's length is calculated
            mu = i - begin;
            //the new one is ready to go
            begin = i;
            //cache the new label to check
            old_val=product[i];
            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];
        }
    }
    //we check the last one, in case it's left hanging
    mu = N - begin;
    H += mu * mylog[mu];
    
    //normalize the result
    H = -H / N + mylog[N];

    double h1 = p1.entropia_shannon,
            h2 = p2.entropia_shannon;
            
    this->dist_fuzzy = 2 * H - h1 - h2;
    


}

