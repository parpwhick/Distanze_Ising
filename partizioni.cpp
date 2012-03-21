#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include "strutture.h"

extern options opts;
extern double *mylog;


void partition::print() {
    // {1,3,4,...}
//    if (opts.verbose > 2) {
//        fprintf(stdout, "Simp: {%d", atom_positions[0]);
//        for (int k = 1; k < n; k++)
//            fprintf(stdout, ",%d", atom_positions[k]);
//        fprintf(stdout, "}\n");
//    }
    if (opts.verbose > 3) {
        // {1,0,0,1,0,1,...}
        fprintf(stdout, "Bin: {%d", binary[0]);
        for (int k = 1; k < N; k++)
            fprintf(stdout, ",%d", binary[k]);
        fprintf(stdout, "}\n");
    }
    fprintf(stdout, "Partitions[n]: %d, Shannon %f, Topological %f\n", n,entropia_shannon,
            entropia_topologica);
}


template <typename T>
void partition::fill(const T* seq, int len) {
    int i;
    
    //Total length of the partition is equal to the sequence
    N = len;
    n = 0;
        
    if(binary==0)
        binary = new int[len];
    
    //first one always start an atom
    binary[0] = 1;

    //checking for a different symbol from the one before
    //when that happens => new atom!
    for (i = 1; i < len; i++){
        binary[i] = seq[i] != seq[i - 1];
        n+= binary[i];
    }
        
    entropia_topologica=mylog[n];
    entropia_shannon=entropy_binary_partition(binary,N);

}
template void partition::fill(const int*,int);