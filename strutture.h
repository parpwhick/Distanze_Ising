/* 
 * File:   strutture.h
 * Author: fake
 *
 * Created on January 13, 2012, 10:51 PM
 */


enum source{
    RANDOM,
    FILES
};

enum algorithm {
    NORMAL,
    SORTED,
    PMATRIX,
    AUTO
};

typedef struct {
    int seq_len;
    int n_seq;
    source from;
    char filename[255];
    bool write;
    bool distance;

    int fuzzy;
    int seed;
    int verbose;
    int n_symbols;
    bool translate;
    
    double beta;
    int threads;
    algorithm alg;

}options;

class partition {

public:
    //array of atom starting positions {1,4,5,...}
    //int *atom_positions;
    //binary atom beginning codes {1,0,0,0,1,1,...}
    int *binary;
    //number of atoms found
    int n;
    //total length of the partition
    int N;
    
    double entropia_shannon;
    double entropia_topologica;
    void print();
    
    template <typename T> void fill(const T*seq, int len);
    partition(int len=0){
        n=0;
        N=len;
        binary=0;
        
    }
    ~partition(){
        if(N && binary)
            delete[] binary;
         n=N=0;
    }
};

class distance{
private:
    void allocate(int n);
public:
    int *common_factor;
    int *reduced1;
    int *reduced2;
    int *product;
    int *product_reduced;

    double dist_s;
    double dist_s_r;
    int N;

    
    void binary_partition(const partition &e1, const partition &e2);
    
    ~distance();
    distance(int n=1000);
    distance(const distance &d1);
    

};


void set_program_options(options &opt, int argc, char**argv) ;
double entropy_binary_partition(int *p, int n);
