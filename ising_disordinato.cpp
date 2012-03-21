/* 
 * Author: Dawid Crivelli
 *
 * Calcolo della distanza media tra sequenze di Ising-1D campionate a temperatura data
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <omp.h>
//#include <vector>
#include "strutture.h"
#include "rand55.h"

extern options opts;
double *mylog;

void print_array(int *array, int len, char *nome) {
    printf("%s {%d", nome, array[0]);
    for (int i = 1; i < len; i++)
        printf(",%d", array[i]);
    printf("}\n\n");
}

void old_create_ising_entries(options opts, int *buffer_sequenze, rand55 &generatore) {
    int L = opts.seq_len;
    int runs = opts.n_seq;
    int i = 0;
    double beta = opts.beta;
    int dH;

    //fprintf(stderr,"Making %d runs long %d, each with %d flips\n\n",runs,L,flips);


    int *chain = new int[L];
    int *nnr = new int[L];
    int *nnl = new int[L];
    int *J = new int[L];
    

    //int entry_count=0;

    //inizializzazione array primi vicini

    for (i = 0; i < L; i++) {
        nnl[i] = i - 1;
    }
    nnl[0] = L - 1;

    for (i = 0; i < L; i++) {
        nnr[i] = i + 1;
    }
    nnr[L - 1] = 0;


    for (int i = 0; i < L; i++) {
        chain[i] = 2 * (generatore() > .5) - 1;
        J[i] = 2 * (generatore() > .5) - 1;
    }

    int metodo = 3;
    double prob = 1 - exp(-2 * beta);
    /* METODO 1, funziona ok per beta<2, termalizza troppo lento dopo
     * aggiornamento di un sito per volta preso a random 
     */
    double E_med = 0;
    double m_med = 0;
    for (i = 0; i < runs; i++) {


        /* METODO 1
         * single spin flip
         */
        if (metodo == 1) {
            int flips = L * 1000;
            for (int k = 0; k < flips; k++) {
                int j = rand() % L;
                dH = 2 * chain[j]*(J[nnl[j]] * chain[nnl[j]] + J[j] * chain[nnr[j]]);
                if (dH < 0 || generatore() < exp(-beta * dH))
                    chain[j] *= -1;
            }
        }            
        /* METODO 2
         * Wolf / Swendsen-Wang
         */
        else if (metodo == 3) {
            int flips = (2 * L)/beta;
            for (int m = 0; m < flips; m++) {
                int k;
                int j = floor(generatore() * L);
                int sign = 2 * (generatore() > .5) - 1;
                int oldvalue = chain[j];
                int propagation;

                chain[j] = sign;

                propagation = 1;
                for (k = j + 1; k < L && (prob > generatore()); k++) {
                    propagation = propagation * J[k - 1];
                    if (chain[k] != propagation * oldvalue)
                        break;
                    chain[k] = sign*propagation;
                }
                //    printf("Changed <===%d ",k-j-1);
                propagation = 1;
                for (k = j - 1; k >= 0 && (prob > generatore()); k--) {
                    propagation = propagation * J[k];
                    if (chain[k] != propagation * oldvalue)
                        break;
                    chain[k] = sign*propagation;
                }
                //    printf("%d ===> signs\n",j-k-1);
            }
        }
        //messa in buffer e calcolo medie
        for (int k = 0; k < opts.seq_len; k++) {
            buffer_sequenze[i * opts.seq_len + k] = chain[k] == 1 ? '+' : '-';
            // printf("%c", chain[k] == 1 ? '+' : '-');
        }
        //printf("\n");
        int E = 0;
        int m = 0;
        for (int k = 0; k < opts.seq_len; k++) {
            int miniE = J[k] * chain[k] * chain[(k + 1) % L];
            if (miniE > 1 || miniE < -1)
                fprintf(stderr, "WRONG! %d = %d * %d * %d\n", miniE, J[k], chain[k], chain[(k + 1) % L]);
            E += miniE;
            m += chain[k];
        }
        E_med += E;
        m_med += m;
        //        fprintf(stderr,"%d:\t %d\t %d\n",i,E,m);
    }

    E_med /= runs;
    m_med /= runs;
    //    fprintf(stderr,"<e>=%f,<m>=%f\n",E_med/L,m_med/L);


}


// J = {-1, 1}
void ising_entries_jflip(options opts, int *buffer_sequenze, rand55 &generatore) {
    int L = opts.seq_len;
    int runs = opts.n_seq;
    double beta = opts.beta;

    int *flipchain = new int [L];
    int *chain = new int[L];
    int *J = new int[L];    
    
    //probabilita di trovare un flip, ovvero -1
    
    double prob = exp(-2 * beta);
    prob = prob / (1 + prob);
       

    for (int i = 0; i < L; i++)         
        J[i] = 2 * (generatore() > .5) - 1;
    
    
    for (int i = 0; i < runs; i++) {
        chain[0] = 2 * (generatore() > .5) - 1;
        
        for (int k = 0; k < L; k++)
            flipchain[k] = (prob  >  generatore() ) ? -1 : 1;
        
        for (int k = 1; k < L; k++) 
            chain[k] = flipchain[k] * chain[k-1] * J[k];
        
        for (int k = 0; k < L; k++) {            
            buffer_sequenze[i * L + k] = chain[k] == 1 ? '+' : '-';
            // printf("%c", chain[k] == 1 ? '+' : '-');
        }
    }

    delete []flipchain;
    delete []chain;
    delete []J;
}


// J = normale(media=0,std=1)
void ising_entries_jnorm(options opts, int *buffer_sequenze, rand55 &generatore) {
    int L = opts.seq_len;
    int runs = opts.n_seq;
    double beta = opts.beta;

    int *flipchain = new int [L];
    int *chain = new int[L];
    int *J = new int[L];     
    double *prob = new double[L];
          
    for (int i = 0; i < L; i++){
        //probabilita di trovare un flip, ovvero -1
		//J gaussiano (positivo)
        prob[i] = exp(-2 * beta * generatore.semi_norm());
        //J uniforme [0,1] (positivo)
		prob[i] = exp(-2 * beta * generatore.rand());
        prob[i] /= (1 + prob[i]);
        // J solo positivi
		J[i] = 1; //2 * (generatore() > .5) - 1;
    }
        
    for (int i = 0; i < runs; i++) {
        chain[0] = 2 * (generatore() > .5) - 1;
        
        for (int k = 0; k < L; k++)
            flipchain[k] = (prob[k]  >  generatore() ) ? -1 : 1;
        
        for (int k = 1; k < L; k++) 
            chain[k] = flipchain[k] * chain[k-1] * J[k];
        
        for (int k = 0; k < L; k++) {            
            buffer_sequenze[i * L + k] = chain[k] ;
        }
    }

    delete []flipchain;
    delete []chain;
    delete []J;
}


template <typename T>
int hamming(T* seq1, T* seq2, int L) {
    int dist = 0;
    for (int i = 0; i < L; i++)
        dist += seq1[i] == seq2[i];
    return (dist);
}

int main(int argc, char** argv) {
   

    FILE *out;

    set_program_options(opts, argc, argv);

    //random number initialization
    out = fopen("/dev/urandom", "r");
    int bytes_read=fread(&opts.seed, sizeof (unsigned int), 1, out);
    if(!bytes_read)
        printf("Can't initialize random number generation\n");
    fclose(out);
    //fprintf(stderr,"Seed initialized to %d\n",opts.seed);
    srand(opts.seed);

    //logarithm lookup table, 6x program speedup
    mylog = new double[3 * opts.seq_len + 10];
    for (int i = 1; i < 3 * opts.seq_len + 10; i++)
        mylog[i] = log(i);
    mylog[0] = 0;

    double media_globale = 0;
    double media_globale_n2 = 0;
    int n_estrazioni = 50;
    int runs = 0;

    
#pragma omp parallel
    {
        partition *partitions = new partition[opts.n_seq];
        int *buf_sequenze = new int[opts.n_seq * opts.seq_len];
        distance d(opts.seq_len);
        rand55 generatore;

        for (int L = 0; L < n_estrazioni; L++) {
            // Generazione di un nuovo vettore J_ij random
            // e di opts.n_seq sequenze che hanno quel J 

            double media_locale = 0;
            double media_locale_n2 = 0;
            ising_entries_jnorm(opts, buf_sequenze, generatore);
            
            //riempi le partizioni, a partire dalle sequenze date
            for (int i = 0; i < opts.n_seq; i++)
                partitions[i].fill(&buf_sequenze[i * opts.seq_len], opts.seq_len);


            //media delle distanze tra le coppie di sequenze generate
            //#pragma omp parallel for firstprivate(d) schedule(dynamic,10) reduction(+: media_n, media_n2)
            for (int i = 0; i < opts.n_seq; i++) {
                for (int j = i + 1; j < opts.n_seq; j++) {
#ifdef RIDUZIONE
                    d.binary_partition(partitions[i], partitions[j]);
                    media_locale += d.dist_s_r;
                    media_locale_n2 += (d.dist_s_r)*(d.dist_s_r);
#else
                    d.binary_partition(partitions[i], partitions[j]);

                    media_locale += d.dist_s;
                    media_locale_n2 += (d.dist_s)*(d.dist_s);
#endif
                }
            }
#pragma omp critical
            {
                media_globale += media_locale;
                media_globale_n2 += media_locale_n2;
                runs += 1;
            }


        }
    }
    


    double varianza_n;
    int Nd = runs * (opts.n_seq * (opts.n_seq - 1)) / 2;
    media_globale /= Nd;
    media_globale_n2 /= Nd;

    varianza_n = media_globale_n2 - media_globale*media_globale;
    printf("%d %f %f\n", opts.seq_len, media_globale, varianza_n);
    fprintf(stderr, "%d %f %f\n", opts.seq_len, media_globale, varianza_n);


    return 0;
}
