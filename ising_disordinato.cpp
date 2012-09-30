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
#include <ctime>
#include <vector>
#include "strutture.h"
#include "rand_mersenne.h"
#include "distance.h"
#include "adj_handler.h"

extern options opts;
double *mylog;

adj_struct topologia;

void print_array(const int *array, int len, const char *nome) {
    printf("%s {%d", nome, array[0]);
    for (int i = 1; i < len; i++)
        printf(",%d", array[i]);
    printf("}\n\n");
}


int *nnl = 0,
        *nnu = 0,
        *nnd = 0,
        *nnr = 0;

#define up(i) (i - (i % lato)+ ((i+lato-1)%lato))
#define down(i) ((i/lato)*lato + ((i+lato+1)%lato))
#define left(i) (i+N-lato)%N
#define right(i) (i+N+lato)%N

void ising_lattice(options opts, RandMT &generatore, general_partition *partitions) {
    int N = opts.seq_len;
    int lato = opts.lato;
    int runs = opts.n_seq;
    int i = 0;
    double beta = opts.beta[0];
    int dH;
    int metodo;
    double prob;

    vector<int> chain(N);

    //int entry_count=0;

    //inizializzazione array primi vicini
    if (nnl == 0)
#pragma omp critical
    {
        nnl = new int[N];
        nnr = new int[N];
        nnu = new int[N];
        nnd = new int[N];
        for (int i = 0; i < opts.seq_len; i++) {
            nnl[i] = left(i);
            nnr[i] = right(i);
            nnu[i] = up(i);
            nnd[i] = down(i);

        }
    }

    if (beta < 0.4)
        prob = 0.5;
    else {
        //genero in base al fit quadratico della probabilita da un esperimento
        prob = log(beta);
        prob = 2.14 - 5.63 * prob - 7.71 * prob*prob;
        prob = exp(prob);
    }


    for (int i = 0; i < N; i++) {
        chain[i] = 2 * (generatore.get_double() > prob) - 1;
    }
    //print_array(chain,50,"chain");

    int flips = 2;
    if (beta < 0.36 || beta > 0.47)
        metodo = 1;
    else {
        metodo = 1;
        flips = 10;
    }

    for (i = -2; i < runs; i++) {
        /* METODO 1
         * single spin flip
         */
        if (metodo == 1) {

            for (int k = 0; k < flips; k++) {

                for (int j = 0; j < N; j += 2) {
                    dH = 2 * chain[j]*(chain[nnl[j]] + chain[nnr[j]] + chain[nnu[j]] + chain[nnd[j]]);
                    if (dH < 0 || generatore.get_double() < exp(-beta * dH))
                        chain[j] = -chain[j];
                }
                for (int j = 1; j < N; j += 2) {
                    dH = 2 * chain[j]*(chain[nnl[j]] + chain[nnr[j]] + chain[nnu[j]] + chain[nnd[j]]);
                    if (dH < 0 || generatore.get_double() < exp(-beta * dH))
                        chain[j] = -chain[j];
                }

            }
        }

        if (i < 0)
            continue;
        partitions[i].from_configuration(chain.data(), topologia);
    }
}

// J = normale(media=0,std=1)

void ising_entries_jnorm(options opts, int *buffer_sequenze, RandMT &generatore) {
    int L = opts.seq_len;
    int runs = opts.n_seq;
    double beta = opts.beta[0];

    vector<int> flipchain(L);
    vector<int> chain(L);
    vector<int> J(L);
    vector<double> prob(L);

    for (int i = 0; i < L; i++) {
        double r;
		//probabilita di trovare un flip, ovvero -1
        //J gaussiano (positivo)
        //r = generatore.semi_norm();
        //J uniforme [0,1] (positivo)
		//r = generatore.rand();
        //J uniforme [0,0.5] (positivo)
		r = generatore.get_double()/2;
		//J cost
		//r=1;

        prob[i] = exp(-2 * beta * r);
        prob[i] /= (1 + prob[i]);

        // J  +- 1
        //J[i] = 2 * (generatore() > .5) - 1;
        // J positivi
        J[i] = 1;
    }

    for (int i = 0; i < runs; i++) {
        chain[0] = 2 * (generatore.get_double() > .5) - 1;

        for (int k = 0; k < L; k++)
            flipchain[k] = (prob[k] > generatore.get_double()) ? -1 : 1;

        for (int k = 1; k < L; k++)
            chain[k] = flipchain[k] * chain[k - 1] * J[k];

        for (int k = 0; k < L; k++) {
            buffer_sequenze[i * L + k] = chain[k];
        }
    }

}

int main(int argc, char** argv) {


    set_program_options(opts, argc, argv);
///Carica l'opportuna struttura di adiacenza, selezionata da linea di comando
    if (opts.topologia == TORO_2D) 
        topologia = adiacenza_toroidal_lattice(opts.lato);        
    else if (opts.topologia == LINEARE){
        opts.partition_type = LINEAR_PARTITION;
        topologia = adiacenza_simple_line(opts.seq_len);
    }
    else {
        printf("Not supported topology\n");
        exit(1);
    }
    opts.seq_len = topologia.N;
    
    //logarithm lookup table, 6x program speedup
    mylog = new double[3 * opts.seq_len + 10];
    for (int i = 1; i < 3 * opts.seq_len + 10; i++)
        mylog[i] = log(i);
    mylog[0] = 0;

    double media_globale = 0;
    double media_globale_n2 = 0;
    int n_estrazioni = 100;
    int runs = 0;

    
    // <editor-fold defaultstate="collapsed" desc="Sequenze monodimensionali">
    if (opts.topologia == LINEARE)
#pragma omp parallel
    {
        linear_partition *partitions = new linear_partition[opts.n_seq];
        int *buf_sequenze = new int[opts.n_seq * opts.seq_len];
        distance d(opts.seq_len);
        RandMT generatore;

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
                    d.dist(partitions[i], partitions[j]);
                    media_locale += d.dist_shan_r;
                    media_locale_n2 += (d.dist_shan_r)*(d.dist_shan_r);
                }
            }
#pragma omp critical
            {
                media_globale += media_locale;
                media_globale_n2 += media_locale_n2;
                runs += 1;
            }


        }
    }// </editor-fold>

    // <editor-fold defaultstate="collapsed" desc="Reticoli bidimensionali">
    if (opts.topologia == TORO_2D) {
        std::clock_t start = std::clock();
        double time_diff;
        double completed_ratio;
#pragma omp parallel num_threads(opts.threads)
        {
            general_partition *partitions = new general_partition[opts.n_seq];
            distance d(opts.seq_len);
            RandMT generatore;


            for (int L = 0; L < n_estrazioni; L++) {
                // Generazione di un nuovo vettore J_ij random
                // e di opts.n_seq sequenze che hanno quel J 

                double media_locale = 0;
                double media_locale_n2 = 0;

                ising_lattice(opts, generatore, partitions);

                //media delle distanze tra le coppie di sequenze generate
                //#pragma omp parallel for firstprivate(d) schedule(dynamic,10) reduction(+: media_n, media_n2)
                for (int i = 0; i < opts.n_seq; i++) {
                    for (int j = i + 1; j < opts.n_seq; j++) {
                        d(partitions[i], partitions[j]);
                        media_locale += d.dist_shan;
                        media_locale_n2 += (d.dist_shan)*(d.dist_shan);
                    }
                }
#pragma omp critical
                {
                    media_globale += media_locale;
                    media_globale_n2 += media_locale_n2;
                    runs += 1;
                }


#ifdef _OPENMP
                int this_thread = omp_get_thread_num();
                if (this_thread)
                    continue;
                double time_ratio = omp_get_num_threads();
#else
                double time_ratio = 1.0;
#endif
                fprintf(stderr, "\r");
                time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC / time_ratio;
                completed_ratio = (L + 1.0) / n_estrazioni;
                fprintf(stderr, "%.1f%% done, ETA %.0fs    ",
                        completed_ratio * 100, ceil(time_diff * (1 / completed_ratio - 1)));
                fflush(stderr);

            }
        }
        time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC;
        fprintf(stderr, "\r100%% done in %.1f seconds of CPU time\n", time_diff);


    }// </editor-fold>


    double varianza_n;
    int Nd = runs * (opts.n_seq * (opts.n_seq - 1)) / 2;
    media_globale /= Nd;
    media_globale_n2 /= Nd;

    varianza_n = media_globale_n2 - media_globale*media_globale;
	int lunghezza;
	if(opts.topologia == TORO_2D)
		lunghezza=opts.lato;
	else
		lunghezza=opts.seq_len;
    printf("%d %f %f\n", lunghezza, media_globale, varianza_n);
    //fprintf(stderr, "%d %f %f\n", opts.seq_len, media_globale, varianza_n);


    return 0;
}
