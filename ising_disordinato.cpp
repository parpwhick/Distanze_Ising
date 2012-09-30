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
double *myexp;

adj_struct topologia;

void print_array(const int *array, int len, const char *nome) {
    printf("%s {%d", nome, array[0]);
    for (int i = 1; i < len; i++)
        printf(",%d", array[i]);
    printf("}\n\n");
}


int *nnl = 0, *nnu = 0, *nnd = 0, *nnr = 0;

inline int up(int i, int lato, int N) {
    return (i - (i % lato)+ ((i + lato - 1) % lato));
}

inline int down(int i, int lato, int N) {
    return ((i / lato) * lato + ((i + lato + 1) % lato));
}

inline int left(int i, int lato, int N) {
    return (i + N - lato) % N;
}

inline int right(int i, int lato, int N) {
    return (i + N + lato) % N;
}

void ising_lattice(options opts, RandMT &generatore, general_partition *partitions) {
    int N = opts.seq_len;
    int lato = opts.lato;
    int runs = opts.n_seq;
    int iteration = 0;
    double beta = opts.beta[0];
    int dH;
    double prob;
    int flips = 2;

    vector<int> chain(N);

    //inizializzazione array primi vicini
    if (nnl == 0)
#pragma omp critical
    {
        nnl = new int[N];
        nnr = new int[N];
        nnu = new int[N];
        nnd = new int[N];
        for (int i = 0; i < opts.seq_len; i++) {
            nnl[i] = left(i,lato,N);
            nnr[i] = right(i,lato,N);
            nnu[i] = up(i,lato,N);
            nnd[i] = down(i,lato,N);

        }
    }

    vector<int> J(2 * N);
    if (opts.ordered) {
        //ordered, J = 1 case, with phase transition
        if (beta < 0.4)
            prob = 0.5;
        else {
            //genero in base al fit quadratico della probabilita da un esperimento
            prob = log(beta);
            prob = 2.14 - 5.63 * prob - 7.71 * prob*prob;
            prob = exp(prob);
        }
        for (int i = 0; i < N; i++)
            chain[i] = 2 * (generatore.get_double() > prob) - 1;

        J.assign(N,1);
    } else {
        //disordered case, with glass-spin transition at T=0, always kinda hot-temperature
        //Random J_ij equally distributed at -1 and 1
        for (int i = 0; i < 2 * N; i++)
            J[i] = 2 * (generatore.get_int() % 2) - 1;
        for (int i = 0; i < N; i++)
            chain[i] = 2 * (generatore.get_int() % 2) - 1;
    }

    if ((beta < 0.36 || beta > 0.47) && opts.ordered)
        flips = 15;
    else
        flips = 30;

    for (iteration = -2; iteration < runs; iteration++) {
        for (int k = 0; k < flips; k++) {
            for (int s = 0; s < N; s += 2) {
                dH = 0;
                //the link DOWN for site s, is J[s]
                dH += J[s] * chain[nnd[s]];
                //the link RIGHT for site s, is J[s+N]
                dH += J[s + N] * chain[nnr[s]];
                //the link UP, is the link down of site nnu[s], that is J[nnu[s]]
                dH += J[nnu[s]] * chain[nnu[s]];
                //the link LEFT, is the link right of site nnl[s], i.e. J[nnl[s]+N]
                dH += J[nnl[s] + N] * chain[nnl[s]];
                dH *= 2 * chain[s];
                if (dH <= 0 || generatore.get_double() < myexp[dH]) // exp(-beta * dH)
                    chain[s] = -chain[s];
            }
            for (int s = 1; s < N; s += 2) {
                dH = 0;
                //the link DOWN for site s, is J[s]
                dH += J[s] * chain[nnd[s]];
                //the link RIGHT for site s, is J[s+N]
                dH += J[s + N] * chain[nnr[s]];
                //the link UP, is the link down of site nnu[s], that is J[nnu[s]]
                dH += J[nnu[s]] * chain[nnu[s]];
                //the link LEFT, is the link right of site nnl[s], i.e. J[nnl[s]+N]
                dH += J[nnl[s] + N] * chain[nnl[s]];
                dH *= 2 * chain[s];
                if (dH <= 0 || generatore.get_double() < myexp[dH]) // exp(-beta * dH)
                    chain[s] = -chain[s];
            }
        }

        if (iteration < 0)
            continue;
        partitions[iteration].from_configuration(chain.data(), topologia);

    }
}

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
    myexp = new double[100];
    for (int i = 0; i < 100; i++)
        myexp[i] = exp(- opts.beta[0] * i);

    double media_globale = 0;
    double media_globale_n2 = 0;
    double media_rid_globale = 0;
    double media_rid_globale_n2 = 0;
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
            double media_rid_locale = 0;
            double media_rid_locale_n2 = 0;
            ising_entries_jnorm(opts, buf_sequenze, generatore);

            //riempi le partizioni, a partire dalle sequenze date
            for (int i = 0; i < opts.n_seq; i++)
                partitions[i].fill(&buf_sequenze[i * opts.seq_len], opts.seq_len);


            //media delle distanze tra le coppie di sequenze generate
            //#pragma omp parallel for firstprivate(d) schedule(dynamic,10) reduction(+: media_n, media_n2)
            for (int i = 0; i < opts.n_seq; i++) {
                for (int j = i + 1; j < opts.n_seq; j++) {
                    d.dist(partitions[i], partitions[j]);
                    media_locale += d.dist_shan;
                    media_rid_locale += d.dist_shan_r;
                    media_locale_n2 += (d.dist_shan)*(d.dist_shan);
                    media_rid_locale_n2 += (d.dist_shan_r)*(d.dist_shan_r);
                }
            }
#pragma omp critical
            {
                media_globale += media_locale;
                media_globale_n2 += media_locale_n2;
                media_rid_globale += media_rid_locale;
                media_rid_globale_n2 += media_rid_locale_n2;
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
            double media_rid_locale = 0;
            double media_rid_locale_n2 = 0;

                ising_lattice(opts, generatore, partitions);

                //media delle distanze tra le coppie di sequenze generate
                //#pragma omp parallel for firstprivate(d) schedule(dynamic,10) reduction(+: media_n, media_n2)
                for (int i = 0; i < opts.n_seq; i++) {
                    for (int j = i + 1; j < opts.n_seq; j++) {
                        d(partitions[i], partitions[j]);
                         media_locale += d.dist_shan;
                    media_rid_locale += d.dist_shan_r;
                    media_locale_n2 += (d.dist_shan)*(d.dist_shan);
                    media_rid_locale_n2 += (d.dist_shan_r)*(d.dist_shan_r);
                    }
                }
#pragma omp critical
                {
                   media_globale += media_locale;
                media_globale_n2 += media_locale_n2;
                media_rid_globale += media_rid_locale;
                media_rid_globale_n2 += media_rid_locale_n2;
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


    double varianza_n, varianza_r;
    int Nd = runs * (opts.n_seq * (opts.n_seq - 1)) / 2;
    media_globale /= Nd;
    media_globale_n2 /= Nd;
    media_rid_globale /= Nd;
    media_rid_globale_n2 /= Nd;

    varianza_n = media_globale_n2 - media_globale*media_globale;
    varianza_r = media_rid_globale_n2 - media_rid_globale*media_rid_globale;
	int lunghezza;
	if(opts.topologia == TORO_2D)
		lunghezza=opts.lato;
	else
		lunghezza=opts.seq_len;
    printf("%d %f %f %f %f\n", lunghezza, media_globale, varianza_n, media_rid_globale, varianza_r);
    //fprintf(stderr, "%d %f %f\n", opts.seq_len, media_globale, varianza_n);


    return 0;
}
