#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>

#include "strutture.h"

double *mylog;
int *colore;

//static rand55 gen;


void general_partition::allocate(label_t len) {
    assert(len != 0);
    N = len;

    try {
        if (!labels){
            labels = new label_t[N];
        }
        if (!prev_site){
            prev_site = new label_t[N];
        }
        if (!atomi)
                atomi = new atom[N];
    }    catch (std::bad_alloc &e) {
        fprintf(stderr, "Error allocating new partition: %s\n", e.what());
        exit(1);
    }
}


general_partition::general_partition(int len) {
    n = 0;
    allocated_n = 0;
    N = len;
    entropia_topologica = 0;
    entropia_shannon = 0;
    prev_site=0;
    atomi=0;
    labels=0;
    NNB=0;

    if (N)
        allocate(N);
}

general_partition::~general_partition(){
    if (N) {
        delete []labels;
        delete []prev_site;
        delete []atomi;
    }    
}


template <typename pointer_t> int findroot(int i,pointer_t *ptr)
{
  if (ptr[i]<0) return i;
  return ptr[i] = findroot(ptr[i],ptr);
}

#define nnu(i) (i - (i % lato)+ ((i+lato-1)%lato))
#define nnd(i) ((i/lato)*lato + ((i+lato+1)%lato))
#define nnl(i) (i+N-lato)%N
#define nnr(i) (i+N+lato)%N

void general_partition::from_square_lattice(int *reticolo, int lato, int lato2){
    label_t s1, s2;
    label_t r1, r2;
    label_t neighbors[2];
    
    N=lato*lato2;
    allocate(N);
    
    for (label_t i = 0; i < N; i++) {
        labels[i] = -1;
    }  
 
    for (s1 = 0; s1 < N; s1++) {
        r1=findroot(s1,labels);
        
        neighbors[0] = (reticolo[s1] == reticolo[nnu(s1)]) ? nnu(s1) : s1;
        neighbors[1] = (reticolo[s1] == reticolo[nnl(s1)]) ? nnl(s1) : s1;
        for (int j = 0; j < 2; j++) {
            s2 = neighbors[j];        
            if (s1 == s2 || s2 < 0)
                continue;
                        
            r2 = findroot(s2, labels);
            // attribution to proper tree root
            if (r1 != r2) {
                if (labels[r1] >= labels[r2]) {
                    labels[r2] += labels[r1];
                    labels[r1] = r2;
                    r1 = r2;
                } else {
                    labels[r1] += labels[r2];
                    labels[r2] = r1;
                }
            }
            //next neightbor
        }
        //next site
    }
    this->relabel();
}

void general_partition::relabel(){
    label_t *new_label=new label_t[N];
    entropia_shannon=0;
    n=0;
    
    // 1-creazione array atomi
    // 2-inizializzazione ogni elemento
    // 3-creazione indice (label atomo) <--> root
    // 4-calcolo entropia a partire dai size nei root
    #define ATOMO atomi[n]
    for(label_t i=0;i<N;i++){
        if(labels[i]<0){
            entropia_shannon+= -labels[i]*mylog[-labels[i]];
            ATOMO.size=-labels[i];
            ATOMO.end=i;
            ATOMO.start=i;
            new_label[i]=n;
            prev_site[i]=i;
            n++;
        } else{
            prev_site[i]=findroot(i,labels);
            new_label[i]=prev_site[i];
        }
    }
    entropia_topologica=mylog[n];
    entropia_shannon= -entropia_shannon/N + mylog[N];    
    // 1-relabeling secondo l'indice dell'atomo, non del sito di appartenenza
    // 2-hashing
    // 3-creazione del collegamento prev_site e atom.end
    for(label_t i=0;i<N;i++){
        int atom_pos=new_label[prev_site[i]];
        labels[i]=atom_pos;
        prev_site[i]=std::min(atomi[atom_pos].end,i);
        atomi[atom_pos].end=i;      
    }
    
    delete []new_label;        
}

using namespace std;
#include <vector>
void general_partition::print_cluster_adjacency(){
    //label_t e' definito in strutture.h come int32
    vector<label_t> riga;
    vector<label_t> colonna;
    int totale=0;
    FILE *vec1 = fopen("vector1.bin", "wb");
    FILE *vec2 = fopen("vector2.bin", "wb");
    //per ogni sito
    for (label_t which = 0; which <N; which++) {
        /* Per ogni sito appartenente al reticolo, recupero le informazioni
           sul cluster di appartenenza (atomo in questa nomenclatura).
           Cio' e' necessario per ottenere tutti i siti (ordinati) del cluster
	   cercato con efficienza massima. */
	const atom &atomo = atomi[labels[which]];
        int quanti=atomo.size;
        riga.reserve(quanti);
        colonna.reserve(quanti);
            
        int sito=0;
        //scorro tutti i siti appartenenti allo stesso atomo, 
        //con l'iteratore ii. Il sito corrispondente e'  *ii
        for (Iter_t ii = this->begin(atomo); ii != this->end(); ii++){
            /* gli elementi nonnulli della matrice di adiacenza A(i,j)
               sono in (which, *ii), salvo i valori delle righe e delle colonne
	       corrispondenti in due vettori */
            
            //salto elemento diagonale
            if(which==*ii)
                continue;
            riga[sito]=which+1;
            colonna[sito]=*ii+1;
            sito+=1;
        }       
        //stampa i vettori cosi costruiti!!!!!!!!!!
        fwrite(&riga[0], sizeof (label_t), sito, vec1);
        fwrite(&colonna[0], sizeof (label_t), sito, vec2);
        totale+=sito;
        
    /* Creazione vettori di adiacenza per matrice sparse in stile Matlab */
    /* La funzione per caricare i dati cosi creati e':
    -------------------------------
    function adiacenza=load_sierpinski()
	indici_riga=fread(fopen('vector1.bin','r'),inf,'int32');
	indici_colonna=fread(fopen('vector2.bin','r'),inf,'int32');
	N=max(max(indici_riga),max(indici_colonna));
	adiacenza=sparse(indici_riga,indici_colonna,1,N,N);
    end
    ******************************/
    }   
    fclose(vec1);
    fclose(vec2);
    printf("Elementi nonnulli della matrice di adiacenza: %d\n",totale);
}

int main(int argc, char **argv){
    int x=100,y=100; //dimensioni a caso
    int N= x*y; //numero elementi
    int *reticolo=new int[N];
    FILE *input=fopen("reticolo.bin","rb");
    if(input==0){
	printf("File non trovato\n");
	return(1);
    }
    fread(reticolo,sizeof(int),N,input);
    
    //preallocazione dei logaritmi, necessaria
    mylog = new double[N+10];
    for (int i = 1; i <  N+10;  i++)
        mylog[i] = log(i);
    mylog[0]=0; 
    
    general_partition partizione(N); //prealloca una partizione per N elementi
    partizione.from_square_lattice(reticolo,x,y); 
    //proprieta utili, popolate durante la funzione from_square_...
    //partizione.n - numero di clusters formati
    //partizione.entropia_shannon
    //partizione.labels[] e' l'array delle etichette dei cluster - puoi usarlo per stamparlo e disegnarlo ecc
    printf("Partizione con %d clusters, entropia: %g\n",partizione.n,partizione.entropia_shannon);
    //crea vettori di adiacenza 
    partizione.print_cluster_adjacency();
    printf("Completato con successo\n");
}
