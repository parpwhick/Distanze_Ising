#include <cstdlib>
#include <cstdio>
#include <string>
#include <cassert>


#include "strutture.h"


extern options opts;
extern double *mylog;
extern int *colore;


FORCE_INLINE void DJBHASH_STEP(u_int32_t &hash, u_int32_t value) {
    hash = (hash << 5) + hash + (0xff000000 & value);
    hash = (hash << 5) + hash + (0x00ff0000 & value);
    hash = (hash << 5) + hash + (0x0000ff00 & value);
    hash = (hash << 5) + hash + (0x000000ff & value);
}

template <typename T>
linear_partition::linear_partition(const T* seq, int len) {
    this->fill(seq,len);
}

template <typename T>
void linear_partition::fill(const T* seq, int len) {
    int i, j;
    
    //Total length of the partition is equal to the sequence
    N = len;
        
    if(!binary)
		binary = new int[len];
    //atom_positions = new int[len];

    //first one always start an atom
    binary[0] = 1;

    //checking for a different symbol from the one before
    //when that happens => new atom!
    for (i = 1; i < len; i++)
        binary[i] = seq[i] != seq[i - 1];

    //now going back and writing the atom positions
    j=0;
    for (i = 0; i < len; i++)
        if (binary[i]) {
            //atom_positions[j++] = i + 1;
            j++;
        }
    //number of atoms found - index++ of last the atom in the array (zero addressing)
    n = j;
   entropia_topologica=mylog[n];
    entropia_shannon=entropy_binary_partition(binary,N);
}
template void linear_partition::fill(const char *, int);
template void linear_partition::fill(const int *, int);


/************************************************
 ************************************************
 *            PARTIZIONI COMPLESSE              *
 ************************************************
 ************************************************
 */



void general_partition::allocate(int len){
    assert(len !=0);
    N=len;
    
    if (!labels){
        labels = new int[N];
		
	}
    if (!prev_site)
        prev_site = new int[N];
//    if (!atomi)
//        atomi = new atom[N];
}

general_partition::general_partition(int len) {
    n = 0;
    N = len;
    entropia_topologica = 0;
    entropia_shannon = 0;

    if (N)
        allocate(N);
    else {
        labels = 0;
        prev_site = 0;
        atomi = 0;
        NNB = 0;
    }    
}

general_partition::~general_partition(){
    if (N) {
        delete []labels;
        delete []prev_site;
        //delete []atomi;
    }    
    if(NNB)
        delete []NNB;
}


template <typename T>
general_partition::general_partition(const T* seq, int len) {
    labels=0;
    N=0;
    prev_site=0;
    this->from_linear_sequence(seq,len);
}


inline int compare (const void * a, const void * b){
  return ( *(int*)a - *(int*)b );
}

void general_partition::sort_entropy(){
    int i;
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;
    
    int *temp=new int32_t[N];
    allocate(N);
    
    for(int i=0; i<N;i++)
        temp[i]=labels[i];
    
    qsort(temp,N,sizeof(temp[0]),compare);
    
    //the first position always starts an atom
    begin = 0;
    label_count=1;
    int old_val=temp[0];
        
    for (i = 1; i < N; i++) {      
        //whenever we find a new atom
        if (temp[i]!=old_val) {
            //a new atom starts
            
            //the closed (old)atom's length is calculated
            mu = i - begin;           
            //the new one is ready to go
            label_count++;
            begin = i;
            //cache the new label to check
            old_val=temp[i];
            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];   
        }
    }
    //the last one, so it's not left hanging
    mu = N - begin;
    H += mu * mylog[mu];
    
    //normalize the result
    H = -H / N + mylog[N];
    
    entropia_topologica=mylog[label_count];
    n=label_count;
    entropia_shannon=H;
    
    delete []temp;
}

bool is_similar_atom(Iter_t from1, Iter_t from2, Iter_t to, int tol=10){               
    int differenze=0;
    
    while (true) {
        // Se il numero delle differenze e' elevato, escludo l'uguaglianza
        if (differenze > tol) {
            return (false);        
        // Se ho raggiunto la fine di entrambi gli atomi senza troppe differenze
        // allora accetto l'uguaglianza
        } else if ((from1 == to) && (from2 == to)) {
            return (true);         
        // Nel caso di atomi con disparita' numerica, continuo a scorrere il piu
        // lungo, per vedere se accumulo abbastanza differenze
        }else if ((from1 == to) != (from2 == to)) {
            differenze++;
            from1++;
            from2++;         
        // Se ho elementi diversi nei due atomi, scorro per riallinearli
        } else if (*from1>*from2) {
            differenze++;
            from1++;
        } else if (*from2>*from1) {
            differenze++;
            from2++;
        // I due atomi hanno lo stesso elemento - scorro in avanti entrambi
        } else {
            from1++;
            from2++;
        }
    }
}

void general_partition::reduce(const general_partition &p1, const general_partition &p2){
    //inizializzazioni
    int fattori_indipendenti=0;
    lato=p1.lato;
    allocate(p1.N);
  
    for(int i=0;i<N;i++)
        labels[i]=1;
    
    //per ogni atomo
    for (int which = 0; which < p1.n; which++) {
        // Considero l'atomo n-esimo del primo
        // Trovo l'atomo che corrisponde nella seconda partizione
        //  attraverso il primo sito in comune
        const atom &atomo1 = p1.atomi[which];
        const atom &atomo2 = p2.find_atom(atomo1);
        //const atom &atomo2 = p2.atomi[p2.labels[atomo1.start]];       
        
        
        // Se gli atomi sono uguali, l'intersezione delle partizioni dicotomiche
        // non e' banale => salto
        if (atomo1 == atomo2)
            continue;
        
        // Creazione degli iteratori, per ottenere tutti i siti in atomo1
        Iter_t ii = p1.begin(atomo1);
        Iter_t end = p1.end();
        
        
        
        // uguaglianza "fuzzy" tra atomi, a meno di 'tol' siti
        // similmente, se sono 'uguali', salto l'atomo nella partizione risultante
//        if(is_similar_atom(ii,p2.begin(atomo2),end))
//            continue;
        
        // altrimenti interseca il fattore dicotomico con i precedenti
        fattori_indipendenti++;
                
        
        // faccio il prodotto rapido - moltiplico i siti di atomo1 per numero
        for (; ii != end; ii++) 
            labels[*ii] *= fattori_indipendenti + 1;
        
//        for (int i = p1.atomi[which].end;; i = p1.prev_site[i]) {
//            labels[i] *= fattori_indipendenti + 1;
//            if (i == p1.prev_site[i])
//                break;
//        }
           
    }   
    //printf("%.1f%% (%d/%d) accettati\n",(fattori_indipendenti+0.0)/p1.n*100.0, fattori_indipendenti,p1.n);
    // Calcolo dell'entropia della partizione cosi generata
    // !!ACHTUNG!! tutti i labels verranno mischiati - TEMP ARRAY NECESSARIO
    this->sort_entropy();
}



int findroot(int i,int *ptr)
{
  if (ptr[i]<0) return i;
  return ptr[i] = findroot(ptr[i],ptr);
}


#define nnu (s1 - (s1 % lato)+ ((s1+lato-1)%lato))
#define nnl (s1+N-lato)%N

template <typename T> 
void general_partition::from_square_lattice(const T *reticolo, int lato,int){
    int s1, s2;
    int i;
    dim=2;
    N=lato*lato;
    allocate(N);
    int neighbors[2];
    
    //presetting labels for this partition to 0
    for (i = 0; i < N; i++) {
        labels[i] = -1;
    }    
    /* FASE DI PERCOLAZIONE
     * E LABELLING
     */
 
    for (s1 = 0; s1 < N; s1++) {
        int r2;
        int r1=findroot(s1,labels);

        
        neighbors[0] = (reticolo[s1] == reticolo[nnu]) ? nnu : s1;
        neighbors[1] = (reticolo[s1] == reticolo[nnl]) ? nnl : s1;
        for (int j = 0; j < dim; j++) {
            s2 = neighbors[j];        
            if (s1 == s2 || s2 < 0)
                continue;
            
            r2 = findroot(s2, labels);
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
        }
    }
    entropia_shannon=0;
    n=0;
    for(i=0;i<N;i++){
        if(labels[i]<0){
            entropia_shannon+= -labels[i]*mylog[-labels[i]];
            prev_site[i]=i;
            n++;
        } else{
            prev_site[i]=findroot(i,labels);
        }
    }
    std::swap(labels,prev_site);
    
    
    
    entropia_topologica=mylog[n];
    entropia_shannon= -entropia_shannon/N+mylog[N];
          
}
template void general_partition::from_square_lattice(const int*, int, int);




void general_partition::linear_intersection(const general_partition &p1, const general_partition &p2){
    int *vicinato[2 * p1.dim];
    lato=p1.lato;
    allocate(p1.N);
    
    
    // Per fare l'intersezione unisco l'insieme dei vicini delle due partizioni
    // ottenendo sostanzialmente una struttura di adiacenza 2x dimensionale
    for (int i = 0; i < p1.dim; i++) {
        vicinato[2 * i] = p1.NNB[i];
        vicinato[2 * i + 1] = p2.NNB[i];
    }
    // Calcolo a partire dall'insieme dei vicini
    from_nnb(vicinato, 2 * p1.dim);

}


void general_partition::from_nnb(int **neighbors, int dimensione){
    int s1, s2;
    int i;
    dim=dimensione;
    
    //presetting labels for this partition to 0
    for (i = 0; i < N; i++) {
        labels[i] = -1;
    }    
    /* FASE DI PERCOLAZIONE
     * E LABELLING
     */
 //   print_array(neighbors[0],N,"neighbors: ");
    for (s1 = 0; s1 < N; s1++) {
        int r2;
        int r1=findroot(s1,labels);
        
        //printf("%d ha come vicini: %d e %d\n",s1,neighbors[0][s1],neighbors[1][s1]);
        for (int j = 0; j < dim; j++) {
            s2=neighbors[j][s1];
                if(s1==s2 || s2 < 0)
                    continue;
		r2=findroot(s2,labels);
//                printf("%d--->%d(%d)\n",s1,s2,r2);
		if(r1!=r2){
			if(labels[r1]>=labels[r2]){ //VERSIONE GIUSTA!
//                           printf("%d=>%d\n",r1,r2);
				labels[r2]+=labels[r1];
				labels[r1]=r2;
				r1=r2;
			}
			else{
//                          printf("%d=>%d\n",r2,r1);
				labels[r1]+=labels[r2];
				labels[r2]=r1;
			}
		}
	}
    }
    int *new_site_label=new int[N];
    entropia_shannon=0;
    int label_count=0;
    
 //   print_array(labels,N,"raw parti: ");
    // 1-creazione array atomi
    // 2-inizializzazione ogni elemento
    // 3-creazione indice (label atomo) <--> root
    // 4-calcolo entropia a partire dai size nei root
    #define ATOMO atomi[label_count]
    for(i=0;i<N;i++){
        if(labels[i]<0){
            entropia_shannon+= -labels[i]*mylog[-labels[i]];
            ATOMO.size=-labels[i];
            ATOMO.hash=5381;
            ATOMO.end=i;
            ATOMO.start=i;
            new_site_label[i]=label_count;
            prev_site[i]=i;
            label_count++;
        } else{
            prev_site[i]=findroot(i,labels);
            new_site_label[i]=prev_site[i];
        }
    }
 //   print_array(prev_site,N,"good labels: ");
    n=label_count;
    entropia_topologica=mylog[n];
    entropia_shannon= -entropia_shannon/N+mylog[N];
    
    // 1-relabeling secondo l'indice dell'atomo, non del sito di appartenenza
    // 2-hashing
    // 3-creazione del collegamento prev_site e atom.end
    #undef HASHVAR
    #define HASHVAR atomi[atom_pos].hash
    for(i=0;i<N;i++){
        int atom_pos=new_site_label[prev_site[i]];
        labels[i]=atom_pos;
        DJBHASH_STEP(HASHVAR,i);
        prev_site[i]=std::min(atomi[atom_pos].end,i);
        //next_site[prev_site[i]]=i;
        atomi[atom_pos].end=i;
      
    }
 //   print_array(labels,N,"best labels: ");
    
    delete []new_site_label;        
}

