/* 
 * File:   strutture.h
 * Author: fake
 *
 * Created on January 13, 2012, 10:51 PM
 */

#ifndef STRUTTURE_H
#define	STRUTTURE_H

#include <string>
#include <stdint.h>

typedef uint32_t u_int32_t;
typedef uint64_t u_int64_t;
#define FORCE_INLINE __attribute__((always_inline))

typedef int32_t label_t;

class basic_partition {
public:
    //number of atoms found
    label_t n;
    //total length of the partition
    label_t N;

    double entropia_shannon;
    double entropia_topologica;


};


class atom {
public:
    label_t size;
    label_t start;
    label_t end;

    atom() {
        size = 0;
        end = 0;
        start = 0;
    }
};

class general_partition : public basic_partition {
private:
    
    void allocate(label_t len);
public:
    label_t **NNB;
    int allocated_n;
    //labels identify generic atoms across the partition
    label_t *labels;
    //nearest neighbors
    label_t *prev_site;

    atom * atomi;
    int lato;

    int dim;

    template <typename T> void from_square_lattice(const T* val, int lato, int dim);
    
    void from_square_lattice(int *reticolo, int lato, int lato2);
    void relabel();

    //void from_atom(int *label, const int which, const int set_to);
    general_partition(int len = 0);
    void print();
    void print_cluster_adjacency();
    atom& find_atom(const atom &atomo1) const {
        return atomi[labels[atomo1.start]];
    }


    ~general_partition();

    class Iterator {
    private:
        label_t _site;
        const label_t *_next;

    public:

        Iterator(int dove, const label_t *vicini) : _site(dove), _next(vicini) {
        };

        int operator*() {
            return (_site);
        }

        bool operator==(Iterator due) {
            return (_site == due._site);
        }

        bool operator!=(Iterator due) {
            return (_site != due._site);
        }

        int operator++() {
            if (_site == _next[_site])
                _site = -1;
            else
                _site = _next[_site];
            return (_site);
        }

        int operator++(int) {
            return (operator++());
        }

    };

    Iterator begin(const int where) const {
        return Iterator(atomi[where].end, prev_site);
    }

    Iterator begin(const atom &where) const {
        return Iterator(where.end, prev_site);
    }

    Iterator end() const {
        return Iterator(-1, prev_site);
    }
};

typedef general_partition::Iterator Iter_t;

#endif
