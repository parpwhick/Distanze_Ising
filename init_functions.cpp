#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "strutture.h"

options opts;


void translate(std::string type,char *testo,int len) ;

void error(const char* message){
    fputs(message,stderr);
    exit(1);
}

void print_help() {
    fprintf(stderr,
            "Usage: distanze [-option1 [arg]] [-option2 [arg]]...\n"
            "\n"
            "Option list for this program:\n"
            "  -random          Turns on random string generation [on]\n"
            "  -file FILENAME   Read sequences from FILENAME [off]\n"
            "  -seqnum N        Limits the number of sequences to N [1000]\n"
            "  -seqlength N     Limits sequence length to N [1000]\n"
            "  -nodistance      Doesn't calculate any distance matrix\n"
            "  -fuzzy N         Set degree of fuzziness in partitioning the sequence [2]\n"
            "  -translate       Translate protein sequence into reduced 5 letter code [off]\n"
            "  -symbols N       Generate random strings with N symbols [2]\n"
            "  -write           Write the distance matrices [off]\n"
            "  -threads N       Use N threads to calculate distance matrix [1]\n"
            "  -seed N          Random number generator seed [37337]\n"
            "  -sorted          Use sorted general partition distance algorithm [auto]\n"
            "  -pmatrix         Use pmatrix general partition distance algorithm [auto]\n"
            "  -standard        Use standard general partition distance algorithm [auto]\n"
            "  -v [-v -v]       Turns on increasingly verbose messages to stderr [off]\n"
            "  -help            Shows this message\n"
            "  -beta N          Set beta to N\n"
            );
}

void set_program_options(options &opts, int argc, char**argv) {
    opts.seq_len = 625;
    opts.lato = 25;
    opts.n_seq = 2550;
    opts.n_symbols = 2;
    opts.from = SEQUENCE | RANDOM;
    opts.seed = 37337;
    opts.translate = false;
    opts.graphics=false;
    opts.verbose = 0;
    opts.write=false;
    opts.distance=true;
    opts.fuzzy=2;
    opts.threads=2;
    opts.alg=AUTO;
    opts.beta=1.0;
    opts.da_calcolare=   SHAN | SHAN_TOP 
                        | RID | RID_TOP 
         //              | GENERAL | GENERAL_TOP 
          //              | GENERAL_RID | GENERAL_RID_TOP
                      //| HAMM                        
            ;
    
    int killswitch=0;

    std::string input;
    if (argc > 1) {
        int read_argvs = 1;
        do {
            input = argv[read_argvs++];
            if (input == "-random") {
                fprintf(stderr, "Specifying random sequence generation\n");
                opts.from |= RANDOM;
                opts.from &= ~FROM_FILE;
            } else if (input == "-file") {
                if (argc - read_argvs < 1)
                    error("Missing filename to read!\n");
                if (argv[read_argvs][0] == '-') 
                    error("Expecting argument, not another option\n");
                
                strncpy(opts.filename, argv[read_argvs++], 255);
                fprintf(stderr, "Reading from filename: %s\n", opts.filename);
                opts.from &= ~RANDOM;
                opts.from |= FROM_FILE;
            } else if (input == "-lattice") {
                fprintf(stderr, "Analysing 2d lattice\n");
                opts.from |= LATTICE;
                opts.from &= ~SEQUENCE;
            } else if (input == "-sequence") {
                fprintf(stderr, "Analysing 1d sequences\n");
                opts.from &= ~LATTICE;
                opts.from |= SEQUENCE;
            } else if (input == "-seqlength") {
                if (argc - read_argvs < 1)
                    error("Need to specify sequence length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");
                
                opts.seq_len = atoi(argv[read_argvs++]);
                fprintf(stderr, "Sequence length limited to %d\n", opts.seq_len);
            } else if (input == "-lato") {
                if (argc - read_argvs < 1)
                    error("Need to specify lattice side length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");
                
                opts.lato = atoi(argv[read_argvs++]);
                fprintf(stderr, "Lattice side set to %d\n", opts.lato);
            } 
            else if (input == "-seqnum") {
                if (argc - read_argvs < 1)
                    error("Missing max number of sequences to read\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.n_seq = atoi(argv[read_argvs++]);
                fprintf(stderr, "Number of sequences limited to: %d\n", opts.n_seq);
            }else if (input == "-fuzzy") {
                if (argc - read_argvs < 1)
                    error("Missing number of partitioning fuzziness\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.fuzzy = atoi(argv[read_argvs++]);
                fprintf(stderr, "Fuzziness degree set to: %d\n", opts.fuzzy);
            } else if (input == "-v") {
                opts.verbose++;
                fprintf(stderr, "Verbosity at %d\n", opts.verbose);
            }  else if (input == "-graphics") {
                opts.graphics = true;
                fprintf(stderr, "Making pretty lattice graphs\n");
            }         
            else if (input == "-translate") {
                opts.translate = true;
                fprintf(stderr, "Simplifying sequence alphabet\n");
            } else if (input == "-sorted") {
                opts.alg=SORTED;
                fprintf(stderr, "Using sorted algorithm\n");
            } else if (input == "-pmatrix") {
                opts.alg=PMATRIX;
                fprintf(stderr, "Using pmatrix algorithm\n");
            } else if (input == "-standard") {
                opts.alg=NORMAL;
                fprintf(stderr, "Using standard algorithm\n");
            } else if (input == "-write") {
                opts.write = true;
                fprintf(stderr, "Writing out the distance matrices\n");
            } else if (input == "-nodistance") {
                opts.distance = false;
                fprintf(stderr, "Not calculating the distance matrix\n");
            } else if (input == "-symbols") {
                if (argc - read_argvs < 1)
                    error("Missing max number of random letters to use\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");
                
                opts.n_symbols = atoi(argv[read_argvs++]);
                fprintf(stderr, "Number of letters limited to: %d\n", opts.n_symbols);
            } else if (input == "-beta") {
                if (argc - read_argvs < 1)
                    error("Missing BETA!!!\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");
                
                opts.beta = atof(argv[read_argvs++]);
            }else if (input == "-seed") {
                if (argc - read_argvs < 1)
                    error("Expecting random number generation seed\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.seed = atoi(argv[read_argvs++]);
                fprintf(stderr, "Seed set to: %d\n", opts.seed);
            } else if (input == "-threads") {
                if (argc - read_argvs < 1)
                    error("Expecting thread number\n");                            
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.threads = atoi(argv[read_argvs++]);
                fprintf(stderr, "Threads limited to: %d\n", opts.threads);
            } else if (input == "-h" || input == "-help" || input == "--help") {
                print_help();
                killswitch=1;
            }            else {
                fprintf(stderr, "Unknown option: %s\n", input.c_str());
                print_help();
                killswitch=1;
            }
        } while (argc - read_argvs > 0);
    }
    else{
        print_help();
        killswitch=1;
    }
    if(killswitch)
        exit(0);
    
    srand(opts.seed);
    
    if(opts.from & LATTICE){
        opts.seq_len = opts.lato * opts.lato;
        opts.da_calcolare &= ~(SHAN | SHAN_TOP | RID | RID_TOP);
    }
    
    fprintf(stderr,"\n");

}