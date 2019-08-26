/* 
*  Amrith Deepak
*  amrithd@andrew.cmu.edu
*  February 26, 2014
*
*  csim.c - This function takes command line arguments of set bits,
*  block bits, associativity, and tracefile, with which it produces
*  output with the hits, misses, and evictions. It's a cache simulator.
*
*/

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "cachelab.h"
#include <getopt.h>
#include <unistd.h>
#define TRUE 1
#define FALSE 0
#define INVALID -1
#define STRSIZE 100
#define ZEROINIT 0
#define ADRSTRSIZE 15
#define HEXSIZE 16
#define FIRSTELEMENT 0

/* Variable Declarations.
   Cache contains valid, tag, set, and bits.
   Cache contains strct Set* which has struct Line *.
*/

typedef long long addr_t; //signed so that we can store -1

typedef struct Line {
    struct Line *prev;
    addr_t tag;
    struct Line *next;
}Line;

typedef struct Set {
    struct Line *line;
}Set;

typedef struct Cache{
    int s, S;
    int E;
    int b, B;
    int v;
    int hits, misses, evictions;
    Set* sets;
}Cache;

// Allocates the correct amount of memory for cache, and initializes it.
Cache init_cache(int s, int E, int b)
{
    int i, j;
	
    Cache cache = {0};
    cache.s = s;
    cache.S = (int) (pow(2, (double) s));
    cache.E = E;
    cache.b = b;
    cache.B = (int) (pow(2, (double) b));
    cache.sets = (Set*) calloc(cache.S, sizeof(Set*));
    if (cache.sets == NULL)    {
        printf("Calloc returned NULL");
        abort();
    }
    for(i = 0; i < cache.S; i++) {
        Line *current_line = cache.sets[i].line
        = (Line *) calloc(1, sizeof(Line)); //Allocate 1st lines
        if(current_line == NULL) {
            printf("Calloc returned NULL");
            abort();
        }
        current_line->tag = -1;
        for(j = 1; j < E; j++) { //Allocate more lines if required
            current_line->next = (Line *) calloc(1, sizeof(Line));
            if(current_line->next == NULL) {
                printf("Calloc returned NULL");
                abort();
            }
            current_line->next->tag = -1;
            current_line->next->prev = current_line;
            current_line = current_line->next;
        }
    }
    return cache;
}

// Frees elements of cache to prevent memory leaks.
void free_cache(Cache *cache)
{
    int i;
    for(i = 0; i < cache->S; i++){
        Line *next;
        do {
            next = cache->sets[i].line->next;
            free(cache->sets[i].line);
            cache->sets[i].line = next;
        }while(next != NULL);
    }
    free(cache->sets);
}

// Move the accessed line to the begining of the set. Needed to implement LRU.
void move_to_front(Line *line)
{
    Line* prev_line = line->prev;
    while(prev_line != NULL) {
    // there is no data in this simulated cache line,
    // so instead of swaping the lines just swap the tags in them
        addr_t temp = prev_line->tag;
        prev_line->tag = line->tag;
        line->tag = temp;
        prev_line = prev_line->prev;
        line = line->prev;
    }
}

// This function performs the lookup and records a hit miss, or eviction.
void lookup_cache(Cache *cache, addr_t addr)
{
    addr_t tag = addr >> (cache->s + cache->b);
    int i, set_index = (addr >> cache->b) % cache->S;
    Line *curr = cache->sets[set_index].line;
  
    for(i = 0; i < cache->E; i++) {
        if(curr->tag == tag) { //hit
            cache->hits++;
            if(cache->v) {
                printf("hit  (set=%d,line=%d) \t", set_index, i);
            }
            break;
	}
	if(curr->next == NULL) { //done searching all lines, so it's a miss
            cache->misses++;
            if(cache->v) {
                printf("miss (set=%d) \t\t", set_index);
            }
            if(curr->tag == -1) { //use unused line, need not evict
                curr->tag = tag;
            }
            else {  //evict
                cache->evictions++;
                curr->tag = tag; //(LRU) evict last line
                if(cache->v) {
                    printf("evic (set=%d,line=%d) \t\t", set_index, i);
                }
            }
        }
        else {
            curr = curr->next;
        }
    }
    move_to_front(curr);
}

// The file is read here. Getopt is used to process the input. All the functions
// are called.
int main(int argc, char **argv)
{
    FILE *fp;
    int s = FALSE;
    int E = FALSE;
    int b = FALSE;
    char *t;
    int h = FALSE;
    int v = FALSE;
    
    int op_count = 1;
    char str[STRSIZE] = {ZEROINIT};
    Cache cache;
    int curr;
	
    while ((curr = getopt(argc, argv, "hvs:E:b:t:")) != INVALID) {
        switch (curr) {
            case 'h':
                h = TRUE;
                break;
            case 'v':
                v = TRUE;
                break;
            case 's':
                s = atoi(optarg);
                break;
            case 'E':
                E = atoi(optarg);
                break;
            case 'b':
                b = atoi(optarg);
                break;
            case 't':
                t = optarg;
                break;
            case '?':
            default:
                h = TRUE;
        }
    }
	
    if(TRUE == h) {
        printf("Usage: ./csim [-hv] -s <num> -E <num> -b <num> -t <file>\n"
               "Options:\n"
               "  -h         Print this help message.\n"
               "  -v         Optional verbose flag.\n"
               "  -s <num>   Number of set index bits.\n"
               "  -E <num>   Number of lines per set.\n"
               "  -b <num>   Number of block offset bits.\n"
               "  -t <file>  Trace file.\n\n"
               "Examples:\n"
               "  linux>  ./csim -s 4 -E 1 -b 4 -t traces/yi.trace\n"
               "  linux>  ./csim -v -s 8 -E 2 -b 4 -t traces/yi.trace\n");
        return 0;
    }
    
    fp = fopen(t, "r");
    if(fp == NULL) {
        printf("Can't read input\n");
        return INVALID;
    }
    
    cache = init_cache(s, E, b);
    cache.v = v;
    
    while(fgets(str, STRSIZE, fp) != NULL) {
        char opn;
        char *bytes_str;
        char addr_str[ADRSTRSIZE] = {ZEROINIT};
        addr_t addr; //64bit addresses
        if(str[FIRSTELEMENT] != ' ') //ignore I
            continue;
        
        sscanf(str, " %c %s", &opn, addr_str);
        addr = strtoull(addr_str, &bytes_str, HEXSIZE); //hexadecimal to decimal
        
        if(v) {
            printf("%d) %c %s \t", op_count++, opn, addr_str);
        }	
        lookup_cache(&cache, addr);
        if(opn == 'M') {
            lookup_cache(&cache, addr);
        }
        if(v) {
            printf("\n");
        }
    }
    printSummary(cache.hits, cache.misses, cache.evictions);
    fclose(fp);
    free_cache(&cache);
    
    return 0;
}

