//
// Created by nilsonneto on 1/14/16.
//

#ifndef PARALLELCLUSTERING_EM_H
#define PARALLELCLUSTERING_EM_H

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "em"
#define DESCRIPTION "expectation maximization example program"
#define VERSION     "version 1.2 (05.01.2002)         " \
                    "(c) 2000-2002   Christian Borgelt"


/*----------------------------------------------------------------------
  Error codes
----------------------------------------------------------------------*/
#define E_NONE        0         /* no error */
#define E_NOMEM     (-1)        /* not enough memory */
#define E_FOPEN     (-2)        /* cannot open file */
#define E_FREAD     (-3)        /* read error on file */
#define E_FWRITE    (-4)        /* write error on file */
#define E_OPTION    (-5)        /* unknown option */
#define E_OPTARG    (-6)        /* missing option argument */
#define E_ARGCNT    (-7)        /* too few/many arguments */
#define E_TABLE     (-8)        /* unknown table */
#define E_ALG       (-9)        /* unknown update algorithm */
#define E_UNKNOWN  (-10)        /* unknown error */


/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- a (reduced) data tuple --- */
    int    A, B, C;               /* values of attributes A, B and C */
    double cnt;                   /* number of occurences */
    double n[2];                  /* frequencies of the values of H */
} TUPLE;                        /* (data tuple) */

typedef struct {                /* --- a probability --- */
    double prob;                  /* current probability */
    double emp;                   /* probability after EM step */
    double prv;                   /* previous EM step */
    double chg;                   /* change value for resilient EM */
} PROB;                         /* (probability) */

typedef struct {                /* --- simple Bayesian network --- */
    PROB pH[2];                   /* prob. dist. of hidden attribute H */
    PROB pAH[2][2];               /* cond. prob. dist. of attribute A */
    PROB pBH[2][2];               /* cond. prob. dist. of attribute B */
    PROB pCH[2][2];               /* cond. prob. dist. of attribute C */
} BAYESNET;                     /* (Bayesian network) */


/*----------------------------------------------------------------------
  Random Number Functions
----------------------------------------------------------------------*/
#ifdef DRAND48                  /* if library for drand48() available */
extern void   srand48 (long seed);
extern double drand48 (void);   /* use drand48 functions */
#define dseed(s) srand48((long)(s))
#define drand    drand48

#else                           /* if only standard rand() available */
#define dseed(s) srand((unsigned)(s))
static double drand (void)
{ return rand()/(RAND_MAX +1.0); }
#endif


#endif //PARALLELCLUSTERING_EM_H
