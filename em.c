/*
 * Created by nilsonneto on 1/14/16.
 *
 * Expectation-Maximization Algorithm
 * How to tackle:
 *
 *
 *
 * Author  : Christian Borgelt
 * History : 24.02.2000 file created
 *           30.06.2001 EM with resilient step width added
 *           05.01.2002 option evaluation added, output improved
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "em.h"


/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
static const char *errmsgs[] = {
        /* E_NONE      0 */  "no error\n",
        /* E_NOMEM    -1 */  "not enough memory\n",
        /* E_FOPEN    -2 */  "cannot open file %s\n",
        /* E_FREAD    -3 */  "read error on file %s\n",
        /* E_FWRITE   -4 */  "write error on file %s\n",
        /* E_OPTION   -5 */  "unknown option -%c\n",
        /* E_OPTARG   -6 */  "missing option argument\n",
        /* E_ARGCNT   -7 */  "wrong number of arguments\n",
        /* E_TABLE    -8 */  "unknown table %d\n",
        /* E_ALG      -9 */  "unknown update algorithm %d\n",
        /* E_UNKNOWN -10 */  "unknown error\n"
};                              /* error message texts */

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
const char *prgname = NULL;     /* program name for error messages */

/* --- data tables --- */
TUPLE table0[8] = { { 0, 0, 0, 14 },
                    { 0, 0, 1, 11 },
                    { 0, 1, 0, 20 },
                    { 0, 1, 1, 20 },
                    { 1, 0, 0,  5 },
                    { 1, 0, 1,  5 },
                    { 1, 1, 0, 11 },
                    { 1, 1, 1, 14 } };
TUPLE table1[8] = { { 0, 0, 0, 38 },
                    { 0, 0, 1, 18 },
                    { 0, 1, 0, 44 },
                    { 0, 1, 1, 36 },
                    { 1, 0, 0, 11 },
                    { 1, 0, 1,  9 },
                    { 1, 1, 0, 17 },
                    { 1, 1, 1, 27 } };

/* --- Bayesian network --- */
BAYESNET bn =   {   { { 0.3 }, { 0.7 } },   /* pH */
                    { { { 0.4 }, { 0.6 } },   /* pAH */
                               { { 0.6 }, { 0.4 } } },
                    { { { 0.7 }, { 0.3 } },   /* pBH */
                               { { 0.8 }, { 0.2 } } },
                    { { { 0.2 }, { 0.8 } },   /* pCH */
                               { { 0.5 }, { 0.5 } } } };

/* --- parameters --- */
double minprob  = 0.01;         /* minimal probability */
double maxprob  = 0.99;         /* maximal probability */
double moment   = 0.95;         /* momentum coefficient */
double growth   = 1.2;          /* growth factor */
double shrink   = 0.5;          /* shrink factor */
double initchg  = 0.1;          /* initial change value */


/*----------------------------------------------------------------------
  Output Functions
----------------------------------------------------------------------*/

static void table (TUPLE *tab)
{                               /* --- show data table */
    int   i;                      /* loop variable */
    TUPLE *tpl;                   /* to traverse the table */

    printf("A B C   cnt   n[0]      n[1]\n");
    for (tpl = tab, i = 8; --i >= 0; tpl++) {
        printf("%d %d %d   ",   tpl->A, tpl->B, tpl->C);
        printf("%2.0f: ",       tpl->cnt);
        printf("%9.6f %9.6f\n", tpl->n[0], tpl->n[1]);
    }                             /* traverse the tuples and */
    printf("\n");                 /* show the expected frequencies */
}  /* table() */

/*--------------------------------------------------------------------*/

static void probs (BAYESNET *bn)
{                               /* --- show probabilities */
    int i, k;                     /* loop variables */

    printf("                   A = 0    A = 1"
                   "       B = 0    B = 1"
                   "       C = 0    C = 1\n");
    for (i = 0; i < 2; i++) {     /* traverse the values of H */
        printf("H = %d: %.6f    ", i, bn->pH[i].prob);
        for (k = 0; k < 2; k++)  printf("%.6f ", bn->pAH[i][k].prob);
        printf("   ");
        for (k = 0; k < 2; k++)  printf("%.6f ", bn->pBH[i][k].prob);
        printf("   ");
        for (k = 0; k < 2; k++)  printf("%.6f ", bn->pCH[i][k].prob);
        printf("\n");               /* print probs. of values of A */
    }
    printf("\n");                 /* add a blank line */
}  /* probs() */

/*----------------------------------------------------------------------
  Expectation Maximization Functions
----------------------------------------------------------------------*/

static void expect (BAYESNET *bn, TUPLE *tab)
{                               /* --- expectation step */
    int    i;                     /* loop variable */
    TUPLE  *tpl;                  /* to traverse the table */
    double s;                     /* sum for normalization */

    for (tpl = tab, i = 8; --i >= 0; tpl++) {
        s  = tpl->n[0] = bn->pH[0].prob * bn->pAH[0][tpl->A].prob
                         * bn->pBH[0][tpl->B].prob
                         * bn->pCH[0][tpl->C].prob;
        s += tpl->n[1] = bn->pH[1].prob * bn->pAH[1][tpl->A].prob
                         * bn->pBH[1][tpl->B].prob
                         * bn->pCH[1][tpl->C].prob;
        tpl->n[0] *= s = tpl->cnt /s;
        tpl->n[1] *= s;             /* compute the expected absolute */
    }                             /* frequencies of the values of H */
}  /* expect() */

/*--------------------------------------------------------------------*/

static void maximize (BAYESNET *bn, TUPLE *tab)
{                               /* --- maximization step */
    int    i, k;                  /* loop variables */
    TUPLE  *tpl;                  /* to traverse the table */
    double s;                     /* sum/factor for normalization */

    for (tpl = tab, i = 8; --i >= 0; tpl++) {
        bn->pH[0].emp          += tpl->n[0];
        bn->pH[1].emp          += tpl->n[1];
        bn->pAH[0][tpl->A].emp += tpl->n[0];
        bn->pAH[1][tpl->A].emp += tpl->n[1];
        bn->pBH[0][tpl->B].emp += tpl->n[0];
        bn->pBH[1][tpl->B].emp += tpl->n[1];
        bn->pCH[0][tpl->C].emp += tpl->n[0];
        bn->pCH[1][tpl->C].emp += tpl->n[1];
    }                             /* compute sufficient statistics */
    for (i = 2; --i >= 0; ) {     /* traverse the values of H */
        s = 1 /bn->pH[i].emp;       /* compute normalization factor */
        for (k = 2; --k >= 0; ) {   /* traverse the values of A, B, C */
            bn->pAH[i][k].emp *= s;   /* compute conditional probabilities */
            bn->pBH[i][k].emp *= s;
            bn->pCH[i][k].emp *= s;
        }
    }
    s = 0;                        /* normalize distribution on H */
    for (i = 2; --i >= 0; ) s += bn->pH[i].emp;
    for (i = 2; --i >= 0; ) bn->pH[i].emp /= s;
}  /* maximize() */

/*----------------------------------------------------------------------
  Update Functions
----------------------------------------------------------------------*/

static double standard (PROB *p)
{                               /* --- standard expect. maximization */
    double t;                     /* temporary buffer */

    t = fabs(p->emp -p->prob);    /* compute the absolute change */
    p->prob = p->emp;             /* set the new probability and */
    p->emp  = 0;                  /* clear the var. for the next step */
    return t;                     /* return the absolute change */
}  /* standard() */

/*--------------------------------------------------------------------*/

static double momentum (PROB *p)
{                               /* --- EM with momentum term */
    double t;                     /* temporary buffer */

    p->chg = p->chg *moment +(p->emp -p->prob);
    p->emp = 0;                   /* compute the next change */
    if      (p->chg < 0) {        /* if the probability is decreased, */
        t = minprob -p->prob;       /* clamp at lower bound */
        if (p->chg > t) p->prob += t = p->chg;
        else            p->prob  = minprob; }
    else if (p->chg > 0) {        /* if the probability is increased, */
        t = maxprob -p->prob;       /* clamp at upper bound */
        if (p->chg < t) p->prob += t = p->chg;
        else            p->prob  = maxprob; }
    else t = 0;                   /* if the probability is unchanged */
    return fabs(t);               /* return the value of the change */
}  /* momentum() */

/*--------------------------------------------------------------------*/

static double selfadapt (PROB *p)
{                               /* --- self adaptive EM */
    double s, t, c;               /* temporary buffers */

    s = p->emp -p->prob; p->emp = 0;
    if      (s < 0) t = -p->prv;
    else if (s > 0) t =  p->prv;
    else            t =  0;
    if      (t < 0) {
        p->chg *= shrink;
        if (p->chg <  1) p->chg =  1;
        p->prv = 0; }
    else if (t > 0) {
        p->chg *= growth;
        if (p->chg > 64) p->chg = 64;
        p->prv = s; }
    else
        p->prv = s;
    c = s *p->chg;
    if      (c < 0) {
        t = minprob -p->prob;
        if (c > t) p->prob += t = c;
        else       p->prob  = minprob; }
    else if (c > 0) {
        t = maxprob -p->prob;
        if (c < t) p->prob += t = c;
        else       p->prob  = maxprob;
    }
    return t;
}  /* selfadapt() */

/*--------------------------------------------------------------------*/

static double resilient (PROB *p)
{                               /* --- resilient expect. maximization */
    double s, t;                  /* temporary buffers */

    s = p->emp -p->prob; p->emp = 0;
    if      (s < 0) t = -p->prv;
    else if (s > 0) t =  p->prv;
    else            t =  0;
    if      (t < 0) { p->chg *= shrink; p->prv = 0; }
    else if (t > 0) { p->chg *= growth; p->prv = s; }
    else            {                   p->prv = s; }
    if      (s < 0) {
        t = p->prob -minprob;
        if (p->chg < t) p->prob -= t = p->chg;
        else            p->prob  = minprob; }
    else if (s > 0) {
        t = maxprob -p->prob;
        if (p->chg < t) p->prob += t = p->chg;
        else            p->prob  = maxprob;
    }
    return t;
}  /* resilient() */

/*--------------------------------------------------------------------*/

static double update (BAYESNET *bn, double updfn (PROB *p))
{                               /* --- update Byesian network */
    int    i, k;                  /* loop variables */
    double max = 0, chg;          /* (maximal) change of a probability */

    for (i = 2; --i >= 0; ) {     /* traverse the values of attribute H */
        chg = updfn(bn->pH +i); if (chg > max) max = chg;
        for (k = 2; --k >= 0; ) {   /* traverse the attribute values */
            chg = updfn(bn->pAH[i] +k); if (chg > max) max = chg;
            chg = updfn(bn->pBH[i] +k); if (chg > max) max = chg;
            chg = updfn(bn->pCH[i] +k); if (chg > max) max = chg;
        }                           /* update the cond. probabilities */
    }                             /* for all attributes and return the */
    return max;                   /* maximal change of a probability */
}  /* update() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/

static void error (int code, ...)
{                               /* --- print an error message */
    va_list    args;              /* list of variable arguments */
    const char *msg;              /* error message */

    assert(prgname);              /* check the program name */
    if (code < E_UNKNOWN) code = E_UNKNOWN;
    if (code < 0) {               /* if to report an error, */
        msg = errmsgs[-code];       /* get the error message */
        if (!msg) msg = errmsgs[-E_UNKNOWN];
        fprintf(stderr, "\n%s: ", prgname);
        va_start(args, code);       /* get variable arguments */
        vfprintf(stderr, msg, args);/* print error message */
        va_end(args);               /* end argument evaluation */
    }
    exit(code);                   /* abort the program */
}  /* error() */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- main function */
    int    i, k = 0;              /* loop variables, buffers */
    char   *s;                    /* to traverse the options */
    char   **optarg = NULL;       /* option argument */
    int    tabno    = 0;          /* number of table to process */
    int    alg      = 0;          /* update algorithm identifier */
    int    verb     = 0;          /* counter for verbose output */
    double epsilon  = 1e-6;       /* change limit for termination */
    int    maxcnt   = -1;         /* maximal number of steps */
    int    random   = 0;          /* flag for random initialization */
    long   seed     = time(NULL); /* seed for random number generator */
    double sum;                   /* sum of probabilities */
    double chg, range;            /* change/range of probabilities */
    int    cnt;                   /* step counter */
    TUPLE  *tab = NULL;           /* table to process */
    double (*updfn)(PROB *p) = 0; /* probability update function */


    prgname = argv[0];            /* get program name for error msgs. */

    /* --- print usage message --- */
    if (argc > 1) {               /* if arguments are given */
        fprintf(stderr, "%s - %s\n", argv[0], DESCRIPTION);
        fprintf(stderr, VERSION); } /* print a startup message */
    else {                        /* if no arguments given */
        printf("usage: %s [options] tabno\n", argv[0]);
        printf("%s\n", DESCRIPTION);
        printf("%s\n", VERSION);
        printf("-e#      epsilon for termination (default: %g)\n", epsilon);
        printf("-t#      maximal number of steps (default: no limit)\n");
        printf("-v#      verbose output (show state every # iterations)\n");
        printf("-a#      update algorithm to use (default: %d)\n", alg);
        printf("         0: standard expectation maximization\n");
        printf("         1: EM with a momentum term\n");
        printf("         2: EM with self-adaptive factor\n");
        printf("         3: EM with resilient step width\n");
        printf("-p#      minimal probability  (default: %g)\n", minprob);
        printf("-P#      maximal probability  (default: %g)\n", maxprob);
        printf("-m#      momentum coefficient (default: %g)\n", moment);
        printf("-g#      growth factor for algorithms 2 and 3 "
                       "(default: %g)\n", growth);
        printf("-k#      shrink factor for algorithms 2 and 3 "
                       "(default: %g)\n", shrink);
        printf("-i#      initial step width for resilient EM  "
                       "(default: %g)\n", initchg);
        printf("-r       initialize probabilities randomly\n");
        printf("-s#      seed for random number generator "
                       "(default: time)\n");
        printf("tabno    identifier of table to process "
                       "(must be 0 or 1)\n");
        return 0;                   /* print a usage message */
    }                             /* and abort the program */

    /* --- evaluate arguments --- */
    for (i = 1; i < argc; i++) {  /* traverse arguments */
        s = argv[i];                /* get option argument */
        if (optarg) { *optarg = s; optarg = NULL; continue; }
        if ((*s == '-') && *++s) {  /* -- if argument is an option */
            while (*s) {              /* traverse options */
                switch (*s++) {         /* evaluate switches */
                    case 'e': epsilon =      strtod(s, &s);    break;
                    case 't': maxcnt  = (int)strtol(s, &s, 0); break;
                    case 'a': alg     = (int)strtol(s, &s, 0); break;
                    case 'v': verb    = (int)strtol(s, &s, 0); break;
                    case 'p': minprob =      strtod(s, &s);    break;
                    case 'P': maxprob =      strtod(s, &s);    break;
                    case 'm': moment  =      strtod(s, &s);    break;
                    case 'g': growth  =      strtod(s, &s);    break;
                    case 'k': shrink  =      strtod(s, &s);    break;
                    case 'i': initchg =      strtod(s, &s);    break;
                    case 'r': random  = 1;                     break;
                    case 's': seed    =      strtol(s, &s, 0); break;
                    default : error(E_OPTION, *--s);           break;
                }                       /* set option variables */
                if (optarg && *s) { *optarg = s; optarg = NULL; break; }
            } }                       /* get option argument */
        else {                      /* -- if argument is no option */
            switch (k++) {            /* evaluate non-options */
                case  0: tabno = atoi(s); break;
                default: error(E_ARGCNT); break;
            }                         /* note filenames */
        }
    }
    if (optarg) error(E_OPTARG);  /* if missing option argument */
    if (k != 1) error(E_ARGCNT);  /* if too few arguments given */
    switch (tabno) {              /* evaluate the table identifier */
        case  0: tab = table0;      break;
        case  1: tab = table1;      break;
        default: error(E_TABLE);    break;
    }                             /* set the table to work on */
    switch (alg) {                /* evaluate the algorithm identifier */
        case  0: updfn = standard;  break;
        case  1: updfn = momentum;  break;
        case  2: updfn = selfadapt; break;
        case  3: updfn = resilient; break;
        default: error(E_ALG);      break;
    }                             /* note the update function */
    epsilon = fabs(epsilon);      /* remove a possible sign */
    fprintf(stderr, "\n\n");      /* terminate the startup message */

    /* --- initialize probabilities --- */
    if (random) {                 /* if to initialize randomly */
        dseed(seed);                /* init. random number generator */
        range = maxprob -minprob;   /* get range of acceptable values */
        for (i = 2; --i >= 0; ) {   /* traverse the values of H */
            bn.pH[i].prob = range *drand() +minprob;
            for (k = 2; --k >= 0; ) { /* traverse the values of A and B */
                bn.pAH[i][k].prob = range *drand() +minprob;
                bn.pBH[i][k].prob = range *drand() +minprob;
                bn.pCH[i][k].prob = range *drand() +minprob;
            }                         /* choose random rel. probabilities */
        }                           /* within the acceptable range */
    }

    /* --- normalize probabilities --- */
    sum = 0;                      /* normalize distribution on H */
    for (i = 2; --i >= 0; ) sum += bn.pH[i].prob;
    for (i = 2; --i >= 0; ) bn.pH[i].prob /= sum;
    for (i = 2; --i >= 0; ) {     /* traverse the values of H */
        sum = 0;                    /* normalize cond. dists. on A */
        for (k = 2; --k >= 0; ) sum += bn.pAH[i][k].prob;
        for (k = 2; --k >= 0; ) bn.pAH[i][k].prob /= sum;
        sum = 0;                    /* normalize cond. dists. on B */
        for (k = 2; --k >= 0; ) sum += bn.pBH[i][k].prob;
        for (k = 2; --k >= 0; ) bn.pBH[i][k].prob /= sum;
        sum = 0;                    /* normalize cond. dists. on C */
        for (k = 2; --k >= 0; ) sum += bn.pCH[i][k].prob;
        for (k = 2; --k >= 0; ) bn.pCH[i][k].prob /= sum;
    }
    printf("initial state:\n");
    probs(&bn);                   /* show the initial probabilities */

    /* --- initialize other values --- */
    if (alg == 2) initchg = 1;    /* adapt the initial change value */
    for (i = 2; --i >= 0; ) {     /* traverse the values of H */
        bn.pH[i].emp = bn.pH[i].prv = 0; bn.pH[i].chg = initchg;
        for (k = 2; --k >= 0; ) {   /* traverse the attribute values */
            bn.pAH[i][k].emp = bn.pBH[i][k].emp = bn.pCH[i][k].emp =
            bn.pAH[i][k].prv = bn.pBH[i][k].prv = bn.pCH[i][k].prv = 0;
            bn.pAH[i][k].chg = bn.pBH[i][k].chg = bn.pCH[i][k].chg = initchg;
        }                           /* clear the EM value and prev. step */
    }                             /* and set the initial change value */

    /* --- expectation maximization --- */
    cnt = 0;                      /* initialize the step counter */
    do {                          /* expectation maximization loop */
        expect  (&bn, tab);         /* expectation  step */
        maximize(&bn, tab);         /* maximization step */
        chg = update(&bn, updfn);   /* update the Bayesian network */
        if ((++cnt <= 1) || ((verb > 0) && (cnt %verb == 0))) {
            printf("after %i iteration%s:\n", cnt, (cnt > 1) ? "s" : "");
            table(tab);               /* show data table */
            probs(&bn);               /* and probabilities */
        }
    } while ((chg > epsilon)      /* while no convergence */
             &&       ((maxcnt < 0) || (cnt < maxcnt)));

    /* --- print result --- */
    printf("result:\n");
    table(tab);                   /* show data table */
    probs(&bn);                   /* and probabilities */
    printf("%d iteration(s)\n", cnt);

    return 0;                     /* return `ok' */
}  /* main() */