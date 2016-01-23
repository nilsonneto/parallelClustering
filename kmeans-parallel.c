/*
 * Created by nilsonneto on 1/14/16.
 *
 * K-Means Algorithm
 * How to tackle:
 *     Starting point for means
 *         Forgy Algorithm
 *         Random Partition Algorithm
 *     Assignment step
 *     Update step
 *
 * Serial source: http://rosettacode.org/wiki/K-means%2B%2B_clustering
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#include "kmeans.h"

int numCores;

/*
 * Generate tiny random double
 */
double randf(double m) {
    return m * rand() / (RAND_MAX - 1.);
}

/*
 * Creates # (count) within the radius (radius)
 * Returns a vector with all points
 */
point gen_xy(long count, double radius, long numMeans) {
    double ang, r;
    point p;
    point pt = malloc(sizeof(point_t) * count);

    /* note: this is not a uniform 2-d distribution */
    for (p = pt + count; p-- > pt;) {
        ang = randf(2 * M_PI);
        r = randf(radius);
        p->x = r * cos(ang);
        p->y = r * sin(ang);
        p->group = (int) numMeans;
    }

    return pt;
}

/*
 * Returns the distance between 2 points
 */
double dist2(point a, point b) {
    double x, y;
    x = a->x - b->x;
    y = a->y - b->y;
    return x*x + y*y;
}

/*
 * Assigns the
 */
int nearest(point pt, point cent, long n_cluster, double *d2) {
    int i, min_i, ret = 0;
    double dist, minDist;

    for (i = 0; i < n_cluster; i++) {
        minDist = HUGE_VAL;
        min_i = pt->group;
        /*
         * For each mean/cluster, find closest to the point
         */
        for (i = 0; i < n_cluster; i++) {
            dist = dist2((cent+i), pt);
            if (minDist > dist) {
                minDist = dist;
                min_i = i;
            }
        }
    }
    if (d2)
        *d2 = minDist;

    if (min_i != pt->group) {
        pt->group = min_i;
        ret = 1;
    }

    return ret;
}

/*
 * K++ solution to initialize the variables differently from Lloyd's random approach
 */
void kpp(point pts, int len, point cent, int n_cent) {
    int i, j;
    int n_cluster;
    double sum, *d = malloc(sizeof(double) * len);

    point p, c;
    cent[0] = pts[ rand() % len ];
    for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
        sum = 0;
        for (j = 0, p = pts; j < len; j++, p++) {
            nearest(p, cent, n_cluster, d + j);
            sum += d[j];
        }
        sum = randf(sum);
        for (j = 0, p = pts; j < len; j++, p++) {
            if ((sum -= d[j]) > 0) continue;
            cent[n_cluster] = pts[j];
            break;
        }
    }
    for (j = 0, p = pts; j < len; j++, p++)
        p->group = nearest(p, cent, n_cluster, 0);

    free(d);
}

/*
 *
 */
void print_eps(point pts, long len, point cent, long n_cluster) {
#	define W 400
#	define H 400
    int i, j;
    double min_x, max_x, min_y, max_y, scale, cx, cy;
    double *colors = malloc(sizeof(double) * n_cluster * 3);

    for (i = 0; i < n_cluster; i++) {
        colors[3*i + 0] = (3 * (i + 1) % 11)/11.;
        colors[3*i + 1] = (7 * i % 11)/11.;
        colors[3*i + 2] = (9 * i % 11)/11.;
    }


    /*
     * Finds the largest and the smallest points
     */
    min_x = min_y = HUGE_VAL;
    max_x = max_y = -(HUGE_VAL);
    for (j = 0; j < len; j++) {
        if (max_x < (pts+j)->x)
            max_x = (pts+j)->x;

        if (min_x > (pts+j)->x)
            min_x = (pts+j)->x;

        if (max_y < (pts+j)->y)
            max_y = (pts+j)->y;

        if (min_y > (pts+j)->y)
            min_y = (pts+j)->y;
    }

    scale = W / (max_x - min_x);
    if (scale > H / (max_y - min_y)) scale = H / (max_y - min_y);
    cx = (max_x + min_x) / 2;
    cy = (max_y + min_y) / 2;

    printf("%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
    printf( "/l {rlineto} def /m {rmoveto} def\n"
                    "/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n"
                    "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
                    "	gsave 1 setgray fill grestore gsave 3 setlinewidth"
                    " 1 setgray stroke grestore 0 setgray stroke }def\n"
    );
    for (i = 0; i < n_cluster; i++) {
        printf("%g %g %g setrgbcolor\n",
               colors[3*i], colors[3*i + 1], colors[3*i + 2]);
        for (j = 0; j < len; j++) {
            if ((pts+j)->group != i) continue;
            printf("%.3f %.3f c\n",
                   ((pts+j)->x - cx) * scale + W / 2,
                   ((pts+j)->y - cy) * scale + H / 2);
        }
        printf("\n0 setgray %g %g s\n",
               ((cent+i)->x - cx) * scale + W / 2,
               ((cent+i)->y - cy) * scale + H / 2);
    }
    printf("\n%%%%EOF");
    free(colors);
}

/*
 * Lloyd's algorithm for resolving the k-means problem
 */
point lloyd(point pts, long len, long n_cluster) {
    int i, j;
    int changed;

    point cent = malloc(sizeof(point_t) * n_cluster);

    /* assign init grouping randomly */
    for (j = 0; j < len; j++) {
        (pts+j)->group = j % (int) n_cluster;
    }

    /* or call k++ init */
    // kpp(pts, len, cent, n_cluster);

    //#pragma omp  num_threads(numCores) default(shared)
    do {
        /* group element for centroids are used as counters */
        for (i = 0; i < n_cluster; i++) {
            (cent+i)->group = 0;
            (cent+i)->x = (cent+i)->y = 0;
        }

        for (j = 0; j < len; j++) {
            (cent+(pts+j)->group)->group++;
            (cent+(pts+j)->group)->x += (pts+j)->x;
            (cent+(pts+j)->group)->y += (pts+j)->y;
        }

        for (i = 0; i < n_cluster; i++) {
            (cent+i)->x /= (cent+i)->group;
            (cent+i)->y /= (cent+i)->group;
        }

        changed = 0;

        /* find closest centroid of each point */
        #pragma omp parallel for private(j) reduction(+:changed)
        for (j = 0; j < len; j++) {
            changed += nearest((pts+j), cent, n_cluster, 0);
        }
    } while (changed > (len >> 10)); /* stop when 99.9% of points are good */

    for (i = 0; i < n_cluster; i++) {
        (cent+i)->group = i;
    }

    return cent;
}

double getTempo(struct timeval start, struct timeval stop){
    double t = (((double)(stop.tv_sec)*1000.0  + (double)(stop.tv_usec / 1000.0)) - \
                   ((double)(start.tv_sec)*1000.0 + (double)(start.tv_usec / 1000.0)));
    return (t);
}

void tempo(struct timeval start, struct timeval stop){
    printf("%g ms\n", getTempo(start, stop));
}

int main(int argc, char *argv[]) {
    struct timeval start, stop;

    long numPoints = strtol(argv[1], NULL, 10);
    long radius = strtol(argv[2], NULL, 10);
    long numMeans = strtol(argv[3], NULL, 10);
    numCores = 2;

    int time = 0;
    int gen = 0;

    gettimeofday(&start, NULL);
    point v = gen_xy(numPoints, radius, numMeans);
    gettimeofday(&stop, NULL);
    if (time && gen) fprintf(stdout, "Generation time: %lf\n", getTempo(start, stop));

    gettimeofday(&start, NULL);
    point c = lloyd(v, numPoints, numMeans);
    gettimeofday(&stop, NULL);
    if (time) fprintf(stdout, "Classification time: %lf\n", getTempo(start, stop));

    if (!time) print_eps(v, numPoints, c, numMeans);

    return 0;
}
