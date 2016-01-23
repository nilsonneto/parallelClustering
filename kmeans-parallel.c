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
#include <sys/time.h>
#include <math.h>
#include <omp.h>

typedef struct _point_t{
    double x, y;
    int group;
} point_t;

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
void gen_xy(point_t *points, long count, double radius, long numMeans) {
    int i;
    double ang, r;

    /* note: this is not a uniform 2-d distribution */
    for (i = 0; i < count; i++) {
        ang = randf(2 * M_PI);
        r = randf(radius);
        points[i].x = r * cos(ang);
        points[i].y = r * sin(ang);
        points[i].group = (int) numMeans;
    }
}

/*
 * Returns the distance between 2 points
 */
double dist2(double ax, double ay, double bx, double by) {
    double x, y;
    x = ax - bx;
    y = ay - by;
    return x*x + y*y;
}

/*
 * Assigns the
 */
int nearest(point_t *points, point_t *means, long n_cluster, double *d2) {
    int i, min_i;
    double dist, minDist = 0;

    for (i = 0; i < n_cluster; i++) {
        minDist = HUGE_VAL;
        min_i = points[i].group;
        /*
         * For each mean/cluster, find closest to the point
         */
        for (i = 0; i < n_cluster; i++) {
            dist = dist2(means[i].x, means[i].y, points[i].x, points[i].y);
            if (minDist > dist) {
                minDist = dist;
                min_i = i;
            }
        }
    }
    if (d2)
        *d2 = minDist;

    return min_i;
}

/*
 * K++ solution to initialize the variables differently from Lloyd's random approach
 *
void kpp(point_t pts, int len, point_t *cent, int n_cent) {
    int i, j;
    int n_cluster;
    double sum, *d = malloc(sizeof(double) * len);

    point_t p, c;
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
void print_eps(point_t *points, long len, point_t *means, long n_cluster) {
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
        if (max_x < points[j].x)
            max_x = points[j].x;

        if (min_x > points[j].x)
            min_x = points[j].x;

        if (max_y < points[j].y)
            max_y = points[j].y;

        if (min_y > points[j].y)
            min_y = points[j].y;
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
            if (points[j].group != i) continue;
            printf("%.3f %.3f c\n",
                   (points[j].x - cx) * scale + W / 2,
                   (points[j].y - cy) * scale + H / 2);
        }
        printf("\n0 setgray %g %g s\n",
               (means[i].x - cx) * scale + W / 2,
               (means[i].y - cy) * scale + H / 2);
    }
    printf("\n%%%%EOF");
    free(colors);
}

/*
 * Lloyd's algorithm for resolving the k-means problem
 */
void lloyd(point_t *points, point_t *means, long len, long numMeans) {
    // struct timeval start, stop;
    int i, j;
    long changed;
    // double count1 = 0, count2 = 0, count3 = 0, count4 = 0;

    /* assign init grouping randomly */
    for (j = 0; j < len; j++) {
        points[j].group = j % (int) numMeans;
    }

    /* or call k++ init */
    // kpp(pts, len, means, numMeans);

    // #pragma omp parallel num_threads(numCores) default(none) private(j, min_i, changed) shared(len, means, pts, numMeans, maxP)
    do {
        /*
         * group element for centroids are used as counters
         * Go through each mean and clear it
         */
        //gettimeofday(&start, NULL);
        for (i = 0; i < numMeans; i++) {
            means[i].group = 0;
            means[i].x = means[i].y = 0;
        }
        //gettimeofday(&stop, NULL);
        //count1 += getTempo(start, stop);

        /*
         * Go through each point, go to it's assigned mean
         * Increase mean group counter
         * Add to x and y of the mean, the x and y of the point
         */
        //gettimeofday(&start, NULL);
        // #pragma omp for
        for (j = 0; j < len; j++) {
            means[points[i].group].group++;
            means[points[i].group].x += points[i].x;
            means[points[i].group].y += points[i].y;
        }
        //gettimeofday(&stop, NULL);
        //count2 += getTempo(start, stop);

        /*
         * Go through each mean and divide the position my the number of assignees
         * This will reposition the mean based on all of their assignees
         */
        //gettimeofday(&start, NULL);
        for (i = 0; i < numMeans; i++) {
            means[i].x /= means[i].group;
            means[i].y /= means[i].group;
        }
        //gettimeofday(&stop, NULL);
        //count3 += getTempo(start, stop);


        /*
         * Go through each point and find its nearest mean
         * If different that previous, change it and increase counter
         */
        //gettimeofday(&start, NULL);
        changed = 0;
        // #pragma omp parallel for num_threads(numCores) default(none) private(j, k, p, c) shared(len, means, pts)
        for (j = 0; j < len; j++) {
            int min_i = nearest(points, means, numMeans, 0);
            if (min_i != points[j].group) {
                changed++;
                points[j].group = min_i;
            }
        }
        //gettimeofday(&stop, NULL);
        //count4 += getTempo(start, stop);
    } while (changed > (len >> 10)); /* stop when 99.9% of points are good */
    //printf("%lf // %lf // %lf // %lf \n", count1, count2, count3, count4);

    for (i = 0; i < numMeans; i++) {
        means[i].group = i;
    }
}

double getTempo(struct timeval start, struct timeval stop){
    double t = (((double) (stop.tv_sec) * 1000.0 + stop.tv_usec / 1000.0) - \
                   ((double)(start.tv_sec)*1000.0 + start.tv_usec / 1000.0));
    return (t);
}

void tempo(struct timeval start, struct timeval stop){
    printf("%g ms\n", getTempo(start, stop));
}

int main(int argc, char *argv[]) {
    struct timeval start, stop;
    // int i;

    long numPoints = strtol(argv[1], NULL, 10);
    long radius = strtol(argv[2], NULL, 10);
    long numMeans = strtol(argv[3], NULL, 10);
    numCores = 4;

    point_t *means;
    means = (point_t *) malloc(numMeans * sizeof(point_t));

    point_t *points;
    points = (point_t *) malloc(numPoints * sizeof(point_t));

    gettimeofday(&start, NULL);
    gen_xy(points, numPoints, radius, numMeans);
    gettimeofday(&stop, NULL);
    //fprintf(stdout, "Generation time: ");
    //tempo(start, stop);

    gettimeofday(&start, NULL);
    lloyd(points, means, numPoints, numMeans);
    gettimeofday(&stop, NULL);
    //fprintf(stdout, "Classification time: ");
    //tempo(start, stop);

    print_eps(points, numPoints, means, numMeans);

    return 0;
}

