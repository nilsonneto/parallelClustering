//
// Created by nilsonneto on 1/14/16.
//

#ifndef PARALLELCLUSTERING_KMEANS_H
#define PARALLELCLUSTERING_KMEANS_H


typedef struct {
    double x, y;
    int group;
} point_t, *point;

double randf(double m);
point gen_xy(long count, double radius, long numMeans);
double dist2(point a, point b);
int nearest(point pt, point cent, long n_cluster, double *d2);
void kpp(point pts, int len, point cent, int n_cent);
point lloyd(point pts, long len, long numMeans);
void print_eps(point pts, long len, point cent, long n_cluster);
double getTempo(struct timeval start, struct timeval stop);
void tempo(struct timeval start, struct timeval stop);

#endif //PARALLELCLUSTERING_KMEANS_H
