//
// Created by nilsonneto on 1/14/16.
//

#ifndef PARALLELCLUSTERING_KMEANS_H
#define PARALLELCLUSTERING_KMEANS_H

#define PTS 100000
#define K 11

typedef struct {
    double x, y;
    int group;
} point_t, *point;

double randf(double m);
point gen_xy(int count, double radius);
double dist2(point a, point b);
int nearest(point pt, point cent, int n_cluster, double *d2);
void kpp(point pts, int len, point cent, int n_cent);
point lloyd(point pts, int len, int n_cluster);
void print_eps(point pts, int len, point cent, int n_cluster);

#endif //PARALLELCLUSTERING_KMEANS_H
