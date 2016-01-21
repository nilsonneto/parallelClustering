//
// Created by nilsonneto on 1/21/16.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "kmeans.h"
void print_eps(point pts, long len, point cent, long n_cluster) {
#	define W 400
#	define H 400
    int i, j;
    point p, c;
    double min_x, max_x, min_y, max_y, scale, cx, cy;
    double *colors = malloc(sizeof(double) * n_cluster * 3);

    for (c = cent, i = 0; i < n_cluster; i++, c++) {
        colors[3*i + 0] = (3 * (i + 1) % 11)/11.;
        colors[3*i + 1] = (7 * i % 11)/11.;
        colors[3*i + 2] = (9 * i % 11)/11.;
    }


    /*
     * Finds the largest and the smallest points
     */
    min_x = min_y = HUGE_VAL;
    max_x = max_y = -(HUGE_VAL);
    for (j = 0, p = pts; j < len; j++, p++) {
        if (max_x < p->x)
            max_x = p->x;

        if (min_x > p->x)
            min_x = p->x;

        if (max_y < p->y)
            max_y = p->y;

        if (min_y > p->y)
            min_y = p->y;
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
    for (c = cent, i = 0; i < n_cluster; i++, c++) {
        printf("%g %g %g setrgbcolor\n",
               colors[3*i], colors[3*i + 1], colors[3*i + 2]);
        for (j = 0, p = pts; j < len; j++, p++) {
            if (p->group != i) continue;
            printf("%.3f %.3f c\n",
                   (p->x - cx) * scale + W / 2,
                   (p->y - cy) * scale + H / 2);
        }
        printf("\n0 setgray %g %g s\n",
               (c->x - cx) * scale + W / 2,
               (c->y - cy) * scale + H / 2);
    }
    printf("\n%%%%EOF");
    free(colors);
}