#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "util.h"
#include "xcms.h"

void ProfBinLin(double *xvals, double *yvals, int *numin,
                double *xstart, double *xend, int *numout, double *out) {

    double dx, xi,
      xpre=-1, ypre=-1,
      xpost, ypost, startx;
    int    i, ipost;

    dx = (*numout != 1) ? (*xend - *xstart)/(*numout - 1) : (*xend - *xstart);
    // why 20 here?
    startx = *xstart - 20*dx;
    //
    FindEqualLess(xvals, numin, &startx, &ipost);
    //ipost = 0;
    // this will be 0 if ipost = 0
    xpost = *xstart + dx*(int)((xvals[ipost] - *xstart)/dx + 0.5);
    ypost = yvals[ipost];

    for (i = 0; i < *numout; i++) {
       xi = *xstart + dx*i;
       if (xi < xvals[0] || xi > xvals[(*numin)-1])
           out[i] = 0;
       else {
           while (xi > xpost && ipost < *numin-1) {
               xpre = xpost;
               ypre = ypost;
               ipost++;
               xpost = *xstart + dx*(int)((xvals[ipost] - *xstart)/dx + 0.5);
               ypost = yvals[ipost];
               while (ipost < *numin-1 && xpost == *xstart + dx*(int)((xvals[ipost+1] - *xstart)/dx + 0.5)) {
                   ipost++;
                   ypost = (ypost > yvals[ipost]) ? ypost : yvals[ipost];
               }
           }
           out[i] = ypre + (xi-xpre)*(ypost-ypre)/(xpost-xpre);
       }
    }
}

void ProfBinLinM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                 double *xstart, double *xend, int *numout, double *out) {

    int i, vectlen;

    for (i = 0; i < *nummi; i++) {
        if (i < *nummi-1)
            vectlen = mindex[i+1] - mindex[i];
        else
            vectlen = *numin - mindex[i];
        ProfBinLin(xvals+mindex[i], yvals+mindex[i], &vectlen, xstart, xend,
                   numout, out+i*(*numout));
    }
}

/*
 * xvals vector of x values; binning is performed on these.
 * yvals vector of y values; these are binned.
 * numin
 */
void ProfBinLinBase(double *xvals, double *yvals, int *numin, double *baselevel, double *basespace,
                    double *xstart, double *xend, int *numout, double *out) {

    double dx, ypre = -1, ypost = -1, startx;
    int    i, ipre = -1, ipost, ix, ibase;

    /* dx: the bin size.
     * ypre: y value.
     * ypost: y value.
     * startx
     * i index
     * ipre: index of last value.
     * ipost: index of next value.
     * ix:
     * ibase: the number of neighboring bins considered for interpolation.
     */
    dx = (*numout != 1) ? (*xend - *xstart)/(*numout - 1) : (*xend - *xstart);
    ibase = floor(*basespace/dx);

    // Initialize the bin after the first interpolation point (post)
    startx = *xstart + 0.5*dx;
    FindEqualLess(xvals, numin, &startx, &ix);
    ipost = round((xvals[ix] - *xstart)/dx);
    ypost = yvals[ix];
    if (ipost <= 0) { // Found one less, search backwards in the bin
        i = ix;
        while (--ix >= 0 && round((xvals[ix] - *xstart)/dx) == ipost)
            if (yvals[ix] > ypost)
                ypost = yvals[ix];
        ix = i + 1;
    }
    else // Found one greater, search forwards in the bin
        while (++ix < *numin && round((xvals[ix] - *xstart)/dx) == ipost)
            if (yvals[ix] > ypost)
                ypost = yvals[ix];

    // loop through the bins.
    for (i = 0; i < *numout; i++) {
        // Move post to pre if less than or equal to interpolating point
        if (ipost <= i && ypost != -1) {
            ypre = ypost;
            ipre = ipost;
            ypost = -1;
        }
        // Find the next post if there are more points in the queue
        if (ypost == -1 && ix < *numin) {
            ipost = round((xvals[ix] - *xstart)/dx);
            ypost = yvals[ix];
            while (++ix < *numin && round((xvals[ix] - *xstart)/dx) == ipost)
                if (yvals[ix] > ypost)
                    ypost = yvals[ix];
        }
	/*
	 * If I understand it correctly, ibase is used to devine whether interpolation
	 * takes place or not. And interpolation takes place if i, for a bin with a missing
	 * value, is close to another bin with a value.
	 */
        if (ipre == i)
	  out[i] = ypre;  // assuming that's the binning...
        else if (ypre != -1 && ypost != -1 && (ipost-ipre <= 2*ibase+1))
	  out[i] = ypre + (ypost-ypre)/(ipost-ipre)*(i-ipre);  // that's interpolation between existing bins
        else if (ypre != -1 && i-ipre <= ibase && (ypost == -1 || ipost-i > ibase))
	  out[i] = ypre + (*baselevel-ypre)/(ibase+1)*(i-ipre);  // interpolating: decreasing to baselevel
        else if ((ypre == -1 || i-ipre > ibase) && ypost != -1 && ipost-i <= ibase)
	  out[i] = *baselevel + (ypost-*baselevel)/(ibase+1)*(i-ipost+ibase+1);  // interpolating: increasing from baselevel
        else
	  out[i] = *baselevel; // base
    }
}

void ProfBinLinBaseM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                     double *baselevel, double *basespace, double *xstart, double *xend,
                     int *numout, double *out) {

    int i, vectlen;

    for (i = 0; i < *nummi; i++) {
        if (i < *nummi-1)
            vectlen = mindex[i+1] - mindex[i];
        else
            vectlen = *numin - mindex[i];
        ProfBinLinBase(xvals+mindex[i], yvals+mindex[i], &vectlen, baselevel, basespace,
                       xstart, xend, numout, out+i*(*numout));
    }
}

void ProfIntLin(double *xvals, double *yvals, int *numin,
                double *xstart, double *xend, int *numout, double *out) {

    // Search for initial j
    int    i, j = 0, thru;
    double dx, x1, x2, xb, xe, yb, ye, totarea, startx;

    dx = (*numout != 1) ? (*xend - *xstart)/(*numout - 1) : (*xend - *xstart);

    startx = *xstart - dx;
    FindEqualLess(xvals, numin, &startx, &j);
    x2 = *xstart - dx*0.5;
    for (i = 0; i < *numout; i++) {
        x1 = x2;
        x2 = *xstart + dx*(i + 0.5);
        totarea = 0;
        if (x2 <= xvals[0] || x1 >= xvals[*numin-1]) {
            out[i] = 0;
            continue;
        }
        thru = 0;
        while ((! thru && j < *numin-1) || (j < *numin-1 && xvals[j+1] <= x2)) {
            if (xvals[j+1] > x2 || j >= *numin-2)
                thru = 1;
            if (xvals[j+1] <= x1) {
                j++;
                continue;
            }

            if (xvals[j] < x1) {
                xb = x1;
                yb = yvals[j] + (yvals[j+1]-yvals[j])*(x1-xvals[j])/(xvals[j+1]-xvals[j]);
            } else {
                xb = xvals[j];
                yb = yvals[j];
            }
            if (xvals[j+1] > x2) {
                xe = x2;
                ye = yvals[j] + (yvals[j+1]-yvals[j])*(x2-xvals[j])/(xvals[j+1]-xvals[j]);
            } else {
                xe = xvals[j+1];
                ye = yvals[j+1];
            }
            totarea += (ye+yb)*(xe-xb)/2;
            if (xvals[j+1] <= x2)
                j++;
        }
        out[i] = totarea/dx;
    }
}

void ProfIntLinM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                 double *xstart, double *xend, int *numout, double *out) {

    int i, vectlen;

    for (i = 0; i < *nummi; i++) {
        if (i < *nummi-1)
            vectlen = mindex[i+1] - mindex[i];
        else
            vectlen = *numin - mindex[i];
        ProfIntLin(xvals+mindex[i], yvals+mindex[i], &vectlen, xstart, xend,
                   numout, out+i*(*numout));
    }
}

void ProfBin(double *xvals, double *yvals, int *numin,
             double *xstart, double *xend, int *numout, double *out) {

    int    i, outi = 0;
    double dx, startx, endx;

    // calculate the step size
    dx = (*numout != 1) ? (*xend - *xstart)/(*numout - 1) : (*xend - *xstart);

    // fill out with 0s.
    for (i = 0; i < *numout; i++)
        out[i] = 0;

    startx = *xstart - dx;
    endx = *xend + dx;
    FindEqualGreater(xvals, numin, &startx, &i);
    for (; i < *numin && xvals[i] < endx; i++) {
      outi = (int)floor((xvals[i] - *xstart)/dx + 0.5);
      if (outi >= 0 && outi < *numout)
	if (out[outi] < yvals[i])
	  out[outi] = yvals[i];
    }
}

void ProfBin_test(double *xvals, double *yvals, int *numin,
		  double *xstart, double *xend, int *numout,
		  double *out, double *the_dx) {

    int    i, outi = 0;
    double dx, startx, endx;

    // calculate the step size
    dx = (*numout != 1) ? (*xend - *xstart)/(*numout - 1) : (*xend - *xstart);

    the_dx[0] = dx;
    // fill out with 0s.
    for (i = 0; i < *numout; i++)
        out[i] = 0;

    startx = *xstart - dx;
    endx = *xend + dx;
    FindEqualGreater(xvals, numin, &startx, &i);
    for (; i < *numin && xvals[i] < endx; i++) {
      outi = (int)floor((xvals[i] - *xstart)/dx + 0.5);
      if (outi >= 0 && outi < *numout)
	if (out[outi] < yvals[i])
	  out[outi] = yvals[i];
    }
}

void ProfBinM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
              double *xstart, double *xend, int *numout, double *out) {

    int i, vectlen;

    for (i = 0; i < *nummi; i++) {
        if (i < *nummi-1)
            vectlen = mindex[i+1] - mindex[i];
        else
            vectlen = *numin - mindex[i];
        ProfBin(xvals+mindex[i], yvals+mindex[i], &vectlen, xstart, xend,
                numout, out+i*(*numout));
    }
}

void ProfMaxIdx(double *xvals, double *yvals, int *numin,
                double *xstart, double *xend, int *numout, int *out) {

    int    i, outi = 0;
    double dx, startx, endx;

    dx = (*numout != 1) ? (*xend - *xstart)/(*numout - 1) : (*xend - *xstart);

    for (i = 0; i < *numout; i++)
        out[i] = INT_MIN;

    startx = *xstart - dx;
    endx = *xend + dx;
    FindEqualGreater(xvals, numin, &startx, &i);
    for (; i < *numin && xvals[i] < endx; i++) {
        outi = (int)floor((xvals[i] - *xstart)/dx + 0.5);
        if (outi >= 0 && outi < *numout)
            if (out[outi] < 0 || yvals[out[outi]] < yvals[i])
                out[outi] = i;
    }
}

void ProfMaxIdxM(double *xvals, double *yvals, int *numin, int *mindex, int *nummi,
                 double *xstart, double *xend, int *numout, int *out) {

    int i, j, vectlen;

    for (i = 0; i < *nummi; i++) {
        if (i < *nummi-1)
            vectlen = mindex[i+1] - mindex[i];
        else
            vectlen = *numin - mindex[i];
        ProfMaxIdx(xvals+mindex[i], yvals+mindex[i], &vectlen, xstart, xend,
                   numout, out+i*(*numout));
        for (j = i*(*numout); j < (i+1)*(*numout); j++)
            if (out[j] >= 0)
                out[j] += mindex[i]+1;
    }
}

void MedianFilter(double *inmat, int *m, int *n, int *mrad, int *nrad, double *outmat) {

    int    i, j, k, l, mmin, mmax, nmin, nmax, bufLen;
    double * sortBuf = (double *) malloc((*mrad*2+1)*(*nrad*2+1)*sizeof(double));

    for (i = 0; i < *m; i++) {
        for (j = 0; j < *n; j++) {
            mmin = (i - *mrad > 0) ? (i - *mrad) : 0;
            mmax = (i + *mrad < *m) ? (i + *mrad) : *m - 1;
            nmin = (j - *nrad > 0) ? (j - *nrad) : 0;
            nmax = (j + *nrad < *n) ? (j + *nrad) : *n - 1;
            bufLen = 0;
            for (k = mmin; k <= mmax; k++) {
                for (l = nmin; l <= nmax; l++) {
                    sortBuf[bufLen] = *(inmat + k + *m*l);
                    bufLen++;
                }
            }
            qsort(sortBuf, bufLen, sizeof(double), CompareDouble);
            if (bufLen % 2 == 1)
                *(outmat + i + *m*j) = sortBuf[(bufLen-1)/2];
            else
                *(outmat + i + *m*j) = (sortBuf[(bufLen-2)/2] + sortBuf[(bufLen)/2])/2;
        }
    }

    free(sortBuf);
}

int CompareDouble(const void *a, const void *b) {

    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}
