/* 
   * distanzen speichern, wenn 2 cluster verschmolzen werden, distanz zu allen
   anderen clustern berechnen
*/

#include <R.h>
#include <Rdefines.h>

void
R_hra_hclust(double *x, int *num, double *d, int *g, double *eppm,
	     double *eabs)
{
	/* anzahl punkte */
	int n = *num;
	int m = n + 2;

	double *d_idx[n - 1];
	double means[n];
	double mrange[2];
	int overlimit = 0;

	int minclust[2];
	double mindst = DBL_MAX;

	/*loop vars */
	int i, j, z;

	int *clust;
	clust = malloc(n * m * sizeof(int));

	/* calculate indices for distance matrix */
	j = n - 1;
	d_idx[0] = d;
	for (i = 1; i < n - 1; i++) {
		d_idx[i] = d_idx[i - 1] + j;
		j = j - 1;
	}

	/* init means */
	for (i = 0; i < n; i++) {
		means[i] = x[i];
	}

	/* init cluster matrix */
	for (z = 0; z < n; z++) {
		clust[z * m] = 1;
		clust[z * m + 1] = z;
	}

	for (z = 0; z < n-1; z++) {
		for (i = 0; i < n; i++) {
			if (clust[i * m] <= 0)
				continue;
			for (j = i + 1; j < n; j++) {
				if (clust[j * m] <= 0)
					continue;
				/*aufsummieren z.b. 5 und 9, index = sum(1:5)+(9-5) */
				if (*(d_idx[i] + (j - i - 1)) < mindst) {
					mindst = *(d_idx[i] + (j - i - 1));
					minclust[0] = i;
					minclust[1] = j;
				}
			}
		}
		/* mindst ist DBL_MAX, keine cluster mehr */
		if (mindst == DBL_MAX)
			break;
		/* berechne intervall fuer abbruchkriterium */
		means[minclust[0]] = (means[minclust[0]] *
				      clust[minclust[0] * m] +
				      means[minclust[1]] *
				      clust[minclust[1] * m]) / (double)(clust[minclust[0] * m] + clust[minclust[1] * m]);
		mrange[0] =
		    means[minclust[0]] - means[minclust[0]] * *eppm -
		    *eabs;
		mrange[1] =
		    means[minclust[0]] + means[minclust[0]] * *eppm +
		    *eabs;

		/* finde ausreisser */
		for (i = 1; i <= clust[minclust[0] * m]; i++) {
			if (x[clust[minclust[0] * m + i]] < mrange[0] ||
			    x[clust[minclust[0] * m + i]] > mrange[1]) {
				overlimit = 1;
				break;
			}
		}
		if (!overlimit)
			for (i = 1; i <= clust[minclust[1] * m]; i++) {
				if (x[clust[minclust[1] * m + i]] < mrange[0]
				    || x[clust[minclust[1] * m + i]] >
				    mrange[1]) {
					overlimit = 1;
					break;
				}
			}

		/* wenn abbruchkriterium erfüllt, cluster deaktivieren */
		if (overlimit) {
			clust[minclust[0] * m] *= -1;
			clust[minclust[1] * m] *= -1;
			overlimit = 0;
		} else {
			/* sonst cluster verschmelzen */
			for (i = 1; i <= clust[minclust[1] * m]; i++) {
				clust[minclust[0] * m +
				      clust[minclust[0] * m] + i] =
				    clust[minclust[1] * m + i];
			}
			/* clusterlaenge aktualisieren */
			clust[minclust[0] * m] += clust[minclust[1] * m];
			/* 2. cluster raus */
			clust[minclust[1] * m] = 0;

			/* distanz aktualisieren */
			/* abstand zwischen cluster a und {c[i],c[j]} ist
			 * max(c[i] zu a, c[j] zu a) */
			for(i = 0; i<n; i++) {
				if(clust[i * m] <= 0 || i == minclust[0])
					continue;
				if(i<minclust[0]) {
					if(i<minclust[1]) {
						if(*(d_idx[i]+(minclust[0]-i-1)) <
								*(d_idx[i]+(minclust[1]-i-1)))
						*(d_idx[i]+(minclust[0]-i-1)) =
								*(d_idx[i]+(minclust[1]-i-1));
					} else {
						if(*(d_idx[i]+(minclust[0]-i-1)) <
								*(d_idx[minclust[1]]+(i-minclust[1]-1)))
						*(d_idx[i]+(minclust[0]-i-1)) =
								*(d_idx[minclust[1]]+(i-minclust[1]-1));
					}
				} else {
					if(i<minclust[1]) {
						if(*(d_idx[minclust[0]]+(i-minclust[0]-1)) <
								*(d_idx[i]+(minclust[1]-i-1)))
						*(d_idx[minclust[0]]+(i-minclust[0]-1)) =
								*(d_idx[i]+(minclust[1]-i-1));
					} else {
						if(*(d_idx[minclust[0]]+(i-minclust[0]-1)) <
								*(d_idx[minclust[1]]+(i-minclust[1]-1)))
						*(d_idx[minclust[0]]+(i-minclust[0]-1)) =
								*(d_idx[minclust[1]]+(i-minclust[1]-1));
					}
				}
			}


		}
		mindst = DBL_MAX;
	}
	int gnum = 1;
	int tmp;
	for (i = 0; i < n; i++) {
		tmp = abs(clust[i * m]);
		if (tmp == 0)
			continue;
		for (j = 1; j <= tmp; j++) {
			g[clust[i * m + j]] = gnum;
		}
		gnum++;
	}
	free(clust);
}
