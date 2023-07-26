#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "algorthm.h"
#include "contdist.h"
#include "correlat.h"
#include "descript.h"
#include "deviate.h"
#include "discdist.h"
#include "eigensys.h"
#include "linalg.h"
#include "matutils.h"
#include "ranking.h"
#include "residual.h"
#include "rootfind.h"
#include "sort.h"
#include "vecutils.h"

#define SIZE 5

int randomValue(int n)
{
	return (int)((double)rand() / ((double)RAND_MAX + 1) * n);
}

void initArray(double* first, double* last)
{
	while(first != last){
		*first = randomValue(SIZE);
		++first;
	}
}

void outputArray(double* array, int n)
{
	int i;

	printf("%.2lf", array[0]);
	for(i = 1; i < n; ++i){
		printf(" %.2lf", array[i]);
	}
	printf("\n");
}

void outputIArray(int* array, int n)
{
	int i;

	printf("%d", array[0]);
	for(i = 1; i < n; ++i){
		printf(" %d", array[i]);
	}
	printf("\n");
}

void outputMatrix(double** mat, int r, int c)
{
	int i;

	for(i = 0; i < r; ++i){
		outputArray(mat[i], c);
	}
}

void copyArray(double* src, int n, double* dest)
{
	int i;

	for(i = 0; i < n; ++i){
		dest[i] = src[i];
	}
}

void copyMatrix(double** src, int r, int c, double** dest)
{
	int i;

	for(i = 0; i < r; ++i){
		copyArray(src[i], c, dest[i]);
	}
}

void testAlgorithm()
{
	double x[SIZE];
	double y[SIZE];
	double* pair;
	double* walsh;

	initArray(x, x + SIZE);
	initArray(y, y + SIZE);

	printf("x: "); outputArray(x, SIZE);
	printf("y: "); outputArray(y, SIZE);

	printf("average(x[0], x[1]): %lf\n", average(x[0], x[1]));

	pair = dvector(0, SIZE * SIZE - 1);
	pairwdiff(x, x + SIZE, y, y + SIZE, pair);
	printf("pairwdiff(x, y): "); outputArray(pair, SIZE * SIZE);
	free_dvector(pair, 0);
	
	walsh = dvector(0, SIZE * (SIZE + 1) / 2);
	walshavg(x, x + SIZE, walsh);
	printf("walshavg(x): "); outputArray(walsh, SIZE * ( SIZE + 1) / 2);
	free_dvector(walsh, 0);
}

void testContDist()
{
	double p;
	double x;

	p = betap(.5, 1, 2); x = betav(p, 1, 2);
	printf("betap(.5, 1, 2): %.2lf     betav(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = cauchyp(.5, 1, 2); x = cauchyv(p, 1, 2);
	printf("cauchyp(.5, 1, 2): %.2lf     cauchyv(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = dblexpp(.5, 1, 2); x = dblexpv(p, 1, 2);
	printf("dblexpp(.5, 1, 2): %.2lf     dblexpv(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = expp(.5, 1); x = expv(p, 1);
	printf("expp(.5, 1, 2): %.2lf     expv(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = gammap(.5, 1, 2); x = gammav(p, 1, 2);
	printf("gammap(.5, 1, 2): %.2lf     gammav(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = logisticp(.5, 1, 2); x = logisticv(p, 1, 2);
	printf("logisticp(.5, 1, 2): %.2lf     logisticv(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = lognormp(.5, 1, 2); x = lognormv(p, 1, 2);
	printf("lognormp(.5, 1, 2): %.2lf     lognormv(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = rayleighp(.5, 1); x = rayleighv(p, 1);
	printf("rayleighp(.5, 1, 2): %.2lf     rayleighv(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = uniformp(1.5, 1, 2); x = uniformv(p, 1, 2);
	printf("uniformp(.5, 1, 2): %.2lf     uniformv(%.2lf, 1, 2): %.2lf\n", p, p, x);

	p = weibullp(.5, 1, 2); x = weibullv(p, 1, 2);
	printf("weibullp(.5, 1, 2): %.2lf     weibullv(%.2lf, 1, 2): %.2lf\n", p, p, x);
}

void testCorrelation()
{
	double x[SIZE];
	double* y;

	y = dvector(0, SIZE - 1);

	initArray(x, x + SIZE);
	initArray(y, y + SIZE);

	printf("x: "); outputArray(x, SIZE);
	printf("y: "); outputArray(y, SIZE);

	printf("covariance(x, y): %lf\n", covariance(x, x + SIZE, y));
	printf("correlation(x, y): %lf\n", correlation(y, y + SIZE, x));
	printf("kendall(x, y): %lf\n", kendall(x, x + SIZE, y));

	free_dvector(y, 0);
}

void testDescript()
{
	double x[SIZE];
	double y[SIZE];
	double* z;

	z = dvector(0, SIZE - 1);

	initArray(x, x + SIZE);
	initArray(y, y + SIZE);
	initArray(z, z + SIZE);

	printf("centmon(x, 3): %lf\n", centmom(x, x + SIZE, 3));
	printf("coeffvar(y): %lf\n", coeffvar(y, y + SIZE));
	printf("geomean(z): %lf\n", geomean(z, z + SIZE));
	printf("harmean(y): %lf\n", harmean(y, y + SIZE));
	printf("iqr(x): %lf\n", iqr(x, x + SIZE));
	printf("kurtosis(y): %lf\n", kurtosis(y, y + SIZE));
	printf("mean(z): %lf\n", mean(z, z + SIZE));
	printf("meandev(y): %lf\n", meandev(y, y + SIZE));
	printf("meddev(x): %lf\n", meddev(x, x + SIZE));
	printf("medain(y): %lf\n", median(y, y + SIZE));
	printf("quantile(z): %lf\n", quantile(z, z + SIZE, .1));
	printf("quartile1(y): %lf\n", quartile1(y, y + SIZE));
	printf("quartile3(x): %lf\n", quartile3(x, x + SIZE));
	printf("range(y): %lf\n", range(y, y + SIZE));
	printf("rms(z): %lf\n", rms(z, z + SIZE));
	printf("skewness(y): %lf\n", skewness(y, y + SIZE));
	printf("stddev(x): %lf\n", stddev(x, x + SIZE));
	printf("stderrmean(y): %lf\n", stderrmean(y, y + SIZE));
	printf("sum(x): %lf\n", sum(x, x + SIZE));
	printf("sum2(y): %lf\n", sum2(y, y + SIZE));
	printf("sump(x): %lf\n", sump(x, x + SIZE, y));
	printf("sumsquare(y): %lf\n", sumsquare(y, y + SIZE));
	printf("variance(z): %lf\n", variance(z, z + SIZE));

	free_dvector(z, 0);
}

void testDeviate()
{
	int n = SIZE;
	int i;
	long seed = 10203040L;

	printf("ran0: ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", ran0(&seed)); } printf("\n");;
	
	printf("ran1: ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", ran1(&seed)); } printf("\n");
	
	printf("ran2: ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", ran2(&seed)); } printf("\n");
	
	printf("ran3: ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", ran3(&seed)); } printf("\n");
	
	printf("berdev(.5): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", berdev(.5, &seed)); } printf("\n");
	
	printf("betadev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", betadev(1, 2, &seed)); } printf("\n");
	
	printf("bnldev(10, .5): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", bnldev(10, .5, &seed)); } printf("\n");
	
	printf("caudev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", caudev(1, 2, &seed)); } printf("\n");
	
	printf("chidev(1): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", chidev(1, &seed)); } printf("\n");
	
	printf("dexpdev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", dexpdev(1, 2, &seed)); } printf("\n");
	
	printf("expdev(1): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", expdev(1, &seed)); } printf("\n");
	
	printf("fdev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", fdev(1, 2, &seed)); } printf("\n");
	
	printf("gamdev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", gamdev(1, 2, &seed)); } printf("\n");
	
	printf("gasdev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", gasdev(1, 2, &seed)); } printf("\n");
	
	printf("logdev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", logdev(1, 2, &seed)); } printf("\n");
	
	printf("nchidev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", nchidev(1, 2, &seed)); } printf("\n");
	
	printf("nfdev(1, 2, 3): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", nfdev(1, 2, 3, &seed)); } printf("\n");
	
	printf("poisdev(2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", poisdev(2, &seed)); } printf("\n");
	
	printf("tdev(1): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", tdev(1, &seed)); } printf("\n");
	
	printf("unifdev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", unifdev(1, 2, &seed)); } printf("\n");
	
	printf("weibdev(1, 2): ");
	for(i = 0; i < SIZE; ++i){ printf("%lf ", weibdev(1, 2, &seed)); } printf("\n");
}

void testDiscDist()
{
	double p;
	int x1, x2;

	p = binomialp(5, 10, .5); binomialv(p, 10, .5, &x1, &x2);
	printf("binomialp(5, 10, .5): %.2lf     binomialv(%.2lf, 10, .5): %d, %d\n",
		p, p, x1, x2);

	p = geomp(5, .5); geomv(p, .5, &x1, &x2);
	printf("geomp(5, 10, .5): %.2lf     geomv(%.2lf, 10, .5): %d, %d\n",
		p, p, x1, x2);

	p = hyperp(5, 10, 20, 30); hyperv(p, 10, 20, 30, &x1, &x2);
	printf("hyperp(5, 10, .5): %.2lf     hyperv(%.2lf, 10, .5): %d, %d\n",
		p, p, x1, x2);

	p = negbnlp(5, 10, .5); negbnlv(p, 10, .5, &x1, &x2);
	printf("negbnlp(5, 10, .5): %.2lf     negbnlv(%.2lf, 10, .5): %d, %d\n",
		p, p, x1, x2);

	p = poissonp(5, 10); poissonv(p, 10, &x1, &x2);
	printf("poissonp(5, 10, .5): %.2lf     poissonv(%.2lf, 10, .5): %d, %d\n",
		p, p, x1, x2);
}

void testEigenSys()
{
	double** a;
	double d[SIZE];
	double* e;
	int i, j;
	double aij;

	a = dmatrix(0, SIZE - 1, 0, SIZE - 1);
	e = dvector(0, SIZE - 1);

	for(i = 0; i < SIZE; ++i){
		for(j = i + 1; j < SIZE; ++j){
			aij = randomValue(SIZE);
			a[i][j] = aij;
			a[j][i] = aij;
		}
		a[i][i] = randomValue(SIZE);
	}

	printf("a:\n"); outputMatrix(a, SIZE, SIZE);
	tridred(a, SIZE, d, e);
	printf("tridred(a):\n");
	printf("Diagonal: "); outputArray(d, SIZE);
	printf("Off Diag: "); outputArray(e, SIZE);

	trideig(d, d + SIZE, e);
	printf("trideig(d, e):\n");
	printf("Eigenvalues: "); outputArray(d, SIZE);

	free_dmatrix(a, 0, SIZE - 1, 0);
	free_dvector(e, 0);
}

void testLinAlg()
{
	double** a;
	double** aa;
	double p[SIZE];
	int i, j;
	double aij;

	a = dmatrix(0, SIZE - 1, 0, SIZE - 1);
	aa = dmatrix(0, SIZE - 1, 0, SIZE - 1);

	for(i = 0; i < SIZE; ++i){
		for(j = i + 1; j < SIZE; ++j){
			aij = randomValue(SIZE);
			aa[i][j] = aij;
			aa[j][i] = aij;
		}
		aa[i][i] = randomValue(SIZE);
		p[i] = randomValue(SIZE);
	}

	copyMatrix(aa, SIZE, SIZE, a);
	printf("a:\n"); outputMatrix(a, SIZE, SIZE);
	printf("determinant(a): %lf", determinant(a, SIZE));

	inverse(a, SIZE);
	printf("\ninverse(a):\n"); outputMatrix(a, SIZE, SIZE);
	
	copyMatrix(aa, SIZE, SIZE, a);
	printf("linsolve(a, b): ");
	printf("\nb: "); outputArray(p, SIZE);
	linsolve(a, p, SIZE, LINSOLVE_LU);
	printf("x: "); outputArray(p, SIZE);

	free_dmatrix(a, 0, SIZE - 1, 0);
	free_dmatrix(aa, 0, SIZE - 1, 0);
}

void testRanking()
{
	double x[SIZE];
	double r[SIZE];
	double* cp;
	int indx[SIZE];

	cp = dvector(0, SIZE - 1);

	initArray(cp, cp + SIZE);

	copyArray(cp, SIZE, x);
	printf("x: "); outputArray(x, SIZE);

	index(x, x + SIZE, indx);
	printf("index(x): "); outputIArray(indx, SIZE);

	copyArray(cp, SIZE, x);
	rank(x, x + SIZE, r, RANK_MEAN);
	printf("rank(x): "); outputArray(r, SIZE);

	copyArray(cp, SIZE, x);
	arank(x, x + SIZE, r, RANK_MEAN);
	printf("arank(x): "); outputArray(r, SIZE);

	copyArray(cp, SIZE, x);
	abrank(x, x + SIZE, r, RANK_MEAN);
	printf("abrank(x): "); outputArray(r, SIZE);
}

void testResidual()
{
	double x[SIZE];
	double* cp;

	cp = dvector(0, SIZE - 1);

	initArray(cp, cp + SIZE);

	copyArray(cp, SIZE, x);
	printf("x: "); outputArray(x, SIZE);

	resid(x, x + SIZE, 3.0);
	printf("resid(x, 3): "); outputArray(x, SIZE);

	copyArray(cp, SIZE, x);
	aresid(x, x + SIZE, 3.0);
	printf("aresid(x, 3): "); outputArray(x, SIZE);

	copyArray(cp, SIZE, x);
	nresid(x, x + SIZE, 5, 2);
	printf("nresid(x, 5, 2): "); outputArray(x, SIZE);

	free_dvector(cp, 0);
}

double function1(double x)
{
	return sin(x) * cos(x) / exp(x);
}

double function2(double x)
{
	return sin(x) - cos(x) / exp(x);
}

double function2Prime(double x)
{
	return cos(x) + (cos(x) + sin(x)) / exp(x);
}

void testRootfind()
{
	double x = bisection(function1, 1.0, 2.0, NUMERICS_MAX_ERROR);
	double y = function1(x);

	printf("bisection(function1, 1, 2): %lf\n", x);
	printf("function1(%lf): %lf\n", x, y);

	x = newton(function2, function2Prime, 3.0, NUMERICS_MAX_ERROR);
	y = function2(x);
	printf("newton(function2, function2Prime, 3): %lf\n", x);
	printf("function2(%lf): %lf\n", x, y);

	x = secant(function1, 3.0, 4.0, NUMERICS_MAX_ERROR);
	printf("secant(function1, 3, 4): %lf\n", x);
	printf("function1(%lf): %lf\n", x, y);
}

void testSort()
{
	double x[SIZE * 10];
	double cpx[SIZE * 10];
	double* y;
	double* cpy;
	int s10 = SIZE * 10;

	y = dvector(0, s10 - 1);
	cpy = dvector(0, s10 - 1);

	initArray(cpx, cpx + s10);
	initArray(cpy, cpy + s10);

	copyArray(cpx, s10, x);
	copyArray(cpy, s10, y);
	printf("x: "); outputArray(x, s10);
	printf("y: "); outputArray(y, s10);

	isort1(x, x + s10);
	printf("isort1(x): "); outputArray(x, s10);

	copyArray(cpx, s10, x);
	isort2(x, x + s10, y);
	printf("isort2(x, y):\n"); outputArray(x, s10); outputArray(y, s10);

	copyArray(cpx, s10, x);
	qsort1(x, x + s10);
	printf("qsort1(x): "); outputArray(x, s10);

	copyArray(cpx, s10, x);
	copyArray(cpy, s10, y);
	qsort2(x, x + s10, y);
	printf("qsort2(x, y):\n"); outputArray(x, s10); outputArray(y, s10);

	copyArray(cpx, s10, x);
	copyArray(cpy, s10, y);
	qsort2key(x, x + s10, y);
	printf("qsort2key(x, y):\n"); outputArray(x, s10); outputArray(y, s10);

	free_dvector(y, 0);
	free_dvector(cpy, 0);
}

int main()
{
	srand(10203040L);

	printf("Testing Algorithms:\n");
	testAlgorithm();

	printf("\nTesting Continuous Distributions:\n");
	testContDist();

	printf("\nTesting Correlation:\n");
	testCorrelation();

	printf("\nTesting Descriptive Statistics:\n");
	testDescript();

	printf("\nTesting Deviate:\n");
	testDeviate();

	printf("\nTesting Discrete Distributions:\n");
	testDiscDist();

	printf("\nTesting Eigen Systems:\n");
	testEigenSys();

	printf("\nTesting Linear Algebra:\n");
	testLinAlg();

	printf("\nTesting Ranking:\n");
	testRanking();

	printf("\nTesting Residuals:\n");
	testResidual();

	printf("\nTesting Rootfind:\n");
	testRootfind();

	printf("\nTesting Sorting:\n");
	testSort();

	return 0;
}
