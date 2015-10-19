//histogram.h
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "square.h"
#include "cells.h"
// #include "MC.h"
#include "hardrods.h"
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <vector>
#include <utility> 
// #include "python.h"
using namespace std;

class Histogram {

	private:
		double binWidth,min,max;
		int binCount;
		int lowerOutlierCount, upperOutlierCount;
		int * counts;
	
	public:
	Histogram(double min, double max, double binWidth);

	void record(double datum);
	int width() const; // Get the width;
	int bins() const;                    // Get the number of bins
	int count(int bin) const;            // Get the number of data points in some bin
	int countLowerOutliers() const;      // Get count of values below the minimum
	int countUpperOutliers() const;      // Get count of values at or above the maximum
	double lowerBound(int bin) const;    // Get the inclusive lower bound of a bin
	double upperBound(int bin) const;    // Get the exclusive upper bound of a bin
    pair<int,double> Minave() const; // get the index of minimum #of data points in some bin
	int Getmax() const; // get the index of maximum #of data points in some bin
	void plot(int a) const; // plot the hisdogram;
	void reset();

    ~Histogram();
};
