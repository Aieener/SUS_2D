//histogram.cpp

#include "histogram.h"


Histogram::Histogram(double Min, double Max, double bW)
{
	min = Min;
	max = Max;
	binWidth = bW;
	binCount = int((max - min)/binWidth);
	lowerOutlierCount=0;
	upperOutlierCount=0;
	counts = new int [binCount];
	for(int i = 0; i < bins(); ++i)
	{
		counts[i] = 0;
	}
}
Histogram::~Histogram()
{
	delete[] counts;
}

int Histogram::width() const // Get the width of bins
{
	return binWidth;
} 

int Histogram::bins() const   // Get the total number of bins
{
	return binCount;	
}                  

int Histogram::count(int bin) const // Get the number of data points in some bin
{
	return counts[bin];
}            

int Histogram::countLowerOutliers() const // Get count of values below the minimum
{
	return lowerOutlierCount;
}
      
int Histogram::countUpperOutliers() const // Get count of values at or above the maximum
{
	return upperOutlierCount;
}
      
double Histogram::lowerBound(int bin) const // Get the inclusive lower bound of a bin
{
	return double(bin) * binWidth;
}   

double Histogram::upperBound(int bin) const // Get the exclusive upper bound of a bin
{
	return double(bin+1) * binWidth;
}

pair<int,double> Histogram::Minave() const
{
	double acc = 0;
	double ave ;
	double a = count(0);
	pair <int,double> minave ;

	for (int i = 0; i< bins(); i++)
	{
		acc = acc + count(i); 
		if (count(i) < a)
		{
			a = count(i);
		}
	}
	ave = acc/bins();
	minave = make_pair(a,ave);

	return minave;
}

int Histogram::Getmax() const
{
	double a = count(0);
	for (int i = 0; i< bins(); i++)
	{
		if (count(i) > a)
		{
			a = count(i);
		}
	}

	return a;
}      

void Histogram::reset()
{
	for(int i = 0; i < bins(); ++i)
	{
		counts[i] = 0;
	}
}

void Histogram::record(double datum) {
	stringstream st;
	string s;
    int bin;
    bin = int((datum - min) / binWidth);

    if (bin < 0) {
        lowerOutlierCount++;
    } else if (bin >= binCount) {
        upperOutlierCount++;
    } else {
        counts[bin]++;
    }
}

void Histogram::plot(int a) const
{
	stringstream sh;
	for(int i = 0; i< binCount; i++)
	{
		sh<< i*binWidth + 1*binWidth<< "  "<< count(i) << endl;
		// sh<< i*binWidth << "  "<< count(i) << endl;
	}
	ofstream myfile ("his.txt");
	string data2 = sh.str();
	myfile << data2;
	myfile.close();
	FILE* gnuplot = popen("gnuplot -persistent","w");
	fprintf(gnuplot, "set terminal png\n set term png truecolor\n set output \"Histo.png\" \n");
	if (a == 0)
	{
		fprintf(gnuplot, "set title \"HistoN\"\n set xlabel \"N\"\n set ylabel \"Freq(N)\"\n");
	}
	else
	{
		fprintf(gnuplot, "set title \"HistoN%d \"\n set xlabel \"N\"\n set ylabel \"Freq(N)\"\n",a);
	}
	fprintf(gnuplot, "plot \'his.txt\' using 1:2 smooth freq with boxes title \"N\" \n");
	fflush(gnuplot);
	pclose(gnuplot);
}
