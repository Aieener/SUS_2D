// cells.cpp
// 2-D Rds
// Author: Yuding Ai
// Date: 2015.06.05


#include "cells.h"
extern const string EXC_INVALID_DESC = "Not a valid description of cells!";

Cells::Cells()
{
	N0 = 1;
	N1 = 1;
	size = 1;
	//initialize my cell
	arr = new Square[1];
	arr[0] = Square();
}


Cells::Cells(int X, int Y)
{
	N0 = X; // length
	N1 = Y; // width
    size = (N0)*(N1);
    arr = new Square[size];

	//initialize my cell
	for(int i = 0; i < size; i++)
	{
		arr[i] = Square();
	}
}

/*
*  Destructor
*/
Cells::~Cells()
{   
	if(arr)
	{
		delete [] arr;
	}
}

// *** Getters *** //
int Cells::getN0() const
{
	return N0;
}
int Cells::getN1() const
{
	return N1;
}

int Cells::getSize() const
{
	return (N0)*(N1);
}


// *** Other Functionality *** //

int Cells::getIdx( int x, int y) const
{
	return  x + N0*y;
}

Square& Cells::getSquare( int x, int y) const
{
	if (x >= N0 || y >= N1)
	{
		throw EXC_INVALID_DESC;
	}

	int idx = getIdx(x,y);
	// cout << idx<<endl;
	return arr[idx];
}
