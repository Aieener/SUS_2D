// cells.h	
// 2-D Rods
// Author: Yuding Ai
// Date: 2015.06.05

#ifndef CELLS_H
#define CELLS_H
#include <vector>
#include "hardrods.h"
#include "square.h"
using namespace std;


class Cells
{
	private:
		int N0;         // Weith; 
		int N1;         // Length; 
		int size;       // the total size of the grid
		Square * arr;   // linearized array of Squares

	public:
		Cells();

		/*
		*  Constructor builds a cells from an input stream
		*  @param      int r - row index
		*  @param      int c - column index
		*  initialize all states to 0;
		*/
		Cells(int x,int y);

		/*
		*  Destructor
		*/
		~Cells();

		// *** Getters *** //
		int getN0() const;
		int getN1() const;
		int getSize() const;
		int getIdx( int x, int y) const;



		// *** Other Functionality *** //
		/*
		*  @function      getSquare
		*  @param         int i - row index, in range 0 through numRows - 1
		*  @param         int j - column index, in range 0 through numCols - 1
		*  @return        Square& - return a reference to the (i,j)-Square stored
		*                 in this cells
		*/
		Square& getSquare(int x, int y) const;

	};


#endif /* CELLS_H */