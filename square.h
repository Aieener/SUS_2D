// square.h
// 2-D lattice gas
// Author: Yuding Ai
// Date: 2015.06.05

#ifndef SQUARE_H
#define SQUARE_H

#include <iostream>
#include <string>

class Square
{
	private:
		// instance variables
		int status;

	public:

		// class constants
		enum STATUS {EMPTY,OCCUPIED};
		// static const int EMPTY = 0; // char for printing : _
		// static const int OCCUPIED = 1;  // char for printing : # 

		/*!
		*  Default constructor
		*/   
	    Square();


		// *** Setters *** //
		void setStatus(int newStatus);

		// *** Getters *** //
		int getStatus() const;

		// ***  Other Functionality ***//

		/*
		*  @function      isEmpty
		*  @return        bool - return true iff this is Empty 
		*/
		bool isEmpty() const;

		/*
		*  @function      isOccupied
		*  @return        bool - return true iff this is Occupied
		*/
		bool isOccupied() const;

		/*
		*  @function      getStatusChar
		*  @return        char - return the char corresponding to the status, for
		*                 printing the cell (Empty = '_', Occupied = '#'.)
		*/
		char getStatusChar() const;

};

#endif /* SQUARE_H */
