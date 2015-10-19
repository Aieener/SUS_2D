// square.cpp
// 2-D lattice gas
// Author: Yuding Ai
// Date: 2015.06.05

#include "square.h"

Square::Square()
{
	status = 0;
}
// *** Setters *** //
/*
*  @function         setStatus
*  @param            int newStatus - an int in the range 0 and 1
*/
void Square::setStatus(int newStatus)
{
	status = newStatus;
}


// *** Getters *** //
int Square::getStatus() const
{
	return status;
}


// *** Other Functionality *** //
/*
*  @function      isEmpty
*  @return        bool - return true iff this is Empty 
*/
bool Square::isEmpty() const
{
	if(status == EMPTY)
	{
		return true;
	}
	return false;
}
/*
*  @function      isOccupied
*  @return        bool - return true iff this is Occupied
*/
bool Square::isOccupied() const
{
	if(status == OCCUPIED)
	{
		return true;
	}
	return false;
}

/*
*  @function      getStatusChar
*  @return        char - return the char corresponding to the status, for
*                 printing the cell (Empty = '_', Occupied = '#'.)
*
*/
char Square::getStatusChar() const
{
	char sout;
	
	switch(status)
	{
		case 0:
		sout = '_';
		break;
		case 1:
		sout = '#';
	}
	return sout;
}

