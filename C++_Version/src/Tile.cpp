/***********************************************************************
 * Tile.cpp
 * Author: jam5pi
 *   Date: Nov 9, 2009
 *
 * <software purpose to be completed>
 ***********************************************************************
 * Revision History
 *
 ***********************************************************************
 * Copyright 2009, Children's Hospital of Cincinnati
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include "Tile.h"

Tile::Tile() {
	// TODO Auto-generated constructor stub

}
Tile::Tile(AMOS::Tile_t aTile)
	: src(aTile.source), off(aTile.offset), clr1(aTile.range.begin), clr2(aTile.range.end)
{

}


Tile::~Tile() {
	// TODO Auto-generated destructor stub
}

string Tile::getAF() {
	stringstream ss;
	if (isForwardRead()) {
		ss << "U " << off;
	} else {
		ss << "C " << off;
	}
	return ss.str();
}
