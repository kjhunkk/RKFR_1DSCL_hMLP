#include "Grid.h"

Grid::Grid(real_t area, real_t sizeX)
{
	_area = area;
	_sizeX = sizeX;
	_num_cell = area / sizeX;
	_num_cell += 2 * GHOST;
	_cell.resize(_num_cell);
	for (int_t icell = 0; icell < GHOST; ++icell)
	{
		_cell[icell] = std::make_shared<Cell>();
		_cell[_num_cell - icell - 1] = std::make_shared<Cell>();
		real_t coord_X = -0.5*_area + sizeX*(0.5 + (icell - GHOST));
		_cell[icell]->initialize(false, sizeX, coord_X);
		_cell[_num_cell - icell - 1]->initialize(false, sizeX, -coord_X);
	}
	for (int_t icell = GHOST; icell < _num_cell - GHOST; ++icell)
	{
		real_t coord_X = -0.5*_area + sizeX*(0.5 + (icell - GHOST));
		_cell[icell] = std::make_shared<Cell>();
		_cell[icell]->initialize(true, sizeX, coord_X);
	}
}

Grid::~Grid()
{

}