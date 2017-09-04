#include "Post.h"

Post::Post(std::shared_ptr<Reader> reader)
{
	_reader = reader;
}

Post::~Post()
{

}

void Post::solution(std::shared_ptr<Zone> zone) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "RKFR_1D_P";
	fileName += zone->getPolyOrder();
	fileName += "_";
	fileName += _reader->getPDE();
	fileName += "_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += ".plt";

	// Cell size
	real_t dx = zone->getGrid()->getSizeX();

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone->getGrid()->getCell();
	std::vector<real_t> solution = zone->getDescSolution();
	std::vector<real_t> X;
	std::vector<real_t> U;
	for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			X.push_back(cell[icell]->getPosX());
			U.push_back(solution[icell]);
		}
	}

	// Write solutions
	write(fileName, "X", "Velocity", X, U);
}

void Post::solution(const std::string& name, std::shared_ptr<Zone> zone) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "RKFR_1D_P";
	fileName += std::to_string(zone->getPolyOrder());
	fileName += "_";
	fileName += _reader->getPDE();
	fileName += "_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += "_";
	fileName += name;
	fileName += ".plt";

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone->getGrid()->getCell();
	std::vector<real_t> solution = zone->getDescSolution();
	std::vector<real_t> X;
	std::vector<real_t> U;
	for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			X.push_back(cell[icell]->getPosX());
			U.push_back(solution[icell]);
		}
	}

	// Write solutions
	write(fileName, "X", "Velocity", X, U);
}

void Post::FRsolution(const std::string& name, std::shared_ptr<Zone> zone) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "RKFR_1D_P";
	fileName += std::to_string(zone->getPolyOrder());
	fileName += "_FRsolution_";
	fileName += _reader->getPDE();
	fileName += "_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += "_";
	fileName += name;
	fileName += ".plt";

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone->getGrid()->getCell();
	std::vector<real_t> X;
	std::vector<real_t> U;
	real_t sizeX = zone->getGrid()->getSizeX();
	real_t dx = zone->getGrid()->getSizeX() / double(POST_GRID_NUM);
	std::shared_ptr<FRbasis> basis = std::make_shared<FRbasis>(zone->getGrid());

	for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			real_t posX = cell[icell]->getPosX();
			for (int_t idegree = 0; idegree <= zone->getPolyOrder(); ++idegree)
			{
				X.push_back(posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone->getPolyOrder()));
				U.push_back(zone->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone->getPolyOrder())));
			}
		}
	}

	// Write solutions
	write(fileName, "X", "Velocity", X, U);
}

void Post::error(std::shared_ptr<Zone> zone) const
{

}

void Post::write(const std::string& fileName, const std::string& VN1, const std::string& VN2, std::vector<real_t> V1, std::vector<real_t> V2) const
{
	std::ofstream file;
	file.open(fileName, std::ios::trunc);
	int_t _num_element = V1.size();
	if (_num_element != V2.size()) ERROR("Number of elements does not match");

	if (file.is_open())
	{
		MESSAGE("Output file open");
		file << "variables = " << VN1 << ", " << VN2 << "\n";
		file << "zone t = \"RKFR 1D\", i=" << _num_element << ", f=point\n";
		for (int_t ielem = 0; ielem < _num_element; ++ielem)
		{
			file << std::to_string(V1[ielem]) << "\t" << std::to_string(V2[ielem]) << "\n";
		}
		file.close();
	}
	else ERROR("cannot open output file");
}