#include "DataType.h"
#include "Reader.h"
#include "Grid.h"
#include "Zone.h"
#include "InitialCondition.h"
#include "TimeIntegEuler.h"
#include "TimeIntegRK.h"
#include "Post.h"
#include "Boundary.h"
#include "OrderTest.h"

// Modified 2017-09-01
// by Juhyeon Kim
// Caution : polynomial order starts from 0

int main()
{
	// Read input file
	std::shared_ptr<Reader> reader = std::make_shared<Reader>();
	if (!reader->readFile("./input.inp")) return 0;

	// Set advection speed
	if (reader->getPDE() == "advection") SET_SPEED(reader->getAdvSpeed());

	// Initializing objects
	std::shared_ptr<Post> post = std::make_shared<Post>(reader);
	std::shared_ptr<Grid> grid = std::make_shared<Grid>(reader->getArea(), reader->getSizeX());
	std::shared_ptr<Zone> zone = std::make_shared<Zone>(grid, reader->getPolyOrder());
	std::shared_ptr<OrderTest> orderTest = std::make_shared<OrderTest>();

	// Initializing solution domain
	zone->initialize(std::make_shared<InitialCondition>(reader->getInitial()));

	// initializing boundary condition
	std::shared_ptr<Boundary> bdry = std::make_shared<Boundary>(reader->getBoundary(), zone);

	// Save exact solution for order test
	std::vector<real_t> exact = orderTest->ZoneToPoly(zone);
	orderTest->setExact(exact);

	// Post initial condition
	post->solution("initial", zone);
	if (reader->getPolyOrder() > 0) post->FRsolution("initial", zone);

	// Initialzing time integrator
	std::shared_ptr<TimeInteg> timeInteg;

	if (reader->getTimeInteg() == "Euler")
		timeInteg = std::make_shared<TimeIntegEuler>
		(reader->getPDE(), reader->getFluxScheme(), reader->getCorrectFunc(), reader->getLimiter(), reader->getCFL(), reader->getTargetT(), zone, bdry);
	else if (reader->getTimeInteg() == "RK3")
		timeInteg = std::make_shared<TimeIntegRK>
		(reader->getPDE(), reader->getFluxScheme(), reader->getCorrectFunc(), reader->getLimiter(), reader->getCFL(), reader->getTargetT(), zone, bdry, 3);
	else ERROR("cannot find time integrator");

	// Print initialized solution
	// zone->print();

	int_t iter = 0;
	// Time marching
	while (timeInteg->march(zone))
	{
		iter++;
		if (iter % 100 == 0) MESSAGE("Iteration = " + std::to_string(iter));
		if (reader->getPolyOrder() > 0) post->FRsolution("result" + std::to_string(iter), zone);
		MESSAGE("current time = " + std::to_string(timeInteg->getTime()) + "     time step    = " + std::to_string(timeInteg->getTimeStep()));
	}

	// Computed solution array
	std::vector<real_t> computed = orderTest->ZoneToPoly(zone);

	// Compute L errors
	real_t L1 = orderTest->L1error(computed);
	real_t L2 = orderTest->L2error(computed);
	real_t Linf = orderTest->Linf_error(computed);

	// Print L errors
	std::cout << "L1 error    = " << L1 << "\n";
	std::cout << "L2 error    = " << L2 << "\n";
	std::cout << "L inf error = " << Linf << "\n";

	// Post solution
	post->solution("result", zone);
	if(reader->getPolyOrder() > 0) post->FRsolution("result", zone);

	return 0;
}