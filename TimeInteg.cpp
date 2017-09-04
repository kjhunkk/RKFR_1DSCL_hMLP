#include "TimeInteg.h"

TimeInteg::TimeInteg(Type PDEtype, Type fluxType, Type correctFunc, Type limiterType, real_t CFL, real_t targetTime, std::shared_ptr<Zone> zone, std::shared_ptr<Boundary> bdry)
{
	// Initializing variables
	_PDEtype = PDEtype; _fluxType = fluxType;
	_correctFunc = correctFunc; _limiterType = limiterType;
	_CFL = CFL; _targetTime = targetTime;
	_currentTime = 0.0; _timeStep = 0.0;

	// Initializing objects
	_zone = zone; _bdry = bdry;
	_basis = std::make_shared<FRbasis>(zone->getGrid());

	// Initializing temporary variables
	_godFlux = std::make_shared<ConvFluxGodunov>(PDEtype, correctFunc, zone);
	_temp_solution.resize(zone->getGrid()->getNumCell());
	_prev_DOF.resize(zone->getPolyOrder() + 1);
	_temp_RHS.resize(zone->getPolyOrder() + 1);
	for (int_t iorder = 0; iorder <= zone->getPolyOrder(); ++iorder)
	{
		_prev_DOF[iorder].resize(zone->getGrid()->getNumCell());
		_temp_RHS[iorder].resize(zone->getGrid()->getNumCell());
	}
}

TimeInteg::~TimeInteg()
{

}

std::vector<std::vector<real_t> > TimeInteg::computeRHS(std::shared_ptr<Zone> zone) const
{
	real_t sizeX = zone->getGrid()->getSizeX();
	real_t inv_sizeX = 1.0/(zone->getGrid()->getSizeX());
	int_t num_cell = zone->getGrid()->getNumCell();
	int_t polyOrder = zone->getPolyOrder();

	// Temporary degree of freedom
	std::vector<std::vector<real_t> > RHS;
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();
	RHS.resize(polyOrder + 1);
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		RHS[iorder].resize(num_cell);

	// Initializing RHS
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		for (int_t icell = 0; icell < num_cell; ++icell)
			RHS[iorder][icell] = 0.0;
	
	// FR flux
	std::vector<std::vector<real_t> > diff_flux;
	diff_flux.resize(polyOrder + 1);
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		diff_flux[iorder].resize(num_cell);
	real_t back_flux, for_flux; // computational flux
	real_t left_face, right_face; // face coordinate
	real_t diff_disc_flux; // derivative of discontinuous flux
	std::shared_ptr<FRbasis> basis = std::make_shared<FRbasis>(_zone->getGrid());

	for (int_t icell = GHOST; icell <= num_cell - GHOST; ++icell)
	{
		for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		{
			// Compute face coordinates
			left_face = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*sizeX;
			right_face = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*sizeX;

			// Computational flux
			back_flux = _godFlux->computeFlux(zone->getPolySolution(icell - 1, left_face), zone->getPolySolution(icell, left_face));
			for_flux = _godFlux->computeFlux(zone->getPolySolution(icell, right_face), zone->getPolySolution(icell + 1, right_face));
			
			// Compute derivative of discontinuous flux
			diff_disc_flux = 0.0;
			for (int_t iorder2 = 0; iorder2 <= polyOrder; ++iorder2)
				diff_disc_flux += ConvPhyFlux::phyFlux(_PDEtype, temp_DOF[iorder2][icell])*basis->diff_LagrangeP(basis->getSolutionPoint(iorder, polyOrder), iorder2, polyOrder);
			
			// Compute derivative of FR flux
			diff_flux[iorder][icell] = diff_disc_flux
				+ (back_flux - ConvPhyFlux::phyFlux(_PDEtype, zone->getPolySolution(icell, left_face)))*_godFlux->DcorrectFunc(false, basis->getSolutionPoint(iorder, polyOrder))
				+ (for_flux - ConvPhyFlux::phyFlux(_PDEtype, zone->getPolySolution(icell, right_face)))*_godFlux->DcorrectFunc(true, basis->getSolutionPoint(iorder, polyOrder));

			// Coordinate conversion from computational to physical
			diff_flux[iorder][icell] *= 2.0*inv_sizeX;
		}
	}

	// Calculate RHS
	for (int_t icell = GHOST; icell < num_cell - GHOST; ++icell)
		for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
			RHS[iorder][icell] = -diff_flux[iorder][icell];
	
	return RHS;
}

void TimeInteg::computeTimeStep(std::shared_ptr<Zone> zone)
{
	if (_PDEtype == "advection")
		_timeStep = _CFL*zone->getGrid()->getSizeX() / abs(GET_SPEED) / double(2 * zone->getPolyOrder() + 1);

	else if (_PDEtype == "burgers")
	{
		real_t temp_sol1;
		real_t temp_sol2;
		std::vector<real_t> shockSpeed(zone->getGrid()->getNumCell() - 1, 0.0);
		for (int_t icell = 0; icell < shockSpeed.size(); ++icell)
		{
			temp_sol1 = zone->getDescSolution()[icell];
			temp_sol2 = zone->getDescSolution()[icell + 1];
			if (temp_sol1 >= temp_sol2) shockSpeed[icell] = 0.5*abs(temp_sol1 + temp_sol2);
			else shockSpeed[icell] = std::max(abs(temp_sol1), abs(temp_sol2));
		}
		_timeStep = _CFL*zone->getGrid()->getSizeX() / *std::max_element(shockSpeed.begin(), shockSpeed.end()) / double(2 * zone->getPolyOrder() + 1);
	}
}

void TimeInteg::print() const
{
	MESSAGE("Marching finished.....");
	MESSAGE("CFL = " + std::to_string(_CFL));
	MESSAGE("Target time = " + std::to_string(_targetTime));
	MESSAGE("Final time step = " + std::to_string(_timeStep));
	MESSAGE("Current time = " + std::to_string(_currentTime));
}