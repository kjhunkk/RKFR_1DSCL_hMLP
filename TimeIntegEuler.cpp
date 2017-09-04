#include "TimeIntegEuler.h"

TimeIntegEuler::TimeIntegEuler(Type PDEtype, Type fluxType, Type correctFunc, Type limiterType, real_t CFL, real_t targetTime, std::shared_ptr<Zone> zone, std::shared_ptr<Boundary> bdry)
	:TimeInteg(PDEtype, fluxType, correctFunc, limiterType, CFL, targetTime, zone, bdry)
{
	_temp_DOF.resize(zone->getPolyOrder() + 1);
	for (int_t iorder = 0; iorder <= zone->getPolyOrder(); ++iorder)
		_temp_DOF[iorder].resize(zone->getGrid()->getNumCell());
}

TimeIntegEuler::~TimeIntegEuler()
{

}

bool TimeIntegEuler::march(std::shared_ptr<Zone> zone)
{
	// Apply boundary condition
	_bdry->apply(zone);

	// Marching starts
	bool procedure = true;
	if (abs(_currentTime) < epsilon) MESSAGE("Marching starts.....");

	// Calculate time step
	if ((_currentTime + _timeStep) > _targetTime)
	{
		_timeStep = _targetTime - _currentTime;
		procedure = false;
	}
	else computeTimeStep(zone);

	// Apply hMLP limiter
	std::shared_ptr<Limiter> limiter = std::make_shared<Limiter>(_limiterType, zone);
	limiter->hMLP_Limiter(zone);

	// Save previous degree of freedom
	_prev_DOF = zone->getDOF();

	// Calculate RHS
	_temp_RHS = computeRHS(zone);

	// Calculate DOF
	for (int_t idegree = 0; idegree <= zone->getPolyOrder(); ++idegree)
	{
		for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
			_temp_DOF[idegree][icell] = _prev_DOF[idegree][icell] + _timeStep*_temp_RHS[idegree][icell];
	}
	zone->setDOF(_temp_DOF);

	// Apply hMLP limiter
	limiter->hMLP_Limiter(zone);

	// Calculate solution
	zone->calSolution();

	// Update current time
	_currentTime += _timeStep;

	// Print finish condition
	if (!procedure) print();

	return procedure;
}