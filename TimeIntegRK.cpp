#include "TimeIntegRK.h"

TimeIntegRK::TimeIntegRK(Type PDEtype, Type fluxType, Type correctFunc, Type limiterType, real_t CFL, real_t targetTime, std::shared_ptr<Zone> zone, std::shared_ptr<Boundary> bdry, int_t RKorder)
	:TimeInteg(PDEtype, fluxType, correctFunc, limiterType, CFL, targetTime, zone, bdry)
{
	_RKorder = RKorder;
	_temp_DOF.resize(RKorder);
	for (int_t iorder = 0; iorder < RKorder; ++iorder)
	{
		_temp_DOF[iorder].resize(zone->getPolyOrder() + 1);
		for (int_t idegree = 0; idegree <= zone->getPolyOrder(); ++idegree)
			_temp_DOF[iorder][idegree].resize(zone->getGrid()->getNumCell());
	}
}

TimeIntegRK::~TimeIntegRK()
{

}

bool TimeIntegRK::march(std::shared_ptr<Zone> zone)
{
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

	// Declare temprorary Zone for TVD Runge-Kutta time integration
	std::shared_ptr<Zone> temp_zone = std::make_shared<Zone>(*zone);

	// Declare local projection limiter object
	std::shared_ptr<Limiter> limiter = std::make_shared<Limiter>(_limiterType, zone);

	// ----------------------First step--------------------------
	// Apply boundary condition
	_bdry->apply(temp_zone);

	// Apply hMLP limiter
	limiter->hMLP_Limiter(temp_zone);

	// Save previous degree of freedom
	_prev_DOF = temp_zone->getDOF();

	// Calculate RHS
	_temp_RHS = computeRHS(temp_zone);

	// Calculate DOF
	for (int_t idegree = 0; idegree <= zone->getPolyOrder(); ++idegree)
	{
		for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
			_temp_DOF[0][idegree][icell] = _prev_DOF[idegree][icell] + _timeStep*_temp_RHS[idegree][icell];
	}

	// Update temporary Zone object
	temp_zone->setDOF(_temp_DOF[0]);
	temp_zone->calSolution();

	// ---------------------Second step---------------------------
	// Apply boundary condition
	_bdry->apply(temp_zone);

	// Apply hMLP limiter
	limiter->hMLP_Limiter(temp_zone);

	// Calculate RHS
	_temp_RHS = computeRHS(temp_zone);

	// Calculate DOF
	for (int_t idegree = 0; idegree <= zone->getPolyOrder(); ++idegree)
	{
		for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
			_temp_DOF[1][idegree][icell] = 0.75*_prev_DOF[idegree][icell] + 0.25*(_temp_DOF[0][idegree][icell] + _timeStep*_temp_RHS[idegree][icell]);
	}

	// Update temporary Zone object
	temp_zone->setDOF(_temp_DOF[1]);
	temp_zone->calSolution();

	// ----------------------Third step---------------------------
	// Apply boundary condition
	_bdry->apply(temp_zone);

	// Apply hMLP limiter
	limiter->hMLP_Limiter(temp_zone);

	// Calculate RHS
	_temp_RHS = computeRHS(temp_zone);

	// Calculate DOF
	for (int_t idegree = 0; idegree <= zone->getPolyOrder(); ++idegree)
	{
		for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
			_temp_DOF[2][idegree][icell] = CONST13*_prev_DOF[idegree][icell] + CONST23*(_temp_DOF[1][idegree][icell] + _timeStep*_temp_RHS[idegree][icell]);
	}

	// Update solution zone
	zone->setDOF(_temp_DOF[2]);

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