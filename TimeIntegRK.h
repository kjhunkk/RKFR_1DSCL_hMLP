#pragma once
#include "DataType.h"
#include "TimeInteg.h"

class TimeIntegRK : public TimeInteg
{
public:
	// Constructor / p.m. Equation type, flux type, correction function type, limiter type, CFL number, target time, Zone(object), Boundary(object), RK order
	TimeIntegRK(Type, Type, Type, Type, real_t, real_t, std::shared_ptr<Zone>, std::shared_ptr<Boundary>, int_t);

	// Destructor
	virtual ~TimeIntegRK();

public:
	// Functions
	// Compute time integration / p.m. Zone(object) / r.t. go/stop
	virtual bool march(std::shared_ptr<Zone>);

protected:
	// Variables
	int_t _RKorder;
	// temporary DOF for TVD-RK / RK order, DG degree, cell index
	std::vector<std::vector<std::vector<real_t> > > _temp_DOF;
};