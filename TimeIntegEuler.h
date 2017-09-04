#pragma once
#include "DataType.h"
#include "TimeInteg.h"

class TimeIntegEuler : public TimeInteg
{
public:
	// Constructor / p.m. Equation type, flux type, correction function type, limiter type, CFL number, target time, Zone(object), Boundary(object)
	TimeIntegEuler(Type, Type, Type, Type, real_t, real_t, std::shared_ptr<Zone>, std::shared_ptr<Boundary>);

	// Destructor
	virtual ~TimeIntegEuler();

public:
	// Functions
	// Compute time integration / p.m. Zone(object) / r.t. go/stop
	virtual bool march(std::shared_ptr<Zone>);	

protected:
	// Variables
	std::vector<std::vector<real_t> > _temp_DOF;
};