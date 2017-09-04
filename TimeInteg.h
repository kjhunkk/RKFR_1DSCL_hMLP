#pragma once
#include "DataType.h"
#include "Zone.h"
#include "ConvFluxGodunov.h"
#include "Boundary.h"
#include "Limiter.h"
#include "FRbasis.h"
#include "ConvPhyFlux.h"

class TimeInteg
{
public:
	// Constructor / p.m. Equation type, Flux type, Correction function type, limiter type, CFL number, target time, Zone(object), Boundary(object)
	TimeInteg(Type, Type, Type, Type, real_t, real_t, std::shared_ptr<Zone>, std::shared_ptr<Boundary>);

	// Destructor
	virtual ~TimeInteg();

public:
	// Functions
	inline void reset() { _currentTime = 0.0; }

	inline void reset(real_t target) { _currentTime = 0.0; _targetTime = target; }

	inline real_t getTime() const { return _currentTime; }

	inline real_t getTimeStep() const { return _timeStep; }

	inline real_t getTargetTime() const { return _targetTime; }

	// Compute time integration / p.m. Zone(object) / r.t. go/stop
	virtual bool march(std::shared_ptr<Zone>) = 0;

protected:
	// Variables
	std::vector<real_t> _temp_solution;
	std::vector<std::vector<real_t> > _prev_DOF;
	std::vector<std::vector<real_t> > _temp_RHS;
	std::shared_ptr<Zone> _zone;
	std::shared_ptr<ConvFluxGodunov> _godFlux;
	std::shared_ptr<Boundary> _bdry;
	std::shared_ptr<FRbasis> _basis;
	Type _PDEtype;
	Type _fluxType;
	Type _correctFunc;
	Type _limiterType;
	real_t _CFL;
	real_t _currentTime;
	real_t _timeStep;
	real_t _targetTime;

protected:
	// Functions
	// Compute right hand side / p.m. Zone to compute
	std::vector<std::vector<real_t> > computeRHS(std::shared_ptr<Zone>) const;

	// Compute time step / p.m. Zone(object)
	void computeTimeStep(std::shared_ptr<Zone>);

	// Print time variables
	void print() const;
};