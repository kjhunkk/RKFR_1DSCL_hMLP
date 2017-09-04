#pragma once
#include "DataType.h"
#include "ConvFlux.h"
#include "ConvPhyFlux.h"

class ConvFluxGodunov : public ConvFlux
{
public:
	// Constructor / p.m. physical flux type, correction function type, Zone(Object)
	ConvFluxGodunov(Type, Type, std::shared_ptr<Zone>);

	// Destructor
	virtual ~ConvFluxGodunov();

public:
	// Functions
	// Compute flux / p.m. begin, end
	virtual real_t computeFlux(real_t, real_t) const;
};