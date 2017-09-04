#pragma once
#include "DataType.h"

class ConvPhyFlux
{
protected:
	ConvPhyFlux();
	~ConvPhyFlux();

public:
	// Set advection speed
	static void setAdvSpeed(real_t advSpeed) { _advSpeed = advSpeed; }

	static real_t getAdvSpeed() { return _advSpeed; }

	// Calculate Physical flux / p.m. physical flux name, solution variable
	static real_t phyFlux(const std::string&, real_t);

	// Calculate derivative of physical flux / p.m. physical flux name, solution variable
	static real_t diff_phyFlux(const std::string&, real_t);

	// Compute characteristic speed / p.m. physical flux name, solution variable
	static real_t phyCharSpeed(const std::string&, real_t);

protected:
	// Physical fluxes
	static real_t burgers(real_t);
	static real_t advection(real_t);

	// Derivative of physical fluxes
	static real_t diff_burgers(real_t);
	static real_t diff_advection(real_t);

	// Physical characteristic speed
	static real_t burgersChar(real_t);

	// Advection speed
	static real_t _advSpeed;

private:
	// ConvPhyFlux variable
	static ConvPhyFlux _phyFlux;
};

// Convective Physical Flux macro
#define PHY_FLUX(name, u) ConvPhyFlux::phyFlux(name, u)

#define PHY_FLUX_BURGERS(u) ConvPhyFlux::phyFlux("burgers", u)

#define PHY_FLUX_ADVEC(u) ConvPhyFlux::phyFlux("advection", u)

#define PHY_CHAR_SPEED(name, u) ConvPhyFlux::phyCharSpeed(name, u)

#define SET_SPEED(a) ConvPhyFlux::setAdvSpeed(a)

#define GET_SPEED ConvPhyFlux::getAdvSpeed()