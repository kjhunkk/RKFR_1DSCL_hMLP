#include "ConvPhyFlux.h"

ConvPhyFlux::ConvPhyFlux()
{
	_advSpeed = 0.0;
}

ConvPhyFlux::~ConvPhyFlux()
{

}

real_t ConvPhyFlux::phyFlux(const std::string& name, real_t u)
{
	if (name == "advection") return advection(u);

	else if (name == "burgers") return burgers(u);

	else
	{
		ERROR("cannot find physical flux");
		return 0.0;
	}
}

real_t ConvPhyFlux::diff_phyFlux(const std::string& name, real_t u)
{
	if (name == "advection") return diff_advection(u);

	else if (name == "burgers") return diff_burgers(u);

	else
	{
		ERROR("cannot find physical flux");
		return 0.0;
	}
}

real_t ConvPhyFlux::phyCharSpeed(const std::string& name, real_t u)
{
	if (name == "advection") return _advSpeed;

	else if (name == "burgers") return burgersChar(u);

	else
	{
		ERROR("cannot find physical flux");
		return 0.0;
	}
}

real_t ConvPhyFlux::burgers(real_t u)
{
	return 0.5*u*u;
}

real_t ConvPhyFlux::advection(real_t u)
{
	return _advSpeed*u;
}

real_t ConvPhyFlux::diff_burgers(real_t u)
{
	return u;
}

real_t ConvPhyFlux::diff_advection(real_t u)
{
	return _advSpeed;
}

real_t ConvPhyFlux::burgersChar(real_t u)
{
	return u;
}

real_t ConvPhyFlux::_advSpeed;
ConvPhyFlux ConvPhyFlux::_phyFlux;