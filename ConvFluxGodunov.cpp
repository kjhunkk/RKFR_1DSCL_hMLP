#include "ConvFluxGodunov.h"

ConvFluxGodunov::ConvFluxGodunov(Type phyFlux, Type correctFunc, std::shared_ptr<Zone> zone)
	: ConvFlux(phyFlux, correctFunc, zone)
{

}

ConvFluxGodunov::~ConvFluxGodunov()
{

}

real_t ConvFluxGodunov::computeFlux(real_t begin, real_t end) const
{
	if (_phyFlux == "advection")
	{
		if (GET_SPEED > 0.0) return PHY_FLUX_ADVEC(begin);
		else return PHY_FLUX_ADVEC(end);
	}

	else if (_phyFlux == "burgers")
	{
		if ((begin - end) > epsilon)
		{
			real_t speed = 0.5*(begin + end);
			if (speed >= 0.0) return PHY_FLUX_BURGERS(begin);
			else return PHY_FLUX_BURGERS(end);
		}

		if ((begin - end) < -epsilon)
		{
			if (begin > 0.0) return PHY_FLUX_BURGERS(begin);
			else if (end < 0.0) return PHY_FLUX_BURGERS(end);
			else return 0.0;
		}

		else return PHY_FLUX_BURGERS(begin);
	}

	else
	{
		ERROR("cannot fine flux");
		return 0.0;
	}
}