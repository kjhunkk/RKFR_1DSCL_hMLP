#include "InitialCondition.h"

InitialCondition::InitialCondition(Type type)
{
	_type = type;
}

InitialCondition::~InitialCondition()
{

}

real_t InitialCondition::initializer(real_t x) const
{
	if (_type == "square")
		return square(x);

	if (_type == "halfdome")
		return halfdome(x);

	if (_type == "gauss")
		return gauss(x);

	if (_type == "shock")
		return shock(x);

	if (_type == "expansion")
		return expansion(x);

	if (_type == "sine")
		return sine(x);

	if (_type == "benchmark1")
		return benchmark1(x);

	if (_type == "constant")
		return 1.0;

	if (_type == "benchmark2")
		return benchmark2(x);

	ERROR("cannot find initial condition");

	return 0.0;
}

real_t InitialCondition::square(real_t x) const
{
	if ((x <= 0.5) && (x >= -0.5)) return 1.0;
	else return 0.0;
}

real_t InitialCondition::halfdome(real_t x) const
{
	if ((x <= sqrt(0.3)) && (x >= -sqrt(0.3))) return sqrt(1 - 10.0*CONST13*x*x);
	else return 0.0;
}

real_t InitialCondition::gauss(real_t x) const
{
	return exp(-300.0*x*x);
}

real_t InitialCondition::shock(real_t x) const
{
	if (x <= -0.3) return 1.2;
	else return 0.4;
}

real_t InitialCondition::expansion(real_t x) const
{
	if (x <= -0.5) return 0.0;
	else return 1.0;
}

real_t InitialCondition::sine(real_t x) const
{
	return (1 + 0.5*sin(M_PI*x));
}

real_t InitialCondition::benchmark1(real_t x) const
{
	return (0.25 + 0.5*sin(M_PI*x));
}

real_t InitialCondition::benchmark2(real_t x) const
{
	return (0.5 + sin(2.0*M_PI*x));
}