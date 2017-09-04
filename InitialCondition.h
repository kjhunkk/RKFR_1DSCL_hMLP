#pragma once
#include "DataType.h"

class InitialCondition
{
public:
	// Constructor / p.m. initial condition type
	InitialCondition(Type);

	// Destructor
	~InitialCondition();

public:
	// Functions
	inline Type getType() { return _type; }

	// calculate initial condition
	real_t initializer(real_t) const;

protected:
	// Variables
	Type _type;

protected:
	// Functions
	// Initial condition functions
	real_t square(real_t) const;
	real_t halfdome(real_t) const;
	real_t gauss(real_t) const;
	real_t shock(real_t) const;
	real_t expansion(real_t) const;
	real_t sine(real_t) const;
	real_t benchmark1(real_t) const;
	real_t benchmark2(real_t) const;
};