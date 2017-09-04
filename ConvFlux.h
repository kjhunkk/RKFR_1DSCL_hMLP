#pragma once
#include "DataType.h"
#include "Zone.h"

class ConvFlux
{
public:
	// Constructor / p.m. physical flux type, correction function type, Zone(Object)
	ConvFlux(Type, Type, std::shared_ptr<Zone>);

	// Destructor
	virtual ~ConvFlux();

public:
	// Functions
	// Set Zone
	void setZone(std::shared_ptr<Zone> zone) { _zone = zone; }

	// Compute flux / p.m. begin, end
	virtual real_t computeFlux(real_t, real_t) const = 0;

	// Compute correction function / p.m. Right(true)/Left(false), computational coordinate
	real_t correctFunc(bool, real_t) const;

	// Compute derivative of correction function / p.m. Right(true)/Left(false), computational coordinate
	real_t DcorrectFunc(bool, real_t) const;

protected:
	// Variables
	std::shared_ptr<Zone> _zone;
	Type _phyFlux;
	Type _correctFunc;

protected:
	// Functions
	// Left boundary DG correction function / p.m. computational coordinate
	real_t correctDG(bool, real_t) const;

	// Left boundary Derivative of DG correction function / p.m. computational coordinate
	real_t DcorrectDG(bool, real_t) const;

	// Left boundary SG correction function / p.m. computational coordinate
	real_t correctSG(bool, real_t) const;

	// Left boundary Derivative of SG correction function / p.m. computational coordinate
	real_t DcorrectSG(bool, real_t) const;

	// Left boundary Lumping for Lobatto points correction function / p.m. computational coordinate
	real_t correct2(bool, real_t) const;

	// Left boundary Derivative of lumping for Lobatto points correction function / p.m. computational coordinate
	real_t Dcorrect2(bool, real_t) const;

	// Left boundary Legendre-Lobatto correction function / p.m. computational coordinate
	real_t correctLo(bool, real_t) const;

	// Left boundary Derivative of Legendre-Lobatto correction function / p.m. computational coordinate
	real_t DcorrectLo(bool, real_t) const;

	// Left boundary Gauss correction function / p.m. computational coordinate
	real_t correctGa(bool, real_t) const;

	// Left boundary Derivative of Gauss correction function / p.m. computational coordinate
	real_t DcorrectGa(bool, real_t) const;
};