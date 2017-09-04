#pragma once
#include "DataType.h"
#include "Zone.h"
#include "Reader.h"
#include "FRbasis.h"

class Post
{
public:
	// Constructor / p.m. Reader(object)
	Post(std::shared_ptr<Reader>);

	// Destructor
	~Post();

public:
	// Functions
	// Export solution file / p.m. Zone(object)
	void solution(std::shared_ptr<Zone>) const;

	// Export solution file / p.m. file name, Zone(object)
	void solution(const std::string&, std::shared_ptr<Zone>) const;

	// Export DG solution file / p.m. file name, Zone(object)
	void FRsolution(const std::string&, std::shared_ptr<Zone>) const;
	
	// Export error log / p.m. Zone(object)
	void error(std::shared_ptr<Zone>) const;

protected:
	// Variables
	std::shared_ptr<Reader> _reader;

protected:
	// Functions
	// Write vector to file / p.m. file name, variable name1, variable name2, variable1, variable2 
	void write(const std::string&, const std::string&, const std::string&, std::vector<real_t>, std::vector<real_t>) const;
};