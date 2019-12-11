//
// Created by alex on 04.12.19.
//

#ifndef MMAPR_SECOND_TRY_METHOD_H
#define MMAPR_SECOND_TRY_METHOD_H

#include "LAE_Solver.h"

class Method
{
public:
	explicit Method(double epsilon);
	void solve(double final_time, double _dt);
private:
	void print_vector(std::vector<double> vector);
	void print_matrix(std::vector<std::vector<double>> matrix);
	double find_max_delta();
	const double eps;
	LAE_Solver solver;
	double dt;
	double current_time;
	std::vector<double> basis;
	std::vector<double> prev_basis;
	std::vector<double> delta;
	unsigned int dimension;
	void find_delta();
};


#endif //MMAPR_SECOND_TRY_METHOD_H
