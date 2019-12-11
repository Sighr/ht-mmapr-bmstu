//
// Created by alex on 04.12.19.
//

#ifndef MMAPR_SECOND_TRY_LAE_SOLVER_H
#define MMAPR_SECOND_TRY_LAE_SOLVER_H


#include <vector>

class LAE_Solver
{
public:
	std::vector<double> solve(std::vector<std::vector<double>> a, std::vector<double> b);
private:
	void straight_gauss();
	void reverse_gauss();
	int find_max_in_A_from(int i);
	void swap_lines(int l, int r);
	std::vector<std::vector<double>> A;
	std::vector<double> B;
	std::vector<double> X;
	unsigned int dimension;
};


#endif //MMAPR_SECOND_TRY_LAE_SOLVER_H
