//
// Created by alex on 04.12.19.
//

#include <cmath>
#include <iostream>
#include <utility>
#include "LAE_Solver.h"

std::vector<double> LAE_Solver::solve(std::vector<std::vector<double>> a, std::vector<double> b)
{
	A = std::move(a);
	B = b;
	X = std::vector<double>(b.size());
	dimension = static_cast<unsigned int>(b.size());
	straight_gauss();
	reverse_gauss();
	return X;
}

void LAE_Solver::straight_gauss()
{
	for(int i = 0; i < dimension; ++i)
	{
		int index = find_max_in_A_from(i);
		swap_lines(i, index);
		for(int j = i + 1; j < dimension; ++j)
		{
			double m = -A[j][i] / A[i][i];
			for(int k = 0; k < dimension; ++k) // sum lines
			{
				A[j][k] += m * A[i][k];
			}
			B[j] += m * B[i];
		}
	}
}

void LAE_Solver::reverse_gauss()
{
	for(int i = dimension - 1; i >= 0; --i)
	{
		double temp = B[i];
		for(int j = i + 1; j < dimension; ++j)
		{
			temp -= X[j] * A[i][j];
		}
		X[i] = temp / A[i][i];
	}
}

int LAE_Solver::find_max_in_A_from(int i)
{
	double max = fabs(A[i][i]);
	int index = i;
	for(int j = 0; j < dimension; ++j)
	{
		double num = fabs(A[j][i]);
		if(num > max)
		{
			max = num;
			index = j;
		}
	}
	return index;
}

void LAE_Solver::swap_lines(int l, int r)
{
	std::swap(A[l], A[r]);
	std::swap(B[l], B[r]);
}

