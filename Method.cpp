//
// Created by alex on 04.12.19.
//

#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "Method.h"
#include "Consts.h"


std::vector<std::vector<double>> make_triangle(std::vector<std::vector<double>> ansambl_K_matrix, std::vector<double> & value_vector ) {
	int size = ansambl_K_matrix.size();
	
	std::vector<std::vector<double>> matrix_result = ansambl_K_matrix;
	
	double value = 0.0;
	for (int i = 0; i < size; ++i) {
		double value_max = 0.0;
		int i_max = i;
		for (int j = i + 1; j < size; ++j) {
			value = matrix_result[j][i];
			value_max = matrix_result[i_max][i];
			if (fabs(value) > fabs(value_max)) {
				i_max = j;
			}
		}
		
		std::swap(value_vector[i], value_vector[i_max]);
		std::swap(matrix_result[i], matrix_result[i_max]);
		for (int j = i + 1; j < size; ++j) {
			double q = 0.0;
			value = matrix_result[j][i];
			value_max = matrix_result[i][i];
			if (value_max) {
				q = -1 * value / value_max;
			}
			value_vector[j] += q * value_vector[i];
			for (int k = i; k < size; ++k) {
				value = matrix_result[j][k];
				value_max = matrix_result[i][k];
				value += q * value_max;
				matrix_result[j][k] = value;
			}
		}
	}
	return matrix_result;
}

std::vector<double> solve_slau(std::vector<std::vector<double>> matrix_not, std::vector<double> value_vector) {
	std::vector<std::vector<double>> matrix = make_triangle(matrix_not, value_vector);
	int size = value_vector.size();
	std::vector<double> vector_result(size, 0);
	
	vector_result[size - 1] = (value_vector[size - 1] / matrix[size -1][size - 1]);
	
	for (int i = size - 2; i != -1; i--) {
		double x = value_vector[i];
		for (int j = i + 1; j != size; j++) {
			x -= vector_result[j] * matrix[i][j];
		}
		x /= matrix[i][i];
		vector_result[i] = x;
	}
	
	return vector_result;
}


Method::Method(double epsilon)
	:
	eps(epsilon)
{}

double Et(double time)
{
	const double d = 10 * sin(time * (2 * M_PI) / 0.0001);
	return d;
}

void Method::solve(double final_time, double _dt)
{
	dimension = 17;
	basis = std::vector<double>(dimension);
	prev_basis = std::vector<double>(dimension);
	delta = std::vector<double>(dimension);
	dt = _dt;
	std::fstream f14;
	f14.open("out14.txt", std::fstream::out);
	std::fstream f12;
	f12.open("out12.txt", std::fstream::out);
	for(current_time = 0; current_time < final_time; current_time += dt)
	{
		std::cout << "======================================" << std::endl;
		std::cout << "Time: " << current_time << std::endl;
		std::cout << "Basis: " << std::endl;
		print_vector(basis);
		do
		{
			find_delta();
//			std::cout << "\t\tDelta: " << std::endl;
//			print_vector(delta);
			for(int i = 0; i < dimension; ++i)
			{
				basis[i] += delta[i];
			}
//			std::cout << "\t\tSum: " << std::endl;
//			print_vector(basis);
			std::cout << "Max: " << find_max_delta() << " Time: " << current_time << std::endl;
		} while(find_max_delta() > eps);
		f14 << current_time << "\t" << basis[14] << std::endl;
		f12 << current_time << "\t" << basis[13] << std::endl;
		prev_basis = basis;
	}
	f14.close();
	f12.close();
}

void Method::find_delta()
{
	const double R1 = Consts::R1;
	const double R2 = Consts::R2;
	const double L = Consts::L;
	const double C1 = Consts::C1;
	const double E = Consts::E;
	const double rb = Consts::rb;
	const double ry = Consts::ry;
	const double Cb = Consts::Cb;
	const double It = Consts::It;
	const double mfit = Consts::mfit;
	const double k = 1;
	std::vector<double> I = {
		basis[0] - 1/dt * (basis[10] - prev_basis[10]),
		basis[1] - 1/dt * (basis[11] - prev_basis[11]),
		basis[2] - 1/dt * (basis[12] - prev_basis[12]),
		basis[3] - 1/dt * (basis[13] - prev_basis[13]),
		basis[4] - 1/dt * (basis[14] - prev_basis[14]),
		basis[5] - (prev_basis[5] + basis[10] * dt),
		basis[6] - (prev_basis[6] + basis[11] * dt),
		basis[7] - (prev_basis[7] + basis[12] * dt),
		basis[8] - (prev_basis[8] + basis[13] * dt),
		basis[9] - (prev_basis[9] + basis[14] * dt),
		basis[15] - basis[16] - 1/R1 * (basis[12] - basis[10]), // 10
		basis[16] - 1/L * (basis[7] - basis[6]), // 11
		1/R1 * (basis[12] - basis[10]) + 1/L * (basis[7] - basis[6]) - 1/rb * (basis[13] - basis[12]), // 12
		1/rb * (basis[13] - basis[12]) + It*(exp((basis[13] - basis[14])/mfit)-1) + 1/ry * (basis[13] - basis[14]) + Cb * (basis[3] - basis[4]), // 13
//		C1 * basis[4] + 1/R2 * basis[14] - It*(exp((basis[13] - basis[14])/mfit)-1) - 1/ry * (basis[13] - basis[14]) - Cb * (basis[3] - basis[4]), // 14
		C1 * basis[4] + k*basis[14]*basis[14] - It*(exp((basis[13] - basis[14])/mfit)-1) - 1/ry * (basis[13] - basis[14]) - Cb * (basis[3] - basis[4]), // 14
		basis[10] - E,
		basis[11] - basis[10] - Et(current_time)
	};
	for(auto& it : I)
	{
		it = -it;
	}
	double dIdt = It/mfit * exp((basis[13] - basis[14])/mfit);
	const double d4d4 = 1/rb + 1/ry + dIdt;
	const double d4d5 = 	 - 1/ry - dIdt;
	const double d5d4 = 	 - 1/ry - dIdt;
	const double d5d5 = 2*k*basis[14] + 1/ry + dIdt;
	std::vector<std::vector<double>> Y = {
		{1, 0, 0, 0, 0,			0, 0, 0, 0, 0,			-1/ dt, 0, 0, 0, 0,				0, 0}, // d/dt
		{0, 1, 0, 0, 0,			0, 0, 0, 0, 0,			0, -1/ dt, 0, 0, 0,				0, 0}, // d/dt
		{0, 0, 1, 0, 0,			0, 0, 0, 0, 0,			0, 0, -1/ dt, 0, 0,				0, 0}, // d/dt
		{0, 0, 0, 1, 0,			0, 0, 0, 0, 0,			0, 0, 0, -1/ dt, 0,				0, 0}, // d/dt
		{0, 0, 0, 0, 1,			0, 0, 0, 0, 0,			0, 0, 0, 0, -1/ dt,				0, 0}, // d/dt

		{0, 0, 0, 0, 0,			1, 0, 0, 0, 0,			-dt, 0, 0, 0, 0,				0, 0}, // integral
		{0, 0, 0, 0, 0,			0, 1, 0, 0, 0,			0, -dt, 0, 0, 0,				0, 0}, // integral
		{0, 0, 0, 0, 0,			0, 0, 1, 0, 0,			0, 0, -dt, 0, 0,				0, 0}, // integral
		{0, 0, 0, 0, 0,			0, 0, 0, 1, 0,			0, 0, 0, -dt, 0,				0, 0}, // integral
		{0, 0, 0, 0, 0,			0, 0, 0, 0, 1,			0, 0, 0, 0, -dt,				0, 0}, // integral

		{0, 0, 0, 0, 0,			0, 0, 0, 0, 0,          1/R1, 0, -1/R1, 0, 0,			1, -1}, // potential
		{0, 0, 0, 0, 0,			0, 1/L, -1/L, 0, 0,  	0, 0, 0, 0, 0,					0, 1}, // potential
		{0, 0, 0, 0, 0,			0, -1/L, 1/L, 0, 0,  	-1/R1, 0, 1/R1+1/rb, -1/rb, 0,	0, 0}, // potential
		{0, 0, 0, Cb, -Cb,		0, 0, 0, 0, 0,          0, 0, -1/rb, d4d4, d4d5,		0, 0}, // potential
		{0, 0, 0, -Cb, Cb+C1,	0, 0, 0, 0, 0,          0, 0, 0, d5d4, d5d5,			0, 0}, // potential

		{0, 0, 0, 0, 0,			0, 0, 0, 0, 0,          1, 0, 0, 0, 0,					0, 0},
		{0, 0, 0, 0, 0,			0, 0, 0, 0, 0,          -1, 1, 0, 0, 0,					0, 0}
	};
//	print_matrix(Y);
//	print_vector(I);
//	delta = solver.solve(Y, I);
	delta = solve_slau(Y, I);
}

void Method::print_vector(std::vector<double> vector)
{
	for(auto& it : vector)
	{
		std::cout << it << std::endl;
	}
}

void Method::print_matrix(std::vector<std::vector<double>> matrix)
{
	for(auto& row : matrix)
	{
		for(auto& it : row)
		{
			std::cout << std::setw(6) << it << " ";
		}
		std::cout << std::endl;
	}
}

double Method::find_max_delta()
{
	double max = fabs(delta[0]);
	for(auto& it : delta)
	{
		double num = fabs(it);
		if(num > max)
		{
			max = num;
		}
	}
	return max;
}
