//
// Created by alex on 05.12.19.
//

#include <iostream>
#include "LAE_Solver.h"

void LAE_test()
{
	LAE_Solver solver;
	auto result = solver.solve(
		{
			{0, 1, 1, 0},
			{1, 0, 1, 0},
			{0, 1, 0, 1},
			{1, 1, 0, 0},
		}, {2, 2, 2, 2});
	bool passed = true;
	std::vector<double> expected({1, 1, 1, 1});
	for(int i = 0; i < expected.size(); ++i)
	{
		if(result[i] != expected[i])
		{
			passed = false;
			break;
		}
	}
	if(passed)
		std::cout << "Success!" << std::endl;
	else
	{
		std::cout << "Fail!" << std::endl;
		std::cout << "Expected:" << std::endl;
		for(auto& it : expected)
		{
			std::cout << it << ' ';
		}
		std::cout << std::endl;
		std::cout << "Got:" << std::endl;
		for(auto& it : result)
		{
			std::cout << it << ' ';
		}
		std::cout << std::endl;
	}
}

int main()
{
	LAE_test();
}