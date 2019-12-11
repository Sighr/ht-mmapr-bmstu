#include <iostream>
#include "Method.h"

int main()
{
	double epsilon = 1e-1;
	Method method(epsilon);
	method.solve(0.001, 1e-9);
}