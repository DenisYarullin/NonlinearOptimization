#include <iostream>
#include "Vertex.h"
#include "NelderMeadMethod.h"
#include "SlidingToleranceMethod.h"


template <size_t N>
struct TestObjectiveFunctionType
{
	typedef double(*fptr)(const PolyhedronVertex<N>& vertex);
};


template <size_t N>
double ObjectiveFunction1(const PolyhedronVertex<N>& vertex)
{
	return 4.0 * (vertex[0] - 5.0) * (vertex[0] - 5.0) + (vertex[1] - 6.0) * (vertex[1] - 6.0);
	//return (1.0 - vertex[0]) * (1.0 - vertex[0]) + 100 * (vertex[1] - vertex[0] * vertex[0]) * (vertex[1] - vertex[0] * vertex[0]); // Функция Розенброка
}


template <size_t N>
double ObjectiveFunction2(const PolyhedronVertex<N>& vertex)
{
	//return (1.0 - vertex[0]) * (1.0 - vertex[0]) + 100 * (vertex[1] - vertex[0] * vertex[0]) * (vertex[1] - vertex[0] * vertex[0]); // Функция Розенброка
}


int main()
{
	PolyhedronVertex<2> initialVertex1;
	initialVertex1[0] = 1.0;
	initialVertex1[1] = 1.0;

	double distanceBetweenVertices1 = 1.0;

	NelderMeadMethod<2, 0, TestObjectiveFunctionType<2>::fptr> nelderMead1(initialVertex1, distanceBetweenVertices1, ObjectiveFunction1);
	nelderMead1.FindMinimum();

	//------------------------------------------------------

	PolyhedronVertex<2> initialVertex2;
	initialVertex2[0] = 5.0;
	initialVertex2[1] = 6.0;

	double distanceBetweenVertices2 = 1.0;

	NelderMeadMethod<2, 0, TestObjectiveFunctionType<2>::fptr> nelderMead2(initialVertex2, distanceBetweenVertices2, ObjectiveFunction2);
	nelderMead2.FindMinimum();

	//----------------------------------------------------

	return 0;
}