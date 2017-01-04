#include <iostream>
#include "NelderMeadMethod.h"
#include "SlidingToleranceMethod.h"

template <size_t N>
struct TestObjectiveFunctionType
{
	typedef double(*fptr)(const PolyhedronVertex<N>& vertex);
};

// test1
double ObjectiveFunction1(const PolyhedronVertex<2>& vertex)
{
	return 4.0 * (vertex[0] - 5.0) * (vertex[0] - 5.0) + (vertex[1] - 6.0) * (vertex[1] - 6.0);
}

// test2 - Функция Розенброка
double ObjectiveFunction2(const PolyhedronVertex<2>& vertex)
{
	return (1.0 - vertex[0]) * (1.0 - vertex[0]) + 100 * (vertex[1] - vertex[0] * vertex[0]) * (vertex[1] - vertex[0] * vertex[0]); 
}

// test3
double ObjectiveFunction3(const PolyhedronVertex<2>& vertex)
{
	return 4.0 * vertex[0] - vertex[1] * vertex[1] - 12.0;
}

double BoundaryCriterionOfEquality3(const PolyhedronVertex<2>& vertex)
{
	return 25.0 - vertex[0] * vertex[0] - vertex[1] * vertex[1];
}

double BoundaryCriterionOfInequality31(const PolyhedronVertex<2>& vertex)
{
	return 10.0 * vertex[0] - vertex[0] * vertex[0] + 10.0 * vertex[1] - vertex[1] * vertex[1] - 34.0;
}

double BoundaryCriterionOfInequality32(const PolyhedronVertex<2>& vertex)
{
	return vertex[0];
}

double BoundaryCriterionOfInequality33(const PolyhedronVertex<2>& vertex)
{
	return vertex[1];
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

	PolyhedronVertex<2> initialVertex3;
	initialVertex3[0] = 1.0;
	initialVertex3[1] = 1.0;

	std::vector<TestObjectiveFunctionType<2>::fptr> equalities;
	equalities.push_back(BoundaryCriterionOfEquality3);

	std::vector<TestObjectiveFunctionType<2>::fptr> inequalities;
	inequalities.push_back(BoundaryCriterionOfInequality31);
	inequalities.push_back(BoundaryCriterionOfInequality32);
	inequalities.push_back(BoundaryCriterionOfInequality33);

	SlidingToleranceMethod<2, 1, TestObjectiveFunctionType<2>::fptr> slidingToleranceMethod(0.3, initialVertex3, ObjectiveFunction3, equalities, inequalities);
	slidingToleranceMethod.FindMinimum();

	return 0;
}