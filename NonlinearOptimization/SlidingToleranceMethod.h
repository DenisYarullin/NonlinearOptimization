#pragma once

#include "NelderMeadMethod.h"

template <size_t N, typename InequalityBoundaryConditionsType>
class SlidingToleranceInterpolation
{
public:
	SlidingToleranceInterpolation(const PolyhedronVertex<N>& previousExternalVertex, const PolyhedronVertex<N>& currentInternalVertex,
		const std::vector<InequalityBoundaryConditionsType>& inequalities);

	void SetBorderVertices(const PolyhedronVertex<N>& previousExternalVertex, const PolyhedronVertex<N>& currentInternalVertex);
	PolyhedronVertex<N> CalculateInterpolationVertex();

private:
	PolyhedronVertex<N> previousExternalVertex_, currentInternalVertex_;
	std::vector<InequalityBoundaryConditionsType> inequalities_;
};


template <size_t N, typename InequalityBoundaryConditionsType>
SlidingToleranceInterpolation<N, InequalityBoundaryConditionsType>::SlidingToleranceInterpolation(const PolyhedronVertex<N>& previousExternalVertex, const PolyhedronVertex<N>& currentInternalVertex,
	const std::vector<InequalityBoundaryConditionsType>& inequalities) :
	previousExternalVertex_(previousExternalVertex), currentInternalVertex_(currentInternalVertex), inequalities_(inequalities)
{
}


template <size_t N, typename InequalityBoundaryConditionsType>
void SlidingToleranceInterpolation<N, InequalityBoundaryConditionsType>::SetBorderVertices(const PolyhedronVertex<N>& previousExternalVertex, const PolyhedronVertex<N>& currentInternalVertex)
{
	previousExternalVertex_ = previousExternalVertex;
	currentInternalVertex_ = currentInternalVertex;
}


template <size_t N, typename InequalityBoundaryConditionsType>
PolyhedronVertex<N> SlidingToleranceInterpolation<N, InequalityBoundaryConditionsType>::CalculateInterpolationVertex()
{
	double distance = CalculateDistanceBetweenVertices(previousExternalVertex_, currentInternalVertex_);
	PolyhedronVertex<N> unitDirectingVector = (1.0 / distance) * (previousExternalVertex_ - currentInternalVertex_);
	PolyhedronVertex<N> middleVertex = currentInternalVertex_ + 0.5 * distance * unitDirectingVector;

	std::vector<PolyhedronVertex<N> > vertices;
	vertices.push_back(currentInternalVertex_);
	vertices.push_back(middleVertex);
	vertices.push_back(previousExternalVertex_);

	std::vector<InequalityBoundaryConditionsType> disturbedInequalities;
	for (auto& inequality: inequalities_)
	{
		if (inequality(previousExternalVertex_) < 0.0)
			disturbedInequalities.push_back(inequality);
	}

	double functionValue = 0.0, squareSumInequalities = 0.0;
	std::vector<double> Z;
	for (auto& vertex: vertices)
	{
		for (auto& disturbedInequality: disturbedInequalities)
		{
			functionValue = disturbedInequality(vertex);
			squareSumInequalities += functionValue * functionValue;
		}

		Z.push_back(squareSumInequalities);
		squareSumInequalities = 0.0;
	}

	double alfa = Z[0] - 2.0 * Z[1] + Z[2];
	double betta = 3.0 * Z[0] - 4.0 * Z[1] + Z[2];

	return currentInternalVertex_ + (betta + sqrt(betta * betta - 8.0 * alfa * Z[0])) / (4.0 * alfa) * distance * unitDirectingVector;
}


template <size_t N, size_t M, typename Fn>
NelderMeadMethod<N, M, Fn> CreateNelderMeadMethodObject(PolyhedronVertex<N> initialVertex, double distanceBetweenVertices, Fn fn)
{
	NelderMeadMethod<N, M, Fn> object(initialVertex, distanceBetweenVertices, fn);
	return object;
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double CalculateA(NelderMeadMethod<N, M, ObjectiveFunctionType>& nelderMeadMethod)
{
	auto& centerOfGravity = nelderMeadMethod.FindCenterOfGravity();
	auto& simplex = nelderMeadMethod.GetSimplex();

	double tempValue = 0.0, sum = 0.0;
	for (size_t i = 0; i < simplex.SimplexSize(); ++i)
	{
		tempValue = nelderMeadMethod.GetValueObjectiveFunction(simplex[i]) - centerOfGravity.second;
		sum += (tempValue * tempValue);
	}
	return sqrt(sum) / (N + 1);
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
class SlidingToleranceMethod
{
public:
	SlidingToleranceMethod(double t, const PolyhedronVertex<N>& initialVertex, ObjectiveFunctionType func, 
		const std::vector<ObjectiveFunctionType>& equalities, const std::vector<ObjectiveFunctionType>& inequalities);
	SlidingToleranceMethod(const PolyhedronVertex<N>& lowBorderOfChange, const PolyhedronVertex<N>& highBorderOfChange, 
		const PolyhedronVertex<N>& initialVertex, ObjectiveFunctionType func, const std::vector<ObjectiveFunctionType>& equalities,
		const std::vector<ObjectiveFunctionType>& inequalities);

	double GetValueObjectiveFunction(const PolyhedronVertex<N>& vertex) const;
	double GetValueEqualityBoundaryCriterion(const PolyhedronVertex<N>& vertex, int index) const;
	double GetValueInequalityBoundaryCriterion(const PolyhedronVertex<N>& vertex, int index) const;

	double CalculateAverageDistanceFromSimplexVerticesToCenterOfGravity();

	double CalculateCriterionToleranceFirstStep() const;
	double CalculateCriterionToleranceKStep();
	double CalculateDistanceBetweenPolyhedronVertices(const PolyhedronVertex<N>& lowBorderOfChange, const PolyhedronVertex<N>& highBorderOfChange) const;
	double CalculateCriterionT(PolyhedronVertex<N> vertex);
	double CalculateCriterionTolerance();

	int FindTolerantVertex(const std::pair<PolyhedronVertex<N>, double>& initialVertex, std::pair<PolyhedronVertex<N>, double>& tolerantVertex, double t);
	int FindMinimum();
	
private:
	double t_, m_, r_;
	PolyhedronVertex<N> initialVertex_;
	ObjectiveFunctionType func_;
	std::vector<ObjectiveFunctionType> equalities_, inequalities_;
	NelderMeadMethod<N, M, ObjectiveFunctionType> mainNelderMead_;
	double criterionTolerancePrev_, criterionToleranceCur_;
	size_t stepK_;
};


template <size_t N, size_t M, typename ObjectiveFunctionType>
SlidingToleranceMethod<N, M, ObjectiveFunctionType>::SlidingToleranceMethod(double t, const PolyhedronVertex<N>& initialVertex, ObjectiveFunctionType func,
	const std::vector<ObjectiveFunctionType>& equalities, const std::vector<ObjectiveFunctionType>& inequalities) :
	t_(t), m_(M), r_(N - m_), initialVertex_(initialVertex), func_(func), equalities_(equalities), inequalities_(inequalities), 
	mainNelderMead_(initialVertex_, t_, func_), criterionTolerancePrev_(-1.0), criterionToleranceCur_(-1.0), stepK_(0)
{
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
SlidingToleranceMethod<N, M, ObjectiveFunctionType>::SlidingToleranceMethod(const PolyhedronVertex<N>& lowBorderOfChange, const PolyhedronVertex<N>& highBorderOfChange,
	const PolyhedronVertex<N>& initialVertex, ObjectiveFunctionType func, const std::vector<ObjectiveFunctionType>& equalities,
	const std::vector<ObjectiveFunctionType>& inequalities) :
	t_(CalculateDistanceBetweenPolyhedronVertices(lowBorderOfChange, highBorderOfChange)), m_(M), r_(N - m_), initialVertex_(initialVertex), func_(func), 
	equalities_(equalities), inequalities_(inequalities), mainNelderMead_(initialVertex_, t_, func_), criterionTolerancePrev_(-1.0), criterionToleranceCur_(-1.0), stepK_(0)
{
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::GetValueObjectiveFunction(const PolyhedronVertex<N>& vertex) const
{
	return func_(vertex);
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::GetValueEqualityBoundaryCriterion(const PolyhedronVertex<N>& vertex, int index) const
{
	return equalities_[index](vertex);
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::GetValueInequalityBoundaryCriterion(const PolyhedronVertex<N>& vertex, int index) const
{
	return inequalities_[index](vertex);
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::CalculateCriterionToleranceFirstStep() const
{
	return 2.0 * (m_ + 1.0) * t_;
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::CalculateDistanceBetweenPolyhedronVertices(const PolyhedronVertex<N>& lowBorderOfChange,
																									   const PolyhedronVertex<N>& highBorderOfChange) const
{
	auto& difference = highBorderOfChange - lowBorderOfChange;
	double sumElem = std::accumulate(difference.begin(), difference.end(), 0.0);
	double minElem = *std::min_element(difference.begin(), difference.end());
	return std::min((0.2 / N) * sumElem, minElem);
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::CalculateAverageDistanceFromSimplexVerticesToCenterOfGravity()
{
	auto& centerOfGravity = mainNelderMead_.FindCenterOfGravity().first;
	auto& simplex = mainNelderMead_.GetSimplex();

	PolyhedronVertex<N> tempVertex;
	double sum = 0.0;
	for (size_t i = 0; i < mainNelderMead_.SimplexSize(); ++i)
	{
		tempVertex = simplex[i] - centerOfGravity;
		for (size_t j = 0; j < N; ++j)
			sum += tempVertex[j] * tempVertex[j];
	}
	return sqrt(sum);
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::CalculateCriterionToleranceKStep()
{
	return (m_ + 1.0) / (r_ + 1.0) * CalculateAverageDistanceFromSimplexVerticesToCenterOfGravity();
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::CalculateCriterionT(PolyhedronVertex<N> vertex)
{
	double functionValue = 0.0, squareSumEqualities = 0.0;
	for (size_t i = 0; i < equalities_.size(); ++i)
	{
		functionValue = GetValueEqualityBoundaryCriterion(vertex, i);
		squareSumEqualities += functionValue * functionValue;
	}

	double heavysideOperator = 0.0, squareSumInequalitiesMultHeavyside = 0.0;
	for (size_t i = 0; i < inequalities_.size(); ++i)
	{
		functionValue = GetValueInequalityBoundaryCriterion(vertex, i);
		heavysideOperator = (functionValue >= 0.0 ? 0.0 : 1.0);
		squareSumInequalitiesMultHeavyside += heavysideOperator * functionValue * functionValue;
	}

	return sqrt(squareSumEqualities + squareSumInequalitiesMultHeavyside);
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
double SlidingToleranceMethod<N, M, ObjectiveFunctionType>::CalculateCriterionTolerance()
{
	if (stepK_ == 0)
		return CalculateCriterionToleranceFirstStep();
	return std::min(criterionTolerancePrev_, CalculateCriterionToleranceKStep());
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
int SlidingToleranceMethod<N, M, ObjectiveFunctionType>::FindTolerantVertex(const std::pair<PolyhedronVertex<N>, double>& initialVertex, 
																			std::pair<PolyhedronVertex<N>, double>& tolerantVertex, double t)
{
	double criterionT = -1.0;
	auto& nelderMeadMethod = CreateNelderMeadMethodObject<N, 0>(initialVertex.first, t, std::bind1st(std::mem_fun(&SlidingToleranceMethod::CalculateCriterionT), this));
	std::pair<PolyhedronVertex<N>, double> prevLowVertexT = initialVertex;

	while (true)
	{
		try
		{
			nelderMeadMethod.MakeOneStepMinimization();
		}
		catch (std::exception& e)
		{
			break;
		}

		std::pair<PolyhedronVertex<N>, double>& curLowVertexT = nelderMeadMethod.GetLowVertex();
		criterionT = curLowVertexT.second;

		if (criterionT == 0.0)
		{
			SlidingToleranceInterpolation<N, ObjectiveFunctionType> interpolation(prevLowVertexT.first, curLowVertexT.first, inequalities_);
			PolyhedronVertex<N>& interpolationVertex = interpolation.CalculateInterpolationVertex();

			while (CalculateCriterionT(interpolationVertex) > criterionToleranceCur_)
			{
				interpolation.SetBorderVertices(curLowVertexT.first, interpolationVertex);
				interpolationVertex = interpolation.CalculateInterpolationVertex();
				curLowVertexT.first = interpolationVertex;
			}
			tolerantVertex = std::make_pair(interpolationVertex, GetValueObjectiveFunction(interpolationVertex));
			break;
		}
		else if (criterionT <= criterionToleranceCur_)
		{
			tolerantVertex = std::make_pair(curLowVertexT.first, GetValueObjectiveFunction(curLowVertexT.first));
			break;
		}
		else
		{
			if (CalculateA(nelderMeadMethod) > 1e-7)
				continue;
			else
				return -1;
		}
		prevLowVertexT = curLowVertexT;
	}
	return 0;
}


template <size_t N, size_t M, typename ObjectiveFunctionType>
int SlidingToleranceMethod<N, M, ObjectiveFunctionType>::FindMinimum()
{
	double criterionT = -1.0;
	while ((criterionToleranceCur_ = CalculateCriterionTolerance()) > std::numeric_limits<double>::epsilon())
	{
		double t = 0.05 * criterionToleranceCur_;

		auto& lowVertex =  mainNelderMead_.GetLowVertex();
		auto& highVertex = mainNelderMead_.GetHighVertex();
		criterionT = CalculateCriterionT(lowVertex.first);

		if (criterionT > criterionToleranceCur_)
		{
			std::pair<PolyhedronVertex<N>, double> newSimplexVertex;
			if (FindTolerantVertex(lowVertex, newSimplexVertex, t) == 0)
				mainNelderMead_.ChangeSimplexVertex(lowVertex, newSimplexVertex);
			else
				return -1;
		}

		auto& centerOfGravity = mainNelderMead_.FindCenterOfGravity();

		if (criterionToleranceCur_ < std::numeric_limits<double>::epsilon())
			return 0;

		auto& reflectionVertex = mainNelderMead_.Reflection();
		criterionT = CalculateCriterionT(reflectionVertex.first);

		if (criterionT > criterionToleranceCur_)
		{
			std::pair<PolyhedronVertex<N>, double> newReflectionVertex;
			if (FindTolerantVertex(reflectionVertex, newReflectionVertex, t) == 0)
				mainNelderMead_.ChangeReflectionVertex(newReflectionVertex);
			else
				return -1;
		}

		if (reflectionVertex.second < lowVertex.second)
		{
			auto& expansionVertex = mainNelderMead_.Expansion();
			criterionT = CalculateCriterionT(expansionVertex.first);

			std::pair<PolyhedronVertex<N>, double> newExpansionVertex;
			if (criterionT > criterionToleranceCur_)
			{
				if (FindTolerantVertex(expansionVertex, newExpansionVertex, t) == 0)
					mainNelderMead_.ChangeExpansionVertex(newExpansionVertex);
				else
					return -1;
			}

			expansionVertex.second < mainNelderMead_.GetLowVertex().second ? mainNelderMead_.ChangeHighVertex(expansionVertex) :
																			 mainNelderMead_.ChangeHighVertex(reflectionVertex);
		}
		else
		{
			if (mainNelderMead_.CompareAllVerticesWithReflectionVertexExceptHighVertex())
			{
				if (reflectionVertex.second <= highVertex.second)
					mainNelderMead_.ChangeHighVertex(reflectionVertex);

				auto& contractionVertex = mainNelderMead_.Contraction();
				criterionT = CalculateCriterionT(contractionVertex.first);

				if (criterionT > criterionToleranceCur_)
				{
					std::pair<PolyhedronVertex<N>, double> newContractionVertex;
					if (FindTolerantVertex(contractionVertex, newContractionVertex, t) == 0)
						mainNelderMead_.ChangeContractionVertex(newContractionVertex);
					else
						return -1;
				}
				
				contractionVertex.second > highVertex.second ? mainNelderMead_.Reduction() : mainNelderMead_.ChangeHighVertex(contractionVertex);
			}
			else
			{
				mainNelderMead_.ChangeHighVertex(reflectionVertex);
			}
		}

		++stepK_;
		lowVertex = mainNelderMead_.GetLowVertex();
		highVertex = mainNelderMead_.GetHighVertex();
		criterionTolerancePrev_ = criterionToleranceCur_;
	}
	return 0;
}

