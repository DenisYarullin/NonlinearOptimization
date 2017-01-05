#pragma once

#include "Simplex.h"
#include <vector>
#include <exception>

#define EPSILON 1e-12

template <size_t N>
struct Less
{
	bool operator() (const std::pair<PolyhedronVertex<N>, double>& lhs, const std::pair<PolyhedronVertex<N>, double>& rhs) const
	{
		return lhs.second < rhs.second;
	}
};


template <size_t N, size_t M, typename Fn>
class NelderMeadMethod
{
public:
	NelderMeadMethod(const PolyhedronVertex<N>& initialVetex, double distanceBetweenVertices, Fn func, double alfa = 1.0, double betta = 0.5, double gamma = 2.0);

	size_t SimplexSize() const;
	double GetValueObjectiveFunction(const PolyhedronVertex<N>& vertex);

	std::pair<PolyhedronVertex<N>, double>& GetLowVertex();
	std::pair<PolyhedronVertex<N>, double>& GetHighVertex();
	std::pair<PolyhedronVertex<N>, double>& FindCenterOfGravity();
	Simplex<N, M> GetSimplex();

	void ChangeHighVertex(const std::pair<PolyhedronVertex<N>, double>& vertex);
	void ChangeSimplexVertex(const std::pair<PolyhedronVertex<N>, double>& destinationVertex, const std::pair<PolyhedronVertex<N>, double>& sourceVertex);
	void ChangeReflectionVertex(const std::pair<PolyhedronVertex<N>, double>& vertex);
	void ChangeExpansionVertex(const std::pair<PolyhedronVertex<N>, double>& vertex);
	void ChangeContractionVertex(const std::pair<PolyhedronVertex<N>, double>& vertex);

	std::pair<PolyhedronVertex<N>, double>& Reflection();
	std::pair<PolyhedronVertex<N>, double>& Expansion();
	std::pair<PolyhedronVertex<N>, double>& Contraction();
	void Reduction();

	bool CompareAllVerticesWithReflectionVertexExceptHighVertex();
	void MakeOneStepMinimization();
	void ExecuteMinimizationFunctions();
	void FindMinimum();

private:
	bool IsAlgorithmCoverged();

private:
	Simplex<N, M> simplex_;
	Fn func_;
	size_t numberSimplexVertices_;
	double alfa_, betta_, gamma_;
	std::vector<std::pair<PolyhedronVertex<N>, double> > verticesAndValuesObjectiveFunction_;
	std::pair<PolyhedronVertex<N>, double> lowVertex_, highVertex_, centerOfGravity_;
	std::pair<PolyhedronVertex<N>, double> reflectionVertex_, expansionVertex_, contractionVertex_;
};


template <size_t N, size_t M, typename Fn>
NelderMeadMethod<N, M, Fn>::NelderMeadMethod(const PolyhedronVertex<N>& initialVetex, double distanceBetweenVertices, Fn func, double alfa, double betta, double gamma) :
	simplex_(initialVetex, distanceBetweenVertices), func_(func), numberSimplexVertices_(simplex_.SimplexSize()), alfa_(alfa), betta_(betta), gamma_(gamma)
{
	simplex_.BuildInitialPolyhedron();
	verticesAndValuesObjectiveFunction_.resize(numberSimplexVertices_);

	for (size_t i = 0; i < numberSimplexVertices_; ++i)
		verticesAndValuesObjectiveFunction_[i] = std::make_pair(simplex_[i], GetValueObjectiveFunction(simplex_[i]));
}


template <size_t N, size_t M, typename Fn>
size_t NelderMeadMethod<N, M, Fn>::SimplexSize() const
{
	return numberSimplexVertices_;
}


template <size_t N, size_t M, typename Fn>
double NelderMeadMethod<N, M, Fn>::GetValueObjectiveFunction(const PolyhedronVertex<N>& vertex)
{
	return func_(vertex);
}


template <size_t N, size_t M, typename Fn>
std::pair<PolyhedronVertex<N>, double>& NelderMeadMethod<N, M, Fn>::GetLowVertex()
{
	lowVertex_ = *std::min_element(verticesAndValuesObjectiveFunction_.begin(), verticesAndValuesObjectiveFunction_.end(), Less<N>());
	return lowVertex_;
}


template <size_t N, size_t M, typename Fn>
std::pair<PolyhedronVertex<N>, double>& NelderMeadMethod<N, M, Fn>::GetHighVertex()
{
	highVertex_ = *std::max_element(verticesAndValuesObjectiveFunction_.begin(), verticesAndValuesObjectiveFunction_.end(), Less<N>());
	return highVertex_;
}


template <size_t N, size_t M, typename Fn>
std::pair<PolyhedronVertex<N>, double>& NelderMeadMethod<N, M, Fn>::FindCenterOfGravity()
{
	centerOfGravity_.first = PolyhedronVertex<N>();

	for (size_t i = 0; i < numberSimplexVertices_; ++i)
		centerOfGravity_.first += simplex_[i];

	centerOfGravity_.first -= highVertex_.first;
	centerOfGravity_.first = (1.0 / N) * centerOfGravity_.first ;
	centerOfGravity_.second = GetValueObjectiveFunction(centerOfGravity_.first);

	return centerOfGravity_;
}


template <size_t N, size_t M, typename Fn>
Simplex<N, M> NelderMeadMethod<N, M, Fn>::GetSimplex()
{
	return simplex_;
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::ChangeHighVertex(const std::pair<PolyhedronVertex<N>, double>& vertex)
{
	simplex_.ChangeVertex(highVertex_.first, vertex.first);
	std::replace(verticesAndValuesObjectiveFunction_.begin(), verticesAndValuesObjectiveFunction_.end(), highVertex_, vertex);
	highVertex_ = vertex;
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::ChangeSimplexVertex(const std::pair<PolyhedronVertex<N>, double>& destinationVertex, const std::pair<PolyhedronVertex<N>, double>& sourceVertex)
{
	simplex_.ChangeVertex(destinationVertex.first, sourceVertex.first);
	std::replace(verticesAndValuesObjectiveFunction_.begin(), verticesAndValuesObjectiveFunction_.end(), destinationVertex, sourceVertex);
	GetLowVertex();
	GetHighVertex();
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::ChangeReflectionVertex(const std::pair<PolyhedronVertex<N>, double>& vertex)
{
	reflectionVertex_ = vertex;
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::ChangeExpansionVertex(const std::pair<PolyhedronVertex<N>, double>& vertex)
{
	expansionVertex_ = vertex;
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::ChangeContractionVertex(const std::pair<PolyhedronVertex<N>, double>& vertex)
{
	contractionVertex_ = vertex;
}


template <size_t N, size_t M, typename Fn>
std::pair<PolyhedronVertex<N>, double>& NelderMeadMethod<N, M, Fn>::Reflection()
{
	reflectionVertex_.first = centerOfGravity_.first + alfa_ * (centerOfGravity_.first - highVertex_.first);
	reflectionVertex_.second = GetValueObjectiveFunction(reflectionVertex_.first);
	return reflectionVertex_;
}


template <size_t N, size_t M, typename Fn>
std::pair<PolyhedronVertex<N>, double>& NelderMeadMethod<N, M, Fn>::Expansion()
{
	expansionVertex_.first = centerOfGravity_.first + gamma_ * (reflectionVertex_.first - centerOfGravity_.first);
	expansionVertex_.second = GetValueObjectiveFunction(expansionVertex_.first);
	return expansionVertex_;
}


template <size_t N, size_t M, typename Fn>
std::pair<PolyhedronVertex<N>, double>& NelderMeadMethod<N, M, Fn>::Contraction()
{
	contractionVertex_.first = centerOfGravity_.first + betta_ * (highVertex_.first - centerOfGravity_.first);
	contractionVertex_.second = GetValueObjectiveFunction(contractionVertex_.first);
	return contractionVertex_;
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::Reduction()
{
	for (size_t i = 0; i < numberSimplexVertices_; ++i)
	{
		simplex_[i] = lowVertex_.first + 0.5 * (simplex_[i] - lowVertex_.first);
		verticesAndValuesObjectiveFunction_[i] = std::make_pair(simplex_[i], GetValueObjectiveFunction(simplex_[i]));
	}
}


template <size_t N, size_t M, typename Fn>
bool NelderMeadMethod<N, M, Fn>::CompareAllVerticesWithReflectionVertexExceptHighVertex()
{
	for (size_t i = 0; i < numberSimplexVertices_; ++i)
	{
		if (verticesAndValuesObjectiveFunction_[i].first != highVertex_.first &&
			verticesAndValuesObjectiveFunction_[i].second >= reflectionVertex_.second)
			return false;
	}
	return true;
}


template <size_t N, size_t M, typename Fn>
bool NelderMeadMethod<N, M, Fn>::IsAlgorithmCoverged()
{
	auto& min_max_values = std::minmax_element(verticesAndValuesObjectiveFunction_.begin(), verticesAndValuesObjectiveFunction_.end(), Less<N>());
	lowVertex_ = *min_max_values.first;
	highVertex_ = *min_max_values.second;

	FindCenterOfGravity();

	double valueConditionExpression = 0.0;
	for (size_t i = 0; i < numberSimplexVertices_; ++i)
		valueConditionExpression += pow(verticesAndValuesObjectiveFunction_[i].second - centerOfGravity_.second, 2.0);

	valueConditionExpression = sqrt((1.0 / numberSimplexVertices_) * valueConditionExpression);

	if (valueConditionExpression < EPSILON)
		return true;
	return false;
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::MakeOneStepMinimization()
{
	if (!IsAlgorithmCoverged())
		ExecuteMinimizationFunctions();
	else
		throw std::exception("Nelder Mead method is coverged!");	
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::ExecuteMinimizationFunctions()
{
	Reflection();

	if (reflectionVertex_.second < lowVertex_.second)
	{
		Expansion();
		expansionVertex_.second < lowVertex_.second ? ChangeHighVertex(expansionVertex_) : ChangeHighVertex(reflectionVertex_);
		return;
	}

	if (CompareAllVerticesWithReflectionVertexExceptHighVertex())
	{
		if (reflectionVertex_.second <= highVertex_.second)
			ChangeHighVertex(reflectionVertex_);
		Contraction();
	}
	else
	{
		ChangeHighVertex(reflectionVertex_);
		return;
	}

	contractionVertex_.second > highVertex_.second ? Reduction() : ChangeHighVertex(contractionVertex_);
}


template <size_t N, size_t M, typename Fn>
void NelderMeadMethod<N, M, Fn>::FindMinimum()
{
	while (!IsAlgorithmCoverged())
	{
		ExecuteMinimizationFunctions();
	}
}