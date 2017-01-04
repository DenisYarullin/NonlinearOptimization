#pragma once

#include "Vertex.h"

template <size_t N, size_t M>
class Simplex
{
public:
	Simplex(const PolyhedronVertex<N>& initialVertex, double distanceBetweenVertices); 

	PolyhedronVertex<N> operator[] (int index) const;
	PolyhedronVertex<N>& operator[] (int index);

	void BuildInitialPolyhedron();
	void ChangeVertex(const PolyhedronVertex<N>& destinationVertex, const PolyhedronVertex<N>& sourceVertex);
	size_t SimplexSize() const;

private:
	static const size_t Q = ((N - M + 1) < 3) ? 3 : N - M + 1;
	PolyhedronVertex<N> initialVertex_;
	double distanceBetweenVertices_;
	array<PolyhedronVertex<N>, Q> simplex_;
};


template <size_t N, size_t M>
Simplex<N, M>::Simplex(const PolyhedronVertex<N>& initialVertex, double distanceBetweenVertices) :
	initialVertex_(initialVertex), distanceBetweenVertices_(distanceBetweenVertices) {}


template <size_t N, size_t M>
PolyhedronVertex<N> Simplex<N, M>::operator[] (int index) const
{
	return simplex_[index];
}


template <size_t N, size_t M>
PolyhedronVertex<N>& Simplex<N, M>::operator[] (int index)
{
	return simplex_[index];
}


template <size_t N, size_t M>
void Simplex<N, M>::BuildInitialPolyhedron()
{
	simplex_.fill(initialVertex_);
	
	double d1 = (distanceBetweenVertices_ / (N * sqrt(2.0))) * (sqrt(N + 1.0) + N - 1.0);
	double d2 = (distanceBetweenVertices_ / (N * sqrt(2.0))) * (sqrt(N + 1.0) - 1.0);

	int indexChangeValue = 0;
	for (auto i = simplex_.begin(); i != simplex_.end(); ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			if (j != indexChangeValue)
				(*i)[j] += d2;
			else
				(*i)[j] += d1;
		}
		++indexChangeValue;
	}
}


template <size_t N, size_t M>
void Simplex<N, M>::ChangeVertex(const PolyhedronVertex<N>& destinationVertex, const PolyhedronVertex<N>& sourceVertex)
{
	std::replace(simplex_.begin(), simplex_.end(), destinationVertex, sourceVertex);
}


template <size_t N, size_t M>
size_t Simplex<N, M>::SimplexSize() const
{
	return Q;
}


template <size_t N, size_t M>
const size_t Simplex<N, M>::Q;