#pragma once

#include <array>
#include <algorithm>
#include <functional>
#include <numeric>

using std::array;

template <size_t N>
class PolyhedronVertex
{
public:
	PolyhedronVertex();
	~PolyhedronVertex() = default;

	double operator[] (int index) const;
	double& operator[] (int index);

	PolyhedronVertex& operator-= (const PolyhedronVertex& other);
	PolyhedronVertex& operator+= (const PolyhedronVertex& other);

	template <size_t N>
	friend PolyhedronVertex operator- (const PolyhedronVertex& lhs, const PolyhedronVertex& rhs);
	template <size_t N>
	friend PolyhedronVertex operator+ (const PolyhedronVertex& lhs, const PolyhedronVertex& rhs);
	template <size_t N>
	friend PolyhedronVertex operator* (double number, PolyhedronVertex& rhs);
	
	template <size_t N>
	friend bool operator==(const PolyhedronVertex& lhs, const PolyhedronVertex& rhs);
	template <size_t N>
	friend bool operator!= (const PolyhedronVertex& lhs, const PolyhedronVertex& rhs);

	typename array<double, N>::iterator begin();
	typename array<double, N>::iterator end();

	friend double CalculateDistanceBetweenVertices(const PolyhedronVertex& firstVertex, const PolyhedronVertex& secondVertex);

private:
	array<double, N> vertex_;
};


template <size_t N>
PolyhedronVertex<N>::PolyhedronVertex()
{
	vertex_.fill(0.0);
}


template <size_t N>
double PolyhedronVertex<N>::operator[] (int index) const
{
	return vertex_[index];
}


template <size_t N>
double& PolyhedronVertex<N>::operator[] (int index)
{
	return vertex_[index];
}


template <size_t N>
PolyhedronVertex<N>& PolyhedronVertex<N>::operator-= (const PolyhedronVertex& other)
{
	for (size_t i = 0; i < N; ++i)
		vertex_[i] -= other.vertex_[i];
	return *this;
}


template <size_t N>
PolyhedronVertex<N>& PolyhedronVertex<N>::operator+= (const PolyhedronVertex& other)
{
	for (size_t i = 0; i < N; ++i)
		vertex_[i] += other.vertex_[i];
	return *this;
}


template <size_t N>
PolyhedronVertex<N> operator- (const PolyhedronVertex<N>& lhs, const PolyhedronVertex<N>& rhs)
{
	PolyhedronVertex<N> tempVertex = lhs;
	tempVertex -= rhs;
	return tempVertex;
}


template <size_t N>
PolyhedronVertex<N> operator+ (const PolyhedronVertex<N>& lhs, const PolyhedronVertex<N>& rhs)
{
	PolyhedronVertex<N> tempVertex = lhs;
	tempVertex += rhs;
	return tempVertex;
}


template <size_t N>
PolyhedronVertex<N> operator* (double number, PolyhedronVertex<N>& rhs)
{
	PolyhedronVertex<N> tempVertex = rhs;
	for (size_t i = 0; i < N; ++i)
		tempVertex[i] *= number;
	return tempVertex;
}


template <size_t N>
bool operator== (const PolyhedronVertex<N>& lhs, const PolyhedronVertex<N>& rhs)
{
	for (size_t i = 0; i < N; ++i)
	{
		if (lhs[i] != rhs[i])
			return false;
	}

	return true;
}


template <size_t N>
bool operator!= (const PolyhedronVertex<N>& lhs, const PolyhedronVertex<N>& rhs)
{
	return !(lhs == rhs);
}


template <size_t N>
typename array<double, N>::iterator PolyhedronVertex<N>::begin()
{
	return vertex_.begin();
}


template <size_t N>
typename array<double, N>::iterator PolyhedronVertex<N>::end()
{
	return vertex_.end();
}


template <size_t N>
double CalculateDistanceBetweenVertices(const PolyhedronVertex<N>& firstVertex, const PolyhedronVertex<N>& secondVertex)
{
	PolyhedronVertex<N>& tempVertex = secondVertex - firstVertex;
	std::transform(tempVertex.begin(), tempVertex.end(), tempVertex.begin(), tempVertex.begin(), std::multiplies<double>());
	return sqrt(std::accumulate(tempVertex.begin(), tempVertex.end()));
}
