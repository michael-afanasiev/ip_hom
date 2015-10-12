#ifndef __homOpt_utilties_hpp__
#define __homOpt_utilties_hpp__

#include <cmath>
#include <iostream>
#include <algorithm>

template <typename T> void
printDebug(T msg)
{
	std::cout << msg << std::endl;
}

template <typename T> void
printError(T msg)
{
	std::cout << msg << std::endl;
	exit(EXIT_FAILURE);
}

template <typename A, typename B> int
getIndForVal(A name, const std::vector<B> &vec)
{
	auto it=std::find(vec.begin(), vec.end(), name);
	if (it == vec.end()) return -1;
	return it-vec.begin();
}

template <typename T> void
printVector(const std::vector<T> &vec)
{
	for_each(vec.begin(), vec.end(), [] (T val)
	{
		std::cout << val << std::endl;
	});
}

template <typename T> void
printVectorRange(const std::vector<T> &vec, int s, int e)
{
	for_each(vec.begin()+s, vec.begin()+e, [] (T val)
	{
		std::cout << val << std::endl;
	});
}

template <typename T> std::vector<T>
getVectorSubset(const std::vector<T> &vec, int s, int e)
{
	std::vector<T> v(vec.begin()+s, vec.begin()+e);
	return v;
}

template <typename T> int
countConvolveHits(const std::vector<T> &vec, const int winLength)
{
	int hits = 0;
	std::vector<double> win(winLength, 1/float(winLength));
	int half = win.size() / 2;
	for (auto i=0; i<vec.size(); i++)
	{
		for (auto j=0; j<win.size(); j++)
		{
			int vecInd = i + (j - half);
			if (vecInd < 0) continue;
			if (vecInd >= vec.size()) continue;
			hits++;
		}
	}
	return hits;
}

template <typename T> std::vector<int>
convolveIndices(const std::vector<T> &vec, const int &winLength, 
			    const int &pInd)
{
	std::vector<int> index(winLength);
	int half = index.size() / 2;
	for (auto j=0; j<index.size(); j++)
	{
		int vecInd = pInd + (j - half);
		index[j] = vecInd;
	}
	return index;
}

template <typename T> std::vector<T>
convolve(const std::vector<T> &vec, const std::vector<T> &win)
{
	int half = win.size() / 2;
	std::vector<T> out(vec.size());
	for (auto i=0; i<vec.size(); i++)
	{
		for (auto j=0; j<win.size(); j++)
		{
			int vecInd = i + (j - half);
			if (vecInd < 0) vecInd = 0;
			if (vecInd >= vec.size()) vecInd = vec.size() - 1;
			out[i] += vec[vecInd] * win[j];
		}
	}
	return out;
}

template <typename T> auto
l2norm(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
	auto arg = 0;
	for (auto i=0; i<vec1.size(); i++)
	{
		arg += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
	}
	return (1.0/2.0) * arg;
}

template <typename T> auto
vecDiff(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
	std::vector<double> diff(vec1.size());
	for (auto i=0; i<vec1.size(); i++)
	{
		diff[i] = vec1[i] - vec2[i];
	}
	return diff;
}

#endif //__homOpt_utilties_hpp__