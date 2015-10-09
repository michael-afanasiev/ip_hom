#ifndef __homOpt_utilties_hpp__
#define __homOpt_utilties_hpp__

#include <iostream>
#include <algorithm>

template <typename T> void
printDebug(T msg)
{
	std::cout << msg << std::endl;
}

template <typename A, typename B> auto
getIndForVal(A name, std::vector<B> vec)
{
	auto it=std::find(vec.begin(), vec.end(), name);
	return it-vec.begin();
}

template <typename T> void
printVector(std::vector<T> vec)
{
	for_each(vec.begin(), vec.end(), [] (T val)
	{
		std::cout << val << std::endl;
	});
}

template <typename T> void
printVectorRange(std::vector<T> vec, int s, int e)
{
	for_each(vec.begin()+s, vec.begin()+e, [] (T val)
	{
		std::cout << val << std::endl;
	});
}

template <typename T> std::vector<T>
getVectorSubset(std::vector<T> vec, int s, int e)
{
	std::vector<T> v(vec.begin()+s, vec.begin()+e);
	return v;
}

template <typename T> int
countConvolveHits(std::vector<T> vec, int winLength)
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

template <typename T> std::vector<T>
convolve(std::vector<T> vec, std::vector<T> win)
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

#endif //__homOpt_utilties_hpp__