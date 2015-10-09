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

#endif //__homOpt_utilties_hpp__