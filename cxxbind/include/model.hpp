#ifndef __homOpt_model_hpp__
#define __homOpt_model_hpp__

#include <iostream>
#include <vector>

class model
{

public:

	model();
	~model();

	static int
	getWrap(const std::vector<double> &x, const std::vector<int> &con);

	static int
	getCtr(const std::vector<double> &x, const std::vector<int> &con);

	static std::vector<double>
	calcA(const std::vector<double> &vph, const std::vector<double> &vpv);

	static std::vector<double>
	calcC(const std::vector<double> &vpv, const std::vector<double> &rho);

	static std::vector<double>
	calcL(const std::vector<double> &vsv, const std::vector<double> &rho);

	static std::vector<double>
	calcF(const std::vector<double> &c, const std::vector<double> &l);

	static std::vector<double>
	calcMu(const std::vector<double> &vsv, const std::vector<double> &rho);

	static std::vector<double>
	calcBigL(const std::vector<double> &mu, const int winLength);

	static std::vector<double>
	calcBigM(const std::vector<double> &mu, const int winLength);

	static std::vector<double>
	calcBigR(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int winLength);

	static std::vector<double>
	calcBigS(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int winLength);

	static std::vector<double>
	calcBigT(const std::vector<double> &theta, const int winLength);

	static std::vector<double>
	calcLambda(const std::vector<double> &vph, 
			   const std::vector<double> &rho, 
			   const std::vector<double> &mu);

	static std::vector<double>
	calcTheta(const std::vector<double> &mu, 
			  const std::vector<double> &lambda);

private:

	static int constexpr nNodePerElem = 4;
	static double constexpr eps = 1e-3;

};

#endif // __homOpt_model_hpp__