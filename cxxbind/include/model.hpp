#ifndef __homOpt_model_hpp__
#define __homOpt_model_hpp__

#include <iostream>
#include <vector>

class model
{

public:

	model();
	~model();

	static bool
	conSen(const int &constraint, const int & var);

	static int
	getWrap(const std::vector<double> &x, const std::vector<int> &con);

	static int
	getCtr(const std::vector<double> &x, const std::vector<int> &con);

	static int
	numJacNonZero(const std::vector<double> &tmp, const int &nCon,
			      const int &nTyp, const int &nLoc, const int &winLength);

	static double
	dLdMu(const std::vector<double> &mu, const std::vector<double> &theta,
		  const int &ind);

	static double
	dMdMu(const std::vector<double> &mu, const std::vector<double> &theta,
		  const int &ind);

	static double
	dRdMu(const std::vector<double> &mu, const std::vector<double> &theta,
		  const int &ind);

	static double
	dSdMu(const std::vector<double> &mu, const std::vector<double> &theta,
		  const int &ind);

	static double
	dRdTheta(const std::vector<double> &mu, const std::vector<double> &theta,
		  const int &ind);

	static double
	dSdTheta(const std::vector<double> &mu, const std::vector<double> &theta,
		  const int &ind);

	static double
	dTdTheta(const std::vector<double> &mu, const std::vector<double> &theta,
		  const int &ind);

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
	dBigLdMu(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int &winLength);

	static std::vector<double>
	dBigMdMu(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int &winLength);

	static std::vector<double>
	dBigRdMu(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int &winLength);

	static double
	dBigSdMu(const std::vector<double> &mu,
			 const std::vector<double> &theta,
			 const std::vector<double> &bigS,
			 const int &ind,
			 const int &winLength);

	static std::vector<double>
	dBigRdTheta(const std::vector<double> &mu, 
			    const std::vector<double> &theta,
			    const int &winLength);

	static double
	dBigSdTheta(const std::vector<double> &mu,
				const std::vector<double> &theta,
				const std::vector<double> &bigS,
				const int &ind,
				const int &winLength);

	static std::vector<double>
	dBigTdTheta(const std::vector<double> &mu,
			 	const std::vector<double> &theta,
			 	const int &winLength);

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