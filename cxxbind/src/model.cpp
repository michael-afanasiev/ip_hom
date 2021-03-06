
#include "includes.hpp"

model::model()
{};

model::~model()
{};

int
model::getWrap(const std::vector<double> &x, const std::vector<int> &con)
{
	double ctrX0 = (x[con[0]-1]+x[con[1]-1]+x[con[2]-1]+x[con[3]-1]) 
		/ float(nNodePerElem);

	size_t itr = 0;
	for (auto i=0; i<con.size()-nNodePerElem; i+=nNodePerElem)
	{	
		double ctrX = (x[con[i+0]-1]+x[con[i+1]-1]
			+x[con[i+2]-1]+x[con[i+3]-1]) / float(nNodePerElem);
		if (std::abs(ctrX - ctrX0) > eps) return itr;
		itr++;
	}
	return -1;
}

int
model::getCtr(const std::vector<double> &x, const std::vector<int> &con)
{
	size_t itr = 0;
	for (auto i=0; i<con.size()-nNodePerElem; i+=nNodePerElem)
	{	
		double ctrX = (x[con[i+0]-1]+x[con[i+1]-1]
			+x[con[i+2]-1]+x[con[i+3]-1]) / float(nNodePerElem);
		if (ctrX < 0) return itr;
		itr++;
	}
	return -1;
}

std::vector<double>
model::calcA(const std::vector<double> &vph, const std::vector<double> &rho)
{
	std::vector<double> a(vph.size());
	for (auto i=0; i<vph.size(); i++) a[i] = rho[i] * vph[i]*vph[i];
	return a;
}

std::vector<double>
model::calcC(const std::vector<double> &vpv, const std::vector<double> &rho)
{
	std::vector<double> c(vpv.size());
	for (auto i=0; i<vpv.size(); i++) c[i] = rho[i] * vpv[i]*vpv[i];
	return c;
}

std::vector<double>
model::calcL(const std::vector<double> &vsv, const std::vector<double> &rho)
{
	std::vector<double> l(vsv.size());
	for (auto i=0; i<vsv.size(); i++) l[i] = rho[i] * vsv[i]*vsv[i];
	return l;
}

std::vector<double>
model::calcF(const std::vector<double> &c, const std::vector<double> &l)
{
	std::vector<double> f(c.size());
	for (auto i=0; i<c.size(); i++) f[i] = c[i] - 2 * l[i];
	return f;
}

std::vector<double>
model::calcMu(const std::vector<double> &vsv, const std::vector<double> &rho)
{
	std::vector<double> mu(vsv.size());
	for (auto i=0; i<vsv.size(); i++) mu[i] = rho[i] * vsv[i]*vsv[i];
	return mu;
}

std::vector<double>
model::calcLambda(const std::vector<double> &vph, 
			      const std::vector<double> &rho, 
			      const std::vector<double> &mu)
{
	std::vector<double> lambda(vph.size());
	for (auto i=0; i<vph.size(); i++)
	{
		lambda[i] = rho[i] * vph[i]*vph[i] - 2 * mu[i];
	}
	return lambda;
}

std::vector<double>
model::calcTheta(const std::vector<double> &mu, 
		         const std::vector<double> &lambda)
{
	std::vector<double> theta(mu.size());
	for (auto i=0; i<mu.size(); i++)
	{
		theta[i] = mu[i] / (lambda[i] + 2 * mu[i]);
	}
	return theta;
}

std::vector<double>
model::calcBigL(const std::vector<double> &mu, const int winLength)
{
	std::vector<double> prep(mu.size());
	for (auto i=0; i<prep.size(); i++)
	{
		prep[i] = 1 / mu[i];
	}
	std::vector<double> win(winLength, 1/float(winLength));
	prep = convolve(prep, win);
	for (auto i=0; i<prep.size(); i++)
	{
		prep[i] = 1 / prep[i];
	}
	return prep;
}

std::vector<double>
model::calcBigM(const std::vector<double> &mu, const int winLength)
{	
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> prep = convolve(mu, win);
	return prep;
}

std::vector<double>
model::calcBigR(const std::vector<double> &mu, 
			    const std::vector<double> &theta, const int winLength)
{
	std::vector<double> prep(mu.size());
	for (auto i=0; i<prep.size(); i++)
	{
		prep[i] = theta[i] / mu[i];
	}
	std::vector<double> win(winLength, 1/float(winLength));
	prep = convolve(prep, win);
	return prep;
}

std::vector<double>
model::calcBigS(const std::vector<double> &mu,
				const std::vector<double> &theta, const int winLength)
{
	std::vector<double> prep(mu.size());
	for (auto i=0; i<prep.size(); i++) {
		prep[i] = theta[i] * mu[i];
	}
	std::vector<double> win(winLength, 1/float(winLength));
	prep = convolve(prep, win);
	return prep;
}

std::vector<double>
model::calcBigT(const std::vector<double> &theta, const int winLength)
{
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> prep = convolve(theta, win);
	return prep;
}

double 
model::dBigLdMu(const std::vector<double> &mu,
				const std::vector<double> &theta, 
				const std::vector<double> &bigL,
				const int &ind,
				const int &winLength)
{
	double out = 0.0;
	int half = winLength / 2;
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> bigLtrial = calcBigL(mu, winLength);
	std::vector<double> prep(mu.size());
	for (auto i=0; i<prep.size(); i++)
	{
		prep[i] = 1 / mu[i];
	}
	prep = convolve(prep, win);
	for (auto i=0; i<prep.size(); i++)
	{
		prep[i] = 1 / (prep[i]*prep[i]);
	}
	for (auto j=0; j<win.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		out += ((bigLtrial[vecInd] - bigL[vecInd]) * prep[vecInd] 
			* (1 / (mu[vecInd]*mu[vecInd]))) * win[j];
	}
	return out;
}

double 
model::dBigMdMu(const std::vector<double> &mu,
				const std::vector<double> &theta, 
				const std::vector<double> &bigM,
				const int &ind,
				const int &winLength)
{
	double out = 0.0;
	int half = winLength / 2;
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> bigMtrial = calcBigM(mu, winLength);
	for (auto j=0; j<win.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		out += (bigMtrial[vecInd] - bigM[vecInd]) * win[j];
	}
	return out;
}

double 
model::dBigRdMu(const std::vector<double> &mu,
				const std::vector<double> &theta, 
				const std::vector<double> &bigR,
				const int &ind,
				const int &winLength)
{
	double out = 0.0;
	int half = winLength / 2;
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> bigRtrial = calcBigR(mu, theta, winLength);
	for (auto j=0; j<win.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		out += (bigRtrial[vecInd] - bigR[vecInd]) * (1 * theta[j]);
	}
	return out;
}

double 
model::dBigSdMu(const std::vector<double> &mu,
				   const std::vector<double> &theta, 
				   const std::vector<double> &bigS,
				   const int &ind,
				   const int &winLength)
{
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> bigStrial = calcBigS(mu, theta, winLength);

	int half = win.size() / 2;
	double out = 0.0;
	for (auto j=0; j<win.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		out += (bigStrial[vecInd] - bigS[vecInd]) * theta[vecInd] * win[j];
	}
	return out;
}

double 
model::dBigRdTheta(const std::vector<double> &mu,
				   const std::vector<double> &theta, 
				   const std::vector<double> &bigR,
				   const int &ind,
				   const int &winLength)
{
	double out = 0.0;
	int half = winLength / 2;
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> bigRtrial = calcBigR(mu, theta, winLength);
	for (auto j=0; j<win.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		out += (bigRtrial[vecInd] - bigR[vecInd]) * (1 / mu[j]);
	}
	return out;
}

double 
model::dBigSdTheta(const std::vector<double> &mu,
				   const std::vector<double> &theta, 
				   const std::vector<double> &bigS,
				   const int &ind,
				   const int &winLength)
{
	double out = 0.0;
	int half = winLength / 2;
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> bigStrial = calcBigS(mu, theta, winLength);
	for (auto j=0; j<win.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		out += (bigStrial[vecInd] - bigS[vecInd]) * mu[vecInd] * win[j];
	}
	return out;
}

double 
model::dBigTdTheta(const std::vector<double> &mu,
				   const std::vector<double> &theta, 
				   const std::vector<double> &bigT,
				   const int &ind,
				   const int &winLength)
{
	double out = 0.0;
	int half = winLength / 2;
	std::vector<double> win(winLength, 1/float(winLength));
	std::vector<double> bigTtrial = calcBigT(theta, winLength);
	for (auto j=0; j<win.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		out += (bigTtrial[vecInd] - bigT[vecInd]) * win[j];
	}
	return out;
}

double 
model::dLdMu(const std::vector<double> &mu,
			 const std::vector<double> &theta, 
			 const int &ind,
			 const int &winLength)
{
	int half = winLength / 2;
	double sum = 0.0;
	std::vector<double> win(winLength, 1/double(winLength));
	std::vector<double> prep(winLength);
	for (auto j=0; j<prep.size(); j++)
	{
		int vecInd = ind + 1 + (j - half);
		if (vecInd < 0) vecInd = 0;
		if (vecInd >= mu.size()) vecInd = mu.size() - 1;
		sum += 1 / mu[vecInd] * win[j];
	}
	sum = 1 / (sum * sum);
	return sum * 1 / (mu[ind]*mu[ind]);
}

double
model::dMdMu(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int &ind)
{
	return 1;
}

double
model::dRdMu(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int &ind)
{
	return (-1) * (theta[ind] / (mu[ind]*mu[ind]));
}

double
model::dSdMu(const std::vector<double> &mu, const std::vector<double> &theta,
			 const int &ind)
{
	return theta[ind];
}

double
model::dRdTheta(const std::vector<double> &mu, const std::vector<double> &theta,
			 	const int &ind)
{
	return 1 / mu[ind];
}

double
model::dSdTheta(const std::vector<double> &mu, const std::vector<double> &theta,
			 	const int &ind)
{
	return mu[ind];
}

double
model::dTdTheta(const std::vector<double> &mu, const std::vector<double> &theta,
			 	const int &ind)
{
	return 1;
}

int
model::numJacNonZero(const std::vector<double> &tmp, const int &nCon, 
					 const int &nTyp, const int &nLoc, const int &winLength)
{
	auto jInd = 0;
	for (auto i=0; i<nCon; i++)
	{
		for (auto j=0; j<nTyp; j++)
		{
			if (! model::conSen(i, j)) continue;
			for (auto k=0; k<nLoc; k++)
			{
				std::vector<int> cInd = convolveIndices(tmp, winLength, k);
				for (auto l : cInd)
				{
					if (l >= 0 && l<nLoc)
					{
						jInd++;
					}
				}
			}
		}
	}
	return jInd;
}

bool 
model::conSen(const int &con, const int &var)
{
	if ((con == 0 || con == 1 || con == 2 || con == 3) && (var == 0)) 
	{
		return true;
	}
	if ((con == 2 || con == 3 || con == 4) && (var == 2))
	{
		return true;
	}
	return false;
}