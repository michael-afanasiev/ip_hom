
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
	for (auto i=0; i<vph.size(); i++)
	{
		a[i] = rho[i] * vph[i]*vph[i];
	}
	return a;
}

std::vector<double>
model::calcC(const std::vector<double> &vpv, const std::vector<double> &rho)
{
	std::vector<double> c(vpv.size());
	for (auto i=0; i<vpv.size(); i++)
	{
		c[i] = rho[i] * vpv[i]*vpv[i];
	}
	return c;
}

std::vector<double>
model::calcL(const std::vector<double> &vsv, const std::vector<double> &rho)
{
	std::vector<double> l(vsv.size());
	for (auto i=0; i<vsv.size(); i++)
	{
		l[i] = rho[i] * vsv[i]*vsv[i];
	}
	return l;
}

std::vector<double>
model::calcF(const std::vector<double> &c, const std::vector<double> &l)
{
	std::vector<double> f(c.size());
	for (auto i=0; i<c.size(); i++)
	{
		f[i] = c[i] - 2 * l[i];
	}
	return f;
}

std::vector<double>
model::calcMu(const std::vector<double> &vsv, const std::vector<double> &rho)
{
	std::vector<double> mu(vsv.size());
	for (auto i=0; i<vsv.size(); i++)
	{
		mu[i] = rho[i] * vsv[i]*vsv[i];
	}
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