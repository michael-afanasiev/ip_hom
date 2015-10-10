/*

	This guy tries to solve the inverse homogenization problem.

	Variable ordering: mu, lambda, theta, where the loop is first over
	co-ordinates.

	Constraint ordering: L, M, R, S, T, where the first loop is also
	over co-ordinates.

*/

#include "includes.hpp"

exoMod1D_Nlp::exoMod1D_Nlp()
{};

exoMod1D_Nlp::~exoMod1D_Nlp()
{};

bool
exoMod1D_Nlp::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
						   Index &nnz_h_lag, IndexStyleEnum &index_style)
{

	// problem has nWrap*3 variables (one whole column, and lambda, mu, theta)
	n = mNvar;

	// each var has 5 constraints (L, M, R, S, T)
	m = mNloc * mNcon;
	
	// number of non-zeros in the jacobian (pre-calculated by 
	// countConvolveHits)
	nnz_jac_g = mNjacNonZero;

	// use C style indexing (0-based)
	index_style = TNLP::C_STYLE;

	return true;
}

bool
exoMod1D_Nlp::get_bounds_info(Index n, Number *x_l, Number *x_u,
						      Index m, Number *g_l, Number *g_u)
{
	// just assert to make sure values make sense
	assert(n == mNvar);
	assert(m == mNloc * mNcon);

	for (Index i=0; i<mNtyp; i++)
	{
		for (Index j=0; j<mNloc; j++)
		{
			error(x_l[j+i*mNloc] = lowerBound(i), "lowerBound");
			error(x_u[j+i*mNloc] = upperBound(i), "upperBound");
		}
	}
	for (Index i=0; i<mNcon; i++)
	{
		for (Index j=0; j<mNloc; j++)
		{
			error(g_l[j+i*mNloc] = constraint(i, j) - 0.10*constraint(i,j), "constraint");
			error(g_u[j+i*mNloc] = constraint(i, j) + 0.10*constraint(i,j), "constraint");
		}
	}
	return true;
}

bool
exoMod1D_Nlp::get_starting_point(Index n, bool init_x, Number* x,
			  		   		     bool init_z, Number* z_L, Number* z_U,
			        		   	 Index m, bool init_lambda, Number* lambda)
{
	std::copy(mMu.begin(), mMu.end(), x);
	std::copy(mLambda.begin(), mLambda.end(), x+mNloc);
	std::copy(mTheta.begin(), mTheta.end(), x+2*mNloc);
	return true;
}

bool
exoMod1D_Nlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	std::vector<double> mu(x, x+mNloc);
	std::vector<double> lambda(x+mNloc, x+2*mNloc);
	std::vector<double> theta(x+2*mNloc, x+3*mNloc);

	obj_value = 0;
	obj_value += l2norm(model::calcBigL(mu, mWinLength), mBigL);
	obj_value += l2norm(model::calcBigM(mu, mWinLength), mBigM);
	obj_value += l2norm(model::calcBigR(mu, theta, mWinLength), mBigR);
	obj_value += l2norm(model::calcBigS(mu, theta, mWinLength), mBigS);
	obj_value += l2norm(model::calcBigT(theta, mWinLength), mBigT);
	return true;
}	

bool
exoMod1D_Nlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	std::vector<double> mu(x, x+mNloc);
	std::vector<double> lambda(x+mNloc, x+2*mNloc);
	std::vector<double> theta(x+2*mNloc, x+3*mNloc);

	for (Index i=0; i<mNtyp; i++)
	{
		std::vector<double> comp = calcGrad(mu, lambda, theta, i);
		std::copy(comp.begin(), comp.end(), grad_f+i*mNloc);
	}

	return true;
}

bool
exoMod1D_Nlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	std::vector<double> mu(x, x+mNloc);
	std::vector<double> lambda(x+mNloc, x+2*mNloc);
	std::vector<double> theta(x+2*mNloc, x+3*mNloc);
	for (Index i=0; i<mNcon; i++)
	{
		std::vector<double> comp = calcConstraint(mu, lambda, theta, i);
		std::copy(comp.begin(), comp.end(), g+i*mNloc);
	}
	return true;

}

bool
exoMod1D_Nlp::eval_jac_g(Index n, const Number* x, bool new_x,
			   		     Index m, Index nele_jac, Index* iRow, Index *jCol,
			   		     Number* values)
{return true;}

bool
exoMod1D_Nlp::eval_h(Index n, const Number* x, bool new_x,
		   		     Number obj_factor, Index m, const Number* lambda,
		   		     bool new_lambda, Index nele_hess, Index* iRow,
		   		     Index* jCol, Number* values)
{return false;}

void exoMod1D_Nlp::finalize_solution(SolverReturn status,
					  			     Index n, const Number* x, const Number *z_L,
					  			     const Number *z_U, Index m, const Number *g,
					  			     const Number *lambda, Number obj_value,
					  			     const IpoptData *ip_data,
					  			     IpoptCalculatedQuantities *ip_cq)
{}

void exoMod1D_Nlp::initFromExodus(exodusFile &exo)
{
	// read from exodus
	mVsv = exo.getVar("vsv");
	mVph = exo.getVar("vph");
	mVpv = exo.getVar("vpv");
	mRho = exo.getVar("rho");
	mX = exo.X();
	mY = exo.Y();

	// calculate love parameters
	mA = model::calcA(mVph, mRho);
	mC = model::calcC(mVpv, mRho);
	mL = model::calcL(mVsv, mRho);
	mF = model::calcF(mC, mL);

	// calculate lame params and theta
	mMu = model::calcMu(mVsv, mRho);
	mLambda = model::calcLambda(mVph, mRho, mMu);
	mTheta = model::calcTheta(mMu, mLambda);

	// get 1-D wrap column
	std::vector<int> con = exo.con();
	int nWrap = model::getWrap(mX, con);
	int nCtr = model::getCtr(mX, con);

	// get 1-D subset vectors
	mA = getVectorSubset(mA, nCtr, nCtr+nWrap);
	mC = getVectorSubset(mC, nCtr, nCtr+nWrap);
	mL = getVectorSubset(mL, nCtr, nCtr+nWrap);
	mF = getVectorSubset(mF, nCtr, nCtr+nWrap);
	mMu = getVectorSubset(mMu, nCtr, nCtr+nWrap);
	mVsv = getVectorSubset(mVsv, nCtr, nCtr+nWrap);
	mVph = getVectorSubset(mVph, nCtr, nCtr+nWrap);
	mVpv = getVectorSubset(mVpv, nCtr, nCtr+nWrap);
	mRho = getVectorSubset(mRho, nCtr, nCtr+nWrap);
	mTheta = getVectorSubset(mTheta, nCtr, nCtr+nWrap);
	mLambda = getVectorSubset(mLambda, nCtr, nCtr+nWrap);

	// calc backus parameters
	int winLength = 20;
	mBigL = model::calcBigL(mMu, winLength);
	mBigM = model::calcBigM(mMu, winLength);
	mBigR = model::calcBigR(mMu, mTheta, winLength);
	mBigS = model::calcBigS(mMu, mTheta, winLength);
	mBigT = model::calcBigT(mTheta, winLength);

	// save some variables
	mNtyp = 3;
	mNcon = 5;
	mNloc = nWrap;
	mNvar = mNloc * mNtyp;
	mWinLength = winLength;
	mNjacNonZero = countConvolveHits(mMu, winLength);
}

Number
exoMod1D_Nlp::lowerBound(const int &i)
{
	if (i == 0)
	{
		return 0;
	}
	else if (i == 1)
	{
		return 0;
	}
	else if (i == 2)		
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

Number
exoMod1D_Nlp::upperBound(const int &i)
{
	if (i == 0)
	{
		return 1000;
	}
	else if (i == 1)
	{
		return 1000;
	}
	else if (i == 2)
	{
		return 3.0/4.0;
	}
	else
	{
		return -1;
	}
}

Number
exoMod1D_Nlp::constraint(const int &i, const int &ind)
{
	if (i == 0)
	{
		return mBigL[ind];
	}
	else if (i == 1)
	{
		return mBigM[ind];
	}
	else if (i == 2)
	{
		return mBigR[ind];
	}
	else if (i == 3)
	{
		return mBigS[ind];
	}
	else if (i == 4)
	{
		return mBigT[ind];
	}
	else
	{
		return -1;
	}
}

std::vector<double>
exoMod1D_Nlp::calcGrad(const std::vector<double> &mu, 
					   const std::vector<double> &lambda,
					   const std::vector<double> &theta,
					   const int &var)
{	
	if (var == 0)
	{
		// initialize vectors that will be needed
		std::vector<double> win(mWinLength, 1/float(mWinLength));
		std::vector<double> t00(mu.size());
		std::vector<double> t01(mu.size());
		std::vector<double> t1(mu.size(), 1.0);
		std::vector<double> t2(mu.size());
		std::vector<double> t3(mu.size());
		for (auto i=0; i<t00.size(); i++)
		{
			t00[i] = 1 / mu[i];
			t01[i] = 1 / (mu[i]*mu[i]);
			t2[i] = theta[i] / (mu[i]*mu[i]);
			t3[i] = theta[i];
		}

		// compute the 'homogenized derivatives'
		t00 = convolve(t00, win);
		t01 = convolve(t01, win);
		t2 = convolve(t2, win);
		t3 = convolve(t3, win);

		// term 1 needs to be squared
		for (auto i=0; i<t00.size(); i++)
		{
			t00[i] = 1 / t00[i]*t00[i];
		}

		// remember these linear terms are left over from the L2 norm
		// TODO check if these need to be convolutions of the difference
		std::vector<double> mul0 = vecDiff(model::calcBigL(mu, 
														   mWinLength), 
							   							   mBigL);
		std::vector<double> mul1 = vecDiff(model::calcBigM(mu, 
														   mWinLength),
														   mBigM);
		std::vector<double> mul2 = vecDiff(model::calcBigR(mu,
														   theta,
														   mWinLength),
														   mBigR);
		std::vector<double> mul3 = vecDiff(model::calcBigS(mu,
														   theta,
														   mWinLength),
														   mBigS);

		// sum up all values and return
		std::vector<double> val(mu.size());
		for (auto i=0; i<t00.size(); i++)
		{
			val[i] = mul0[i] * t00[i]*t01[i] + mul1[i] - mul2[i] * t2[i] +
				mul3[i] * t3[i];
		}

		return val;
	}
	else if (var == 1)
	{
		std::vector<double> val(mu.size(), 0.0);
		return val;
	}
	else if (var == 2)
	{
		// initialize vectors that will be needed
		std::vector<double> win(mWinLength, 1/float(mWinLength));
		std::vector<double> t0(mu.size());
		std::vector<double> t1(mu.size());
		for (auto i=0; i<t0.size(); i++)
		{
			t0[i] = 1 / mu[i];
			t1[i] = mu[i];
		}

		// compute the 'homogenized derivatives'
		t0 = convolve(t0, win);
		t1 = convolve(t1, win);

		// remember these linear terms are left over from the L2 norm
		// TODO check if these need to be convolutions of the difference
		std::vector<double> mul0 = vecDiff(model::calcBigR(mu,
														   theta,
														   mWinLength),
		 												   mBigR);
		std::vector<double> mul1 = vecDiff(model::calcBigS(mu,
														   theta,
														   mWinLength),
														   mBigS);
		std::vector<double> mul2 = vecDiff(model::calcBigT(theta,
														   mWinLength),
														   mBigT);
		// sum up all values and return
		std::vector<double> val(mu.size());
		for (auto i=0; i<t0.size(); i++)
		{
			val[i] = mul0[i] * t0[i] + mul1[i] * t1[i] + mul2[i];
		}
		return val;
	}
	std::vector<double> err(mu.size(), -1);
	return err;
}

std::vector<double>
exoMod1D_Nlp::calcConstraint(const std::vector<double> &mu,
							 const std::vector<double> &lambda,
							 const std::vector<double> &theta,
							 const int &var)
{
	if (var == 0)
	{
		return model::calcBigL(mu, mWinLength);
	}
	else if (var == 1)
	{
		return model::calcBigM(mu, mWinLength);
	}
	else if (var == 2)
	{
		return model::calcBigR(mu, theta, mWinLength);
	}
	else if (var == 3)
	{
		return model::calcBigS(mu, theta, mWinLength);
	}
	else if (var == 4)
	{
		return model::calcBigT(theta, mWinLength);
	}
	return std::vector<double> (mu.size(), -1);

}

void
exoMod1D_Nlp::error(const int &e, const std::string &func)
{
	if (e < 0)
	{
		std::cout << "Failure in: " << func << std::endl;
		exit(EXIT_FAILURE);
	}
}