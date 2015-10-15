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

	// 741 in total
	// problem has nWrap*3 variables (one whole column, and lambda, mu, theta)
	n = mNvar;

	// each location has 5 constraints (L, M, R, S, T)
	m = mNloc * mNcon;
	
	// number of non-zeros in the jacobian (pre-calculated by 
	// countConvolveHits)
	nnz_jac_g = mNjacNonZero;

	nnz_h_lag = 0;

	// use C style indexing (0-based)
	index_style = TNLP::C_STYLE;

	// iteration count
	mItr = 0;

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
			error(g_l[j+i*mNloc] = constraint(i, j) - 0.001*constraint(i, j), "constraint");
			error(g_u[j+i*mNloc] = constraint(i, j) + 0.001*constraint(i, j), "constraint");
		}
	}
	return true;
}

bool
exoMod1D_Nlp::get_starting_point(Index n, bool init_x, Number* x,
			  		   		     bool init_z, Number* z_L, Number* z_U,
			        		   	 Index m, bool init_lambda, Number* lambda)
{
	// sMu = vecPerturb(mMu, 0.1);
	// sTheta = vecPerturb(mTheta, 0.1);
	sMu = vecAvg(mMu);
	sTheta = vecAvg(mTheta);
	// std::vector<double> tmp(mMu.size(), 0.0);
	// std::vector<double> tmp1(mTheta.size(), 0.0);
	// sMu = tmp;
	// sTheta = tmp1;
	std::copy(sMu.begin(), sMu.end(), x);
	std::copy(mLambda.begin(), mLambda.end(), x+mNloc);
	std::copy(sTheta.begin(), sTheta.end(), x+2*mNloc);

	writeParam(mMu, "./dump/mu_true.txt");
	writeParam(mTheta, "./dump/theta_true.txt");
	return true;
}

bool
exoMod1D_Nlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	std::vector<double> mu(x, x+mNloc);
	std::vector<double> lambda(x+mNloc, x+2*mNloc);
	std::vector<double> theta(x+2*mNloc, x+3*mNloc);
	obj_value = 0;
	// obj_value += l2norm(model::calcBigL(mu, mWinLength), mBigL);
	obj_value += l2norm(model::calcBigM(mu, mWinLength), mBigM);
	obj_value += l2norm(model::calcBigR(mu, theta, mWinLength), mBigR);
	obj_value += l2norm(model::calcBigS(mu, theta, mWinLength), mBigS);
	obj_value += l2norm(model::calcBigT(theta, mWinLength), mBigT);

	writeParam(mu, "./dump/mu_" + numToStr(mItr) + ".txt");
	writeParam(theta, "./dump/theta_" + numToStr(mItr) + ".txt");
	mItr++;
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
		for (Index j=0; j<mNloc; j++)
		{
			grad_f[j+i*mNloc] = calcGrad(mu, lambda, theta, i, j);
		}
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
{
	if (values == NULL)
	{
		Index jInd = 0;
		for (auto i=0; i<mNcon; i++)
		{
			for (auto j=0; j<mNtyp; j++)
			{
				if (! model::conSen(i, j)) continue;
				for (auto k=0; k<mNloc; k++)
				{
					std::vector<int> cInd = convolveIndices(mMu, 
														    mWinLength, 
														    k);
					for (auto l : cInd)
					{
						if (l >= 0 && l < mNloc)
						{
							iRow[jInd] = i * (mNloc) + k;
							jCol[jInd] = j * mNloc + l;
							mIrow[jInd] = i * (mNloc) + k;
							mJcol[jInd] = j * mNloc + l;
							jInd++;
						}
					}
				}
			}
		}
	}
	else
	{
		std::vector<double> mu(x, x+mNloc);
		std::vector<double> lambda(x+mNloc, x+2*mNloc);
		std::vector<double> theta(x+2*mNloc, x+3*mNloc);
		for (Index i=0; i<mNjacNonZero; i++)
		{
			values[i] = (1 / float(mWinLength)) * 
				calcConstraintDeriv(mu, lambda, theta, mIrow[i], mJcol[i]);
		}
	}
	return true;
}

bool
exoMod1D_Nlp::eval_h(Index n, const Number* x, bool new_x,
		   		     Number obj_factor, Index m, const Number* lambda,
		   		     bool new_lambda, Index nele_hess, Index* iRow,
		   		     Index* jCol, Number* values)
{
	return false;
}

void exoMod1D_Nlp::finalize_solution(SolverReturn status,
					  			     Index n, const Number* x, const Number *z_L,
					  			     const Number *z_U, Index m, const Number *g,
					  			     const Number *lambda, Number obj_value,
					  			     const IpoptData *ip_data,
					  			     IpoptCalculatedQuantities *ip_cq)
{
	std::cout << std::endl << std::endl 
	<< "Solution of the primal variables, x" << std::endl;
	for (Index i=0; i<n; i++)
	{
		std::cout << "x[" <<i << "] = " << x[i] << std::endl;
	}

	std::cout << std::endl << std::endl
	<< "Solution of the bound multipliers, z_L and z_U" << std::endl;
	for (Index i=0; i<n; i++)
	{
		std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
	}
	for (Index i=0; i<n; i++)
	{
		std::cout << "z_U[" << i << "] = " << z_L[i] << std::endl;
	}

	std::cout << std::endl << std::endl
	<< "Objective Value" << std::endl;
	std::cout << "f(x*) = " << obj_value << std::endl;

	std::cout << std::endl 
	<< "Final value of the constraints:" << std::endl;
	for (Index i=0; i<m; i++)
	{
		std::cout << "g(" << i << ") = " << g[i] << std::endl;
	}

	std::vector<double> mu(x, x+mNloc);
	std::vector<double> theta(x+2*mNloc, x+3*mNloc);

	std::ofstream file;
	file.open("true_mu.txt");
	for (auto i=0; i<mMu.size(); i++)
	{
		file << mMu[i] << "\n";
	}
	file.close();

	file.open("start_mu.txt");
	for (auto i=0; i<sMu.size(); i++)
	{
		file << sMu[i] << "\n";
	}
	file.close();

	file.open("solution_mu.txt");
	for (auto i=0; i<mu.size(); i++)
	{
		file << mu[i] << "\n";
	}
	file.close();

	file.open("true_theta.txt");
	for (auto i=0; i<mTheta.size(); i++)
	{
		file << mTheta[i] << "\n";
	}
	file.close();

	file.open("start_theta.txt");
	for (auto i=0; i<sTheta.size(); i++)
	{
		file << sTheta[i] << "\n";
	}
	file.close();

	file.open("solution_theta.txt");
	for (auto i=0; i<theta.size(); i++)
	{
		file << theta[i] << "\n";
	}
	file.close();
}

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
	mNjacNonZero = model::numJacNonZero(mMu, mNcon, mNtyp, mNloc, mWinLength);
	mJcol.resize(mNjacNonZero);
	mIrow.resize(mNjacNonZero);
}

Number
exoMod1D_Nlp::lowerBound(const int &i)
{
	if (i == 0)
	{
		return 10;
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
		return 20;
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

Number
exoMod1D_Nlp::calcGrad(const std::vector<double> &mu, 
					   const std::vector<double> &lambda,
					   const std::vector<double> &theta,
					   const int &var, const int &ind)
{	
	if (var == 0)
	{
		double dMdM = model::dBigMdMu(mu, theta, mBigM, ind, mWinLength);
		double dRdM = model::dBigRdMu(mu, theta, mBigR, ind, mWinLength);
		double dSdM = model::dBigSdMu(mu, theta, mBigS, ind, mWinLength);
		return dSdM + dRdM + dMdM;
	}
	else if (var == 1)
	{
		return 0.0;
	}
	else if (var == 2)
	{
		double dRdT = model::dBigRdTheta(mu, theta, mBigR, ind, mWinLength);
		double dTdT = model::dBigTdTheta(mu, theta, mBigT, ind, mWinLength);
		double dSdT = model::dBigSdTheta(mu, theta, mBigS, ind, mWinLength);
		return dSdT+ dTdT + dRdT;
	}
	return - 1;
}

std::vector<double>
exoMod1D_Nlp::calcConstraint(const std::vector<double> &mu,
							 const std::vector<double> &lambda,
							 const std::vector<double> &theta,
							 const int &var)
{
	if (var == 0)
	{
		std::vector<double> res(mu.size(), 1.0);
		return model::calcBigL(mu, mWinLength);
	}
	else if (var == 1)
	{
		std::vector<double> res(mu.size(), 1.0);
		return model::calcBigM(mu, mWinLength);
	}
	else if (var == 2)
	{
			std::vector<double> res(mu.size(), 1.0);

		return model::calcBigR(mu, theta, mWinLength);
	}
	else if (var == 3)
	{
		std::vector<double> res(mu.size(), 1.0);
		return model::calcBigS(mu, theta, mWinLength);
	}
	else if (var == 4)
	{
		std::vector<double> res(mu.size(), 1.0);
		return model::calcBigT(theta, mWinLength);
	}
	return std::vector<double> (mu.size(), -1);

}

double
exoMod1D_Nlp::calcConstraintDeriv(const std::vector<double> &mu,
							 	  const std::vector<double> &lambda,
							 	  const std::vector<double> &theta,
							 	  const int &con, const int &ind)
{
	int var;
	int locCon;
	for (auto i=0; i<mNtyp; i++)
	{
		if ((ind >= i*mNloc) && (ind < (i+1)*mNloc)) var = i;
	}
	for (auto i=0; i<mNcon; i++)
	{
		if ((con >= i*mNloc) && (i < (i+1)*mNloc)) locCon = i;
	}
	int locInd = ind % mNloc;
	if ((locCon == 0) && (var == 0))
	{
		return model::dLdMu(mu, theta, locInd);
	}
	else if ((locCon == 1) && (var == 0))
	{
		return model::dMdMu(mu, theta, locInd);
	}
	else if ((locCon == 2) && (var == 0))
	{
		return model::dRdMu(mu, theta, locInd);
	}
	else if ((locCon == 2) && (var == 2))
	{
		return model::dRdTheta(mu, theta, locInd);
	}
	else if ((locCon == 3) && (var == 0))
	{
		return model::dSdMu(mu, theta, locInd);
	}
	else if ((locCon == 3) && (var == 2))
	{
		return model::dSdTheta(mu, theta, locInd);
	}
	else if ((locCon == 4) && (var == 2))
	{
		return model::dTdTheta(mu, theta, locInd);
	}
	else
	{
		printDebug("ERROR!");
		exit(0);
	}
	return -1;
}

void
exoMod1D_Nlp::writeParam(const std::vector<double> &par, const std::string fname)
{
	std::ofstream file;
	file.open(fname);
	for (auto i=0; i<par.size(); i++)
	{
		file << par[i] << "\n";
	}
	file.close();
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