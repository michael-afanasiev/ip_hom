/*

	This guy tries to solve the inverse homogenization problem.

	Variable ordering: mu, lambda, theta, where the loop is first over
	co-ordinates.

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
	m = mNvar * mNcon;
	
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
	assert(m == mNvar * mNcon);


	for (Index i=0; i<mNtyp; i++)
	{
		for (Index j=0; j<mNloc; j++)
		{
			error(x_l[j+i*mNtyp] = lowerBound(i), "lowerBound");
			error(x_u[j+i*mNtyp] = upperBound(i), "upperBound");
		}
	}

	return true;
}

bool
exoMod1D_Nlp::get_starting_point(Index n, bool init_x, Number* x,
			  		   		     bool init_z, Number* z_L, Number* z_U,
			        		   	 Index m, bool init_lambda, Number* lambda)
{return true;}

bool
exoMod1D_Nlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{return true;}

bool
exoMod1D_Nlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{return true;}

bool
exoMod1D_Nlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{return true;}

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

void
exoMod1D_Nlp::error(const int &e, const std::string &func)
{
	if (e < 0)
	{
		std::cout << "Failure in: " << func << std::endl;
		exit(EXIT_FAILURE);
	}
}