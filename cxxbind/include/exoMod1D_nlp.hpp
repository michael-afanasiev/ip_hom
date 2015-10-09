#ifndef __homOpt_exoMod1D_nlp_hpp__
#define __homOpt_exoMod1D_nlp_hpp__

#include "IpTNLP.hpp"

#include <limits>

class exodusFile;

using namespace Ipopt;
class exoMod1D_Nlp: public TNLP
{

public:

	exoMod1D_Nlp();

	virtual
	~exoMod1D_Nlp();

	virtual bool
	get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
				 Index &nnz_h_lag, IndexStyleEnum &index_style);

	virtual bool
	get_bounds_info(Index n, Number *x_l, Number *x_u,
					Index m, Number *g_l, Number *g_u);

	virtual bool 
  	get_starting_point(Index n, bool init_x, Number* x,
					   bool init_z, Number* z_L, Number* z_U,
					   Index m, bool init_lambda, Number* lambda);
	virtual bool 
	eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	virtual bool 
	eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	virtual bool 
	eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

	virtual bool 
	eval_jac_g(Index n, const Number* x, bool new_x,
			   Index m, Index nele_jac, Index* iRow, Index *jCol,
			   Number* values);


	virtual bool 
	eval_h(Index n, const Number* x, bool new_x,
		   Number obj_factor, Index m, const Number* lambda,
		   bool new_lambda, Index nele_hess, Index* iRow,
		   Index* jCol, Number* values);

	virtual void
	finalize_solution(SolverReturn status,
					  Index n, const Number* x, const Number *z_L, 
					  const Number *z_U, Index m, const Number *g,
					  const Number *lambda, Number obj_value,
					  const IpoptData *ip_data,
					  IpoptCalculatedQuantities *ip_cq);

	void
	initFromExodus(exodusFile &exo);

private:

	static constexpr Number zero = 0;
	static constexpr Number infty = std::numeric_limits<double>::max();

	exoMod1D_Nlp(const exoMod1D_Nlp&);
	exoMod1D_Nlp& operator=(const exoMod1D_Nlp&);

	int mNvar;					// number of variables
	int mNtyp;					// number of var types
	int mNcon;					// number of constraints
	int mNloc;					// number of unique locations
	int mWinLength;				// window length of convolution
	int mNjacNonZero;			// number of non-zero in jacobian

	std::vector<double> mX;
	std::vector<double> mY;

	std::vector<double> mA;
	std::vector<double> mC;
	std::vector<double> mL;
	std::vector<double> mF;
	std::vector<double> mMu;
	std::vector<double> mVph;
	std::vector<double> mVpv;
	std::vector<double> mVsv;
	std::vector<double> mRho;
	std::vector<double> mTheta;
	std::vector<double> mLambda;

	std::vector<double> mBigL;	
	std::vector<double> mBigM;
	std::vector<double> mBigR;
	std::vector<double> mBigS;
	std::vector<double> mBigT;

	void
	error(const int &e, const std::string &func);

	Number
	lowerBound(const int &i);

	Number
	upperBound(const int &i);

};

#endif // __homOpt_exoMod1D_nlp_hpp__