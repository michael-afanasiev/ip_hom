#include "hs071_nlp.hpp"

#include <cassert>
#include <iostream>

using namespace Ipopt;

// constructor
HS071_NLP::HS071_NLP()
{}

// destructor
HS071_NLP::~HS071_NLP()
{}

bool
HS071_NLP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
						Index &nnz_h_lag, IndexStyleEnum &index_style)
{
	// Problem has 4 variables (x[0] through x[3]).
	n = 4;

	// one equality constraint and one inequality constraint.
	m = 2;

	// jacobian is dense and contains 8 nonzeros.
	nnz_jac_g = 8;

	// hessian is dense and has 16 nonzeros, but we only need lower left
	// corner (it is symmetric).
	nnz_h_lag = 10;

	// use C style indexing (0-based).
	index_style = TNLP::C_STYLE;

	return true;
}

bool
HS071_NLP::get_bounds_info(Index n, Number *x_l, Number *x_u,
						   Index m, Number *g_l, Number *g_u)
{
	// just assert to make sure values make sense.
	assert(n == 4);
	assert(m == 2);

	// variables have a lower bound of 1.
	for (Index i=0; i<4; i++)
	{
		x_l[i] = 1.0;
	}

	// variables have an upper bound of 5.
	for (Index i=0; i<4; i++)
	{
		x_u[i] = 5.0;
	}

	// first constraint has a lower bound of 25, upper of \infty.
	g_l[0] = 25;
	g_u[0] = infty;

	// second constraint is equality constraint so set u and l to same.
	g_l[1] = g_u[1] = 40.0;

	return true;
}

bool
HS071_NLP::get_starting_point(Index n, bool init_x, Number* x,
					   		  bool init_z, Number* z_L, Number* z_U,
					   	      Index m, bool init_lambda, Number* lambda)
{
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	// initialize to some starting point.
	x[0] = 1.0;
	x[1] = 4.743;
	x[2] = 3.82115;
	x[3] = 1.37941;

	return true;
}

bool
HS071_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	assert(n == 4);

	std::cout << x[1] << std::endl;
	// evaulate the objective function
	obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
	return true;
}

bool
HS071_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	assert(n == 4);
	
	// evaluate grad_{x} f(x)
	grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
	grad_f[1] = x[0] * x[3];
	grad_f[2] = x[0] * x[3] + 1;
	grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
	return true;
}

bool
HS071_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	assert(n == 4);
	assert(m == 2);

	g[0] = x[0] * x[1] * x[2] * x[3];
	g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
	return true;
}

bool
HS071_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
			   		  Index m, Index nele_jac, Index* iRow, Index *jCol,
			   		  Number* values)
{
	if (values == NULL)
	{
		// dense sparse matrix structure.
		iRow[0] = 0;
		jCol[0] = 0;
		iRow[1] = 0;
		jCol[1] = 1;
		iRow[2] = 0;
		jCol[2] = 2;
		iRow[3] = 0;
		jCol[3] = 3;
		iRow[4] = 1;
		jCol[4] = 0;
		iRow[5] = 1;
		jCol[5] = 1;
		iRow[6] = 1;
		jCol[6] = 2;
		iRow[7] = 1;
		jCol[7] = 3;
		std::cout << jCol[0] << std::endl;
	} 
	else 
	{
		// values of jacobian.
		values[0] = x[1]*x[2]*x[3]; // 0,0
		values[1] = x[0]*x[2]*x[3]; // 0,1
		values[2] = x[0]*x[1]*x[3]; // 0,2
		values[3] = x[0]*x[1]*x[2]; // 0,3
		values[4] = 2*x[0]; // 1,0
		values[5] = 2*x[1]; // 1,1
		values[6] = 2*x[2]; // 1,2
		values[7] = 2*x[3]; // 1,3
	}

	return true;
}

bool
HS071_NLP::eval_h(Index n, const Number* x, bool new_x,
		   		  Number obj_factor, Index m, const Number* lambda,
		   		  bool new_lambda, Index nele_hess, Index* iRow,
		   		  Index* jCol, Number* values)
{
	return false;
}

void HS071_NLP::finalize_solution(SolverReturn status,
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
}