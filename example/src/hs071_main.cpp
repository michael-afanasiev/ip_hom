

#include "IpIpoptApplication.hpp"
#include "hs071_nlp.hpp"

#include <iostream>

int
main(int argv, char* argc[])
{

	SmartPtr<TNLP> mynlp = new HS071_NLP();
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	
	app->RethrowNonIpoptException(true);
	app->Options()->SetNumericValue("tol", 1e-15);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	// app->Options()->SetStringValue("derivative_test", "first-order");
	// app->Options()->SetNumericValue("derivative_test_perturbation", 1 * 10e-15);

	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded)
	{
		std::cout << std::endl << std::endl 
		<< "*** Error during initialization!" << std::endl;
		return (int) status;
	}

	status = app->OptimizeTNLP(mynlp);

}