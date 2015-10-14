
#include "includes.hpp"

int main(int argc, char const *argv[])
{

	if (argc < 2)
	{
		printError("Usage: HomTest [exodus_file]");
	}

	exodusFile exoFile;
	exoFile.initFromFile(argv[1]);
	exoFile.readVariables();
	exoFile.readCoordinates();
	exoFile.readConnectivity();

	SmartPtr<exoMod1D_Nlp> mynlp = new exoMod1D_Nlp();
	mynlp->initFromExodus(exoFile);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);
	app->Options()->SetNumericValue("tol", 1e-2);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	// app->Options()->SetNumericValue("derivative_test_perturbation", 1e-6);
	// app->Options()->SetStringValue("derivative_test", "first-order");
	// app->Options()->SetNumericValue("point_perturbation_radius", 1e-5);
	// app->Options()->SetIntegerValue("derivative_test_first_index", 494);
	// app->Options()->SetNumericValue("derivative_test_tol", 1e-8);
	// app->Options()->SetStringValue("derivative_test_print_all", "yes");

	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded)
	{
		std::cout << std::endl << std::endl
		<< "*** Error during initialization!" << std::endl;
		return (int) status;
	}

	status = app->OptimizeTNLP(mynlp);

	return 0;
}