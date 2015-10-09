
#include "includes.hpp"

int main(int argc, char const *argv[])
{

	exodusFile exoFile;
	exoFile.initFromFile("/Users/michaelafanasiev/Desktop/effective_media/layered_iteration_00.ex2");
	exoFile.readVariables();
	exoFile.readCoordinates();
	exoFile.readConnectivity();

	SmartPtr<exoMod1D_Nlp> mynlp = new exoMod1D_Nlp();
	mynlp->initFromExodus(exoFile);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);
	app->Options()->SetNumericValue("tol", 1e-7);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

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