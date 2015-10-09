
#include "includes.hpp"

int main(int argc, char const *argv[])
{
	
	std::cout << "Hello world." << std::endl;

	exodusFile exoFile;
	exoFile.initFromFile("/Users/michaelafanasiev/Desktop/effective_media/model_iteration_00.ex2");
	exoFile.readVariables();
	exoFile.readCoordinates();

	SmartPtr<TNLP> mynlp = new exoMod1D_Nlp();
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();


	return 0;
}