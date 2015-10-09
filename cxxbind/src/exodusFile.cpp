

#include "includes.hpp"

exodusFile::exodusFile()
{mInit = false;}

exodusFile::~exodusFile()
{ex_close(mIdExo);}

void
exodusFile::initFromFile(std::string fname)
{
	mCompWs = 8;
	mIoWs = 0;
	mFname = fname;
	openError(mIdExo = ex_open(mFname.c_str(), EX_READ, &mCompWs, &mIoWs, 
						   	   &mVers));
	error(ex_get_init(mIdExo, mTitle, &mNdim, &mNnode, &mNelem, &mNelemBlk, 
			          &mNsideSet, &mNnodeSet), "ex_get_init");
	error(ex_get_var_param(mIdExo, "e", &mNumElemVars), "ex_get_var_param");
	mInit = true;
}

void
exodusFile::readCoordinates()
{
	mX.resize(mNnode);
	mY.resize(mNnode);
	ex_get_coord(mIdExo, mX.data(), mY.data(), NULL);
}

void
exodusFile::readVariables()
{
	char *nm[mNumElemVars];
	for (size_t i=0; i<mNumElemVars; i++) 
	{
		nm[i] = (char *) calloc((MAX_STR_LENGTH+1), sizeof(char));
	}
	error(ex_get_var_names(mIdExo, "E", mNumElemVars, nm), "ex_get_var_names");
	for (size_t i=0; i<mNumElemVars; i++) 
	{
		mVarNames.push_back(std::string(nm[i]));
		free(nm[i]);
		std::vector<double> val(mNelem);
		error(ex_get_elem_var(mIdExo, 1, i+1, 1, mNelem, val.data()),
			  "ex_get_elem_var");
		mVarValues.push_back(val);
	}
}

void
exodusFile::error(int e, std::string function)
{
	if (e)
	{
		std::cout << "Error in " << function << std::endl;
		exit(EXIT_FAILURE);
	}
}

void 
exodusFile::openError(int e)
{
	if (e < 0)
	{
		std::cout << "Error opening EXODUS file." << std::endl;
		exit(EXIT_FAILURE);
	}
}

std::vector<double> 
exodusFile::getVar(std::string name)
{return mVarValues[getIndForVal(name, mVarNames)];}