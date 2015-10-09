#ifndef __homOpt_exodusFile_hpp__
#define __homOpt_exodusFile_hpp__

#include <string>
#include <vector>
extern "C"
{
	#include "exodusII.h"
}

class exodusFile
{

	bool mInit;

	int mIoWs;
	int	mNdim;
	int	mNnode;
	int	mNelem;
	int	mIdExo;
	int	mCompWs;
	int	mNelemBlk;		
	int mNsideSet;
	int mNnodeSet;
	int mNumElemVars;

	char mTitle[MAX_LINE_LENGTH+1];

	float mVers;

	std::string mFname;
	std::vector<int> mCon;
	std::vector<double> mX;
	std::vector<double> mY;
	std::vector<std::string> mVarNames;
	std::vector<std::vector<double>> mVarValues;

	void
	error(int e, std::string function);

	void
	openError(int e);

public:

	exodusFile();
	~exodusFile();

	// modifiers.
	void
	initFromFile(std::string fname);
	void
	readVariables();
	void
	readCoordinates();
	void
	readConnectivity();

	// accessors.
	std::vector<int> con() {return mCon;};
	std::vector<double> getVar(std::string name);
	std::vector<double> X() {return mX;};
	std::vector<double> Y() {return mY;};

};


#endif // __homOpt_exodusFile_hpp__