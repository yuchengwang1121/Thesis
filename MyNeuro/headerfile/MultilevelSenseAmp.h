#ifndef MULTILEVELSENSEAMP_H_
#define MULTILEVELSENSEAMP_H_

#include <vector>
#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"
#include "CurrentSenseAmp.h"

using namespace std;

class MultilevelSenseAmp: public FunctionUnit {
public:
	MultilevelSenseAmp(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~MultilevelSenseAmp() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numCol, int _levelOutput, double _clkFreq, int _numReadCellPerOperationNeuro, bool _parallel, bool _currentMode);
	void CalculateArea(double heightArray, double widthArray, AreaModify _option);
	void CalculateLatency(const vector<double> &columnResistance, double numColMuxed, double numRead);
	void CalculatePower(const vector<double> &columnResistance, double numRead);
	double GetColumnLatency(double columnRes);
	double GetColumnPower(double columnRes);

	/* Properties */
	bool initialized;		/* Initialization flag */
	int numCol;				/* Number of columns */
	
    int levelOutput;
	double widthNmos, widthPmos;
	
	bool FPGA;
	bool parallel;
	bool neuro;
	bool currentMode;
	double clkFreq;
	int numReadCellPerOperationNeuro;
	vector<double> Rref;

	CurrentSenseAmp currentSenseAmp;
};

#endif /* MULTILEVELSENSEAMP_H_ */


