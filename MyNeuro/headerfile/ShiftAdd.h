#ifndef SHIFTADD_H_
#define SHIFTADD_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"
#include "Adder.h"
#include "DFF.h"

class ShiftAdd: public FunctionUnit {
public:
	ShiftAdd(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~ShiftAdd() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(int _numUnit, int _numAdderBit, double _clkFreq, SpikingMode _spikingMode, int _numReadPulse);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double numRead);
	void CalculatePower(double numRead);
	void CalculateUnitArea();

	/* Properties */
	bool initialized;	/* Initialization flag */
	int numUnit;
	double layoutWidth;
	double widthInvN, widthInvP, widthNandN, widthNandP;
	int numInv, numNand, numAdder, numDff;
	int numAdderBit, numBitPerDff;
	double rampInput, rampOutput;
	double clkFreq;
	SpikingMode spikingMode;
	int numReadPulse;

	Adder adder;
	DFF dff;
};

#endif /* SHIFTADD_H_ */
