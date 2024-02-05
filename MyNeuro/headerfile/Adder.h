#ifndef ADDER_H_
#define ADDER_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class Adder: public FunctionUnit {
public:
	Adder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~Adder() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(int _numBit, int _numAdder);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double _capLoad, double numRead);
	void CalculatePower(double numRead, int numAdderPerOperation);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double capLoad;
	double capNandInput, capNandOutput;
	int numBit;
	int numAdder;
	double widthNandN, widthNandP;
	double hAdder, wAdder;
	double rampInput, rampOutput;

};

#endif /* ADDER_H_ */
