#ifndef MultilevelSAEncoder_H_
#define MultilevelSAEncoder_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class MultilevelSAEncoder: public FunctionUnit {
public:
	MultilevelSAEncoder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~MultilevelSAEncoder() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numLevel, int _numEncoder);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double numRead);
	void CalculatePower(double numRead);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double capNandInput, capNandOutput, capNandLgInput, capNandLgOutput, capInvInput, capInvOutput;
	double widthInvN, widthInvP, widthNandN, widthNandP;
	double rampInput, rampOutput;
	int numEncoder;     // number of encoder needed
	int numLevel;       // number of levels from MultilevelSA
	int numInput;       // number of NAND gate in encoder
	int numGate;        // number of NAND gate in encoder 

};

#endif /* MultilevelSAEncoder_H_ */