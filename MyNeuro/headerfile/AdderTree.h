#ifndef ADDERTREE_H_
#define ADDERTREE_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"
#include "Adder.h"

class AdderTree: public FunctionUnit {
public:
	AdderTree(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~AdderTree() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numSubcoreRow, int _numAdderBit, int _numAdderTree);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double numRead, int numUnitAdd, double _capLoad);
	void CalculatePower(double numRead, int numUnitAdd);

	/* Properties */
	bool initialized;	/* Initialization flag */
    int numSubcoreRow;                    // # of row of subcore in the synaptic core
	int numStage;
	int numTotalAdder;
	int numAdderBit;                      // # of input bits of the Adder
	int numAdderTree;                     // # of Adder Tree
	int numReadPulse;

	Adder adder;
};

#endif /* ADDERTREE_H_ */