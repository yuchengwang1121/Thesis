#ifndef DecoderDriver_H_
#define DecoderDriver_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class DecoderDriver: public FunctionUnit {
public:
	DecoderDriver(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~DecoderDriver() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(int _mode, int _numOutput, int numLoad);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double capLoad1, double capLoad2, double _resLoad, double numRead, double numWrite);
	void CalculatePower(double numReadCellPerOp, double numWriteCellPerOp, double numRead, double numWrite);
	//DecoderDriver & operator=(const DecoderDriver &);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double capLoad1, capLoad2;	/* Output capacitance, unit: F */
	double resLoad;	/* Output resistance, unit: ohm */
	int numOutput;
	double widthInvN, widthInvP, widthTgN, widthTgP;
	double capInvInput, capInvOutput, capTgDrain, capTgGateN, capTgGateP;
	double resTg;
	double rampInput, rampOutput;
	int mode;
	int numRowTg, numColTg;
	double TgHeight, TgWidth;

};

#endif /* DecoderDriver_H_ */
