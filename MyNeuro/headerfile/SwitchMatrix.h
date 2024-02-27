#ifndef SWITCHMATRIX_H_
#define SWITCHMATRIX_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"
#include "DFF.h"

class SwitchMatrix: public FunctionUnit {
public:
	SwitchMatrix(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~SwitchMatrix() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(int _mode, int _numOutput, double _resTg, bool _neuro, bool _parallelWrite, double _activityRowRead, double _activityColWrite, int _numWriteCellPerOperationMemory, int _numWriteCellPerOperationNeuro, double _numWritePulse, double _clkFreq);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite);
	void CalculatePower(double numRead, double numWrite, double activityRowRead, double activityColWrite);

	/* Properties */
	bool initialized;	/* Initialization flag */
	int numOutput;
	double capLoad;
	double resLoad;
	double widthTgN, widthTgP;
	double TgHeight, TgWidth;
	double capTgDrain, capTgGateN, capTgGateP;
	double resTg;
	double rampInput, rampOutput;
	int mode;
	int numRowTgPair, numColTgPair;

	bool neuro;
	bool parallelWrite;
	double activityRowRead;
	double activityColWrite;
	int numWriteCellPerOperationMemory;
	int numWriteCellPerOperationNeuro;
	double numWritePulse;
	double clkFreq;
	DFF dff;
};

#endif /* SWITCHMATRIX_H_ */
