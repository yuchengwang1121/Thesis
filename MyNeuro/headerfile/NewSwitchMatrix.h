#ifndef NEWSWITCHMATRIX_H_
#define NEWSWITCHMATRIX_H_

#include "FunctionUnit.h"
#include "constant.h"
#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "DFF.h"

class NewSwitchMatrix: public FunctionUnit {
public:
	NewSwitchMatrix(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~NewSwitchMatrix();
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(int _numOutput, double _activityRowRead, double _clkFreq);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite);
	void CalculatePower(double numRead, double numWrite, double activityRowRead);
	//Mux & operator=(const Mux &);

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
	int numRowTg, numColTg;

	bool neuro;
	bool parallelWrite;
	double activityRowRead;
	double activityColWrite;
	int numWriteCellPerOperationNeuro;
	double numWritePulse;
	double clkFreq;
	DFF dff;

};

#endif /* NEWSWITCHMATRIX_H_ */
