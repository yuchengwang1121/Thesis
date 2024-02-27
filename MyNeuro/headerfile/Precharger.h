#ifndef PRECHARGER_H_
#define PRECHARGER_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class Precharger: public FunctionUnit {
public:
	Precharger(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~Precharger() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numCol, double _resLoad, double _activityColWrite, int _numReadCellPerOperationNeuro, int _numWriteCellPerOperationNeuro);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double _capLoad, double numRead, double numWrite);
	void CalculatePower(double numRead, double numWrite);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double capLoad, resLoad;
	double capOutputBitlinePrecharger;
	double capWireLoadPerColumn, resWireLoadPerColumn;
	double enableLatency;
	int numCol;			/* Number of columns */
	double widthPMOSBitlinePrecharger, widthPMOSBitlineEqual;
	double capLoadPerColumn;
	double rampInput, rampOutput;
	double activityColWrite;
	int numReadCellPerOperationNeuro;
	int numWriteCellPerOperationNeuro;
};

#endif /* PRECHARGER_H_ */
