#ifndef SARADC_H_
#define SARADC_H_

#include <vector>
#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

using namespace std;

class SarADC: public FunctionUnit {
public:
	SarADC(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~SarADC() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numCol, int _levelOutput, double _clkFreq, int _numReadCellPerOperationNeuro);
	void CalculateUnitArea();
	void CalculateArea(double heightArray, double widthArray, AreaModify _option);
	void CalculateLatency(double numRead);
	void CalculatePower(const vector<double> &columnResistance, double numRead);
	double GetColumnPower(double columnRes);

	/* Properties */
	bool initialized;		/* Initialization flag */
	int numCol;				/* Number of columns */
	double widthNmos, widthPmos;
    int levelOutput;
	
	bool FPGA;
	bool neuro;
	double clkFreq, areaUnit;
	int numReadCellPerOperationNeuro;
	vector<double> Rref;

};

#endif /* SARADC_H_ */


