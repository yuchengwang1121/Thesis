#ifndef CURRENTSENSEAMP_H_
#define CURRENTSENSEAMP_H_

#include <vector>
#include "FunctionUnit.h"
#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"

using namespace std;

class CurrentSenseAmp: public FunctionUnit {
public:
	CurrentSenseAmp(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~CurrentSenseAmp();
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numCol, bool _parallel, bool _rowbyrow, double _clkFreq, int _numReadCellPerOperationNeuro);
	void CalculateArea(double _widthCurrentSenseAmp);
	void CalculateLatency(const vector<double> &columnResistance, double numColMuxed, double numRead);
	void CalculatePower(const vector<double> &columnResistance, double numRead);
	void CalculateUnitArea();
	double GetColumnLatency(double columnRes);
	double GetColumnPower(double columnRes);


	/* Properties */
	bool initialized;	/* Initialization flag */
	bool invalid;		/* Indicate that the current configuration is not valid */
	int numCol;		/* Number of columns */
	double widthNmos, widthPmos;
	double hNmosL, wNmosL, hNmosS, wNmosS, hNmosM, wNmosM;
	double areaUnit;
	double widthArray;
	bool parallel;
	bool rowbyrow;
	double clkFreq, Rref;
	int numReadCellPerOperationNeuro;
};

#endif /* CURRENTSENSEAMP_H_ */