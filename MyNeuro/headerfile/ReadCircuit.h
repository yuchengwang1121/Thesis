#ifndef READCIRCUIT_H_
#define READCIRCUIT_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class ReadCircuit: public FunctionUnit {
public:
	ReadCircuit(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~ReadCircuit() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(ReadCircuitMode _mode, int _numReadCol, int _maxNumIntBit, SpikingMode _spikingMode, double _clkFreq);
	void CalculateArea(double _newWidth);
	void CalculateLatency(double numRead);
	void CalculatePower(double numof1, double numof2, double numof3, double numof4, double numof5, double numof6, double numof7, double numof8, double numof9, double numof10, double numof20, double numof30, double numof40, double numof50, double numof60, double numof70, double numof80, double numof90, double numof100, double numRead);
	void CalculateUnitArea();

	/* Properties */
	bool initialized;	/* Initialization flag */
	int numReadCol;
	double widthDffTgN, widthDffTgP, widthDffInvN, widthDffInvP, widthDffNorN, widthDffNorP, widthTgN, widthTgP, widthInvN, widthInvP;
	double widthNmos1, widthPmos1, widthNmos2, widthNmos3, widthPmos3, widthNmos4, widthPmos4, widthNmos5, widthPmos5, widthNmos6, widthNmos7, widthNmos8, widthPmos8;
	double rampInput, rampOutput;
	double capInput;
	double capDffTgGateN, capDffTgGateP, capDffTgDrain, capDffInvInput, capDffInvOutput, capNorInput, capNorOutput;
	double capTgGateN, capTgGateP, capTgDrain, capNmosGate, capNmosDrain, capPmosGate, capPmosDrain, capInvInput, capInvOutput;
	double areaUnit, hUnit, wUnit, areaReadBody, hReadBody, wReadBody, areaDff, hDff, wDff;
	double maxNumIntPerCycle;
	double voltageIntThreshold;
	int numDff;
	ReadCircuitMode mode;
	double Vhold, Vth, Vrow, Vcol, R_OSC_OFF;
	int maxNumIntBit;
	SpikingMode spikingMode;
	double clkFreq;
	int numUnitPerRow, numRowUnit;
};

#endif /* READCIRCUIT_H_ */
