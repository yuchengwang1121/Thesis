#ifndef BUS_H_
#define BUS_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class Bus: public FunctionUnit {
public:
	Bus(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~Bus() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(BusMode _mode, int _numRow, int _numCol, double _delaytolerance, double _busWidth, double _unitHeight, double _unitWidth);
	void CalculateArea(double foldedratio, bool overLap);
	void CalculateLatency(double numRead);
	void CalculatePower(double numBitAccess, double numRead);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double widthInvN, widthInvP, wInv, hInv, capInvInput, capInvOutput;
	double widthMinInvN, widthMinInvP, wMinInv, hMinInv, capMinInvInput, capMinInvOutput, wRep, hRep, capRepInput, capRepOutput;
	double AR, Rho, unitLengthWireResistance, minDist, minDelay, resOnRep;
	int numRow, numCol, numRepeater, repeaterSize;
	double unitHeight, unitWidth, wireWidth;
	double busWidth, delaytolerance, unitLengthWireCap, wireLength;
	double unitLatencyRep, unitLatencyWire, unitLengthLeakage, unitLengthEnergyRep, unitLengthEnergyWire;
	BusMode mode;
};

#endif /* BUS_H_ */
