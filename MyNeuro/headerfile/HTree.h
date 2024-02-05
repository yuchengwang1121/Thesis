#ifndef HTREE_H_
#define HTREE_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class HTree: public FunctionUnit {
public:
	HTree(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~HTree() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(int _numRow, int _numCol, double _delaytolerance, double _busWidth);
	void CalculateArea(double unitHeight, double unitWidth, double foldedratio);
	void CalculateLatency(int x_init, int y_init, int x_end, int y_end, double unitHeight, double unitWidth, double numRead);
	void CalculatePower(int x_init, int y_init, int x_end, int y_end, double unitHeight, double unitWidth, double numBitAccess, double numRead);
	double GetUnitLengthRes(double wireLength);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double widthInvN, widthInvP, wInv, hInv, capInvInput, capInvOutput;
	double widthMinInvN, widthMinInvP, wMinInv, hMinInv, capMinInvInput, capMinInvOutput, wRep, hRep, capRepInput, capRepOutput;
	double numStage, numTree, AR, Rho, unitLengthWireResistance, minDist, minDelay, resOnRep;
	int numRow, numCol, numRepeater, numTotalRepeater, repeaterSize;
	double unitHeight, unitWidth;
	double numRep_vertical, numRep_horizontal;
	double busWidth, delaytolerance, unitLengthWireCap, totalWireLength;
	double unitLatencyRep, unitLatencyWire, unitLengthLeakage, unitLengthEnergyRep, unitLengthEnergyWire;
	double find_stage;
	int x_center, y_center, hit, skipVer;

};

#endif /* HTREE_H_ */
