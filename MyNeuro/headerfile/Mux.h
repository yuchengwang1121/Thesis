#ifndef MUX_H_
#define MUX_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class Mux: public FunctionUnit {
public:
	Mux(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~Mux() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(int _numInput, int _numSelection, double _resTg, bool _FPGA);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double _capLoad, double numRead);
	void CalculatePower(double numRead);

	/* Properties */
	bool initialized;	/* Initialization flag */
	int numInput;
	int numSelection;		/* Number of Selections */
	bool FPGA;
	double capLoad;
	double minDriverCurrent;
	double widthNandN, widthNandP, widthEnInvN, widthEnInvP, widthTgN, widthTgP, widthMuxInvN, widthMuxInvP;
	double capNandInput, capNandOutput, capEnInvInput, capEnInvOutput, capMuxInvInput, capMuxInvOutput, capTgGateN, capTgGateP, capTgDrain;
	double resTg;
	double widthTgShared;
	double rampInput, rampOutput;
	int numRowTg;
	double TgWidth;
};

#endif /* MUX_H_ */
