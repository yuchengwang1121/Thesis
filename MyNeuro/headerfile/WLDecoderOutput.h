#ifndef WLDecoderOutput_H_
#define WLDecoderOutput_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"

class WLDecoderOutput: public FunctionUnit {
public:
	WLDecoderOutput(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);
	virtual ~WLDecoderOutput() {}
	const InputParameter& inputParameter;
	const Technology& tech;
	const MemCell& cell;

	/* Functions */
	void PrintProperty(const char* str);
	void SaveOutput(const char* str);
	void Initialize(int _numWLRow, bool _multifunctional, bool _neuro);
	void CalculateArea(double _newHeight, double _newWidth, AreaModify _option);
	void CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite);
	void CalculatePower(double numRead, double numWrite);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double capLoad;	/* Output capacitance, unit: F */
	double resLoad;	/* Output resistance, unit: ohm */
	int numWLRow;
	double widthNorN, widthNorP, widthInvN, widthInvP, widthTgN, widthTgP, widthNmos;
	double capNorInput, capNorOutput, capInvInput, capInvOutput, capTgGateN, capTgGateP, capTgDrain, capNmosGate, capNmosDrain;
	double resTg;
	double rampInput, rampOutput;
	bool multifunctional;
	bool neuro;
};

#endif /* WLDecoderOutput_H_ */
