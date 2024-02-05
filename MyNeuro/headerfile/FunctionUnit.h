#ifndef FUNCTIONUNIT_H_
#define FUNCTIONUNIT_H_

class FunctionUnit {
public:
	FunctionUnit();
	virtual ~FunctionUnit() {}

	/* Functions */
	virtual void PrintProperty(const char* str);
	virtual void SaveOutput(const char* str);
	virtual void MagicLayout();
	virtual void OverrideLayout();

	/* Properties */
	double height;		/* Unit: m */
	double width;		/* Unit: m */
	double area;		/* Unit: m^2 */
	double emptyArea;		/* Unit: m^2 */
	double usedArea;		/* Unit: m^2 */
	double totalArea;		/* Unit: m^2 */
	double readLatency, writeLatency;		/* Unit: s */
	double readDynamicEnergy, writeDynamicEnergy;	/* Unit: J */
	double leakage;		/* Unit: W */
	double newWidth, newHeight;
	double readPower, writePower;
};

#endif /* FUNCTIONUNIT_H_ */
