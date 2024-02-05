#include <cmath>
#include <iostream>
#include "../headerfile/constant.h"
#include "../headerfile/formula.h"
#include "../headerfile/DFF.h"

using namespace std;

DFF::DFF(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void DFF::Initialize(int _numDff, double _clkFreq){
	if (initialized)
		cout << "[DFF] Warning: Already initialized!" << endl;
	
	numDff = _numDff;
	clkFreq = _clkFreq;
	
	widthTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void DFF::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[DFF] Error: Require initialization first!" << endl;
	} else {
		double hDffInv, wDffInv, hDff, wDff;
		area = 0;
		height = 0;
		width = 0;
		// Assume DFF size is 12 minimum-size standard cells put together
		CalculateGateArea(INV, 1, MIN_NMOS_SIZE*tech.featureSize, tech.pnSizeRatio*MIN_NMOS_SIZE*tech.featureSize, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hDffInv, &wDffInv);
		hDff = hDffInv;
		wDff = wDffInv * 12;
		
		if (_newHeight && _option==NONE) {	// DFF in multiple columns given the total height
			// Calculate the number of DFF per column
			int numDFFPerCol = (int)(_newHeight/hDff);
			if (numDFFPerCol > numDff) {
				numDFFPerCol = numDff;
			} 
			if (numDFFPerCol == 0) {
				cout << "ERROR: Height of subArray is even smaller than a single DFF !!! " << endl;
			}
			int numColDFF = (int)ceil((double)numDff / numDFFPerCol);
			height = _newHeight;
			width = wDff * numColDFF;

		} else if (_newWidth && _option==NONE) {	// DFF in multiple rows given the total width
			// Calculate the number of DFF per row
			int numDFFPerRow = (int)(_newWidth/wDff);
			if (numDFFPerRow > numDff) {
				numDFFPerRow = numDff;
			} 
			if (numDFFPerRow == 0) {
				cout << "ERROR: Width of subArray is even smaller than a single DFF !!! " << endl;
			}
			int numRowDFF = (int)ceil((double)numDff / numDFFPerRow);
			width = _newWidth;
			height = hDff * numRowDFF;
			
		} else {	// Assume one row of DFF by default
			width = wDff * numDff;
			height = hDff;
		}
		area = height * width;
		
		// Modify layout
		newHeight = _newHeight;
		newWidth = _newWidth;
		switch (_option) {
			case MAGIC:
				MagicLayout();
				break;
			case OVERRIDE:
				OverrideLayout();
				break;
			default:    // NONE
				break;
		}
		
		// Capacitance
		// INV
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hDffInv, tech, &capInvInput, &capInvOutput);
		// TG
		capTgGateN = CalculateGateCap(widthTgN, tech);
		capTgGateP = CalculateGateCap(widthTgP, tech);
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hDffInv, tech, NULL, &capTgDrain);
	}
}

void DFF::CalculateLatency(double _rampInput, double numRead){
	if (!initialized) {
		cout << "[DFF] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		rampInput = _rampInput;
		
		readLatency += 1/clkFreq/2;
		readLatency *= numRead;
	}
}

void DFF::CalculatePower(double numRead, double numDffPerOperation) {
	if (!initialized) {
		cout << "[DFF] Error: Require initialization first!" << endl;
	} else {
		readDynamicEnergy = 0;
		/* Leakage power */
		leakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 8 * numDff;
		
		// Assume input D=1 and the energy of CLK INV and CLK TG are for 1 clock cycles
		// CLK INV (all DFFs have energy consumption)
		readDynamicEnergy += (capInvInput + capInvOutput) * tech.vdd * tech.vdd * 4 * numDff;
		// CLK TG (all DFFs have energy consumption)
		readDynamicEnergy += capTgGateN * tech.vdd * tech.vdd * 2 * numDff;
		readDynamicEnergy += capTgGateP * tech.vdd * tech.vdd * 2 * numDff;
		// D to Q path (only selected DFFs have energy consumption)
		readDynamicEnergy += (capTgDrain * 3 + capInvInput) * tech.vdd * tech.vdd * MIN(numDffPerOperation, numDff);	// D input side
		readDynamicEnergy += (capTgDrain  + capInvOutput) * tech.vdd * tech.vdd * MIN(numDffPerOperation, numDff);	// D feedback side
		readDynamicEnergy += (capInvInput + capInvOutput) * tech.vdd * tech.vdd * MIN(numDffPerOperation, numDff);	// Q output side

		readDynamicEnergy *= numRead;
	}
}

void DFF::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void DFF::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}


