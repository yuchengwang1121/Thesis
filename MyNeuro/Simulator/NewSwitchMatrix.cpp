//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// This NewSwitchMatrix is used for new BNN Parallel RRAM mode, only for WL and BL... Row connected //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <iostream>
#include "../headerfile/constant.h"
#include "../headerfile/formula.h"
#include "../headerfile/NewSwitchMatrix.h"

using namespace std;

NewSwitchMatrix::NewSwitchMatrix(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), dff(_inputParameter, _tech, _cell), FunctionUnit() {
	// TODO Auto-generated constructor stub
	initialized = false;
}

NewSwitchMatrix::~NewSwitchMatrix() {
	// TODO Auto-generated destructor stub
}

void NewSwitchMatrix::Initialize(int _numOutput, double _activityRowRead, double _clkFreq){
	if (initialized)
		cout << "[NewSwitchMatrix] Warning: Already initialized!" << endl;
	
	numOutput = _numOutput;
	activityRowRead = _activityRowRead;
	clkFreq = _clkFreq;
    
	// DFF
	dff.Initialize(numOutput, clkFreq); 
	widthTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	resTg = CalculateOnResistance(widthTgN, NMOS, inputParameter.temperature, tech) * LINEAR_REGION_RATIO;
	
	initialized = true;
}

void NewSwitchMatrix::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[NewSwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		double hTg, wTg;
		area = 0;
		height = 0;
		width = 0;
		double minCellHeight = MAX_TRANSISTOR_HEIGHT * tech.featureSize;   // min standard cell height for 1 Tg 
		if (_newHeight && _option==NONE) {
			if (_newHeight < minCellHeight) {
				cout << "[NewSwitchMatrix] Error: pass gate height is even larger than the array height" << endl;
			}
			int numTgPerCol = (int)(_newHeight / minCellHeight);	// Get max # Tg per column (this is not the final # Tg per column because the last column may have less # Tg)
			numColTg = (int)ceil((double)numOutput / numTgPerCol);	// Get min # columns based on this max # Tg per column
			numTgPerCol = (int)ceil((double)numOutput / numColTg);		// Get # Tg per column based on this min # columns
			TgHeight = _newHeight / numTgPerCol;        // release TG height
			CalculateGateArea(INV, 1, widthTgN, widthTgP, TgHeight, tech, &hTg, &wTg);         // calculate released TG layout height and width
			// DFF
			dff.CalculateArea(_newHeight, NULL, NONE);	
			height = _newHeight;
			width = (wTg * 4) * numColTg + dff.width;      // add switch matrix and dff together
		} else {       // MAGIC or OVERRIDE ...
			CalculateGateArea(INV, 1, widthTgN, widthTgP, minCellHeight, tech, &hTg, &wTg); // Pass gate with folding
			height = hTg * numOutput;
			dff.CalculateArea(height, NULL, NONE);	// Need to give the height information, otherwise by default the area calculation of DFF is in column mode
			width = (wTg * 4) + dff.width;
		}
	    area = height * width;

	    // Modify layout
	    newHeight = _newHeight;
	    newWidth = _newWidth;
	    switch (_option) {
		    case MAGIC:
			    MagicLayout();       // if MAGIC, call Magiclayout() in FunctionUnit.cpp
			    break;
		    case OVERRIDE:
			    OverrideLayout();    // if OVERRIDE, call Overridelayout() in FunctionUnit.cpp
			    break;
		    default:    // NONE
			    break;
		}

	// Capacitance
	// TG
	capTgGateN = CalculateGateCap(widthTgN, tech);
	capTgGateP = CalculateGateCap(widthTgP, tech);
	CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hTg, tech, NULL, &capTgDrain);
	}
}


void NewSwitchMatrix::CalculateLatency(double _rampInput, double _capLoad, double _resLoad, double numRead, double numWrite) {	// For simplicity, assume shift register is ideal
	if (!initialized) {
		cout << "[NewSwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		rampInput = _rampInput;
		capLoad = _capLoad;
		resLoad = _resLoad;
		double capOutput;
		double tr;  /* time constant */
		readLatency = 0;

		// DFF
		dff.CalculateLatency(1e20, numRead);

		// TG
		capOutput = capTgDrain * 5;         // pass 2 TG, 5 loading drain capacitance
		tr = resTg * (capOutput + capLoad) + resLoad * capLoad / 2;     // elmore delay model
		readLatency += horowitz(tr, 0, rampInput, &rampOutput);	// get from chargeLatency in the original SubArray.cpp
		
		readLatency *= numRead;
		readLatency += dff.readLatency;

		writeLatency = horowitz(tr, 0, rampInput, &rampOutput);     // write latency determined by write pulse width
		writeLatency *= numWrite;
		if (numWrite != 0)
		    writeLatency += dff.readLatency;	// Use DFF read latency here because no write in the DFF module
	}
}

void NewSwitchMatrix::CalculatePower(double numRead, double numWrite, double activityRowRead) {      
	if (!initialized) {
		cout << "[NewSwitchMatrix] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		
		// DFF
		dff.CalculatePower(numRead, numOutput);	// Use numOutput since every DFF will pass signal (either 0 or 1)

		// Leakage power
		leakage += dff.leakage;	// Only DFF has leakage, assuming TG do not have leakage

		// Read dynamic energy
		readDynamicEnergy += (capTgDrain * 2) * cell.accessVoltage * cell.accessVoltage * numOutput * activityRowRead;   // 1 TG pass Vaccess to CMOS gate to select the row
		readDynamicEnergy += (capTgDrain * 5) * cell.readVoltage * cell.readVoltage * numOutput * activityRowRead;    // 2 TG pass Vread to BL, total loading is 5 Tg Drain capacitance
		readDynamicEnergy += (capTgGateN + capTgGateP) * 3 * tech.vdd * tech.vdd * numOutput * activityRowRead;    // open 3 TG when selected

		readDynamicEnergy *= numRead;
		readDynamicEnergy += dff.readDynamicEnergy;
		
		// Write dynamic energy (2-step write and average case half SET and half RESET)
		// 1T1R
		// connect to rows, when writing, pass GND to BL, no transmission energy acrossing BL
		writeDynamicEnergy += (capTgDrain * 2) * cell.accessVoltage * cell.accessVoltage;    // 1 TG pass Vaccess to CMOS gate to select the row
		writeDynamicEnergy += (capTgGateN + capTgGateP) * 2 * tech.vdd * tech.vdd * 2;    // open 2 TG when Q selected, and *2 means switching from one selected row to another
        writeDynamicEnergy += (capTgGateN + capTgGateP) * tech.vdd * tech.vdd;    // always open one TG when writing	
		writeDynamicEnergy *= numWrite;
		if (numWrite != 0)
			writeDynamicEnergy += dff.readDynamicEnergy;	// Use DFF read energy here because no write in the DFF module
	}
}


void NewSwitchMatrix::PrintProperty(const char* str) {
	//cout << "NewSwitchMatrix Properties:" << endl;
	FunctionUnit::PrintProperty(str);
}

void NewSwitchMatrix::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}

