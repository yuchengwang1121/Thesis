#include <cmath>
#include <iostream>
#include "../headerfile/constant.h"
#include "../headerfile/formula.h"
#include "../headerfile/MultilevelSAEncoder.h"

using namespace std;

MultilevelSAEncoder::MultilevelSAEncoder(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): inputParameter(_inputParameter), tech(_tech), cell(_cell), FunctionUnit() {
	initialized = false;
}

void MultilevelSAEncoder::Initialize(int _numLevel, int _numEncoder){
	if (initialized)
		cout << "[MultilevelSAEncoder] Warning: Already initialized!" << endl;
	
	numEncoder = _numEncoder;      // number of encoder needed
	numLevel= _numLevel;           // number of levels from MultilevelSA
	numInput = ceil(numLevel/2);       // number of NAND gate in encoder
	numGate = ceil(log2(numLevel));      // number of NAND gate in encoder 
	
	widthInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	widthNandN = 2 * MIN_NMOS_SIZE * tech.featureSize;
	widthNandP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;

	initialized = true;
}

void MultilevelSAEncoder::CalculateArea(double _newHeight, double _newWidth, AreaModify _option) {
	if (!initialized) {
		cout << "[MultilevelSAEncoder] Error: Require initialization first!" << endl;
	} else {
        	double wEncoder, hEncoder, wNand, hNand, wNandLg, hNandLg, wInv, hInv;
		area = 0;
		height = 0;
		width = 0;
		// NAND2
		CalculateGateArea(NAND, 2, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		// INV
		CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		// Large NAND in Encoder
		CalculateGateArea(NAND, numInput, widthNandN, widthNandP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hNandLg, &wNandLg);
		
		wEncoder = 2*wInv + wNand + wNandLg;
		hEncoder = max( (numLevel-1)*hInv, (numLevel-1)*hNand );
	    
		if (_newWidth && _option==NONE) {
			int numEncoderPerRow = (int)ceil(_newWidth/wEncoder);
			if (numEncoderPerRow > numEncoder) {
				numEncoderPerRow = numEncoder;
			}
			int numRowEncoder = (int)ceil((double)numEncoder / numEncoderPerRow);
			width = MAX(_newWidth, wEncoder);
			height = hEncoder * numRowEncoder;
		} else if (_newHeight && _option==NONE) {
			int numEncoderPerColumn = (int) ceil(_newHeight/hEncoder);
			if (numEncoderPerColumn > numEncoder) {
				numEncoderPerColumn = numEncoder;
			}
			int numColEncoder = (int)ceil((double)numEncoder / numEncoderPerColumn);
			height = MAX(_newHeight, hEncoder);
			width = wEncoder*numColEncoder;
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
		CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);
		// NAND2
		CalculateGateCapacitance(NAND, 2, widthNandN, widthNandP, hNand, tech, &capNandInput, &capNandOutput);
		// Large NAND in Encoder
		CalculateGateCapacitance(NAND, numInput, widthNandN, widthNandP, hNandLg, tech, &capNandLgInput, &capNandLgOutput);
	}
}

void MultilevelSAEncoder::CalculateLatency(double _rampInput, double numRead){
	if (!initialized) {
		cout << "[MultilevelSAEncoder] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		rampInput = _rampInput;
		double tr;		/* time constant */
		double gm;		/* transconductance */
		double beta;	/* for horowitz calculation */
		double resPullUp, resPullDown;
		double readLatencyIntermediate = 0;
		double ramp[10];
		
		ramp[0] = rampInput;

		// 1st INV to NAND2
		resPullDown = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capInvOutput + capNandInput * 2);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatency += horowitz(tr, beta, ramp[0], &ramp[1]);
		
		// 2nd NAND2 to Large NAND
		resPullUp = CalculateOnResistance(widthNandP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * (capNandOutput + capNandLgInput * numInput);
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatency += horowitz(tr, beta, ramp[1], &ramp[2]);
		
		// 3rd large NAND to INV
		resPullDown = CalculateOnResistance(widthNandN, NMOS, inputParameter.temperature, tech) * 2;
		tr = resPullDown * (capNandLgOutput + capInvInput);
		gm = CalculateTransconductance(widthNandN, NMOS, tech);
		beta = 1 / (resPullDown * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[2], &ramp[3]);

		// 4th INV
		resPullUp = CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		tr = resPullUp * capInvOutput;
		gm = CalculateTransconductance(widthNandP, PMOS, tech);
		beta = 1 / (resPullUp * gm);
		readLatencyIntermediate += horowitz(tr, beta, ramp[3], &ramp[4]);
		
		readLatency *= numRead;
		rampOutput = ramp[4];
	}
}

void MultilevelSAEncoder::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[MultilevelSAEncoder] Error: Require initialization first!" << endl;
	} else {
		readDynamicEnergy = 0;
		leakage = 0;

		leakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * (numLevel+numGate) * numEncoder
		          + CalculateGateLeakage(NAND, 2, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * (numLevel+numGate) * numEncoder
				  + CalculateGateLeakage(NAND, numInput, widthNandN, widthNandP, inputParameter.temperature, tech) * tech.vdd * numGate * numEncoder;
		
		readDynamicEnergy += (capInvInput + capInvOutput) * tech.vdd * tech.vdd * (numLevel+numGate) * numEncoder;
		readDynamicEnergy += (capNandInput + capNandOutput) * tech.vdd * tech.vdd * (numLevel+numGate) * numEncoder;
		readDynamicEnergy += (capNandLgInput + capNandLgOutput) * tech.vdd * tech.vdd * numGate * numEncoder;
		readDynamicEnergy *= numRead;
		
		if (!readLatency) {
			//cout << "[MultilevelSenseAmp] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}

void MultilevelSAEncoder::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}


