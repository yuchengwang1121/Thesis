#include <cmath>
#include <iostream>
#include <vector>
#include "../headerfile/constant.h"
#include "../headerfile/formula.h"
#include "../headerfile/SubArray.h"


using namespace std;

SubArray::SubArray(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell):
						inputParameter(_inputParameter), tech(_tech), cell(_cell),
						wlDecoder(_inputParameter, _tech, _cell),
						wlNewDecoderDriver(_inputParameter, _tech, _cell),
						wlNewSwitchMatrix(_inputParameter, _tech, _cell),
						slSwitchMatrix(_inputParameter, _tech, _cell),
						mux(_inputParameter, _tech, _cell),
						muxDecoder(_inputParameter, _tech, _cell),
						precharger(_inputParameter, _tech, _cell),
						sramWriteDriver(_inputParameter, _tech, _cell),
						adder(_inputParameter, _tech, _cell),
						dff(_inputParameter, _tech, _cell),
						multilevelSenseAmp(_inputParameter, _tech, _cell),
						multilevelSAEncoder(_inputParameter, _tech, _cell),
						sarADC(_inputParameter, _tech, _cell),
						shiftAddInput(_inputParameter, _tech, _cell),
						shiftAddWeight(_inputParameter, _tech, _cell){
	initialized = false;
	readDynamicEnergyArray = writeDynamicEnergyArray = 0;
}

void SubArray::Initialize(int _numRow, int _numCol, double _unitWireRes){  //initialization module 
	
	numRow = _numRow;    //import parameters
	numCol = _numCol;
	unitWireRes = _unitWireRes;
	
	double MIN_CELL_HEIGHT = MAX_TRANSISTOR_HEIGHT;  //set real layout cell height
	double MIN_CELL_WIDTH = (MIN_GAP_BET_GATE_POLY + POLY_WIDTH) * 2;  //set real layout cell width
	// else if (cell.memCellType == Type::RRAM ||  cell.memCellType == Type::FeFET) {  //if array is RRAM
	double cellHeight = cell.heightInFeatureSize; 
	double cellWidth = cell.widthInFeatureSize;  // 1T1R
	if (relaxArrayCellWidth) {
		lengthRow = (double)numCol * MAX(cellWidth, MIN_CELL_WIDTH*2) * tech.featureSize;	// Width*2 because generally switch matrix has 2 pass gates per column, even the SL/BL driver has 2 pass gates per column in traditional 1T1R memory
	} else {
		lengthRow = (double)numCol * cellWidth * tech.featureSize;
	}
	if (relaxArrayCellHeight) {
		lengthCol = (double)numRow * MAX(cellHeight, MIN_CELL_HEIGHT) * tech.featureSize;
	} else {
		lengthCol = (double)numRow * cellHeight * tech.featureSize;
	}
	// }      //finish setting array size
	
	capRow1 = lengthRow * 0.2e-15/1e-6;	// BL for 1T1R, WL for Cross-point and SRAM
	capRow2 = lengthRow * 0.2e-15/1e-6;	// WL for 1T1R
	capCol = lengthCol * 0.2e-15/1e-6;
	
	resRow = lengthRow * unitWireRes; 
	resCol = lengthCol * unitWireRes;
	
	//start to initializing the subarray modules
	// else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
	// if (cell.accessType == CMOS_access) {	// 1T1R
	cell.resCellAccess = cell.resistanceOn * IR_DROP_TOLERANCE;    //calculate access CMOS resistance
	cell.widthAccessCMOS = CalculateOnResistance(tech.featureSize, NMOS, inputParameter.temperature, tech) * LINEAR_REGION_RATIO / cell.resCellAccess;   //get access CMOS width
	if (cell.widthAccessCMOS > cell.widthInFeatureSize) {	// Place transistor vertically
		printf("Transistor width of 1T1R=%.2fF is larger than the assigned cell width=%.2fF in layout\n", cell.widthAccessCMOS, cell.widthInFeatureSize);
		exit(-1);
	}

	cell.resMemCellOn = cell.resCellAccess + cell.resistanceOn;        //calculate single memory cell resistance_ON
	cell.resMemCellOff = cell.resCellAccess + cell.resistanceOff;      //calculate single memory cell resistance_OFF
	cell.resMemCellAvg = cell.resCellAccess + cell.resistanceAvg;      //calculate single memory cell resistance_AVG

	capRow2 += CalculateGateCap(cell.widthAccessCMOS * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap
	capCol += CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech) * numRow;	// If capCol is found to be too large, increase cell.widthInFeatureSize to relax the limit
	// }
		
	if (conventionalSequential) {  
		double capBL = lengthCol * 0.2e-15/1e-6;
		int numAdder = (int)ceil(numCol/numColMuxed);   // numCol is divisible by numCellPerSynapse
		int numInput = numAdder;        //XXX input number of MUX, 
		double resTg = cell.resMemCellOn;     //transmission gate resistance
		int adderBit = (int)ceil(log2(numRow)) + avgWeightBit;  
		
		wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);          
		// if (cell.accessType == CMOS_access) {
		wlNewDecoderDriver.Initialize(numRow);          
		// }
		slSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, numWritePulseAVG, clkFreq);     //SL use switch matrix
		if (numColMuxed>1) {
			mux.Initialize(numInput, numColMuxed, resTg, FPGA);     
			muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
		}
		
		if (SARADC) {
			sarADC.Initialize(numCol/numColMuxed, pow(2, avgWeightBit), clkFreq, numReadCellPerOperationNeuro);
		} else {
			multilevelSenseAmp.Initialize(numCol/numColMuxed, pow(2, avgWeightBit), clkFreq, numReadCellPerOperationNeuro, false, currentMode);
			if (avgWeightBit > 1) {
				multilevelSAEncoder.Initialize(pow(2, avgWeightBit), numCol/numColMuxed);
			}
		}

		dff.Initialize((adderBit+1)*numAdder, clkFreq); 
		adder.Initialize(adderBit, numAdder);
		if (numCellPerSynapse > 1) {
			shiftAddWeight.Initialize(numAdder, adderBit, clkFreq, spikingMode, numCellPerSynapse);
		}
		if (numReadPulse > 1) {
			shiftAddInput.Initialize(numAdder, adderBit+numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
		}
		
	} else if (conventionalParallel) { 
		double resTg = cell.resMemCellOn / numRow;
		
		// if (cell.accessType == CMOS_access) {
		wlNewSwitchMatrix.Initialize(numRow, activityRowRead, clkFreq);         
		// }
		slSwitchMatrix.Initialize(COL_MODE, numCol, resTg*numRow, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, numWritePulseAVG, clkFreq);     
		if (numColMuxed>1) {
			mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);       
			muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
		}
		
		if (SARADC) {
			sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
		} else {
			multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
			multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
		}
		if (numCellPerSynapse > 1) {
			shiftAddWeight.Initialize(ceil(numCol/numColMuxed), log2(levelOutput), clkFreq, spikingMode, numCellPerSynapse);
		}
		if (numReadPulse > 1) {
			shiftAddInput.Initialize(ceil(numCol/numColMuxed), log2(levelOutput)+numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
		}
	} else {
			double resTg = cell.resMemCellOn / numRow;
			
			// if (cell.accessType == CMOS_access) {
			wlNewSwitchMatrix.Initialize(numRow, activityRowRead, clkFreq);         
			// }
			slSwitchMatrix.Initialize(COL_MODE, numCol, resTg*numRow, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, numWritePulseAVG, clkFreq);     
			if (numColMuxed>1) {
				mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);      
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}
			if (SARADC) {
				sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
			} else {
				multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
				multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
			}
			if (numCellPerSynapse > 1) {
				shiftAddWeight.Initialize(ceil(numCol/numColMuxed), log2(levelOutput), clkFreq, spikingMode, numCellPerSynapse);
			}
			if (numReadPulse > 1) {
				shiftAddInput.Initialize(ceil(numCol/numColMuxed), log2(levelOutput)+numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
			}
		}
	// } 
	initialized = true;  //finish initialization
}



void SubArray::CalculateArea() {  //calculate layout area for total design
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;  //ensure initialization first
	} else {  //if initialized, start to do calculation
		area = 0;
		usedArea = 0;
		// else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		// Array only
		heightArray = lengthCol;
		widthArray = lengthRow;
		areaArray = heightArray * widthArray;
		
		if (conventionalSequential) {  
			wlDecoder.CalculateArea(heightArray, NULL, NONE);
			// if (cell.accessType == CMOS_access) {
			wlNewDecoderDriver.CalculateArea(heightArray, NULL, NONE);
			// } 
			slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
			
			if (numColMuxed > 1) {
				mux.CalculateArea(NULL, widthArray, NONE);
				muxDecoder.CalculateArea(NULL, NULL, NONE);
				double minMuxHeight = MAX(muxDecoder.height, mux.height);
				mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
			}
			if (SARADC) {
				sarADC.CalculateUnitArea();
				sarADC.CalculateArea(NULL, widthArray, NONE);
			} else {
				multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
				if (avgWeightBit > 1) {
					multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
				}
			}

			dff.CalculateArea(NULL, widthArray, NONE);
			adder.CalculateArea(NULL, widthArray, NONE);
			if (numReadPulse > 1) {
				shiftAddInput.CalculateArea(NULL, widthArray, NONE);
			}
			if (numCellPerSynapse > 1) {
				shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
			}
			height = slSwitchMatrix.height + heightArray + mux.height + multilevelSenseAmp.height + multilevelSAEncoder.height + adder.height + dff.height + shiftAddInput.height + shiftAddWeight.height + sarADC.height;
			width = MAX(wlDecoder.width + wlNewDecoderDriver.width, muxDecoder.width) + widthArray;
			usedArea = areaArray + wlDecoder.area + wlNewDecoderDriver.area + slSwitchMatrix.area + mux.area + multilevelSenseAmp.area + multilevelSAEncoder.area + muxDecoder.area + adder.area + dff.area + shiftAddInput.area + shiftAddWeight.area + sarADC.area;
			
			areaADC = multilevelSenseAmp.area + multilevelSAEncoder.area + sarADC.area;
			areaAccum = adder.area + dff.area + shiftAddInput.area + shiftAddWeight.area;
			areaOther = wlDecoder.area + wlNewDecoderDriver.area + slSwitchMatrix.area + mux.area + muxDecoder.area;
			

			area = height * width;
			emptyArea = area - usedArea;
			
		} else if (conventionalParallel) {
			// if (cell.accessType == CMOS_access) {
			wlNewSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
			// }
			slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
			if (numColMuxed > 1) {
				mux.CalculateArea(NULL, widthArray, NONE);
				muxDecoder.CalculateArea(NULL, NULL, NONE);
				double minMuxHeight = MAX(muxDecoder.height, mux.height);
				mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
			}
			if (SARADC) {
				sarADC.CalculateUnitArea();
				sarADC.CalculateArea(NULL, widthArray, NONE);
			} else {
				multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
				multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
			}
			if (numReadPulse > 1) {
				shiftAddInput.CalculateArea(NULL, widthArray, NONE);
			}
			if (numCellPerSynapse > 1) {
				shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
			}
			height = slSwitchMatrix.height + heightArray + mux.height + multilevelSenseAmp.height + multilevelSAEncoder.height + shiftAddInput.height + shiftAddWeight.height + sarADC.height;
			width = MAX(wlNewSwitchMatrix.width , muxDecoder.width) + widthArray;
			usedArea = areaArray  + wlNewSwitchMatrix.area + slSwitchMatrix.area + mux.area + multilevelSenseAmp.area + muxDecoder.area + multilevelSAEncoder.area + shiftAddInput.area + shiftAddWeight.area + sarADC.area;
			
			areaADC = multilevelSenseAmp.area + multilevelSAEncoder.area + sarADC.area;
			areaAccum = shiftAddInput.area + shiftAddWeight.area;
			areaOther = wlNewSwitchMatrix.area + slSwitchMatrix.area + mux.area + muxDecoder.area;
			area = height * width;
			emptyArea = area - usedArea;
		}else {   
			// if (cell.accessType == CMOS_access) {
			wlNewSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
			// }
			slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
			if (numColMuxed > 1) {
				mux.CalculateArea(NULL, widthArray, NONE);
				muxDecoder.CalculateArea(NULL, NULL, NONE);
				double minMuxHeight = MAX(muxDecoder.height, mux.height);
				mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
			}
			if (SARADC) {
				sarADC.CalculateUnitArea();
				sarADC.CalculateArea(NULL, widthArray, NONE);
			} else {
				multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
				multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
			}
			
			height = slSwitchMatrix.height + heightArray + mux.height + multilevelSenseAmp.height + multilevelSAEncoder.height + sarADC.height;
			width = MAX(wlNewSwitchMatrix.width, muxDecoder.width) + widthArray;
			area = height * width;
			usedArea = areaArray + wlNewSwitchMatrix.area + slSwitchMatrix.area + mux.area + multilevelSenseAmp.area + muxDecoder.area + multilevelSAEncoder.area + sarADC.area;
			emptyArea = area - usedArea;
		}
	}
}

void SubArray::CalculateLatency(double columnRes, const vector<double> &columnResistance, const vector<double> &rowResistance) {   //calculate latency for different mode 
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		
		readLatency = 0;
		writeLatency = 0;

		if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (conventionalSequential) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				double colRamp = 0;
				double tau = (capCol)*(cell.resMemCellAvg);
				colDelay = horowitz(tau, 0, 1e20, &colRamp);	// Just to generate colRamp
				colDelay = tau * 0.2 * numColMuxed;  // assume the 15~20% voltage drop is enough for sensing
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				
				wlDecoder.CalculateLatency(1e20, capRow2, NULL, numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				// if (cell.accessType == CMOS_access) {
				wlNewDecoderDriver.CalculateLatency(wlDecoder.rampOutput, capRow2, resRow, numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);	
				// }
				slSwitchMatrix.CalculateLatency(1e20, capCol, resCol, 0, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				if (numColMuxed > 1) {
					mux.CalculateLatency(colRamp, 0, numColMuxed);
					muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), numColMuxed, 0);
				}
				if (SARADC) {
					sarADC.CalculateLatency(numColMuxed*numRow*activityRowRead);
				} else {
					multilevelSenseAmp.CalculateLatency(columnResistance, numColMuxed, numRow*activityRowRead);
					if (avgWeightBit > 1) {
						multilevelSAEncoder.CalculateLatency(1e20, numColMuxed*numRow*activityRowRead);
					}
				}
				
				adder.CalculateLatency(1e20, dff.capTgDrain, numColMuxed*numRow*activityRowRead);
				dff.CalculateLatency(1e20, numColMuxed*numRow*activityRowRead);
				if (numCellPerSynapse > 1) {				 
					shiftAddWeight.CalculateLatency(numColMuxed);							// There are numReadPulse times of shift-and-add
				}																								
				if (numReadPulse > 1) {
					shiftAddInput.CalculateLatency(ceil(numColMuxed/numCellPerSynapse));	// There are numReadPulse times of shift-and-add
				}
				
				// Read
				readLatency = 0;
				readLatency += MAX(wlDecoder.readLatency + wlNewDecoderDriver.readLatency , ( ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) )/numReadPulse);
				readLatency += multilevelSenseAmp.readLatency;
				readLatency += multilevelSAEncoder.readLatency;
				readLatency += adder.readLatency;
				readLatency += dff.readLatency;
				readLatency += shiftAddInput.readLatency + shiftAddWeight.readLatency;
				readLatency += colDelay/numReadPulse;
				readLatency += sarADC.readLatency;
				
				readLatencyADC = multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency;
				readLatencyAccum = adder.readLatency + dff.readLatency + shiftAddInput.readLatency + shiftAddWeight.readLatency;
				readLatencyOther = MAX(wlDecoder.readLatency + wlNewDecoderDriver.readLatency , ( ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) )/numReadPulse) + colDelay/numReadPulse;
				
				// Write
				writeLatency = 0;
				writeLatencyArray = 0;
				writeLatencyArray += totalNumWritePulse * cell.writePulseWidth;
				writeLatency += MAX(wlDecoder.writeLatency + wlNewDecoderDriver.writeLatency , slSwitchMatrix.writeLatency);
				writeLatency += writeLatencyArray;
				
			} else if (conventionalParallel) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double colRamp = 0;
				double tau = (capCol)*(cell.resMemCellAvg/(numRow/2));
				colDelay = horowitz(tau, 0, 1e20, &colRamp);
				colDelay = tau * 0.2 * numColMuxed;  // assume the 15~20% voltage drop is enough for sensing
				
				// if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.CalculateLatency(1e20, capRow2, resRow, numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				// } 
				slSwitchMatrix.CalculateLatency(1e20, capCol, resCol, 0, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				if (numColMuxed>1) {
					mux.CalculateLatency(colRamp, 0, numColMuxed);
					muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), numColMuxed, 0);
				}
				if (SARADC) {
					sarADC.CalculateLatency(numColMuxed);
				} else {
					multilevelSenseAmp.CalculateLatency(columnResistance, numColMuxed, 1);
					multilevelSAEncoder.CalculateLatency(1e20, numColMuxed);
				}
				
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculateLatency(numColMuxed);	
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculateLatency(ceil(numColMuxed/numCellPerSynapse));		
				}
				
				// Read
				readLatency = 0;
				readLatency += MAX(wlNewSwitchMatrix.readLatency , ( ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) )/numReadPulse);
				readLatency += multilevelSenseAmp.readLatency;
				readLatency += multilevelSAEncoder.readLatency;
				readLatency += shiftAddInput.readLatency + shiftAddWeight.readLatency;
				readLatency += colDelay/numReadPulse;
				readLatency += sarADC.readLatency;
				
				readLatencyADC = multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency;
				readLatencyAccum = shiftAddInput.readLatency + shiftAddWeight.readLatency;
				readLatencyOther = MAX(wlNewSwitchMatrix.readLatency , ( ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) )/numReadPulse) + colDelay/numReadPulse;

				// Write
				writeLatency = 0;
				writeLatencyArray = 0;
				writeLatencyArray += totalNumWritePulse * cell.writePulseWidth;
				writeLatency += MAX(wlNewSwitchMatrix.writeLatency , slSwitchMatrix.writeLatency);
				writeLatency += writeLatencyArray;
				
			} else {
				double capBL = lengthCol * 0.2e-15/1e-6;
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double colRamp = 0;
				double tau = (capCol)*(cell.resMemCellAvg/(numRow/2));
				colDelay = horowitz(tau, 0, 1e20, &colRamp);
				colDelay = tau * 0.2 * numColMuxed;  // assume the 15~20% voltage drop is enough for sensing
				
				// if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.CalculateLatency(1e20, capRow2, resRow, numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				// }
				slSwitchMatrix.CalculateLatency(1e20, capCol, resCol, 0, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				if (numColMuxed > 1) {
					mux.CalculateLatency(colRamp, 0, numColMuxed);
					muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), numColMuxed, 0);
				}
				if (SARADC) {
					sarADC.CalculateLatency(numColMuxed);
				} else {
					multilevelSenseAmp.CalculateLatency(columnResistance, numColMuxed, 1);
					multilevelSAEncoder.CalculateLatency(1e20, numColMuxed);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculateLatency(numColMuxed);	
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculateLatency(ceil(numColMuxed/numCellPerSynapse));		
				}

				// Read
				readLatency = 0;
				readLatency += MAX(wlNewSwitchMatrix.readLatency, ( ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) )/numReadPulse);
				readLatency += multilevelSenseAmp.readLatency;
				readLatency += multilevelSAEncoder.readLatency;
				readLatency += shiftAddInput.readLatency + shiftAddWeight.readLatency;
				readLatency += colDelay/numReadPulse;
				readLatency += sarADC.readLatency;
				// Write
				
				writeLatency = 0;
				writeLatencyArray = 0;
				writeLatencyArray += totalNumWritePulse * cell.writePulseWidth;
				writeLatency += MAX(wlNewSwitchMatrix.writeLatency , slSwitchMatrix.writeLatency);
				writeLatency += writeLatencyArray;
				
			}
		}
	}
}

void SubArray::CalculatePower(const vector<double> &columnResistance, const vector<double> &rowResistance) {
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		readDynamicEnergyArray = 0;
		
		double numReadOperationPerRow;   // average value (can be non-integer for energy calculation)
		if (numCol > numReadCellPerOperationNeuro)
			numReadOperationPerRow = numCol / numReadCellPerOperationNeuro;
		else
			numReadOperationPerRow = 1;

		double numWriteOperationPerRow;   // average value (can be non-integer for energy calculation)
		if (numCol * activityColWrite > numWriteCellPerOperationNeuro)
			numWriteOperationPerRow = numCol * activityColWrite / numWriteCellPerOperationNeuro;
		else
			numWriteOperationPerRow = 1;

		// if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		leakage = 0;
		if (conventionalSequential) {
			double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
			double numWriteCells = (int)ceil((double)numCol/*numWriteColMuxed*/); 
			int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
			double capBL = lengthCol * 0.2e-15/1e-6;
			
			wlDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
			// if (cell.accessType == CMOS_access) {
			wlNewDecoderDriver.CalculatePower(numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
			// } 
			slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
			if (numColMuxed > 1) {
				mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
				muxDecoder.CalculatePower(numColMuxed, 1);
			}
			
			if (SARADC) {
				sarADC.CalculatePower(columnResistance, numRow*activityRowRead);
			} else {
				multilevelSenseAmp.CalculatePower(columnResistance, numRow*activityRowRead);
				if (avgWeightBit > 1) {
					multilevelSAEncoder.CalculatePower(numRow*activityRowRead*numColMuxed);
				}
			}

			adder.CalculatePower(numColMuxed*numRow*activityRowRead, numReadCells);
			dff.CalculatePower(numColMuxed*numRow*activityRowRead, numReadCells*(adder.numBit+1)); 
			if (numCellPerSynapse > 1) {
				shiftAddWeight.CalculatePower(numColMuxed);	
			}
			if (numReadPulse > 1) {
				shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));		
			}
			// Read
			readDynamicEnergyArray = 0;
			readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells; // Selected BLs activityColWrite
			readDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd; // Selected WL
			readDynamicEnergyArray *= numRow * activityRowRead * numColMuxed;
			readDynamicEnergy = 0;
			readDynamicEnergy += wlDecoder.readDynamicEnergy;
			readDynamicEnergy += wlNewDecoderDriver.readDynamicEnergy;
			readDynamicEnergy +=  ( ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0) )/numReadPulse;
			readDynamicEnergy += adder.readDynamicEnergy;
			readDynamicEnergy += dff.readDynamicEnergy;
			readDynamicEnergy += shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
			readDynamicEnergy += readDynamicEnergyArray;
			readDynamicEnergy += sarADC.readDynamicEnergy;
			
			readDynamicEnergyADC = readDynamicEnergyArray + multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy + sarADC.readDynamicEnergy;
			readDynamicEnergyAccum = adder.readDynamicEnergy + dff.readDynamicEnergy + shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
			readDynamicEnergyOther = wlDecoder.readDynamicEnergy + wlNewDecoderDriver.readDynamicEnergy + ( ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0) )/numReadPulse;

			// Write
			writeDynamicEnergyArray = writeDynamicEnergyArray;
			writeDynamicEnergy = 0;
			writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
			writeDynamicEnergy += wlNewDecoderDriver.writeDynamicEnergy;
			writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
			writeDynamicEnergy += writeDynamicEnergyArray;
			
			// Leakage
			
			leakage += wlDecoder.leakage;
			leakage += wlNewDecoderDriver.leakage;
			leakage += slSwitchMatrix.leakage;
			leakage += mux.leakage;
			leakage += muxDecoder.leakage;
			leakage += multilevelSenseAmp.leakage;
			leakage += multilevelSAEncoder.leakage;
			leakage += dff.leakage;
			leakage += adder.leakage;
			leakage += shiftAddInput.leakage + shiftAddWeight.leakage;
				
		} else if (conventionalParallel) {
			double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
			int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
			double capBL = lengthCol * 0.2e-15/1e-6;
		
			// if (cell.accessType == CMOS_access) {
			wlNewSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
			// }
			slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
			if (numColMuxed > 1) {
				mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
				muxDecoder.CalculatePower(numColMuxed, 1);
			}
			
			if (SARADC) {
				sarADC.CalculatePower(columnResistance, 1);
			} else {
				multilevelSenseAmp.CalculatePower(columnResistance, 1);
				multilevelSAEncoder.CalculatePower(numColMuxed);
			}
			if (numCellPerSynapse > 1) {
				shiftAddWeight.CalculatePower(numColMuxed);	
			}
			if (numReadPulse > 1) {
				shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));		
			}

			// Read
			readDynamicEnergyArray = 0;
			readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells; // Selected BLs activityColWrite
			readDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd * numRow * activityRowRead; // Selected WL
			readDynamicEnergyArray *= numColMuxed;
			
			readDynamicEnergy = 0;
			readDynamicEnergy += wlNewSwitchMatrix.readDynamicEnergy;
			readDynamicEnergy += ( ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0) )/numReadPulse;
			readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
			readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
			readDynamicEnergy += shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
			readDynamicEnergy += readDynamicEnergyArray;
			if (!isnan(sarADC.readDynamicEnergy))
			{
				readDynamicEnergy += sarADC.readDynamicEnergy;
			}
			
			
			readDynamicEnergyADC = readDynamicEnergyArray + multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy + sarADC.readDynamicEnergy;
			readDynamicEnergyAccum = shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
			readDynamicEnergyOther = wlNewSwitchMatrix.readDynamicEnergy + ( ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0) )/numReadPulse;
			
			// Write
			writeDynamicEnergyArray = writeDynamicEnergyArray;
			writeDynamicEnergy = 0;
			writeDynamicEnergy += wlNewSwitchMatrix.writeDynamicEnergy;
			writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
			writeDynamicEnergy += writeDynamicEnergyArray;
			
			// Leakage
			leakage += wlNewSwitchMatrix.leakage;
			leakage += slSwitchMatrix.leakage;
			leakage += mux.leakage;
			leakage += muxDecoder.leakage;
			leakage += multilevelSenseAmp.leakage;
			leakage += multilevelSAEncoder.leakage;
			leakage += shiftAddWeight.leakage + shiftAddInput.leakage;
		} else {
			double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
			int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
			double capBL = lengthCol * 0.2e-15/1e-6;
		
			// if (cell.accessType == CMOS_access) {
			wlNewSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
			// }
			slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
			if (numColMuxed > 1) {
				mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
				muxDecoder.CalculatePower(numColMuxed, 1);
			}
			if (SARADC) {
				sarADC.CalculatePower(columnResistance, numColMuxed);
			} else {
				multilevelSenseAmp.CalculatePower(columnResistance, numColMuxed);
				multilevelSAEncoder.CalculatePower(numColMuxed);
			}
			if (numCellPerSynapse > 1) {
				shiftAddWeight.CalculatePower(numColMuxed);	
			}
			if (numReadPulse > 1) {
				shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));		
			}
			// Read
			readDynamicEnergyArray = 0;
			readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells; // Selected BLs activityColWrite
			readDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd * numRow * activityRowRead; // Selected WL
			readDynamicEnergyArray *= numReadPulse * numColMuxed;
			
			readDynamicEnergy = 0;
			readDynamicEnergy += wlNewSwitchMatrix.readDynamicEnergy;
			readDynamicEnergy += ( ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0) )/numReadPulse;
			readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
			readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
			readDynamicEnergy += shiftAddInput.readDynamicEnergy + shiftAddWeight.readDynamicEnergy;
			readDynamicEnergy += readDynamicEnergyArray;
			readDynamicEnergy += sarADC.readDynamicEnergy;
			
			// Write
			
			writeDynamicEnergyArray = writeDynamicEnergyArray;
			writeDynamicEnergy = 0;
			writeDynamicEnergy += wlNewSwitchMatrix.writeDynamicEnergy;
			writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
			writeDynamicEnergy += writeDynamicEnergyArray;
			
			
			// Leakage
			leakage = 0;
			leakage += wlNewSwitchMatrix.leakage;
			leakage += slSwitchMatrix.leakage;
			leakage += mux.leakage;
			leakage += muxDecoder.leakage;
			leakage += multilevelSenseAmp.leakage;
			leakage += multilevelSAEncoder.leakage;
			leakage += shiftAddInput.leakage + shiftAddWeight.leakage;
		}
		// } 
	}
}

void SubArray::PrintProperty() {

	if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		
		cout << endl << endl;
	    cout << "Array:" << endl;
	    cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
	    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
	    //cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;
		//cout << "Write Latency = " << writeLatencyArray*1e9 << "ns" << endl;

		if (conventionalSequential) {
			wlDecoder.PrintProperty("wlDecoder");
			// if (cell.accessType == CMOS_access) {
			wlNewDecoderDriver.PrintProperty("wlNewDecoderDriver");
			// }
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			multilevelSenseAmp.PrintProperty("multilevelSenseAmp or single-bit SenseAmp");
			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			adder.PrintProperty("adder");
			dff.PrintProperty("dff");
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput"); 
			}
		} else if (conventionalParallel) {
			// if (cell.accessType == CMOS_access) {
			wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
			// }
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");
			}
		} else {
			// if (cell.accessType == CMOS_access) {
			wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
			// } 
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																					   
			}
		}
	} 
	FunctionUnit::PrintProperty("SubArray");
	cout << "Used Area = " << usedArea*1e12 << "um^2" << endl;
	cout << "Empty Area = " << emptyArea*1e12 << "um^2" << endl;
}

