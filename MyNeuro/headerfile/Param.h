#ifndef PARAM_H_
#define PARAM_H_

class Param {
public:
	Param();

	int operationmode, operationmodeBack, memcelltype, accesstype, transistortype, deviceroadmap;      		
	
	double heightInFeatureSizeSRAM, widthInFeatureSizeSRAM, widthSRAMCellNMOS, widthSRAMCellPMOS, widthAccessCMOS, minSenseVoltage;
	
	double heightInFeatureSize1T1R, widthInFeatureSize1T1R, heightInFeatureSizeCrossbar, widthInFeatureSizeCrossbar;
	
	int relaxArrayCellHeight, relaxArrayCellWidth;
	
	bool globalBufferType, tileBufferType, peBufferType, chipActivation, reLu, novelMapping, pipeline, trainingEstimation, parallelBP, nonlinearIV, SARADC, currentMode;
	int globalBufferCoreSizeRow, globalBufferCoreSizeCol, tileBufferCoreSizeRow, tileBufferCoreSizeCol;
	
	double clkFreq, featuresize, readNoise, resistanceOn, resistanceOff, maxConductance, minConductance, gateCapFeFET, polarization;
	int temp, technode, wireWidth, multipleCells;
	double maxNumLevelLTP, maxNumLevelLTD, readVoltage, readPulseWidth, writeVoltage;
	double accessVoltage, resistanceAccess;
	double nonlinearity;
	double writePulseWidth, numWritePulse;
	double globalBusDelayTolerance, localBusDelayTolerance;
	double treeFoldedRatio, maxGlobalBusWidth;
	double algoWeightMax, algoWeightMin;
	double activityRowReadWG, activityRowWriteWG, activityColWriteWG;
	double bufferOverHeadConstraint;
	
	int neuro, multifunctional, parallelWrite, parallelRead;
	int numlut, numColMuxed, numWriteColMuxed, levelOutput, avgWeightBit, numBitInput, numRowMuxedAG, levelOutputAG, numRowMuxedWG, levelOutputWG;
	int numRowSubArray, numColSubArray, numRowSubArrayWG, numColSubArrayWG;
	int cellBit, synapseBit;
	int speedUpDegree, dramType, batchSize, numIteration;
	
	int conventionalParallel, conventionalSequential; 
	int numRowPerSynapse, numColPerSynapse;
	double AR, Rho, wireLengthRow, wireLengthCol, unitLengthWireResistance, wireResistanceRow, wireResistanceCol;
};

#endif