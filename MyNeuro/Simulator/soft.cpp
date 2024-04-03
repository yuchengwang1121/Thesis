#include <cstdio>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <string>
#include "../headerfile/function.h"
#include "../headerfile/constant.h"
#include "../headerfile/formula.h"
#include "../headerfile/Param.h"
#include "../headerfile/SubArray.h"
#include "../headerfile/AdderTree.h"
#include "../headerfile/Definition.h"
#include "../headerfile/Bus.h"
#include "../headerfile/DFF.h"

using namespace std;
SubArray *CAM;
SubArray *LUT;

int main(int argc, char *argv[])
{
    gen.seed(0);

    vector<vector<double>> netStructure;
	// cout << "get net "<<argv[2] << endl;
	bool STAR = (atoi(argv[2]) == 1);;	// <=== Modify When Change Structure ===>
    netStructure = getNetStructure(argv[3]);                                        // get file from trace.command.sh
	// cout << "get bit "<<argv[3] << endl;
    // define weight/input/memory precision from wrapper
    param->synapseBit = atoi(argv[4]);  // precision of synapse weight
    param->numBitInput = atoi(argv[5]); // precision of input neural activation

	string softtype = STAR ? "STAR" : "TransSeg";
	cout << "------------------------------ Softmax Type --------------------------------" << endl;
	cout << "==> Type : "<< softtype << endl;

    if (param->cellBit > param->synapseBit)
    {
        cout << "ERROR!: Memory precision is even higher than synapse precision, please modify 'cellBit' in Param.cpp!" << endl;
        param->cellBit = param->synapseBit;
    }

    /*** initialize operationMode as default ***/
    param->conventionalParallel = 0;
    param->conventionalSequential = 0;
    switch (param->operationmode)
    {
    case 2:
        param->conventionalParallel = 1;
        break;
    case 1:
        param->conventionalSequential = 1;
        break;
    case -1:
        break;
    default:
        exit(-1);
    }
    param->numRowPerSynapse = 1;
    param->numColPerSynapse = ceil((double)param->synapseBit / (double)param->cellBit);
    inputParameter.transistorType = conventional;
    inputParameter.deviceRoadmap = HP;

    /* Create SubArray object and link the required global objects (not initialization) */
    inputParameter.temperature = param->temp;     // Temperature (K)
    inputParameter.processNode = param->technode; // Technology node
    //Initial the Technology "tech" defined at ../headerfile/Definition.h
    //Give the currentOn/OffNmos, currentOn/OffPmos values
    tech.Initialize(inputParameter.processNode, inputParameter.deviceRoadmap, inputParameter.transistorType);
	
    /******************************************************** Initialize ********************************************************/
    cout << "------------------------------ FloorPlan --------------------------------" << endl;
    cout << endl;
    
	/*** circuit level parameters ***/
	cell.memCellType = Type::RRAM;
	cell.accessType = CMOS_access;
	
	switch(param->deviceroadmap) {
		case 2:	    inputParameter.deviceRoadmap = LSTP;  break;
		case 1:	    inputParameter.deviceRoadmap = HP;    break;
		case -1:	break;
		default:	exit(-1);
	}
	
    /*** build object from each class ***/
	CAM 			= new SubArray(inputParameter, tech, cell);
	LUT 			= new SubArray(inputParameter, tech, cell);
	
    /*** RRAM cell's property ***/ 
	cell.resistanceOn = param->resistanceOn;	                                // Ron resistance at Vr in the reported measurement data (need to recalculate below if considering the nonlinearity)
	cell.resistanceOff = param->resistanceOff;	                                // Roff resistance at Vr in the reported measurement dat (need to recalculate below if considering the nonlinearity)
	cell.resistanceAvg = (cell.resistanceOn + cell.resistanceOff)/2;            // Average resistance (for energy estimation)
	cell.readVoltage = param->readVoltage;	                                    // On-chip read voltage for memory cell
	cell.readPulseWidth = param->readPulseWidth;
	cell.accessVoltage = param->accessVoltage;                                  // Gate voltage for the transistor in 1T1R
	cell.resistanceAccess = param->resistanceAccess;
	cell.featureSize = param->featuresize; 
	cell.maxNumLevelLTP = param->maxNumLevelLTP;	                            // Maximum number of conductance states during LTP or weight increase
	cell.maxNumLevelLTD = param->maxNumLevelLTD;	                            // Maximum number of conductance states during LTD or weight decrease
	double writeVoltageLTP = param->writeVoltage;
	double writeVoltageLTD = param->writeVoltage;
	cell.writeVoltage = sqrt(writeVoltageLTP * writeVoltageLTP + writeVoltageLTD * writeVoltageLTD);    // Use an average value of write voltage for NeuroSim
	double writePulseWidthLTP = param->writePulseWidth;
	double writePulseWidthLTD = param->writePulseWidth;
	cell.writePulseWidth = (writePulseWidthLTP + writePulseWidthLTD) / 2;
	cell.nonlinearIV = param->nonlinearIV; 										 // This option is to consider I-V nonlinearity in cross-point array or not
	cell.nonlinearity = param->nonlinearity; 									 // This is the nonlinearity for the current ratio at Vw and Vw/2
    cell.heightInFeatureSize = param->heightInFeatureSize1T1R;         // Cell height in feature size
	cell.widthInFeatureSize = param->widthInFeatureSize1T1R;

    /*** CAM's property ***/
	int CAMnumRow = 64;		//each vector has 64 bit
	int CAMnumCol = 64;		//hit vector from -20 ~ 43
	int CAMclk    = (STAR)? 20e9 : 15e9;

	CAM->conventionalParallel = param->conventionalParallel;                  
	CAM->conventionalSequential = param->conventionalSequential;   
	CAM->parallelBP = param->parallelBP;	
	CAM->numRow = CAMnumRow;
	CAM->numCol = CAMnumCol;
	CAM->levelOutput = param->levelOutput;
	CAM->levelOutputBP = param->levelOutputAG;
	CAM->numColMuxed = param->numColMuxed;               // How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
	CAM->numRowMuxedBP = param->numRowMuxedAG;
    CAM->clkFreq = CAMclk;                       // Clock frequency
	CAM->relaxArrayCellHeight = param->relaxArrayCellHeight;
	CAM->relaxArrayCellWidth = param->relaxArrayCellWidth;
	CAM->numReadPulse = param->numBitInput;
	CAM->avgWeightBit = param->cellBit;
	CAM->numCellPerSynapse = param->numColPerSynapse;
	CAM->numReadPulseBP = 8;
	CAM->activityBPColRead = 0.5;
	CAM->SARADC = param->SARADC;
	CAM->currentMode = param->currentMode;
	CAM->spikingMode = NONSPIKING;
	
	
	if (CAM->numColMuxed > CAMnumCol) {                      // Set the upperbound of numColMuxed
		CAM->numColMuxed = CAMnumCol;
	}

	CAM->numReadCellPerOperationFPGA = CAMnumCol;	           // Not relevant for IMEC
	CAM->numWriteCellPerOperationFPGA = CAMnumCol;	       // Not relevant for IMEC
	CAM->numReadCellPerOperationMemory = CAMnumCol;          // Define # of SRAM read cells in memory mode because SRAM does not have S/A sharing (not relevant for IMEC)
	CAM->numWriteCellPerOperationMemory = CAMnumCol/8;       // # of write cells per operation in SRAM memory or the memory mode of multifunctional memory (not relevant for IMEC)
	CAM->numReadCellPerOperationNeuro = CAMnumCol;           // # of SRAM read cells in neuromorphic mode
	CAM->numWriteCellPerOperationNeuro = CAMnumCol;	       // For SRAM or analog RRAM in neuro mode
    CAM->maxNumWritePulse = MAX(cell.maxNumLevelLTP, cell.maxNumLevelLTD);

	/*** LUT's property ***/
	int LUTnumRow = (STAR)? 64 : 16;		//from -15 ~ 0
	int LUTnumCol = 32;		//use 32bit ti reprecent precision 8 LUT
	int LUTclk    = (STAR)? 20e9 : 15e9;

	LUT->conventionalParallel = param->conventionalParallel;                  
	LUT->conventionalSequential = param->conventionalSequential;   
	LUT->parallelBP = param->parallelBP;	
	LUT->numRow = LUTnumRow;
	LUT->numCol = LUTnumCol;
	LUT->levelOutput = param->levelOutput;
	LUT->levelOutputBP = param->levelOutputAG;
	LUT->numColMuxed = param->numColMuxed;               // How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
	LUT->numRowMuxedBP = param->numRowMuxedAG;
    LUT->clkFreq = LUTclk;                       // Clock frequency
	LUT->relaxArrayCellHeight = param->relaxArrayCellHeight;
	LUT->relaxArrayCellWidth = param->relaxArrayCellWidth;
	LUT->numReadPulse = param->numBitInput;
	LUT->avgWeightBit = param->cellBit;
	LUT->numCellPerSynapse = param->numColPerSynapse;
	LUT->numReadPulseBP = 8;
	LUT->activityBPColRead = 0.5;
	LUT->SARADC = param->SARADC;
	LUT->currentMode = param->currentMode;
	LUT->spikingMode = NONSPIKING;
	
	
	if (LUT->numColMuxed > LUTnumCol) {                      // Set the upperbound of numColMuxed
		LUT->numColMuxed = LUTnumCol;
	}

	LUT->numReadCellPerOperationFPGA = LUTnumCol;	           // Not relevant for IMEC
	LUT->numWriteCellPerOperationFPGA = LUTnumCol;	       // Not relevant for IMEC
	LUT->numReadCellPerOperationMemory = LUTnumCol;          // Define # of SRAM read cells in memory mode because SRAM does not have S/A sharing (not relevant for IMEC)
	LUT->numWriteCellPerOperationMemory = LUTnumCol/8;       // # of write cells per operation in SRAM memory or the memory mode of multifunctional memory (not relevant for IMEC)
	LUT->numReadCellPerOperationNeuro = LUTnumCol;           // # of SRAM read cells in neuromorphic mode
	LUT->numWriteCellPerOperationNeuro = LUTnumCol;	       // For SRAM or analog RRAM in neuro mode
    LUT->maxNumWritePulse = MAX(cell.maxNumLevelLTP, cell.maxNumLevelLTD);

	/*** User defined num of elements ***/
	int CAMnumSubArrayRow = 2;		// The number of subarray's row
	int CAMnumSubArrayCol = 1; 	// The number of subarray's col
	int LUTnumSubArrayRow = 2;		// The number of subarray's row
	int LUTnumSubArrayCol = 1; 	// The number of subarray's col
	/*** initialize modules ***/
	CAM->Initialize(CAMnumRow, CAMnumCol, param->unitLengthWireResistance);        // initialize CAM
	LUT->Initialize(LUTnumRow, LUTnumCol, param->unitLengthWireResistance);        // initialize LUT

	cout << "number of Segment is " << atoi(argv[6]) << " with Input size is " << 16/atoi(argv[6]) << endl;
    cout << "number of CAM subarray's row is " << CAMnumSubArrayRow << " with number of subarray's col is " << CAMnumSubArrayCol << endl;
	cout << "number of LUT subarray's row is " << LUTnumSubArrayRow << " with number of subarray's col is " << LUTnumSubArrayCol << endl;
    cout << endl;
    cout << "---------------------------- FloorPlan Done ------------------------------" << endl;
    cout << endl;

    /******************************************************** Initialize ********************************************************/

    /****************************************************** CalculateArea *******************************************************/
    double CAMArea, CAMAreaArray, CAMAreaADC, CAMAreaAccum, CAMAreaOther;
	double LUTArea, LUTAreaArray, LUTAreaADC, LUTAreaAccum, LUTAreaOther;
	double heightCAM, heightLUT = 0;
	double widthCAM, widthLUT= 0;
	double busarea = 0;
	double bufferarea = 0;
	double TotalareaCAM, TotalareaLUT = 0;
	double widthArray = 0;
	vector<double> CAMareaResults, LUTareaResults;
	
	CAM->CalculateArea();
    LUT->CalculateArea();
    TotalareaCAM = CAM->usedArea * (CAMnumSubArrayRow*CAMnumSubArrayCol);
	TotalareaLUT = LUT->usedArea * (LUTnumSubArrayRow*LUTnumSubArrayCol);

    
    heightCAM = sqrt(TotalareaCAM);
    widthCAM = TotalareaCAM/(heightCAM);
	heightLUT = sqrt(TotalareaLUT);
    widthLUT = TotalareaCAM/(heightLUT);

    CAMareaResults.push_back(TotalareaCAM);
    CAMareaResults.push_back(CAM->areaADC*(CAMnumSubArrayRow*CAMnumSubArrayCol));
    CAMareaResults.push_back(CAM->areaAccum*(CAMnumSubArrayRow*CAMnumSubArrayCol));
    CAMareaResults.push_back(CAM->areaOther*(CAMnumSubArrayRow*CAMnumSubArrayCol));
    CAMareaResults.push_back(CAM->areaArray*(CAMnumSubArrayRow*CAMnumSubArrayCol));
	
    CAMArea        = CAMareaResults[0];
    CAMAreaArray   = CAMareaResults[1];
    CAMAreaADC     = CAMareaResults[2];
    CAMAreaAccum   = CAMareaResults[3];
    CAMAreaOther   = CAMareaResults[4];

	LUTareaResults.push_back(TotalareaLUT);
    LUTareaResults.push_back(LUT->areaADC*(LUTnumSubArrayRow*LUTnumSubArrayCol));
    LUTareaResults.push_back(LUT->areaAccum*(LUTnumSubArrayRow*LUTnumSubArrayCol));
    LUTareaResults.push_back(LUT->areaOther*(LUTnumSubArrayRow*LUTnumSubArrayCol));
    LUTareaResults.push_back(LUT->areaArray*(LUTnumSubArrayRow*LUTnumSubArrayCol));
	
    LUTArea        = LUTareaResults[0];
    LUTAreaArray   = LUTareaResults[1];
    LUTAreaADC     = LUTareaResults[2];
    LUTAreaAccum   = LUTareaResults[3];
    LUTAreaOther   = LUTareaResults[4];

    cout << "-------------------------------------- Hardware Performance --------------------------------------" << endl;
	
    // save breakdown results of each layer to csv files
    ofstream breakdownfile;
    string breakdownfile_name = "./NeuroSim_Breakdown";
    breakdownfile_name.append(".csv");
    breakdownfile.open(breakdownfile_name);
    if (breakdownfile.is_open()){
        // firstly save the area results to file
        breakdownfile << "CAM Area(m^2), ADC Area(m^2), Accumulation Area(m^2), Other Logic&Storage Area(m^2), CAM CIM (FW+AG) Area (m^2)" << endl;
        breakdownfile << CAMArea << "," << CAMAreaADC << "," << CAMAreaAccum << "," << CAMAreaOther << "," << CAMAreaArray << endl;
		breakdownfile << "LUT Area(m^2), ADC Area(m^2), Accumulation Area(m^2), Other Logic&Storage Area(m^2), LUT CIM (FW+AG) Area (m^2)" << endl;
        breakdownfile << LUTArea << "," << LUTAreaADC << "," << LUTAreaAccum << "," << LUTAreaOther << "," << LUTAreaArray << endl;
		breakdownfile << "Total Area(m^2), ADC Area(m^2), Accumulation Area(m^2), Other Logic&Storage Area(m^2), Total CIM (FW+AG) Area (m^2)" << endl;
        breakdownfile << CAMArea + LUTArea << "," << CAMAreaADC + LUTAreaADC << "," 
					  << CAMAreaAccum + LUTAreaAccum << "," << CAMAreaOther + LUTAreaOther << ","
					  << CAMAreaArray + LUTAreaArray << endl;
        breakdownfile << endl;
    }
    else
    {
        cout << "Error: the breakdown file cannot be opened!" << endl;
    }
    /****************************************************** CalculateArea *******************************************************/

    /*************************************************** CalculatePerformance ***************************************************/

    double numComputation = 0;
	numComputation = 2*(netStructure[0][0] * netStructure[0][1] * netStructure[0][2] * netStructure[0][3] * netStructure[0][4] * netStructure[0][5]);
    /*** define how many subArray are used to map the whole layer ***/
	double ICReadLatency, ICReadDynamicEnergy, ICLeakage, ICLatencyADC, ICLatencyAccum, ICLatencyOther, ICEnergyADC, ICEnergyAccum, ICEnergyOther;
	double SLReadLatency, SLReadDynamicEnergy, SLLeakage, SLLatencyADC, SLLatencyAccum, SLLatencyOther, SLEnergyADC, SLEnergyAccum, SLEnergyOther;
	double CVReadLatency, CVReadDynamicEnergy, CVLeakage, CVLatencyADC, CVLatencyAccum, CVLatencyOther, CVEnergyADC, CVEnergyAccum, CVEnergyOther;
	double ReadLatency, ReadDynamicEnergy, Leakage, LatencyADC, LatencyOther, LatencyAccum, EnergyADC, EnergyAccum, EnergyOther;
	double BufferLatency, LeakageEnergy, BufferDynamicEnergy, BusLatency, BusDynamicEnergy;

	//Input * CAM
	ICReadLatency		=0;
	ICReadDynamicEnergy	=0;
	ICLeakage			=0;
	ICLatencyADC		=0;
	ICLatencyAccum		=0;
	ICLatencyOther		=0;
	ICEnergyADC			=0;
	ICEnergyAccum		=0;
	ICEnergyOther		=0;
	//Sub * LUT
	SLReadLatency		=0;
	SLReadDynamicEnergy	=0;
	SLLeakage			=0;
	SLLatencyADC		=0;
	SLLatencyAccum		=0;
	SLLatencyOther		=0;
	SLEnergyADC			=0;
	SLEnergyAccum		=0;
	SLEnergyOther		=0;
	//Counter * LUT
	CVReadLatency		=0;
	CVReadDynamicEnergy	=0;
	CVLeakage			=0;
	CVLatencyADC		=0;
	CVLatencyAccum		=0;
	CVLatencyOther		=0;
	CVEnergyADC			=0;
	CVEnergyAccum		=0;
	CVEnergyOther		=0;
	//Other
	BufferLatency 		= 0;
	BufferDynamicEnergy = 0;
	BusLatency 			= 0;
	BusDynamicEnergy 	= 0;
	LeakageEnergy 		= 0;
	Leakage				= 0;
	ReadLatency			= 0;
	ReadDynamicEnergy	= 0;
	LatencyADC 			= 0;
	LatencyAccum		= 0;
	LatencyOther		= 0;
	EnergyADC 			= 0;
	EnergyAccum			= 0;
	EnergyOther			= 0;

    // weight matrix is further partitioned inside PE (among subArray) --> no duplicated
	int numInVector, numSeg;
	numSeg = atoi(argv[6]);
	numInVector = 16/numSeg;

	// load in whole file 
	vector<vector<double>> CAMMemory,LUTMemory;
	vector<vector<double>> CAMinputVector, LUTinputVector;
	CAMMemory = LoadInWeightData(argv[7], 1, 1, param->maxConductance, param->minConductance);
	LUTMemory = LoadInWeightData(argv[8], 1, 1, param->maxConductance, param->minConductance);
	CAMinputVector = LoadInInputData(argv[9]);
	LUTinputVector = LoadInInputData(argv[10]);

    /*** assign weight and input to specific subArray ***/
    vector<vector<double>> subCAMMemory, subLUTMemory;
	vector<vector<double>> CAMArrayInput, VMMArrayInput;

    subCAMMemory = CopySubArray(CAMMemory, 0, 0, CAMnumRow, CAMnumCol);		//64*64
	subLUTMemory = CopySubArray(LUTMemory, 0, 0, LUTnumRow, LUTnumCol);		//16*32
	CAMArrayInput = CopySubInput(CAMinputVector, 0, numInVector, CAMnumRow);//16*64
	VMMArrayInput = CopySubInput(LUTinputVector, 0, 1, LUTnumRow);			//1*16
	// cout << "CopySubInput : " << endl;
	
    for (int k=0; k<numInVector*3; k++) {	// calculate three times to FindMax & FindSub & Find SUBMV
        // cout << "k is : " << k << " with numInVector " << numInVector << endl;
		double CAMactivityRowRead = 0;
        vector<double> input;
        input = GetInputVector(CAMArrayInput, k, &CAMactivityRowRead);
		CAM->activityRowRead = CAMactivityRowRead;
        int cellRange = pow(2, param->cellBit);
        if (param->parallelRead) {
            CAM->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
        } else {
            CAM->levelOutput = cellRange;
        }
        
        vector<double> columnResistance, rowResistance;
		// cout << "columnResistance" << endl;
        columnResistance = GetColumnResistance(input, subCAMMemory, cell, param->parallelRead, CAM->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subCAMMemory, cell, param->parallelBP, CAM->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        CAM->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		CAM->CalculatePower(columnResistance, rowResistance);
        
        ICReadLatency += CAM->readLatency;
        ICReadDynamicEnergy += CAM->readDynamicEnergy;
        ICLeakage += CAM->leakage;
        
        ICLatencyADC += CAM->readLatencyADC;
        ICLatencyAccum += CAM->readLatencyAccum;
        ICLatencyOther += CAM->readLatencyOther;
        
        ICEnergyADC += CAM->readDynamicEnergyADC;
        ICEnergyAccum += CAM->readDynamicEnergyAccum;
        ICEnergyOther += CAM->readDynamicEnergyOther;
    }

	for (int k=0; k<numInVector; k++) {	// calculate MV get exponential x_i-x_max
        // cout << "k is : " << k << " with numInVector " << numInVector << endl;
		double LUTactivityRowRead = 0;
        vector<double> input;
        input = GetInputVector(CAMArrayInput, k, &LUTactivityRowRead);
		LUT->activityRowRead = LUTactivityRowRead;
        int cellRange = pow(2, param->cellBit);
        if (param->parallelRead) {
            LUT->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
        } else {
            LUT->levelOutput = cellRange;
        }
        
        vector<double> columnResistance, rowResistance;
		// cout << "columnResistance" << endl;
        columnResistance = GetColumnResistance(input, subLUTMemory, cell, param->parallelRead, LUT->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subLUTMemory, cell, param->parallelBP, LUT->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        LUT->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		LUT->CalculatePower(columnResistance, rowResistance);
        
        SLReadLatency += LUT->readLatency;
        SLReadDynamicEnergy += LUT->readDynamicEnergy;
        SLLeakage += LUT->leakage;
        
        SLLatencyADC += LUT->readLatencyADC;
        SLLatencyAccum += LUT->readLatencyAccum;
        SLLatencyOther += LUT->readLatencyOther;
        
        SLEnergyADC += LUT->readDynamicEnergyADC;
        SLEnergyAccum += LUT->readDynamicEnergyAccum;
        SLEnergyOther += LUT->readDynamicEnergyOther;
    }

	for (int k=0; k<1; k++) {	// calculate one counter vector * LUT get sum of exponential x_i-x_max
        // cout << "k is : " << k << " with numInVector " << numInVector << endl;
		double VMMactivityRowRead = 0;
        vector<double> input;
        input = GetInputVector(VMMArrayInput, k, &VMMactivityRowRead);
		LUT->activityRowRead = VMMactivityRowRead;
        int cellRange = pow(2, param->cellBit);
        if (param->parallelRead) {
            LUT->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
        } else {
            LUT->levelOutput = cellRange;
        }
        
        vector<double> columnResistance, rowResistance;
		// cout << "columnResistance" << endl;
        columnResistance = GetColumnResistance(input, subLUTMemory, cell, param->parallelRead, LUT->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subLUTMemory, cell, param->parallelBP, LUT->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        LUT->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		LUT->CalculatePower(columnResistance, rowResistance);
        
        CVReadLatency += LUT->readLatency;
        CVReadDynamicEnergy += LUT->readDynamicEnergy;
        CVLeakage += LUT->leakage;
        
        CVLatencyADC += LUT->readLatencyADC;
        CVLatencyAccum += LUT->readLatencyAccum;
        CVLatencyOther += LUT->readLatencyOther;
        
        CVEnergyADC += LUT->readDynamicEnergyADC;
        CVEnergyAccum += LUT->readDynamicEnergyAccum;
        CVEnergyOther += LUT->readDynamicEnergyOther;
    }


    Leakage = numInVector*numSeg*(ICLeakage + SLLeakage + CVLeakage); //Weight Q,K,V & transpose I & I
	
	LatencyADC = numInVector*numSeg*(ICLatencyADC + SLLatencyADC + CVLatencyADC);
	LatencyOther = numInVector*numSeg*(ICLatencyOther + SLLatencyOther + CVLatencyOther);
	LatencyAccum = numInVector*numSeg*(ICLatencyAccum + SLLatencyAccum + CVLatencyAccum);

	EnergyADC = numInVector*numSeg*(ICEnergyADC + SLEnergyADC + CVEnergyADC);
	EnergyOther = numInVector*numSeg*(ICEnergyOther + SLEnergyOther + CVEnergyOther);
	EnergyAccum = numInVector*numSeg*(ICEnergyAccum + SLEnergyAccum + CVEnergyAccum);
    ReadLatency = numInVector*numSeg*(ICReadLatency + SLReadLatency + CVReadLatency);
    ReadDynamicEnergy = numInVector*numSeg*(ICReadDynamicEnergy + SLReadDynamicEnergy + CVReadDynamicEnergy);
    LeakageEnergy = Leakage * ReadLatency;

    cout << "------------------------------ Summary --------------------------------" << endl;
    cout << endl;
	cout << "---------- Area of CAM ----------" << endl;
    cout << " ===========>> Area : " << CAMArea * 1e12 << "um^2" << endl;
    cout << "Total CIM (Forward+Activation Gradient) array : " << CAMAreaArray * 1e12 << "um^2" << endl;
    cout << "Total ADC (or S/As and precharger for SRAM) Area on chip : " << CAMAreaADC * 1e12 << "um^2" << endl;
    cout << "Total Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) on chip : " << CAMAreaAccum * 1e12 << "um^2" << endl;
    cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, pooling and activation units) : " << CAMAreaOther * 1e12 << "um^2" << endl;
    cout << endl;
	cout << "---------- Area of LUT ----------" << endl;
    cout << "Area : " << LUTArea * 1e12 << "um^2" << endl;
    cout << "Total CIM (Forward+Activation Gradient) array : " << LUTAreaArray * 1e12 << "um^2" << endl;
    cout << "Total ADC (or S/As and precharger for SRAM) Area on chip : " << LUTAreaADC * 1e12 << "um^2" << endl;
    cout << "Total Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) on chip : " << LUTAreaAccum * 1e12 << "um^2" << endl;
    cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, pooling and activation units) : " << LUTAreaOther * 1e12 << "um^2" << endl;
    cout << endl;
	cout << "---------- Area Summary ----------" << endl;
    cout << "Area : " << (CAMArea + LUTArea) * 1e12 << "um^2" << endl;
    cout << "Total CIM (Forward+Activation Gradient) array : " << (CAMAreaArray + LUTAreaArray) * 1e12 << "um^2" << endl;
    cout << "Total ADC (or S/As and precharger for SRAM) Area on chip : " << (CAMAreaADC + LUTAreaADC) * 1e12 << "um^2" << endl;
    cout << "Total Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) on chip : " << (CAMAreaAccum + LUTAreaAccum) * 1e12 << "um^2" << endl;
    cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, pooling and activation units) : " << (CAMAreaOther + LUTAreaOther) * 1e12 << "um^2" << endl;
    cout << endl;
    cout << "-----------------------------------Chip layer-by-layer Estimation---------------------------------" << endl;

    cout << "readLatency  is: " << ReadLatency * 1e9 << "ns" << endl;
    cout << "readDynamicEnergy  is: " << ReadDynamicEnergy * 1e12 << "pJ" << endl;
    cout << " ===========>> leakage Energy (Leakage * ReadLatency) is: " << LeakageEnergy * 1e12 << "pJ" << endl;
    cout << " ===========>> leakage Power (Leakage) is: " << Leakage * 1e6 << "uW" << endl;
    cout << endl;
    cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
    cout << endl;
    cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << LatencyADC * 1e9 << "ns" << endl;
    cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << LatencyAccum * 1e9 << "ns" << endl;
    cout << "----------- Synaptic Array w/o ADC (Forward + Activate Gradient) readLatency is : " << LatencyOther * 1e9 << "ns" << endl;
    cout << "----------- Buffer readLatency is: " << BufferLatency * 1e9 << "ns" << endl;
	cout << "----------- Bus readLatency is: " << BusLatency * 1e9 << "ns" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << EnergyADC * 1e12 << "pJ" << endl;
    cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << EnergyAccum * 1e12 << "pJ" << endl;
    cout << "----------- Synaptic Array w/o ADC (Forward + Activate Gradient) readDynamicEnergy is : " << EnergyOther * 1e12 << "pJ" << endl;
    cout << "----------- Buffer readDynamicEnergy is: " << BufferDynamicEnergy * 1e12 << "pJ" << endl;
	cout << "----------- Bus readDynamicEnergy is: " << BusDynamicEnergy * 1e12 << "pJ" << endl;
    cout << endl;
    cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
    cout << endl;
    cout << endl;
    cout << "-----------------------------------Chip layer-by-layer Performance---------------------------------" << endl;

    cout << "Energy Efficiency TOPS/W: " << numComputation / ((ReadDynamicEnergy + LeakageEnergy) * 1e12) << endl;
    cout << "Throughput TOPS: " << numComputation / (ReadLatency) * 1e-12 << endl;
    cout << "Throughput FPS: " << 1 / (ReadLatency) << endl;

    cout << "-------------------------------------- Hardware Performance Done --------------------------------------" << endl;
    cout << endl;
    
    // save results to top level csv file (only total results)
    ofstream outfile;
    outfile.open("NeuroSim_Output.csv");
    if (outfile.is_open())
    {
        outfile << ReadLatency << ",";
        outfile << ReadDynamicEnergy;
        outfile << numComputation / ((ReadDynamicEnergy + LeakageEnergy) * 1e12) << ",";
        outfile << numComputation / (ReadLatency) * 1e-12 << ",";
    }
    else
    {
        cout << "Error: the output file cannot be opened!" << endl;
    }
    outfile.close();

    return 0;
    /*************************************************** CalculatePerformance ***************************************************/
}

    /*************************************************** Function ***************************************************/
vector<vector<double>> getNetStructure(const string &inputfile)
{
    ifstream infile(inputfile.c_str());
    string inputline;
    string inputval;
	cout << inputfile.c_str() << endl;
    int ROWin = 0, COLin = 0;
    if (!infile.good())
    {
        cerr << "Error: the input file cannot be opened!" << endl;
        exit(1);
    }
    else // Get the network's total ROW & COL
    {
        while (getline(infile, inputline, '\n'))
        {
            ROWin++;
        }
        infile.clear();
        infile.seekg(0, ios::beg);
        if (getline(infile, inputline, '\n'))
        {
            istringstream iss(inputline);
            while (getline(iss, inputval, ','))
            {
                COLin++;
            }
        }
    }
    infile.clear();
    infile.seekg(0, ios::beg);

    vector<vector<double>> netStructure;
    for (int row = 0; row < ROWin; row++)
    {
        vector<double> netStructurerow;
        getline(infile, inputline, '\n');
        istringstream iss;
        iss.str(inputline);
        for (int col = 0; col < COLin; col++)
        {
            while (getline(iss, inputval, ','))
            {
                istringstream fs;
                fs.str(inputval);
                double f = 0;
                fs >> f;
                netStructurerow.push_back(f);
            }
        }
        netStructure.push_back(netStructurerow);
    }
    infile.close();

    return netStructure;
    netStructure.clear();
}

vector<vector<double>> LoadInWeightData(const string &weightfile, int numRowPerSynapse, int numColPerSynapse, double maxConductance, double minConductance) {
	
	ifstream fileone(weightfile.c_str());                           
	string lineone;
	string valone;
	
	int ROW = 0;
	int COL = 0;
	
	if (!fileone.good()) {                                       
		cerr << "Error: the fileone cannot be opened!" << endl;
		exit(1);
	}else{
		while (getline(fileone, lineone, '\n')) {                   
			ROW++;                                             
		}
		fileone.clear();
		fileone.seekg(0, ios::beg);                               
		if (getline(fileone, lineone, '\n')) {                      
			istringstream iss (lineone);                         
			while (getline(iss, valone, ',')) {                   
				COL++;
			}
		}	
	}
	fileone.clear();
	fileone.seekg(0, ios::beg);                   
	
	double NormalizedMin = 0;
	double NormalizedMax = pow(2, param->synapseBit);
	
	double RealMax = param->algoWeightMax;
	double RealMin = param->algoWeightMin;
	
	vector<vector<double> > weight;            
	// load the data into a weight matrix ...
	for (int row=0; row<ROW; row++) {	
		vector<double> weightrow;
		vector<double> weightrowb;
		getline(fileone, lineone, '\n');              
		istringstream iss;
		iss.str(lineone);
		for (int col=0; col<COL; col++) {       
			while(getline(iss, valone, ',')){	
				istringstream fs;
				fs.str(valone);
				double f=0;
				fs >> f;
				// training version: linear mapping
				weightrow.push_back((f+1)/2*(maxConductance-minConductance)+minConductance);
			}
		}
		weight.push_back(weightrow);
		weightrow.clear();
	}
	fileone.close();
	
	return weight;
	weight.clear();
}

vector<vector<double>> LoadInInputData(const string &inputfile) {
	
	ifstream infile(inputfile.c_str());     
	string inputline;
	string inputval;
	
	int ROWin=0, COLin=0;      
	if (!infile.good()) {       
		cerr << "Error: the input file cannot be opened!" << endl;
		exit(1);
	}else{
		while (getline(infile, inputline, '\n')) {      
			ROWin++;                               
		}
		infile.clear();
		infile.seekg(0, ios::beg);    
		if (getline(infile, inputline, '\n')) {        
			istringstream iss (inputline);      
			while (getline(iss, inputval, ',')) {       
				COLin++;
			}
		}	
	}
	infile.clear();
	infile.seekg(0, ios::beg);    
	// cout << "ROWin, COLin" << ROWin << " , " << COLin << endl;      
	
	vector<vector<double> > inputvector;              
	// load the data into inputvector ...
	for (int row=0; row<ROWin; row++) {	
		vector<double> inputvectorrow;
		vector<double> inputvectorrowb;
		getline(infile, inputline, '\n');             
		istringstream iss;
		iss.str(inputline);
		for (int col=0; col<COLin; col++) {
			// cout << "at ROW, COL" << row << " , " << col << endl;
			while(getline(iss, inputval, ',')){	
				istringstream fs;
				fs.str(inputval);
				double f=0;
				fs >> f;
				inputvectorrow.push_back(f);
			}
		}
		inputvector.push_back(inputvectorrow);
		inputvectorrow.clear();
	}
	// close the input file ...
	infile.close();
	
	return inputvector;
	inputvector.clear();
}

vector<vector<double>> CopySubArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol) {
	vector<vector<double>> copy;
	for (int i=0; i<numRow; i++) {
		vector<double> copyRow;
		for (int j=0; j<numCol; j++) {
			copyRow.push_back(orginal[positionRow+i][positionCol+j]);
		}
		copy.push_back(copyRow);
		copyRow.clear();
	}
	return copy;
	copy.clear();
} 

vector<vector<double>> CopySubInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow) {
	vector<vector<double> > copy;
	// cout << "nR , nIV" << numRow << " " << numInputVector << endl;
	// cout << "CopySubInput" << orginal[positionRow][0] << endl;
	for (int i=0; i<numInputVector; i++) {
		vector<double> copyRow;
		for (int j=0; j<numRow; j++) {
			// cout << "at i, j" << positionRow+i << " , " << j << endl;
			copyRow.push_back(orginal[positionRow+i][j]);
		}
		copy.push_back(copyRow);
		copyRow.clear();
	}
	// cout << "fin CSI" << endl;
	return copy;
	copy.clear();
}

vector<double> GetInputVector(const vector<vector<double> > &input, int numInput, double *activityRowRead) {
	vector<double> copy;
	for (int i=0; i<input.size(); i++) {
		double x = input[i][numInput];
		copy.push_back(x);   
	}  
	double numofreadrow = 0;  // initialize readrowactivity parameters
	for (int i=0; i<input.size(); i++) {
		if (copy[i] != 0) {
			numofreadrow += 1;
		}else {
			numofreadrow += 0;
		}
	}
	double totalnumRow = input.size();
	*(activityRowRead) = numofreadrow/totalnumRow;
	return copy;
	copy.clear();
} 

vector<double> GetColumnResistance(const vector<double> &input, const vector<vector<double> > &weight, MemCell& cell, bool parallelRead, double resCellAccess) {
	vector<double> resistance;
	vector<double> conductance;
	double columnG = 0;
	for (int j=0; j<weight[0].size(); j++) {
		int activatedRow = 0;
		columnG = 0;
		for (int i=0; i<weight.size(); i++) {
			// eNVM
			double totalWireResistance;
			// if (cell.accessType == CMOS_access) {
			totalWireResistance = (double) 1.0/weight[i][j] + (j + 1) * param->wireResistanceRow + (weight.size() - i) * param->wireResistanceCol + cell.resistanceAccess;
			// } 
			if ((int) input[i] == 1) {
				columnG += (double) 1.0/totalWireResistance;
				activatedRow += 1 ;
			} else {
				columnG += 0;
			}
		}
		
		// if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		if (!parallelRead) {  
			conductance.push_back((double) columnG/activatedRow);
		} else {
			conductance.push_back(columnG);
		}
	}
	// covert conductance to resistance
	for (int i=0; i<weight[0].size(); i++) {
		resistance.push_back((double) 1.0/conductance[i]);
	}
		
	return resistance;
	resistance.clear();
} 

vector<double> GetRowResistance(const vector<double> &input, const vector<vector<double> > &weight, MemCell& cell, bool parallelRead, double resCellAccess) {
	vector<double> resistance;
	vector<double> conductance;
	double rowG = 0; 
	double totalWireResistance;
	
	for (int i=0; i<weight.size(); i++) {
		int activatedCol = ceil(weight[0].size()/2);  // assume 50% of the input vector is 1
		rowG = 0;
		for (int j=0; j<weight[0].size(); j++) {
			// eNVM
			// if (cell.accessType == CMOS_access) {
			totalWireResistance = (double) 1.0/weight[i][j] + (i + 1) * param->wireResistanceRow + (weight[0].size() - j) * param->wireResistanceCol + cell.resistanceAccess;
			// }
		}
		rowG = (double) 1.0/totalWireResistance * activatedCol;
		
		// if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		if (!parallelRead) {  
			conductance.push_back((double) rowG/activatedCol);
		} else {
			conductance.push_back(rowG);
		}
	}
	// covert conductance to resistance
	for (int i=0; i<weight.size(); i++) {
		resistance.push_back((double) 1.0/conductance[i]);
		
	}
		
	return resistance;
	resistance.clear();
} 
    /*************************************************** Function ***************************************************/