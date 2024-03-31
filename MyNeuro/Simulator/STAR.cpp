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
SubArray *subArray;
Bus *busInputCM;
Bus *busOutputCM;
Bus *busToQ;
Bus *busToR;
Bus *busToP;
Bus *busToMT1;
Bus *busToMT2;
Bus *busToMT3;
Bus *busToMT4;
Bus *busToSoft;
DFF *bufferInputCM;
DFF *bufferOutputCM;
DFF *bufferSoft;
DFF *bufferP;
DFF *bufferQ;
DFF *bufferR;

int main(int argc, char *argv[])
{
    gen.seed(0);

    vector<vector<double>> netStructure;
	// cout << "get net "<<argv[2] << endl;
    netStructure = getNetStructure(argv[2]);                                        // get file from trace.command.sh
	// cout << "get bit "<<argv[3] << endl;
    // define weight/input/memory precision from wrapper
    param->synapseBit = atoi(argv[3]);  // precision of synapse weight
    param->numBitInput = atoi(argv[4]); // precision of input neural activation

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
    cout << "User-defined SubArray Size: " << param->numRowSubArray << "x" << param->numColSubArray << endl;
    
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
	subArray 			= new SubArray(inputParameter, tech, cell);
    busInputCM          = new Bus(inputParameter, tech, cell);
    busOutputCM         = new Bus(inputParameter, tech, cell);
    busToQ              = new Bus(inputParameter, tech, cell);
    busToR              = new Bus(inputParameter, tech, cell);
    busToP              = new Bus(inputParameter, tech, cell);
    busToMT1            = new Bus(inputParameter, tech, cell);
    busToMT2            = new Bus(inputParameter, tech, cell);
    busToMT3            = new Bus(inputParameter, tech, cell);
    busToMT4            = new Bus(inputParameter, tech, cell);
    busToSoft           = new Bus(inputParameter, tech, cell);
    bufferInputCM      = new DFF(inputParameter, tech, cell);
    bufferOutputCM     = new DFF(inputParameter, tech, cell);
    bufferSoft         = new DFF(inputParameter, tech, cell);
    bufferP            = new DFF(inputParameter, tech, cell);
    bufferQ            = new DFF(inputParameter, tech, cell);
    bufferR            = new DFF(inputParameter, tech, cell);
	
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

    /*** subArray's property ***/         
	subArray->conventionalParallel = param->conventionalParallel;                  
	subArray->conventionalSequential = param->conventionalSequential;   
	subArray->parallelBP = param->parallelBP;	
	subArray->numRow = param->numRowSubArray;
	subArray->numCol = param->numRowSubArray;
	subArray->levelOutput = param->levelOutput;
	subArray->levelOutputBP = param->levelOutputAG;
	subArray->numColMuxed = param->numColMuxed;               // How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
	subArray->numRowMuxedBP = param->numRowMuxedAG;
    subArray->clkFreq = 20e9;                       // Clock frequency
	subArray->relaxArrayCellHeight = param->relaxArrayCellHeight;
	subArray->relaxArrayCellWidth = param->relaxArrayCellWidth;
	subArray->numReadPulse = param->numBitInput;
	subArray->avgWeightBit = param->cellBit;
	subArray->numCellPerSynapse = param->numColPerSynapse;
	subArray->numReadPulseBP = 8;
	subArray->activityBPColRead = 0.5;
	subArray->SARADC = param->SARADC;
	subArray->currentMode = param->currentMode;
	subArray->spikingMode = NONSPIKING;
	
	int numRow = param->numRowSubArray;
	int numCol = param->numColSubArray;
	
	if (subArray->numColMuxed > numCol) {                      // Set the upperbound of numColMuxed
		subArray->numColMuxed = numCol;
	}

	subArray->numReadCellPerOperationFPGA = numCol;	           // Not relevant for IMEC
	subArray->numWriteCellPerOperationFPGA = numCol;	       // Not relevant for IMEC
	subArray->numReadCellPerOperationMemory = numCol;          // Define # of SRAM read cells in memory mode because SRAM does not have S/A sharing (not relevant for IMEC)
	subArray->numWriteCellPerOperationMemory = numCol/8;       // # of write cells per operation in SRAM memory or the memory mode of multifunctional memory (not relevant for IMEC)
	subArray->numReadCellPerOperationNeuro = numCol;           // # of SRAM read cells in neuromorphic mode
	subArray->numWriteCellPerOperationNeuro = numCol;	       // For SRAM or analog RRAM in neuro mode
    subArray->maxNumWritePulse = MAX(cell.maxNumLevelLTP, cell.maxNumLevelLTD);

	/*** User defined num of elements ***/
	int numSubArrayRow = 1;		// The number of subarray's row
	int numSubArrayCol = 5; 	// The number of subarray's col
	int numbusRow = 1;			// The number of bus's in row
	int numbusCol = 1; 			// The number of bus's in col
	/*** initialize modules ***/
	subArray->Initialize(numRow, numCol, param->unitLengthWireResistance);        // initialize subArray

	// For Buffer
	bufferInputCM->Initialize(param->numBitInput*numRow, 20e9);
	bufferP->Initialize(param->numBitInput*numRow, 20e9);
	bufferQ->Initialize(param->numBitInput*numRow, 20e9);
    bufferR->Initialize(param->numBitInput*numRow, 20e9);
    bufferSoft->Initialize(param->numBitInput*numRow, 20e9);
	if (param->parallelRead) {
		bufferOutputCM->Initialize((numCol/param->numColMuxed)*(log2((double)param->levelOutput)+param->numBitInput+param->numColPerSynapse), 20e9);
	} else {
		bufferOutputCM->Initialize((numCol/param->numColMuxed)*((log2((double)numRow)+param->cellBit-1)+param->numBitInput+param->numColPerSynapse), 20e9);
	}
	// For Bus
	busInputCM->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
	busToMT1->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
	busToMT2->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);    
	busToMT3->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
	busToMT4->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
	busToSoft->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
	busToQ->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
	busToR->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
    busToP->Initialize(HORIZONTAL, numbusRow, numbusCol, 0, numRow, subArray->height, subArray->width);
	busOutputCM->Initialize(VERTICAL, numbusRow, numbusCol, 0, numCol, subArray->height, subArray->width);

    cout << " ----------- number of subarray's row is " << numSubArrayRow << " with number of subarray's col is " << numSubArrayCol << endl;
	cout << " ----------- number of bus's row is " << numbusRow << " with number of bus's col is " << numbusCol << endl;
    cout << endl;
    cout << "---------------------------- FloorPlan Done ------------------------------" << endl;
    cout << endl;

    /******************************************************** Initialize ********************************************************/

    /****************************************************** CalculateArea *******************************************************/
    double Area, AreaArray, AreaADC, AreaAccum, AreaOther, OverallArea;
	double height = 0;
	double width = 0;
	double busarea = 0;
	double bufferarea = 0;
	double Totalarea = 0;
	double widthArray = 0;
	vector<double> areaResults;
	
	subArray->CalculateArea();
	//buffer
	bufferInputCM->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferP->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferQ->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferR->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferSoft->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferOutputCM->CalculateArea(NULL, numCol*subArray->width, NONE);
	//bus
    busInputCM->CalculateArea(1, true); 
    busOutputCM->CalculateArea(1, true);
	busToMT1->CalculateArea(1, true); 
    busToMT2->CalculateArea(1, true);
	busToMT3->CalculateArea(1, true); 
    busToMT4->CalculateArea(1, true);
	busToSoft->CalculateArea(1, true); 
    busToP->CalculateArea(1, true);
	busToQ->CalculateArea(1, true);
	busToR->CalculateArea(1, true);
	
	busarea = busInputCM->area + busOutputCM->area + 
			  busToMT1->area + busToMT2->area + busToMT3->area + busToMT4->area + 
			  busToP->area + busToQ->area + busToSoft->area + busToR->area;
	bufferarea = bufferInputCM->area + bufferOutputCM->area + 
				 bufferP->area + bufferQ->area + bufferR->area +
				 bufferSoft->area;
	
	// cout << " ----------- Buffer Area -----------" << endl;
	// cout << " ----------- bufferInputCM is : " << bufferInputCM->area * 1e12 << endl;
	// cout << " ----------- bufferOutputCM is : " << bufferOutputCM->area * 1e12 << endl;
	// cout << " ----------- bufferP is : " << bufferP->area * 1e12 << endl;
	// cout << " ----------- bufferQ is : " << bufferQ->area * 1e12  << endl;
	// cout << " ----------- bufferR is : " << bufferR->area * 1e12 << endl;
	// cout << " ----------- bufferSoft is : " << bufferSoft->area  * 1e12 << endl;
	

    Totalarea = subArray->usedArea * (numSubArrayRow*numSubArrayCol) + bufferarea + busarea;
	cout << " ----------- The subArray->usedArea is : " << subArray->usedArea  * 1e12 << endl;
	cout << " ----------- The subArray Total num is : " << numSubArrayRow*numSubArrayCol << endl;
	cout << " ----------- The bufferarea is : \t " << bufferarea  * 1e12 << endl;
	cout << " ----------- The busarea is : \t\t " << busarea  * 1e12 << endl;
    
    height = sqrt(Totalarea);
    width = Totalarea/(height);
    
    areaResults.push_back(Totalarea);
    areaResults.push_back(subArray->areaADC*(numRow*numCol));
    areaResults.push_back(subArray->areaAccum*(numRow*numCol));
    areaResults.push_back(subArray->areaOther*(numRow*numCol)+ bufferInputCM->area + bufferOutputCM->area);
    areaResults.push_back(subArray->areaArray*(numRow*numCol));
	
    Area        = areaResults[0];
    AreaArray   = areaResults[1];
    AreaADC     = areaResults[2];
    AreaAccum   = areaResults[3];
    AreaOther   = areaResults[4];
	for (double area : areaResults) {
        OverallArea += area;
    }

    cout << "-------------------------------------- Hardware Performance --------------------------------------" << endl;
	
    // save breakdown results of each layer to csv files
    ofstream breakdownfile;
    string breakdownfile_name = "./NeuroSim_Breakdown";
    breakdownfile_name.append(".csv");
    breakdownfile.open(breakdownfile_name);
    if (breakdownfile.is_open()){
        // firstly save the area results to file
        breakdownfile << "Total Area(m^2), ADC Area(m^2), Accumulation Area(m^2), Other Logic&Storage Area(m^2), Total CIM (FW+AG) Area (m^2)" << endl;
        breakdownfile << Area << "," << AreaADC << "," << AreaAccum << "," << AreaOther << "," << AreaArray << endl;
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
	double SWReadLatency, SWReadDynamicEnergy, SWLeakage, SWLatencyADC, SWLatencyAccum, SWLatencyOther, SWEnergyADC, SWEnergyAccum, SWEnergyOther;
	double STReadLatency, STReadDynamicEnergy, STLeakage, STLatencyADC, STLatencyAccum, STLatencyOther, STEnergyADC, STEnergyAccum, STEnergyOther;
	double SSReadLatency, SSReadDynamicEnergy, SSLeakage, SSLatencyADC, SSLatencyAccum, SSLatencyOther, SSEnergyADC, SSEnergyAccum, SSEnergyOther;
	double ReadLatency, ReadDynamicEnergy, Leakage, LatencyADC, LatencyOther, LatencyAccum, EnergyADC, EnergyAccum, EnergyOther;
	double BufferLatency, LeakageEnergy, BufferDynamicEnergy, BusLatency, BusDynamicEnergy;

	// Segment * Weight's P&L  
	SWReadLatency		=0;
	SWReadDynamicEnergy	=0;
	SWLeakage			=0;
	SWLatencyADC		=0;
	SWLatencyAccum		=0;
	SWLatencyOther		=0;
	SWEnergyADC			=0;
	SWEnergyAccum		=0;
	SWEnergyOther		=0;
	// Segment * Transepose's P&L
    STReadLatency 		= 0;
	STReadDynamicEnergy = 0;
	STLeakage 			= 0;
	STEnergyADC 		= 0;
	STEnergyAccum 		= 0;
	STEnergyOther 		= 0;
	STLatencyADC 		= 0;
	STLatencyAccum 		= 0;
	STLatencyOther 		= 0;
	// Soft * Segment's P&L
    SSReadLatency 		= 0;
	SSReadDynamicEnergy = 0;
	SSLeakage 			= 0;
	SSEnergyADC 		= 0;
	SSEnergyAccum 		= 0;
	SSEnergyOther 		= 0;
	SSLatencyADC 		= 0;
	SSLatencyAccum 		= 0;
	SSLatencyOther 		= 0;
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

	/*** get weight matrix file Size ***/
    // I = W*L*D  = 1*16*128  = netStructure[0][0]*netStructure[0][1]*netStructure[0][2]
    // K = K*K'*D = 1*1*128 = netStructure[0][3]*netStructure[0][4]*netStructure[0][5]
	int weightMatrixRow = netStructure[0][2]*netStructure[0][3]*netStructure[0][4];
	int weightMatrixCol = netStructure[0][5];

    // weight matrix is further partitioned inside PE (among subArray) --> no duplicated
    int numRowMatrix = min(param->numRowSubArray, weightMatrixRow);
    int numColMatrix = min(param->numColSubArray, weightMatrixCol);
	int numInVector = (netStructure[0][0]-netStructure[0][3]+1)/netStructure[0][7]*(netStructure[0][1]-netStructure[0][4]+1)/netStructure[0][7];
	
	// load in whole file 
	vector<vector<double>> weightMemory, transMemory, inputMemory;
	vector<vector<double>> inputVector, softVector;
	
	weightMemory = LoadInWeightData(argv[5], 1, 1, param->maxConductance, param->minConductance);
	transMemory = LoadInWeightData(argv[6], 1, 1, param->maxConductance, param->minConductance);
	inputMemory = LoadInWeightData(argv[7], 1, 1, param->maxConductance, param->minConductance);
	inputVector = LoadInInputData(argv[8]);
	softVector = LoadInInputData(argv[9]);

    /*** assign weight and input to specific subArray ***/
    vector<vector<double>> subArrayWeight,subArrayTrans,subArrayIn;
	vector<vector<double>> subArrayInput, subArraySoft;

    subArrayWeight 	= CopySubArray(weightMemory, 0, 0, numRowMatrix, numColMatrix);	//128*128
	subArrayTrans 	= CopySubArray(transMemory, 0, 0, numRowMatrix, numInVector);	//128*16
	subArrayIn  	= CopySubArray(inputMemory, 0, 0, numInVector, numRowMatrix);	//16*128
	subArrayInput 	= CopySubInput(inputVector, 0, numInVector, numRowMatrix);		//16*128
	subArraySoft 	= CopySubInput(softVector, 0, numInVector, numInVector);		//16*16
	// cout << "CopySubInput : " << endl;
	
    for (int k=0; k<numInVector; k++) {                 // calculate WeightsubArray through the total Segment vectors
        // cout << "k is : " << k << " with numInVector " << numInVector << endl;
		double activityRowRead = 0;
        vector<double> input;
        input = GetInputVector(subArrayInput, k, &activityRowRead);
		subArray->activityRowRead = activityRowRead;
        int cellRange = pow(2, param->cellBit);
        if (param->parallelRead) {
            subArray->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
        } else {
            subArray->levelOutput = cellRange;
        }
        
        vector<double> columnResistance, rowResistance;
		// cout << "columnResistance" << endl;
        columnResistance = GetColumnResistance(input, subArrayWeight, cell, param->parallelRead, subArray->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subArrayWeight, cell, param->parallelBP, subArray->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        subArray->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		subArray->CalculatePower(columnResistance, rowResistance);
        
        SWReadLatency += subArray->readLatency;
        SWReadDynamicEnergy += subArray->readDynamicEnergy;
        SWLeakage += subArray->leakage;
        
        SWLatencyADC += subArray->readLatencyADC;
        SWLatencyAccum += subArray->readLatencyAccum;
        SWLatencyOther += subArray->readLatencyOther;
        
        SWEnergyADC += subArray->readDynamicEnergyADC;
        SWEnergyAccum += subArray->readDynamicEnergyAccum;
        SWEnergyOther += subArray->readDynamicEnergyOther;
    }

	
    for (int k=0; k<numInVector; k++) {                 // calculate transposesubArray through the total segment vectors
        // cout << "k is : " << k << " with numInVector " << numInVector << endl;
		double activityRowRead = 0;
        vector<double> input;
        input = GetInputVector(subArrayInput, k, &activityRowRead);
		subArray->activityRowRead = activityRowRead;
        int cellRange = pow(2, param->cellBit);
        if (param->parallelRead) {
            subArray->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
        } else {
            subArray->levelOutput = cellRange;
        }
        
        vector<double> columnResistance, rowResistance;
		// cout << "columnResistance" << endl;
        columnResistance = GetColumnResistance(input, subArrayTrans, cell, param->parallelRead, subArray->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subArrayTrans, cell, param->parallelBP, subArray->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        subArray->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		subArray->CalculatePower(columnResistance, rowResistance);
        
        STReadLatency += subArray->readLatency;
        STReadDynamicEnergy += subArray->readDynamicEnergy;
        STLeakage += subArray->leakage;
        
        STLatencyADC += subArray->readLatencyADC;
        STLatencyAccum += subArray->readLatencyAccum;
        STLatencyOther += subArray->readLatencyOther;
        
        STEnergyADC += subArray->readDynamicEnergyADC;
        STEnergyAccum += subArray->readDynamicEnergyAccum;
        STEnergyOther += subArray->readDynamicEnergyOther;
    }

	for (int k=0; k<numInVector; k++) {                 // calculate segmentsubArray through the total softmax vectors
        // cout << "k is : " << k << " with numInVector " << numInVector << endl;
		double activityRowRead = 0;
        vector<double> input;
        input = GetInputVector(subArraySoft, k, &activityRowRead);
		subArray->activityRowRead = activityRowRead;
        int cellRange = pow(2, param->cellBit);
        if (param->parallelRead) {
            subArray->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
        } else {
            subArray->levelOutput = cellRange;
        }
        
        vector<double> columnResistance, rowResistance;
		// cout << "columnResistance" << endl;
        columnResistance = GetColumnResistance(input, subArrayIn, cell, param->parallelRead, subArray->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subArrayIn, cell, param->parallelBP, subArray->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        subArray->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		subArray->CalculatePower(columnResistance, rowResistance);
        
        SSReadLatency += subArray->readLatency;
        SSReadDynamicEnergy += subArray->readDynamicEnergy;
        SSLeakage += subArray->leakage;
        
        SSLatencyADC += subArray->readLatencyADC;
        SSLatencyAccum += subArray->readLatencyAccum;
        SSLatencyOther += subArray->readLatencyOther;
        
        SSEnergyADC += subArray->readDynamicEnergyADC;
        SSEnergyAccum += subArray->readDynamicEnergyAccum;
        SSEnergyOther += subArray->readDynamicEnergyOther;
    }

	//considering buffer activation: no matter speedup or not, the total number of data transferred is fixed
	// input buffer: total num of data loaded in = weightMatrixRow*numInVector
	// output buffer: total num of data transferred = weightMatrixRow*numInVector/param->numBitInput (total num of IFM in the PE) *adderTree->numAdderTree*adderTree->numAdderBit (bit precision of OFMs) 
    bufferInputCM->CalculateLatency(0, numInVector*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
    bufferOutputCM->CalculateLatency(0, numInVector/param->numBitInput);
	bufferP->CalculateLatency(0, numInVector*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
	bufferQ->CalculateLatency(0, numInVector*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
	bufferR->CalculateLatency(0, numInVector*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
	bufferSoft->CalculateLatency(0, numInVector*ceil((double) weightMatrixRow/(double) param->numRowSubArray));

    bufferInputCM->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInVector);
    bufferOutputCM->CalculatePower(weightMatrixCol/param->numColPerSynapse, numInVector/param->numBitInput);
	bufferP->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInVector);
	bufferQ->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInVector);
	bufferR->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInVector);
	bufferSoft->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInVector);
    
    busInputCM->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth));
	busToMT1->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT1->busWidth));
	busToMT2->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT2->busWidth));
	busToMT3->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT3->busWidth));
	busToMT4->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT4->busWidth));
	busToP->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToP->busWidth));
	busToQ->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToQ->busWidth));
	busToR->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToR->busWidth));
	busToSoft->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busToSoft->busWidth));

    busInputCM->CalculatePower(busInputCM->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth));
	busToMT1->CalculatePower(busToMT1->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT1->busWidth));
	busToMT2->CalculatePower(busToMT2->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT2->busWidth));
	busToMT3->CalculatePower(busToMT1->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT3->busWidth));
	busToMT4->CalculatePower(busToMT2->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToMT4->busWidth));
	busToP->CalculatePower(busToP->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToP->busWidth));
	busToQ->CalculatePower(busToQ->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToQ->busWidth));
	busToR->CalculatePower(busToR->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToR->busWidth));
	busToSoft->CalculatePower(busToSoft->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busToSoft->busWidth));
    
    if (param->parallelRead) {
        busOutputCM->CalculateLatency((weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
        busOutputCM->CalculatePower(busOutputCM->numRow*busOutputCM->busWidth, (weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
	} else {
        busOutputCM->CalculateLatency((weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
        busOutputCM->CalculatePower(busOutputCM->numRow*busOutputCM->busWidth, (weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
    }

    Leakage = 3*SWLeakage + STLeakage +	SSLeakage + //Weight Q,K,V & transpose I & I
			  bufferInputCM->leakage + bufferOutputCM->leakage + bufferSoft->leakage +	//Buffer
			  bufferP->leakage + bufferQ->leakage + bufferR->leakage;
	cout << "~~~~~~~~ !!!!!!!! HERE !!!!!!!! ~~~~~~~~ " << endl;
	cout << "The IW is " << STLeakage << endl;
	cout << "The ST is" << STLeakage << endl;
	cout << "The SS is " << SSLeakage << endl;
	cout << "The bufferInputCM->leakage is " << bufferInputCM->leakage << endl;
	cout << "The bufferOutputCM->leakage is " << bufferOutputCM->leakage << endl;
	cout << "The bufferSoft->leakage is " << bufferSoft->leakage << endl;
	cout << "The bufferP->leakage is " << bufferP->leakage << endl;
	cout << "The bufferQ->leakage is " << bufferQ->leakage << endl;
	cout << "The bufferQ->leakage is " << bufferR->leakage << endl;


	LatencyADC = 3*SWLatencyADC + STLatencyADC + SSLatencyADC;
	LatencyOther = 3*SWLatencyOther + STLatencyOther + SSLatencyOther;
	LatencyAccum = 3*SWLatencyAccum + STLatencyAccum + SSLatencyAccum;
	BufferLatency =  bufferP -> readLatency + bufferQ->readLatency + bufferR->readLatency +
					 bufferInputCM->readLatency + bufferOutputCM->readLatency + bufferSoft->readLatency;
    BufferDynamicEnergy = bufferP -> readDynamicEnergy + bufferQ->readDynamicEnergy +	bufferR->readDynamicEnergy +
						  bufferInputCM->readDynamicEnergy + bufferOutputCM->readDynamicEnergy + bufferSoft->readDynamicEnergy;
	BusLatency = busInputCM->readLatency + busOutputCM->readLatency +					//bus
				 busToMT1->readLatency + busToMT2->readLatency + busToMT3->readLatency + busToMT4->readLatency +
				 busToP->readLatency + busToQ->readLatency + busToR->readLatency + busToSoft->readLatency;
	BusDynamicEnergy = busInputCM->readDynamicEnergy + busOutputCM->readDynamicEnergy +		//bus
				  		busToMT1->readDynamicEnergy + busToMT2->readDynamicEnergy + busToMT3->readDynamicEnergy + busToMT4->readDynamicEnergy+
				  		busToP->readDynamicEnergy + busToQ->readDynamicEnergy + busToR->readLatency + busToSoft->readDynamicEnergy;

	EnergyADC = SWEnergyADC + STEnergyADC + SSEnergyADC;
	EnergyOther = SWEnergyOther + STEnergyOther + SSEnergyOther;
	EnergyAccum = SWEnergyAccum + STEnergyAccum + SSEnergyAccum;
    ReadLatency = 3*SWReadLatency + STReadLatency + SSReadLatency + BufferLatency + BusLatency;
    ReadDynamicEnergy = 3*SWReadDynamicEnergy + STReadDynamicEnergy + SSReadDynamicEnergy + BufferDynamicEnergy + BusDynamicEnergy;
    LeakageEnergy = Leakage * ReadLatency;

    cout << "------------------------------ Summary --------------------------------" << endl;
    cout << endl;
    cout << " ===========>> Area : " << Area * 1e12 << "um^2" << endl;
    cout << " ----------- Total CIM (Forward+Activation Gradient) array : " << AreaArray * 1e12 << "um^2" << endl;
    cout << " ----------- Total ADC (or S/As and precharger for SRAM) Area on chip : " << AreaADC * 1e12 << "um^2" << endl;
    cout << " ----------- Total Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) on chip : " << AreaAccum * 1e12 << "um^2" << endl;
    cout << " ----------- Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, pooling and activation units) : " << AreaOther * 1e12 << "um^2" << endl;
    cout << " ----------- OverallArea : " << OverallArea * 1e12 << "um^2" << endl;
	cout << endl;
    cout << "-----------------------------------Chip layer-by-layer Estimation---------------------------------" << endl;

    cout << " ----------- readLatency  is: " << ReadLatency * 1e9 << "ns" << endl;
    cout << " ----------- readDynamicEnergy  is: " << ReadDynamicEnergy * 1e12 << "pJ" << endl;
    cout << " ===========>> leakage Energy (Leakage * ReadLatency) is: " << LeakageEnergy * 1e12 << "pJ" << endl;
    cout << " ===========>> leakage Power (Leakage) is: " << Leakage * 1e6 << "uW" << endl;
    cout << endl;
    cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
    cout << endl;
    cout << " ----------- ADC (or S/As and precharger for SRAM) readLatency is : " << LatencyADC * 1e9 << "ns" << endl;
    cout << " ----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << LatencyAccum * 1e9 << "ns" << endl;
    cout << " ----------- Synaptic Array w/o ADC (Forward + Activate Gradient) readLatency is : " << LatencyOther * 1e9 << "ns" << endl;
    cout << " ----------- Buffer readLatency is: " << BufferLatency * 1e9 << "ns" << endl;
	cout << " ----------- Bus readLatency is: " << BusLatency * 1e9 << "ns" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << " ----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << EnergyADC * 1e12 << "pJ" << endl;
    cout << " ----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << EnergyAccum * 1e12 << "pJ" << endl;
    cout << " ----------- Synaptic Array w/o ADC (Forward + Activate Gradient) readDynamicEnergy is : " << EnergyOther * 1e12 << "pJ" << endl;
    cout << " ----------- Buffer readDynamicEnergy is: " << BufferDynamicEnergy * 1e12 << "pJ" << endl;
	cout << " ----------- Bus readDynamicEnergy is: " << BusDynamicEnergy * 1e12 << "pJ" << endl;
    cout << endl;
    cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
    cout << endl;
    cout << endl;
    cout << "-----------------------------------Chip layer-by-layer Performance---------------------------------" << endl;

    cout << " ----------- Energy Efficiency TOPS/W: " << numComputation / ((ReadDynamicEnergy + LeakageEnergy) * 1e12) << endl;
    cout << " ----------- Throughput TOPS: " << numComputation / (ReadLatency) * 1e-12 << endl;
    cout << " ----------- Throughput FPS: " << 1 / (ReadLatency) << endl;

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