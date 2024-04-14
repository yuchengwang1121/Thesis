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
#include <iomanip>
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
AdderTree *adderTreeCM;
Bus *busInputCM;
Bus *busOutputCM;
Bus *busToReg1;
Bus *busToMT1;
Bus *busToReg2;
Bus *busToMT2;
Bus *busToSoft;
Bus *busSoftToInput;
Bus *busRingBroadcast;
DFF *bufferInputCM;
DFF *bufferOutputCM;
DFF *bufferSoft;
DFF *bufferReg1;
DFF *bufferReg2;
Mux *mux;

int main(int argc, char *argv[])
{
    gen.seed(0);

    vector<vector<double>> netStructure;
	// cout << "get net "<<argv[2] << endl;
	int Segnum = atoi(argv[2]);
	double clk = 11e9;
    netStructure = getNetStructure(argv[3]);                                        // get file from trace.command.sh
	// cout << "get bit "<<argv[3] << endl;
    // define weight/input/memory precision from wrapper
    param->synapseBit = atoi(argv[4]);  // precision of synapse weight
    param->numBitInput = atoi(argv[5]); // precision of input neural activation

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
	adderTreeCM 		= new AdderTree(inputParameter, tech, cell);
	busInputCM 			= new Bus(inputParameter, tech, cell);
	busOutputCM 		= new Bus(inputParameter, tech, cell);
	busToReg1			= new Bus(inputParameter, tech, cell);
	busToMT1			= new Bus(inputParameter, tech, cell);
	busToReg2			= new Bus(inputParameter, tech, cell);
	busToMT2			= new Bus(inputParameter, tech, cell);
	busToSoft			= new Bus(inputParameter, tech, cell);
	busSoftToInput		= new Bus(inputParameter, tech, cell);
	busRingBroadcast	= new Bus(inputParameter, tech, cell);
	bufferInputCM 		= new DFF(inputParameter, tech, cell);
	bufferOutputCM 		= new DFF(inputParameter, tech, cell);
	bufferReg1 			= new DFF(inputParameter, tech, cell);
	bufferReg2	 		= new DFF(inputParameter, tech, cell);
	bufferSoft 			= new DFF(inputParameter, tech, cell);
	mux 				= new Mux(inputParameter, tech, cell);
	
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
    subArray->clkFreq = clk;                       // Clock frequency
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
	int numSubArrayRow = Segnum;		// The number of subarray's row
	int numSubArrayCol = 1; 			// The number of subarray's col
	int nummuxin = 1;					// The number of output before mux
	int nummuxout = 2;					// The number of output after mux
	double resTg = cell.resistanceOn * IR_DROP_TOLERANCE + cell.resistanceOn;
	
	/*** initialize modules ***/
	subArray->Initialize(numRow, numCol, param->unitLengthWireResistance);        // initialize subArray
	if (param->parallelRead) {
		adderTreeCM->Initialize(Segnum, log2((double)param->levelOutput)+param->numBitInput+param->numColPerSynapse+1, ceil((double)numCol/(double)param->numColMuxed));
	} else {
		adderTreeCM->Initialize(Segnum, (log2((double)numRow)+param->cellBit-1)+param->numBitInput+param->numColPerSynapse+1, ceil((double)numCol/(double)param->numColMuxed));
	}
	
	// For Buffer
	bufferInputCM->Initialize(param->numBitInput*numRow, clk);
	bufferReg1->Initialize(param->numBitInput*numRow, clk);
	bufferReg2->Initialize(param->numBitInput*numRow/Segnum, clk);
	bufferSoft->Initialize(param->numBitInput*numRow, clk);
	if (param->parallelRead) {
		bufferOutputCM->Initialize((numCol/param->numColMuxed)*(log2((double)param->levelOutput)+param->numBitInput+param->numColPerSynapse+adderTreeCM->numStage), clk);
	} else {
		bufferOutputCM->Initialize((numCol/param->numColMuxed)*((log2((double)numRow)+param->cellBit-1)+param->numBitInput+param->numColPerSynapse+adderTreeCM->numStage), clk);
	}
	
	// For Bus
	busInputCM->Initialize(HORIZONTAL, 1, 1, 0, numRow, subArray->height, subArray->width);
	busToReg1->Initialize(HORIZONTAL, 1, 1, 0, numRow, subArray->height, subArray->width);
	busToMT1->Initialize(HORIZONTAL, 1, 1, 0, numRow, subArray->height, subArray->width);
	busToReg2->Initialize(HORIZONTAL, Segnum, 1, 0, numRow/Segnum, subArray->height, subArray->width);
	busToMT2->Initialize(HORIZONTAL, Segnum, 1, 0, numRow/Segnum, subArray->height, subArray->width);
	busToSoft->Initialize(HORIZONTAL, Segnum, 1, 0, numRow/Segnum, subArray->height, subArray->width);
	busSoftToInput->Initialize(HORIZONTAL, 1, 1, 0, numRow, subArray->height, subArray->width);
	busRingBroadcast->Initialize(VERTICAL, 1, 1, 0, numCol, subArray->height, subArray->width);
	busOutputCM->Initialize(VERTICAL, 1, 1, 0, numCol, subArray->height, subArray->width);

	//For mux
	mux->Initialize(nummuxout, nummuxin, resTg, 0);


    cout << " ----------- Ring : number of subarray's row is " << numSubArrayRow << " with number of subarray's col is " << numSubArrayCol << endl;
	cout << " ----------- number of bus's row is " << Segnum << " with number of bus's col is " << 1 << endl;
    cout << endl;
    cout << "---------------------------- FloorPlan Done ------------------------------" << endl;
	cout << endl;
    /******************************************************** Initialize ********************************************************/

    /****************************************************** CalculateArea *******************************************************/
    double Area, AreaArray, AreaADC, AreaAccum, AreaOther,OverallArea;
	double height = 0;
	double width = 0;
	double busarea = 0;
	double bufferarea = 0;
	double Totalarea = 0;
	double widthArray = 0;
	double totalnumofsubArray = 0;
	vector<double> areaResults;
	
	subArray->CalculateArea();
	adderTreeCM->CalculateArea(NULL, subArray->width, NONE);
	//buffer
	bufferInputCM->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferReg1->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferReg2->CalculateArea((numRow*subArray->height)/Segnum, NULL, NONE);
	bufferSoft->CalculateArea(numRow*subArray->height, NULL, NONE);
	bufferOutputCM->CalculateArea(NULL, numCol*subArray->width, NONE);
	//bus
    busInputCM->CalculateArea(1, true); 
    busOutputCM->CalculateArea(1, true);
	busToMT1->CalculateArea(1, true); 
    busToMT2->CalculateArea(1, true);
	busToReg1->CalculateArea(1, true); 
    busToReg2->CalculateArea(1, true);
	busToSoft->CalculateArea(1, true); 
    busSoftToInput->CalculateArea(1, true);
	busRingBroadcast->CalculateArea(1, true);
	//mux
	widthArray = (double)numCol * cell.widthInFeatureSize * tech.featureSize;
	totalnumofsubArray = numSubArrayRow*numSubArrayCol + 2;
	mux->CalculateArea(NULL, widthArray, NONE);
	
	busarea = busInputCM->area + busOutputCM->area + busToMT1->area + busToMT2->area + 
			  busToReg1->area + busToReg2->area + busToSoft->area + busSoftToInput->area + busRingBroadcast->area;
	bufferarea = bufferInputCM->area + bufferOutputCM->area + bufferReg1->area +
				 numSubArrayRow *bufferReg2->area + bufferSoft->area;

	// cout << " ----------- Buffer Area -----------" << endl;
	// cout << " ----------- bufferInputCM is : " << bufferInputCM->area * 1e12 << endl;
	// cout << " ----------- bufferOutputCM is : " << bufferOutputCM->area * 1e12 << endl;
	// cout << " ----------- bufferReg1 is : " << bufferReg1->area * 1e12 << endl;
	// cout << " ----------- bufferReg2 is : " << bufferReg2->area * 1e12  << endl;
	// cout << " ----------- bufferSoft is : " << bufferSoft->area * 1e12 << endl;
	// cout << " ----------- numSubArrayRow is : " << numSubArrayRow << endl;

    Totalarea = subArray->usedArea * totalnumofsubArray + adderTreeCM->area + bufferarea + busarea + numSubArrayRow * mux->area;
    cout << " ----------- The subArray->usedArea is : " << subArray->usedArea  * 1e12 << endl;
	cout << " ----------- The subArray Total num is : " << totalnumofsubArray << endl;
	cout << " ----------- The bufferarea is : \t " << bufferarea  * 1e12 << endl;
	cout << " ----------- The busarea is : \t\t " << busarea  * 1e12 << endl;
	cout << " ----------- numSubArrayRow * mux->area : " << numSubArrayRow * mux->area  * 1e12 << endl;
	cout << endl;


    height = sqrt(Totalarea);
    width = Totalarea/(height);
    
    areaResults.push_back(Totalarea);
    areaResults.push_back(subArray->areaADC*(numRow*numCol));
    areaResults.push_back(subArray->areaAccum*(numRow*numCol)+adderTreeCM->area);
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
	double IWReadLatency, IWReadDynamicEnergy, IWLeakage, IWLatencyADC, IWLatencyAccum, IWLatencyOther, IWEnergyADC, IWEnergyAccum, IWEnergyOther;
	double STReadLatency, STReadDynamicEnergy, STLeakage, STLatencyADC, STLatencyAccum, STLatencyOther, STEnergyADC, STEnergyAccum, STEnergyOther;
	double SSReadLatency, SSReadDynamicEnergy, SSLeakage, SSLatencyADC, SSLatencyAccum, SSLatencyOther, SSEnergyADC, SSEnergyAccum, SSEnergyOther;
	double ReadLatency, ReadDynamicEnergy, Leakage, LatencyADC, LatencyOther, LatencyAccum, EnergyADC, EnergyAccum, EnergyOther;
	double BufferLatency, LeakageEnergy, BufferDynamicEnergy, BusLatency, BusDynamicEnergy;

	// Input * Weight's P&L  
	IWReadLatency		=0;
	IWReadDynamicEnergy	=0;
	IWLeakage			=0;
	IWLatencyADC		=0;
	IWLatencyOther		=0;
	IWEnergyADC			=0;
	IWEnergyAccum		=0;
	IWEnergyOther		=0;
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
	BusLatency			= 0;
	BusDynamicEnergy	= 0;
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
	int numInSeg = numInVector/Segnum;
	
	// load in whole file 
	vector<vector<double>> weightMemory, transMemory, inputMemory;
	vector<vector<double>> inputVector, softVector, segVector;
	weightMemory = LoadInWeightData(argv[6], 1, 1, param->maxConductance, param->minConductance);
	transMemory  = LoadInWeightData(argv[7], 1, 1, param->maxConductance, param->minConductance);
	inputMemory  = LoadInWeightData(argv[8], 1, 1, param->maxConductance, param->minConductance);
	inputVector  = LoadInInputData(argv[8]);
	segVector    = LoadInInputData(argv[9]);
	softVector   = LoadInInputData(argv[10]);

    /*** assign weight and input to specific subArray ***/
    vector<vector<double>> subArrayWeight,subArrayTransWeight,subArrayInWeight;
	vector<vector<double>> subArraySeg, subArraySoft, subArrayInput;

    subArrayWeight 		= CopySubArray(weightMemory, 0, 0, numRowMatrix, numColMatrix);	//128*128
	subArrayTransWeight = CopySubArray(transMemory, 0, 0, numRowMatrix, numInSeg);		//128*4
	subArrayInWeight  	= CopySubArray(inputMemory, 0, 0, numInVector, numRowMatrix);	//16*128
	subArrayInput 		= CopySubInput(inputVector, 0, numInVector, numRowMatrix);		//16*128
	subArraySeg 		= CopySubInput(segVector, 0, numInSeg, numRowMatrix);			//4*128
	subArraySoft 		= CopySubInput(softVector, 0, numInSeg, numInSeg);				//4*4
	// cout << "subArrayWeight : " << subArrayWeight.size() << endl;
	// cout << "subArrayTransWeight : " << subArrayTransWeight.size() << endl;
	// cout << "subArrayInWeight : " << subArrayInWeight.size() << endl;
	// cout << "subArrayInputt : " << subArrayInput.size() << endl;
	// cout << "subArraySeg : " << subArraySeg.size() << endl;
	// cout << "subArraySoft : " << subArraySoft.size() << endl;
	
	for (int k=0; k<numInVector; k++) {                 // Calculate Q,k 
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
        
        IWReadLatency += subArray->readLatency;
        IWReadDynamicEnergy += subArray->readDynamicEnergy;
        IWLeakage += subArray->leakage;
		
        
        IWLatencyADC += subArray->readLatencyADC;
        IWLatencyAccum += subArray->readLatencyAccum;
        IWLatencyOther += subArray->readLatencyOther;
        
        IWEnergyADC += subArray->readDynamicEnergyADC;
        IWEnergyAccum += subArray->readDynamicEnergyAccum;
        IWEnergyOther += subArray->readDynamicEnergyOther;
    }

    for (int k=0; k<numInSeg; k++) {                 // calculate transposesubArray through the total segment vectors
        // cout << "k is : " << k << " with numInVector " << numInVector << endl;
		double activityRowRead = 0;
        vector<double> input;
        input = GetInputVector(subArraySeg, k, &activityRowRead);
		subArray->activityRowRead = activityRowRead;
        int cellRange = pow(2, param->cellBit);
        if (param->parallelRead) {
            subArray->levelOutput = param->levelOutput;               // # of levels of the multilevelSenseAmp output
        } else {
            subArray->levelOutput = cellRange;
        }
        
        vector<double> columnResistance, rowResistance;
		// cout << "columnResistance" << endl;
        columnResistance = GetColumnResistance(input, subArrayTransWeight, cell, param->parallelRead, subArray->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subArrayTransWeight, cell, param->parallelBP, subArray->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        subArray->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		subArray->CalculatePower(columnResistance, rowResistance);
        
        STReadLatency += subArray->readLatency;
        STReadDynamicEnergy += subArray->readDynamicEnergy;
        STLeakage += subArray->leakage;
		// cout << "at i "<< k <<"STLeakage is"<<STLeakage<<endl;
        
        STLatencyADC += subArray->readLatencyADC;
        STLatencyAccum += subArray->readLatencyAccum;
        STLatencyOther += subArray->readLatencyOther;
        
        STEnergyADC += subArray->readDynamicEnergyADC;
        STEnergyAccum += subArray->readDynamicEnergyAccum;
        STEnergyOther += subArray->readDynamicEnergyOther;
    }

	for (int k=0; k<Segnum; k++) {                 // calculate segmentsubArray through the total softmax vectors
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
        columnResistance = GetColumnResistance(input, subArrayInWeight, cell, param->parallelRead, subArray->resCellAccess);
        // cout << "rowResistance" << endl;
		rowResistance = GetRowResistance(input, subArrayInWeight, cell, param->parallelBP, subArray->resCellAccess);
        
	    // cout << "CalculateLatency" << endl;
        subArray->CalculateLatency(1e20, columnResistance, rowResistance);
        // cout << "CalculatePower" << endl;
		subArray->CalculatePower(columnResistance, rowResistance);
        
        SSReadLatency += subArray->readLatency;
        SSReadDynamicEnergy += subArray->readDynamicEnergy;
        SSLeakage += subArray->leakage;
		// cout << "at i "<< k <<"SSLeakage is"<<SSLeakage<<endl;
        
        SSLatencyADC += subArray->readLatencyADC;
        SSLatencyAccum += subArray->readLatencyAccum;
        SSLatencyOther += subArray->readLatencyOther;
        
        SSEnergyADC += subArray->readDynamicEnergyADC;
        SSEnergyAccum += subArray->readDynamicEnergyAccum;
        SSEnergyOther += subArray->readDynamicEnergyOther;
    }


    adderTreeCM->CalculateLatency((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray), 0);
    adderTreeCM->CalculatePower((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray));
    // LatencyAccum = adderTreeCM->readLatency;
    // EnergyAccum += adderTreeCM->readDynamicEnergy;
	//considering buffer activation: no matter speedup or not, the total number of data transferred is fixed
	// input buffer: total num of data loaded in = weightMatrixRow*numInVector
	// output buffer: total num of data transferred = weightMatrixRow*numInVector/param->numBitInput (total num of IFM in the PE) *adderTree->numAdderTree*adderTree->numAdderBit (bit precision of OFMs) 
    bufferInputCM->CalculateLatency(0, numInVector*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
    bufferOutputCM->CalculateLatency(0, numInVector/param->numBitInput);
	bufferReg1->CalculateLatency(0, numInSeg*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
	bufferReg2->CalculateLatency(0, numInSeg*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
	bufferSoft->CalculateLatency(0, numInSeg*ceil((double) weightMatrixRow/(double) param->numRowSubArray));

    bufferInputCM->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInVector);
    bufferOutputCM->CalculatePower(weightMatrixCol/param->numColPerSynapse*adderTreeCM->numAdderBit, numInVector/param->numBitInput);
	bufferReg1->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInSeg*Segnum);
	bufferReg2->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInSeg*Segnum);
	bufferSoft->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInSeg*Segnum);
    
    busInputCM->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth));
	busToMT1->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToMT1->busWidth));
	busToMT2->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToMT2->busWidth));
	busToReg1->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToReg1->busWidth));
	busToReg2->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToReg2->busWidth));
	busSoftToInput->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInSeg/(busSoftToInput->busWidth));
	busToSoft->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToSoft->busWidth));

    busInputCM->CalculatePower(busInputCM->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth));
	busToMT1->CalculatePower(busToMT1->busWidth, weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToMT1->busWidth));
	busToMT2->CalculatePower(busToMT2->busWidth, weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToMT2->busWidth));
	busToReg1->CalculatePower(busToReg1->busWidth, weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToReg1->busWidth));
	busToReg2->CalculatePower(busToReg2->busWidth, weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToReg2->busWidth));
	busSoftToInput->CalculatePower(busSoftToInput->busWidth, weightMatrixRow/param->numRowPerSynapse*numInSeg/(busSoftToInput->busWidth));
	busToSoft->CalculatePower(busToSoft->busWidth, weightMatrixRow/param->numRowPerSynapse*numInSeg/(busToSoft->busWidth));
    
    if (param->parallelRead) {
        busOutputCM->CalculateLatency((weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
        busOutputCM->CalculatePower(busOutputCM->numRow*busOutputCM->busWidth, (weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
		busRingBroadcast->CalculateLatency((weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInSeg/param->numBitInput)/(busRingBroadcast->numRow*busRingBroadcast->busWidth));
        busRingBroadcast->CalculatePower(busRingBroadcast->numRow*busRingBroadcast->busWidth, (weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInSeg/param->numBitInput)/(busRingBroadcast->numRow*busRingBroadcast->busWidth));
	} else {
        busOutputCM->CalculateLatency((weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
        busOutputCM->CalculatePower(busOutputCM->numRow*busOutputCM->busWidth, (weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
		busRingBroadcast->CalculateLatency((weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInSeg/param->numBitInput)/(busRingBroadcast->numRow*busRingBroadcast->busWidth));
        busRingBroadcast->CalculatePower(busRingBroadcast->numRow*busRingBroadcast->busWidth, (weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInSeg/param->numBitInput)/(busRingBroadcast->numRow*busRingBroadcast->busWidth));
    }

    Leakage = 3*IWLeakage +  Segnum*STLeakage +  SSLeakage + 							//Weight Q,K,V & transpose I & I
			  adderTreeCM->leakage +							 						//AdderTree
			  bufferInputCM->leakage + bufferOutputCM->leakage + bufferSoft->leakage +	//Buffer
			  bufferReg1->leakage + Segnum*(bufferReg2->leakage);
	
	LatencyADC = 3*IWLatencyADC + STLatencyADC + SSLatencyADC;
	LatencyOther = 3*IWLatencyOther + STLatencyOther + SSLatencyOther;
	LatencyAccum = 3*IWLatencyAccum + STLatencyAccum + SSLatencyAccum;
	BufferLatency =  bufferReg1 -> readLatency + bufferReg2->readLatency + 
					 bufferInputCM->readLatency + bufferOutputCM->readLatency + bufferSoft->readLatency;
    BufferDynamicEnergy = bufferReg1 -> readDynamicEnergy + bufferReg2->readDynamicEnergy + 
						  bufferInputCM->readDynamicEnergy + bufferOutputCM->readDynamicEnergy + bufferSoft->readDynamicEnergy;
	BusLatency	= busInputCM->readLatency + busOutputCM->readLatency +		//bus
				  busToMT1->readLatency + busToMT2->readLatency + busToReg1->readLatency + busToReg2->readLatency +
				  busRingBroadcast->readLatency + busSoftToInput->readLatency + busToSoft->readLatency +
				  adderTreeCM->readLatency;
	BusDynamicEnergy =  busInputCM->readDynamicEnergy + busOutputCM->readDynamicEnergy +		//bus
				  		busToMT1->readDynamicEnergy + busToMT2->readDynamicEnergy + busToReg1->readDynamicEnergy + busToReg2->readDynamicEnergy+
				  		busRingBroadcast->readDynamicEnergy + busSoftToInput->readDynamicEnergy + busToSoft->readDynamicEnergy +
						adderTreeCM->readDynamicEnergy;  
	
	EnergyADC = IWEnergyADC + pow(Segnum,2)*STEnergyADC + SSEnergyADC;			//num of seg per time * do all seg per total
	EnergyOther = IWEnergyOther + pow(Segnum,2)*STEnergyOther + SSEnergyOther;
	EnergyAccum = IWEnergyAccum + pow(Segnum,2)*STEnergyAccum + SSEnergyAccum;
    ReadLatency = 3*IWReadLatency + STReadLatency + SSReadLatency + BufferLatency + BusLatency;
    ReadDynamicEnergy = IWReadDynamicEnergy + pow(Segnum,2)*STReadDynamicEnergy + SSReadDynamicEnergy + BufferDynamicEnergy + BusDynamicEnergy;
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

    cout << " ===========>> readLatency  is: " << ReadLatency * 1e9 << "ns" << endl;
	cout << " IWReadLatency is: " << IWReadLatency * 1e9 << "ns" << endl;
	cout << " STReadLatency is: " << STReadLatency * 1e9 << "ns" << endl;
	cout << " SSReadLatency is: " << SSReadLatency * 1e9 << "ns" << endl;
	cout << " BufferLatency is: " << BufferLatency * 1e9 << "ns" << endl;
    cout << " ----------- readDynamicEnergy  is: " << ReadDynamicEnergy * 1e12 << "pJ" << endl;
    cout << " ----------- leakage Energy (Leakage * ReadLatency) is: " << LeakageEnergy * 1e12 << "pJ" << endl;
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
	cout << " ----------- Bus readDynamicEnergy is: " << BusDynamicEnergy * 1e9 << "ns" << endl;
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