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
#include "../headerfile/constant.h"
#include "../headerfile/formula.h"
#include "../headerfile/Param.h"
#include "../headerfile/SubArray.h"
#include "../headerfile/Buffer.h"
#include "../headerfile/HTree.h"
#include "../headerfile/AdderTree.h"
#include "../headerfile/Definition.h"
#include "../headerfile/Bus.h"
#include "../headerfile/DFF.h"

using namespace std;
AdderTree *adderTreeCM;
Bus *busInputCM;
Bus *busOutputCM;
DFF *bufferInputCM;
DFF *bufferOutputCM;

vector<vector<double>> getNetStructure(const string &inputfile);

int main(int argc, char *argv[])
{

    auto start = chrono::high_resolution_clock::now();

    gen.seed(0);

    vector<vector<double>> netStructure;
    netStructure = getNetStructure(argv[2]);                                        // get file from trace.command.sh

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
    cout << endl;
    
	/*** circuit level parameters ***/
	cell.memCellType = Type::RRAM;
	cell.accessType = CMOS_access;
	inputParameter.transistorType = conventional;
	
	switch(param->deviceroadmap) {
		case 2:	    inputParameter.deviceRoadmap = LSTP;  break;
		case 1:	    inputParameter.deviceRoadmap = HP;    break;
		case -1:	break;
		default:	exit(-1);
	}
	
    /*** build object from each class ***/
	subArray = new SubArray(inputParameter, tech, cell);
	adderTreeCM = new AdderTree(inputParameter, tech, cell);
	busInputCM = new Bus(inputParameter, tech, cell);
	busOutputCM = new Bus(inputParameter, tech, cell);
	bufferInputCM = new DFF(inputParameter, tech, cell);
	bufferOutputCM = new DFF(inputParameter, tech, cell);
	
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
    cell.heightInFeatureSize = (cell.accessType==CMOS_access)? param->heightInFeatureSize1T1R : param->heightInFeatureSizeCrossbar;         // Cell height in feature size
	cell.widthInFeatureSize = (cell.accessType==CMOS_access)? param->widthInFeatureSize1T1R : param->widthInFeatureSizeCrossbar;

    /*** subArray's property ***/
	subArray->trainingEstimation = param->trainingEstimation;            
	subArray->conventionalParallel = param->conventionalParallel;                  
	subArray->conventionalSequential = param->conventionalSequential;   
	subArray->parallelBP = param->parallelBP;	
	subArray->numRow = param->numRowSubArray;
	subArray->numCol = param->numRowSubArray;
	subArray->levelOutput = param->levelOutput;
	subArray->levelOutputBP = param->levelOutputAG;
	subArray->numColMuxed = param->numColMuxed;               // How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
	subArray->numRowMuxedBP = param->numRowMuxedAG;
    subArray->clkFreq = param->clkFreq;                       // Clock frequency
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

	int numSubArrayRow = 4;	// The number of subarray's row
	int numSubArrayCol = 3; // The number of subarray's col
	
	/*** initialize modules ***/
	subArray->Initialize(numRow, numCol, param->unitLengthWireResistance);        // initialize subArray
	if (param->parallelRead) {
		adderTreeCM->Initialize(numSubArrayRow, log2((double)param->levelOutput)+param->numBitInput+param->numColPerSynapse+1, ceil((double)numSubArrayColCM*(double)numCol/(double)param->numColMuxed));
	} else {
		adderTreeCM->Initialize(numSubArrayRow, (log2((double)numRow)+param->cellBit-1)+param->numBitInput+param->numColPerSynapse+1, ceil((double)numSubArrayColCM*(double)numCol/(double)param->numColMuxed));
	}
	
	bufferInputCM->Initialize(param->numBitInput*numRow, param->clkFreq);
	if (param->parallelRead) {
		bufferOutputCM->Initialize((numCol/param->numColMuxed)*(log2((double)param->levelOutput)+param->numBitInput+param->numColPerSynapse+adderTreeCM->numStage), param->clkFreq);
	} else {
		bufferOutputCM->Initialize((numCol/param->numColMuxed)*((log2((double)numRow)+param->cellBit-1)+param->numBitInput+param->numColPerSynapse+adderTreeCM->numStage), param->clkFreq);
	}
	
	busInputCM->Initialize(HORIZONTAL, numSubArrayRow, numSubArrayColCM, 0, numRow, subArray->height, subArray->width);
	busOutputCM->Initialize(VERTICAL, numSubArrayRow, numSubArrayColCM, 0, numCol, subArray->height, subArray->width);
    cout << endl;
    cout << "---------------------------- FloorPlan Done ------------------------------" << endl;
    cout << endl;

    /******************************************************** Initialize ********************************************************/

    /****************************************************** CalculateArea *******************************************************/
    double Area, AreaArray, AreaADC, AreaAccum, AreaOther;
    vector<double> areaResults;
	*height = 0;
	*width = 0;
	double area = 0;
	
	subArray->CalculateArea();
    adderTreeCM->CalculateArea(NULL, subArray->width, NONE);
    bufferInputCM->CalculateArea(numSubArrayRow*subArray->height, NULL, NONE);
    bufferOutputCM->CalculateArea(NULL, numSubArrayCol*subArray->width, NONE);
    
    busInputCM->CalculateArea(1, true); 
    busOutputCM->CalculateArea(1, true);	
    area += subArray->usedArea * (numSubArrayRow*numSubArrayCol) + adderTreeCM->area + bufferInputCM->area + bufferOutputCM->area;
    
    *height = sqrt(area);
    *width = area/(*height);
    
    areaResults.push_back(area);
    areaResults.push_back(subArray->areaADC*(numSubArrayRow*numSubArrayCol));
    areaResults.push_back(subArray->areaAccum*(numSubArrayRow*numSubArrayCol)+adderTreeCM->area);
    areaResults.push_back(subArray->areaOther*(numSubArrayRow*numSubArrayCol)+ bufferInputCM->area + bufferOutputCM->area);
    areaResults.push_back(subArray->areaArray*(numSubArrayRow*numSubArrayCol));

    Area        = areaResults[0];
    AreaArray   = areaResults[1];
    AreaADC     = areaResults[2];
    AreaAccum   = areaResults[3];
    AreaOther   = areaResults[4];

    cout << "-------------------------------------- Hardware Performance --------------------------------------" << endl;

    // save breakdown results of each layer to csv files
    ofstream breakdownfile;
    string breakdownfile_name = "./NeuroSim_Breakdown";
    breakdownfile_name.append(".csv");
    breakdownfile.open(breakdownfile_name, ios::app);
    if (breakdownfile.is_open())
    {
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

    /*************************************************** CalculatePerformance {}***************************************************/
    double numComputation = 0;
	numComputation += 2*(netStructure[0][0] * netStructure[0][1] * netStructure[0][2] * netStructure[0][3] * netStructure[0][4] * netStructure[0][5]);


    /*** define how many subArray are used to map the whole layer ***/
	double ReadLatency, ReadDynamicEnergy, Leakage, LatencyADC, LatencyAccum, LatencyOther, BufferLatency;
    double EnergyADC, EnergyAccum, EnergyOther, LeakageEnergy, BufferDynamicEnergy;
    ReadLatency = 0;
	ReadDynamicEnergy = 0;
	Leakage = 0;
	BufferLatency = 0;
	BufferDynamicEnergy = 0;
	EnergyADC = 0;
	EnergyAccum = 0;
	EnergyOther = 0;
	LatencyADC = 0;
	LatencyAccum = 0;
	LatencyOther = 0;

	/*** get weight matrix file Size ***/
    // I = W*L*D  = 16*128*1  = netStructure[0][0]*netStructure[0][1]*netStructure[0][2]
    // K = K*K'*D = 128*1*128 = netStructure[0][3]*netStructure[0][4]*netStructure[0][5]
	int weightMatrixRow = netStructure[0][2]*netStructure[0][3]*netStructure[0][4];
	int weightMatrixCol = netStructure[0][5];

    // weight matrix is further partitioned inside PE (among subArray) --> no duplicated
    int numRowMatrix = min(param->numRowSubArray, weightMatrixRow-param->numRowSubArray);
    int numColMatrix = min(param->numColSubArray, weightMatrixCol-param->numColSubArray);
	int numInVector;

	// load in whole file 
	vector<vector<double>> oldMemory;
	vector<vector<double>> newMemory;
	vector<vector<double>> inputVector;
	
	newMemory = LoadInWeightData(argv[4*i+5], 1, 1, param->maxConductance, param->minConductance);
	oldMemory = LoadInWeightData(argv[4*i+6], 1, 1, param->maxConductance, param->minConductance);
	inputVector = LoadInInputData(argv[4*i+7]);

    /*** assign weight and input to specific subArray ***/
	vector<vector<double>> subArrayMemoryOld;
    vector<vector<double>> subArrayMemory;
	vector<vector<double>> subArrayInput;

    subArrayMemoryOld = CopySubArray(oldMemory, param->numRowSubArray, param->numColSubArray, numRowMatrix, numColMatrix);
    subArrayMemory = CopySubArray(newMemory, param->numRowSubArray, param->numColSubArray, numRowMatrix, numColMatrix);
    subArrayInput = CopySubInput(inputVector, param->numRowSubArray, numInVector, numRowMatrix);

    for (int k=0; k<numInVector; k++) {                 // calculate single subArray through the total input vectors
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
        columnResistance = GetColumnResistance(input, subArrayMemory, cell, param->parallelRead, subArray->resCellAccess);
        rowResistance = GetRowResistance(input, subArrayMemory, cell, param->parallelBP, subArray->resCellAccess);
        
        subArray->CalculateLatency(1e20, columnResistance, rowResistance);
        subArray->CalculatePower(columnResistance, rowResistance);
        
        ReadLatency += subArray->readLatency;
        ReadDynamicEnergy += subArray->readDynamicEnergy;
        Leakage = subArray->leakage;
        
        LatencyADC += subArray->readLatencyADC;
        LatencyAccum += subArray->readLatencyAccum;
        LatencyOther += subArray->readLatencyOther;
        
        EnergyADC += subArray->readDynamicEnergyADC;
        EnergyAccum += subArray->readDynamicEnergyAccum;
        EnergyOther += subArray->readDynamicEnergyOther;
    }

    adderTreeCM->CalculateLatency((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray), 0);
    adderTreeCM->CalculatePower((int)(numInVector/param->numBitInput)*ceil(param->numColMuxed/param->numColPerSynapse), ceil((double) weightMatrixRow/(double) param->numRowSubArray));
    ReadLatency += adderTreeCM->readLatency;
    LatencyAccum += adderTreeCM->readLatency*((param->trainingEstimation)&&(layerNumber!=0)==true? 2:1);
    ReadDynamicEnergy += adderTreeCM->readDynamicEnergy;
    EnergyAccum += adderTreeCM->readDynamicEnergy*((param->trainingEstimation)&&(layerNumber!=0)==true? 2:1);
	
	//considering buffer activation: no matter speedup or not, the total number of data transferred is fixed
	// input buffer: total num of data loaded in = weightMatrixRow*numInVector
	// output buffer: total num of data transferred = weightMatrixRow*numInVector/param->numBitInput (total num of IFM in the PE) *adderTree->numAdderTree*adderTree->numAdderBit (bit precision of OFMs) 
    bufferInputCM->CalculateLatency(0, numInVector*ceil((double) weightMatrixRow/(double) param->numRowSubArray));
    bufferOutputCM->CalculateLatency(0, numInVector/param->numBitInput);
    bufferInputCM->CalculatePower(weightMatrixRow/param->numRowPerSynapse, numInVector);
    bufferOutputCM->CalculatePower(weightMatrixCol/param->numColPerSynapse*adderTreeCM->numAdderBit, numInVector/param->numBitInput);
    
    busInputCM->CalculateLatency(weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth)); 
    busInputCM->CalculatePower(busInputCM->busWidth, weightMatrixRow/param->numRowPerSynapse*numInVector/(busInputCM->busWidth));
    
    if (param->parallelRead) {
        busOutputCM->CalculateLatency((weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
        busOutputCM->CalculatePower(busOutputCM->numRow*busOutputCM->busWidth, (weightMatrixCol/param->numColPerSynapse*log2((double)param->levelOutput)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
    } else {
        busOutputCM->CalculateLatency((weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
        busOutputCM->CalculatePower(busOutputCM->numRow*busOutputCM->busWidth, (weightMatrixCol/param->numColPerSynapse*(log2((double)param->numRowSubArray)+param->cellBit-1)*numInVector/param->numBitInput)/(busOutputCM->numRow*busOutputCM->busWidth));
    }
    Leakage = Leakage*numSubArrayRow*numSubArrayCol + adderTreeCM->leakage + bufferInputCM->leakage + bufferOutputCM->leakage;
    
    ReadLatency += (bufferInputCM->readLatency + bufferOutputCM->readLatency + busInputCM->readLatency + busOutputCM->readLatency);
    ReadDynamicEnergy += (bufferInputCM->readDynamicEnergy + bufferOutputCM->readDynamicEnergy + busInputCM->readDynamicEnergy + busOutputCM->readDynamicEnergy);
    
    LeakageEnergy = Leakage * ReadLatency;
    BufferLatency = (bufferInputCM->readLatency + bufferOutputCM->readLatency)*((param->trainingEstimation)&&(layerNumber!=0)==true? 2:1);
    BufferDynamicEnergy = (bufferInputCM->readDynamicEnergy + bufferOutputCM->readDynamicEnergy)*((param->trainingEstimation)&&(layerNumber!=0)==true? 2:1);

    cout << "------------------------------ Summary --------------------------------" << endl;
    cout << endl;
    cout << "Area : " << Area * 1e12 << "um^2" << endl;
    cout << "Total CIM (Forward+Activation Gradient) array : " << AreaArray * 1e12 << "um^2" << endl;
    cout << "Total ADC (or S/As and precharger for SRAM) Area on chip : " << AreaADC * 1e12 << "um^2" << endl;
    cout << "Total Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) on chip : " << AreaAccum * 1e12 << "um^2" << endl;
    cout << "Other Peripheries (e.g. decoders, mux, switchmatrix, buffers, pooling and activation units) : " << AreaOther * 1e12 << "um^2" << endl;
    cout << endl;
    cout << "-----------------------------------Chip layer-by-layer Estimation---------------------------------" << endl;

    cout << "readLatency of Forward (per epoch) is: " << ReadLatency * 1e9 << "ns" << endl;
    cout << "readDynamicEnergy of Forward (per epoch) is: " << ReadDynamicEnergy * 1e12 << "pJ" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "leakage Energy is: " << LeakageEnergy * 1e12 << "pJ" << endl;
    cout << "leakage Power is: " << Leakage * 1e6 << "uW" << endl;
    cout << endl;
    cout << "************************ Breakdown of Latency and Dynamic Energy *************************" << endl;
    cout << endl;
    cout << "----------- ADC (or S/As and precharger for SRAM) readLatency is : " << LatencyADC * 1e9 << "ns" << endl;
    cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readLatency is : " << chipLatencyAccum * 1e9 << "ns" << endl;
    cout << "----------- Synaptic Array w/o ADC (Forward + Activate Gradient) readLatency is : " << LatencyOther * 1e9 << "ns" << endl;
    cout << "----------- Buffer readLatency is: " << BufferLatency * 1e9 << "ns" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "----------- ADC (or S/As and precharger for SRAM) readDynamicEnergy is : " << EnergyADC * 1e12 << "pJ" << endl;
    cout << "----------- Accumulation Circuits (subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units) readDynamicEnergy is : " << chipEnergyAccum * 1e12 << "pJ" << endl;
    cout << "----------- Synaptic Array w/o ADC (Forward + Activate Gradient) readDynamicEnergy is : " << EnergyOther * 1e12 << "pJ" << endl;
    cout << "----------- Buffer readDynamicEnergy is: " << BufferDynamicEnergy * 1e12 << "pJ" << endl;
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
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "------------------------------ Simulation Performance --------------------------------" << endl;
    cout << "Total Run-time of NeuroSim: " << duration.count() << " seconds" << endl;
    cout << "------------------------------ Simulation Performance --------------------------------" << endl;

    // save results to top level csv file (only total results)
    ofstream outfile;
    outfile.open("NeuroSim_Output.csv", ios::app);
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

vector<vector<double> > CopySubArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol) {
	vector<vector<double> > copy;
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

vector<vector<double> > CopySubInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow) {
	vector<vector<double> > copy;
	for (int i=0; i<numRow; i++) {
		vector<double> copyRow;
		for (int j=0; j<numInputVector; j++) {
			copyRow.push_back(orginal[positionRow+i][j]);
		}
		copy.push_back(copyRow);
		copyRow.clear();
	}
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
			if (cell.memCellType == Type::RRAM) {	// eNVM
				double totalWireResistance;
				if (cell.accessType == CMOS_access) {
					totalWireResistance = (double) 1.0/weight[i][j] + (j + 1) * param->wireResistanceRow + (weight.size() - i) * param->wireResistanceCol + cell.resistanceAccess;
				} else {
					totalWireResistance = (double) 1.0/weight[i][j] + (j + 1) * param->wireResistanceRow + (weight.size() - i) * param->wireResistanceCol;
				}
				if ((int) input[i] == 1) {
					columnG += (double) 1.0/totalWireResistance;
					activatedRow += 1 ;
				} else {
					columnG += 0;
				}
			} else if (cell.memCellType == Type::FeFET) {
				double totalWireResistance;
				totalWireResistance = (double) 1.0/weight[i][j] + (j + 1) * param->wireResistanceRow + (weight.size() - i) * param->wireResistanceCol;
				if ((int) input[i] == 1) {
					columnG += (double) 1.0/totalWireResistance;
					activatedRow += 1 ;
				} else {
					columnG += 0;
				}
				
			} else if (cell.memCellType == Type::SRAM) {	
				// SRAM: weight value do not affect sense energy --> read energy calculated in subArray.cpp (based on wireRes wireCap etc)
				double totalWireResistance = (double) (resCellAccess + param->wireResistanceCol);
				if ((int) input[i] == 1) {
					columnG += (double) 1.0/totalWireResistance;
					activatedRow += 1 ;
				} else {
					columnG += 0;
				}
			}
		}
		
		if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (!parallelRead) {  
				conductance.push_back((double) columnG/activatedRow);
			} else {
				conductance.push_back(columnG);
			}
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
			if (cell.memCellType == Type::RRAM) {	// eNVM
				if (cell.accessType == CMOS_access) {
					totalWireResistance = (double) 1.0/weight[i][j] + (i + 1) * param->wireResistanceRow + (weight[0].size() - j) * param->wireResistanceCol + cell.resistanceAccess;
				} else {
					totalWireResistance = (double) 1.0/weight[i][j] + (i + 1) * param->wireResistanceRow + (weight[0].size() - j) * param->wireResistanceCol;
				}
			} else if (cell.memCellType == Type::FeFET) {
				totalWireResistance = (double) 1.0/weight[i][j] + (i + 1) * param->wireResistanceRow + (weight[0].size() - j) * param->wireResistanceCol;
			} else if (cell.memCellType == Type::SRAM) {	
				// SRAM: weight value do not affect sense energy --> read energy calculated in subArray.cpp (based on wireRes wireCap etc)
				totalWireResistance = (double) (resCellAccess + param->wireResistanceCol);
			}
		}
		rowG = (double) 1.0/totalWireResistance * activatedCol;
		
		if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (!parallelRead) {  
				conductance.push_back((double) rowG/activatedCol);
			} else {
				conductance.push_back(rowG);
			}
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

vector<vector<double> > LoadInWeightData(const string &weightfile, int numRowPerSynapse, int numColPerSynapse, double maxConductance, double minConductance) {
	
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
				if ((param->memcelltype != 1)&&(param->synapseBit == param->cellBit)) { // training version: linear mapping
					weightrow.push_back((f+1)/2*(maxConductance-minConductance)+minConductance);
				} else {
					//normalize weight to integer
					double newdata = ((NormalizedMax-NormalizedMin)/(RealMax-RealMin)*(f-RealMax)+NormalizedMax);
					if (newdata >= 0) {
						newdata += 0.5;
					}else {
						newdata -= 0.5;
					}
					// map and expend the weight in memory array
					int cellrange = pow(2, param->cellBit);
					vector<int> synapsevector(numColPerSynapse);       
					int value = newdata; 
					if (param->BNNparallelMode) {
						if (value == 1) {
							weightrow.push_back(maxConductance);
							weightrow.push_back(minConductance);
						} else {
							weightrow.push_back(minConductance);
							weightrow.push_back(maxConductance);
						}
					} else if (param->XNORparallelMode || param->XNORsequentialMode) {
						if (value == 1) {
							weightrow.push_back(maxConductance);
							weightrowb.push_back(minConductance);
						} else {
							weightrow.push_back(minConductance);
							weightrowb.push_back(maxConductance);
						}
					} else {
						int remainder;   
						for (int z=0; z<numColPerSynapse; z++) {   
							remainder = (int) value%cellrange;
							value = (int) value/cellrange;
							synapsevector.insert(synapsevector.begin(), value/*remainder*/);
						}
						for (int u=0; u<numColPerSynapse; u++) {
							int cellvalue = synapsevector[u];
							double conductance = cellvalue/(cellrange-1) * (maxConductance-minConductance) + minConductance;
							weightrow.push_back(conductance);
						}
					}
				}
			}
		}
		weight.push_back(weightrow);
			weightrow.clear();
	}
	fileone.close();
	
	return weight;
	weight.clear();
}

vector<vector<double> > LoadInInputData(const string &inputfile) {
	
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
	
	vector<vector<double> > inputvector;              
	// load the data into inputvector ...
	for (int row=0; row<ROWin; row++) {	
		vector<double> inputvectorrow;
		vector<double> inputvectorrowb;
		getline(infile, inputline, '\n');             
		istringstream iss;
		iss.str(inputline);
		for (int col=0; col<COLin; col++) {
			while(getline(iss, inputval, ',')){	
				istringstream fs;
				fs.str(inputval);
				double f=0;
				fs >> f;
				
				if (param->BNNparallelMode) {
					if (f == 1) {
						inputvectorrow.push_back(1);
					} else {
						inputvectorrow.push_back(0);
					}
				} else if (param->XNORparallelMode || param->XNORsequentialMode) {
					if (f == 1) {
						inputvectorrow.push_back(1);
						inputvectorrowb.push_back(0);
					} else {
						inputvectorrow.push_back(0);
						inputvectorrowb.push_back(1);
					}
				} else {
					inputvectorrow.push_back(f);
				}
			}
		}
		if (param->XNORparallelMode || param->XNORsequentialMode) {
			inputvector.push_back(inputvectorrow);
			inputvectorrow.clear();
			inputvector.push_back(inputvectorrowb);
			inputvectorrowb.clear();
		} else {
			inputvector.push_back(inputvectorrow);
			inputvectorrow.clear();
		}
	}
	// close the input file ...
	infile.close();
	
	return inputvector;
	inputvector.clear();
}

    /*************************************************** Function ***************************************************/