#ifndef NEUROSIM_FUNCTION_H_
#define NEUROSIM_FUNCTION_H_

using namespace std;
#include <stdlib.h>
#include <vector>
#include <string>
#include "MemCell.h"

vector<vector<double>> getNetStructure(const string &inputfile);
vector<vector<double>> LoadInWeightData(const string &weightfile, int numRowPerSynapse, int numColPerSynapse, double maxConductance, double minConductance);
vector<vector<double>> LoadInInputData(const string &inputfile);
vector<vector<double>> CopySubArray(const vector<vector<double> > &orginal, int positionRow, int positionCol, int numRow, int numCol);
vector<vector<double>> CopySubInput(const vector<vector<double> > &orginal, int positionRow, int numInputVector, int numRow);
vector<double> GetInputVector(const vector<vector<double> > &input, int numInput, double *activityRowRead);
vector<double> GetColumnResistance(const vector<double> &input, const vector<vector<double> > &weight, MemCell& cell, bool parallelRead, double resCellAccess);
vector<double> GetRowResistance(const vector<double> &input, const vector<vector<double> > &weight, MemCell& cell, bool parallelRead, double resCellAccess);

#endif /* FUNCTION_H_ */