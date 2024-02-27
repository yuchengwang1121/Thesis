#ifndef NEUROSIM_FORMULA_H_
#define NEUROSIM_FORMULA_H_

#include "Technology.h"

#define MAX(a,b) (((a)> (b))?(a):(b))
#define MIN(a,b) (((a)< (b))?(a):(b))

/* Calculate MOSFET gate capacitance */
double CalculateGateCap(double width, Technology tech);

double CalculateGateArea(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *height, double *width);

/* Calculate the capacitance of a logic gate */
void CalculateGateCapacitance(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *capInput, double *capOutput);

double CalculateDrainCap(
		double width, int type,
		double heightTransistorRegion, Technology tech);

double CalculateGateLeakage(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double temperature, Technology tech);

double CalculateOnResistance(double width, int type, double temperature, Technology tech);

double CalculateTransconductance(double width, int type, Technology tech);

double horowitz(double tr, double beta, double rampInput, double *rampOutput);

double CalculatePassGateArea(double widthNMOS, double widthPMOS, Technology tech, int numFold, double *height, double *width);

double NonlinearResistance(double R, double NL, double Vw, double Vr, double V);

#endif /* FORMULA_H_ */
