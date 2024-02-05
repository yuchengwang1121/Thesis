#ifndef TYPEDEF_H_
#define TYPEDEF_H_

namespace Type {	// To prevent name collision
	enum MemCellType {
		SRAM,
		RRAM,
		FeFET,
	};
}
enum CellAccessType
{
	CMOS_access,
	BJT_access,
	diode_access,
	none_access
};

enum DeviceRoadmap
{
	HP,		/* High performance */
	LSTP	/* Low standby power */
};

enum TransistorType
{
	conventional,	/* conventional CMOS */
	FET_2D,			/* 2D FET */
	TFET
};

enum AreaModify
{
	NONE,		/* No action, just use the original area calculation */
	MAGIC,		/* Use magic folding based on the original area */
	OVERRIDE	/* directly modify the height and width and calculate new area */
};

enum DecoderMode
{
	REGULAR_ROW,	/* Regular row mode */
	REGULAR_COL,	/* Regular column mode */
};

enum ReadCircuitMode
{
	CMOS,		/* Normal read circuit */
	OSCILLATION	/* NbO2 */
};

enum SpikingMode
{
	NONSPIKING,	/* Binary format */
	SPIKING
};

enum BusMode
{
	HORIZONTAL,	/* horizontal bus */
	VERTICAL,	/* vertical bus */
};

#endif /* TYPEDEF_H_ */

