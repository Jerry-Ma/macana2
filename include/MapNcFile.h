#ifndef _MAPNCFILE_H_
#define _MAPNCFILE_H_

#include <netcdfcpp.h>
#include <nr3.h>

class MapNcFile{
	private:
		NcFile *filep;
		VecDoub raCoords;
		VecDoub decCoords;
		MatDoub signalMap;
		MatDoub weightMap;
		bool hasData;

		double raMax;
		double raMin;
		double decMax;
		double decMin;
		void init();
		void setLimits();
	public:
		MapNcFile(string filename);
		VecDoub mapSignal(VecDoub rapoint, VecDoub decpoint);
		VecDoub fastMapSignal (VecDoub rapoint, VecDoub decpoint);
		~MapNcFile();
};

#endif
