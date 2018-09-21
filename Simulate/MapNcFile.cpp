#include <netcdfcpp.h>

#include "vector_utilities.h"
#include "MapNcFile.h"

using namespace std;

MapNcFile::MapNcFile(string filename){
	this->init();
	this->filep = new NcFile (filename.c_str(), NcFile::ReadOnly);

	if (filep->is_valid()){
		NcVar *xcoords = this->filep->get_var("rowCoordsPhys");
		NcVar *ycoords = this->filep->get_var("colCoordsPhys");
		NcVar *mapvar = this->filep->get_var("filteredSignal");
		NcVar *weightvar = this->filep->get_var("filteredWeight");
		size_t nx=0;
		size_t ny=0;
		long *dimtmp;
		if (xcoords->is_valid()){
			dimtmp = xcoords->edges();
			nx = dimtmp[0];
			this->raCoords.resize(nx);
			if (!xcoords->get(this->raCoords.getData(),nx)){
				cerr<<"Cannot retrieve RA coordinates from filename: "<<filename<<". Imploding"<<endl;
				exit(-1);
			}
			delete [] dimtmp;
		} else{
			cerr<<"Cannot retrieve RA coordinates from filename: "<<filename<<". Imploding"<<endl;
			exit(-1);
		}


		if (ycoords->is_valid()){
			dimtmp = ycoords->edges();
			ny = dimtmp[0];
			this->decCoords.resize(ny);
			if (!ycoords->get(this->decCoords.getData(),ny)){
				cerr<<"Cannot retrieve Dec coordinates from filename: "<<filename<<". Imploding"<<endl;
				exit(-1);
			}
			delete [] dimtmp;
		}else{
			cerr<<"Cannot retrieve Dec coordinates from filename: "<<filename<<". Imploding"<<endl;
			exit(-1);
		}
//		for (int ii=0; ii < ny; ii++){
//			cerr<<decCoords[ii]<<",";
//		}
//		exit(-1);
		if (mapvar->is_valid() && weightvar->is_valid()){
			this->signalMap.resize(nx,ny);
			this->weightMap.resize(nx,ny);
			dimtmp = mapvar->edges();
			if ((size_t)dimtmp[0] != nx || (size_t)dimtmp[1] != ny){
				cerr<<"Something is wrong with simulated map signals. Please check simulated nc file"<<endl;
				exit(-1);
			}
			delete[] dimtmp;
			double *tmpArray = new double [nx*ny];
			if (!mapvar->get(tmpArray,nx,ny,0,0,0)){
				cerr<<"Cannot retrieve Map Signals from filename: "<<filename<<". Imploding"<<endl;
				exit(-1);
			}
			for (size_t ix = 0; ix < nx; ix++)
				for (size_t jx = 0; jx <ny; jx++)
					this->signalMap[ix][jx] = tmpArray[ix*ny +jx];
			if (!weightvar->get(tmpArray,nx,ny,0,0,0)){
							cerr<<"Cannot retrieve Weight Signals from filename: "<<filename<<". Imploding"<<endl;
							exit(-1);
						}
						for (size_t ix = 0; ix < nx; ix++)
							for (size_t jx = 0; jx <ny; jx++)
								this->weightMap[ix][jx] = tmpArray[ix*ny +jx];
			delete [] tmpArray;
		}else{
			cerr<<"Cannot retrieve Map Signals from filename: "<<filename<<". Imploding"<<endl;
			exit(-1);
		}

		this->filep->close();
		this->hasData = true;
		cout<<"Map read. Size is ("<<nx<<","<<ny<<")"<<endl;
		this->setLimits();
//		char buff[200];
//		sprintf(buff, "rowCoords%d.txt",0);
//		writeVecOut(buff,this->raCoords.getData(), this->raCoords.size() );
//		sprintf(buff, "colCoords%d.txt",0);
//		writeVecOut(buff,this->decCoords.getData(), this->decCoords.size() );
//
//		exit(-1);
	}else{
		cerr<<"MapFile()::Cannot access file: "<<filename<<endl;
		exit(-1);
	}


}

MapNcFile::~MapNcFile(){
	if (this->filep)
		delete filep;
	this->init();
}

void MapNcFile::init(){
	this->hasData = false;
	this->filep=NULL;
	this->decCoords.resize(0);
	this->raCoords.resize(0);
	this->signalMap.resize(0,0);
	this->weightMap.resize(0,0);
}

void MapNcFile::setLimits(){
	maxmin(this->raCoords, &(this->raMax), &(this->raMin));
	maxmin(this->decCoords, &(this->decMax), &(this->decMin));
	cout<<"RaCoords range: ("<<this->raMin<<","<<this->raMax<<")"<<endl;
	cout<<"DecCoords range: ("<<this->decMin<<","<<this->decMax<<")"<<endl;
}

VecDoub MapNcFile::mapSignal(VecDoub rapoint, VecDoub decpoint){
	size_t npoints = rapoint.size();
	VecDoub outSignal (npoints,0.0);

	size_t xposb;
	size_t yposb;
	size_t xposb1;
	size_t yposb1;

	double intx;
	double intx1;
	double inty;

	xposb=0;
	yposb=0;
	xposb1=0;
	yposb1=0;

	for (size_t ip = 0; ip<npoints; ip ++){


		if (rapoint[ip]>= raMax || rapoint[ip] <raMin || decpoint[ip]>=decMax || decpoint[ip]< decMin){
			outSignal[ip] = 0.0;
			continue;
		}

		for (size_t ira = 0; ira<this->raCoords.size(); ira++){
			if (rapoint[ip] < this->raCoords[ira]){
				xposb=ira-1;
				break;
			}
		}
		for (size_t idec = 0; idec<this->decCoords.size(); idec++){
			if (decpoint[ip] < this->decCoords[idec]){
				yposb=idec-1;
				break;
			}
		}
		xposb1 = xposb+1;
		yposb1 = yposb+1;
		if (xposb1 >= this->raCoords.size() || xposb1 == 0){
			xposb1 = xposb;
		}
		if (yposb1 >= this->decCoords.size() || yposb1 == 0){
			yposb1 = yposb;
		}

		intx = (this->signalMap[xposb1][yposb]-this->signalMap[xposb][yposb])*(rapoint[ip]-raCoords[xposb]);
		intx1 = (this->signalMap[xposb1][yposb1]-this->signalMap[xposb][yposb1])*(rapoint[ip]-raCoords[xposb]);
		if (xposb == xposb1){
			intx = this->signalMap[xposb][yposb];
			intx1 = this->signalMap[xposb][yposb1];
		}else{
			intx = intx/(raCoords[xposb1]-raCoords[xposb]) + this->signalMap[xposb][yposb];
			intx1 = intx1/(raCoords[xposb1]-raCoords[xposb]) + this->signalMap[xposb][yposb1];
		}

		inty = (intx1-intx)*(decpoint[ip]- decCoords[yposb]);

		if (yposb == yposb1){
			inty = intx;
		}else{
			inty =inty / (decCoords[xposb1]-decCoords[xposb]) + intx;
		}
		outSignal[ip] = inty;
	}

	return outSignal;
}


VecDoub MapNcFile::fastMapSignal (VecDoub rapoint, VecDoub decpoint){
	size_t npoints = rapoint.size();
		VecDoub outSignal (npoints,0.0);

		size_t xposb;
		size_t yposb;
		size_t xposb1;
		size_t yposb1;

		double intx;
		double intx1;
		double inty;

		xposb=0;
		yposb=0;
		xposb1=0;
		yposb1=0;

		int dx = 1;
		int dy = 1;

		for (size_t ip = 0; ip<npoints; ip ++){

			if (rapoint[ip]>= raMax || rapoint[ip] <raMin || decpoint[ip]>=decMax || decpoint[ip]< decMin){
				outSignal[ip] = 0.0;
			}else {
				for (size_t ira = xposb; ira<this->raCoords.size()-1; ira=ira+dx){
					if (rapoint[ip] >= this->raCoords[ira] && rapoint[ip] < this->raCoords[ira+1]){
						xposb=ira;
						xposb1=ira+1;
						break;
					}
					if (ira == 0 && dx == -1){
						xposb = 0;
						xposb1 = xposb;
						break;
					}
					if (ira == raCoords.size()-2 && dx == 1){
						xposb = raCoords.size()-1;
						xposb1 = xposb;
						break;
					}
				}
				for (size_t idec = yposb; idec<this->decCoords.size()-1; idec = idec +dy){
					if (decpoint[ip] >= this->decCoords[idec] && decpoint[ip] < this->decCoords[idec+1]){
						yposb=idec;
						yposb1=idec+1;
						break;
					}
					if (idec == 0 && dy == -1){
						yposb = 0;
						yposb1 = yposb;
						break;
					}
					if (idec == decCoords.size()-2 && dy == 1){
						yposb = decCoords.size()-1;
						yposb1 = yposb;
						break;
					}
				}
				//cout<< "Ra[i]"<<rapoint[ip] <<"XPos:"<< xposb<< " YPos: "<<yposb<<"dx: "<<dx << " dy: "<<dy<<endl;

				intx = (this->signalMap[xposb1][yposb]-this->signalMap[xposb][yposb])*(rapoint[ip]-raCoords[xposb]);
				intx1 = (this->signalMap[xposb1][yposb1]-this->signalMap[xposb][yposb1])*(rapoint[ip]-raCoords[xposb]);
				if (xposb == xposb1){
					intx = this->signalMap[xposb][yposb];
					intx1 = this->signalMap[xposb][yposb1];
				}else{
					intx = intx/(raCoords[xposb1]-raCoords[xposb]) + this->signalMap[xposb][yposb];
					intx1 = intx1/(raCoords[xposb1]-raCoords[xposb]) + this->signalMap[xposb][yposb1];
				}

				inty = (intx1-intx)*(decpoint[ip]- decCoords[yposb]);

				if (yposb == yposb1 ){
					inty = intx;
				}else{
					inty =inty / (decCoords[yposb1]-decCoords[yposb]) + intx;
				}
				if (!finite (inty)){
					cerr<<" Not a valid floating point detected: "<<endl;
					cerr << "intx: "<< intx<< endl;
					cerr << "intx1: "<< intx1<< endl;
					cerr << "Dec coords: "<<(decCoords[yposb1]-decCoords[yposb]) <<endl;
				}
				outSignal[ip] = inty;
			}

			if (ip == npoints-1)
				break;
			else{
				if (rapoint[ip+1]<rapoint[ip])
					dx = -1;
				else
					dx = 1;
				if (decpoint[ip+1]<decpoint[ip])
					dy = -1;
				else
					dy = 1;
			}
		}

		return outSignal;

}
