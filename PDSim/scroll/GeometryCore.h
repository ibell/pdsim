//File: GeometryCore.h

#ifndef SCROLLGEOMETRYCORE_H
#define SCROLLGEOMETRYCORE_H

	// ***************************************************************
	//                         Constants/Macros 
	// ***************************************************************
	#define PI 3.141592653589793

	struct phiVals{
		double phi_fi0;
		double phi_fis;
		double phi_fie;
		double phi_fo0;
		double phi_fos;
		double phi_foe;
		double phi_oi0;
		double phi_ois;
		double phi_oie;
		double phi_oo0;
		double phi_oos;
		double phi_ooe;
	};

	struct discVals{
		double x0;
		double y0;
		double R;
		double xa_arc1;
		double ya_arc1;
		double ra_arc1;
		double t1_arc1;
		double t2_arc1;
		double m_line;
		double b_line;
		double t1_line;
		double t2_line;
		double xa_arc2;
		double ya_arc2;
		double ra_arc2;
		double t1_arc2;
		double t2_arc2;
		char Type[100];
	};

	struct wallVals{
		double x0;
		double y0;
		double r;
	};

	struct geoVals{
		struct phiVals phi;
		struct discVals disc;
		struct wallVals wall;
		double rb;
		double ro;
		double t;
		double hs;
	};

	// ***************************************************************
	//                     Function Prototypes
	// ***************************************************************

	void coords_inv(double phi, double theta,struct geoVals * geo, char *code, double *x, double *y);
	void coords_norm(double phi, char *code, double *nx, double *ny);
	void sortAnglesCW(double *t1, double *t2);
	void sortAnglesCCW(double *t1, double *t2);
	double Vs_Wang(struct geoVals *geo ,double theta);
	double phi_s_sa(struct geoVals * geo,double theta);
	double phi_d_dd(struct geoVals * geo,double theta);
	double polyArea(double *x,double*y,int N);
	double Vdisp(struct geoVals *geo);
	double Vratio(struct geoVals *geo);
	int nC(struct geoVals *geo,double theta);
	int nC_Max(struct geoVals *geo);
	double theta_d(struct geoVals *geo);
   
	//~ int buildFlowVec(struct flowVecVals *flowVec, struct geoVals *geo,double theta,int useDDD,int LeftDischarge);
	//~ int flushFlowVec(struct flowVecVals *flowVec);
	//~ void addRadialLeak(struct flowVecVals *flowVec,struct geoVals *geo, int CV1, int CV2, double phi_max,double phi_min, int radialFlowModel, int *count);
	//~ void addFlankLeak(struct flowVecVals *flowVec,struct geoVals *geo, int CV1, int CV2, int flankFlowModel, int *count);
	//~ void addPrimaryFlow(struct flowVecVals *flowVec,struct geoVals *geo, int CV1, int CV2, double theta,char *path,int FlowModel, int *count);
	//~ void GeometryModel(struct geoVals *geo,double theta, int Istep, int useDDD, int LeftDischarge, double *V, double *dV, struct flowVecVals * flowVec);

	//~ double A_disc(struct geoVals *geo,double theta);
	double Vc(struct geoVals *geo,double theta,int alpha);

	// New functions
	void Vs1_calcs(double theta,struct geoVals *geo,double*V, double *dV,double *cx,double*cy);
	void Vc1_calcs(double theta,struct geoVals *geo,int alpha,double*V, double *dV,double *cx,double*cy);
	void Vd1_calcs(double theta,struct geoVals *geo,double*V, double *dV,double *cx,double*cy);
	void Vdd_calcs(double theta,struct geoVals *geo,double*V, double *dV,double *cx,double*cy);
	void setDiscGeo(struct geoVals *geo,char *Type);
    void addPrimaryFlow(struct geoVals *geo, double theta, char *pathKey, double *A);
#endif


