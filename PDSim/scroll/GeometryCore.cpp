//GeometryCore.cpp
// -------------------------------
// To make sense of all the geometry equations, refer to 
//  Bell, I., "Theoretical and Experimental Analysis of 
//  Liquid Flooded Compression in Scroll Compressors", PhD Thesis, 
//  Purdue University, May 2011

#ifndef __GNUC__
#define _CRT_SECURE_NO_WARNINGS
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include <stdio.h>
#include <string.h>

#include "gpc.h"
#include "math.h"
#include "GeometryCore.h"

int Isa,Is1,Is2,Ic1[20],Ic2[20],Id1,Id2,Idd,Iddd,Isuction,Idischarge,Iinjection[20],NCV;
struct geoVals geo_;

//***************************************************
//       Angles and Basic Geometry Functions 
//***************************************************
void coords_inv(double phi, double theta, struct geoVals * geo, char *code, double *x, double *y)
{
	// Function calculates the involute coordinates of a point given the involute angle and the scroll geometry
	// outputs are x and y, all other parameters are inputs

	double theta_m;
	double rb,phi_fie,ro;
	rb=geo->rb;
	ro=geo->ro;
	phi_fie=geo->phi.phi_fie;

	if (strcmp(code,"fi")==0)
	{
    	//case 'fi'  %fixed inner involute
		*x=rb*(cos(phi)+(phi-geo->phi.phi_fi0)*sin(phi));
		*y=rb*(sin(phi)-(phi-geo->phi.phi_fi0)*cos(phi));
	}
	else if (strcmp(code,"fo")==0)
	{
    	//case 'fo' %fixed outer involute
        *x=rb*(cos(phi)+(phi-geo->phi.phi_fo0)*sin(phi));
		*y=rb*(sin(phi)-(phi-geo->phi.phi_fo0)*cos(phi));
	}
    else if (strcmp(code,"oi")==0)
    {
	    //case 'oi' %Orbiting inner involutesc
        theta_m=-PI/2.0+phi_fie-theta;
        *x=-rb*(cos(phi)+(phi-geo->phi.phi_fi0)*sin(phi))+ro*cos(theta_m);
        *y=-rb*(sin(phi)-(phi-geo->phi.phi_fi0)*cos(phi))+ro*sin(theta_m);
    }
    else if (strcmp(code,"oo")==0)
    {
	    //case 'oo' %Orbiting outer involute
        theta_m=-PI/2.0+phi_fie-theta;
        *x=-rb*(cos(phi)+(phi-geo->phi.phi_fo0)*sin(phi))+ro*cos(theta_m);
        *y=-rb*(sin(phi)-(phi-geo->phi.phi_fo0)*cos(phi))+ro*sin(theta_m);        
    }
    else 
    {
	    // Didn't succeed
	    printf("Uh oh. Case incorrect in coords_inv");
    }
}
void coords_norm(double phi, char *code, double *nx, double *ny)
{
	int done=0;
	if (strcmp(code,"fi")==0)
	{
    	//case 'fi'  %fixed inner involute
		*nx=sin(phi);
		*ny=-cos(phi);
		done=1;
	}
	if (strcmp(code,"fo")==0)
	{
    	//case 'fo'  %fixed outer involute
		*nx=-sin(phi);
		*ny=cos(phi);
		done=1;
	}
	if (strcmp(code,"oi")==0)
	{
    	//case 'oi'  %orbiting inner involute
		*nx=-sin(phi);
		*ny=cos(phi);
		done=1;
	}
	if (strcmp(code,"oo")==0)
	{
    	//case 'oo'  %orbiting outer involute
		*nx=sin(phi);
		*ny=-cos(phi);
		done=1;
	}

	if (done==0)
    {
	    // Didn't succeed
	    printf("Uh oh. Case incorrect in coords_inv");
    }
}

void LoadCVIndices(struct geoVals * geo)
{
	int N,i;
	Isa=0;
	Is1=1;
	Id1=2;
	Is2=3;
	Id2=4;
	Idd=5;
	Iddd=6;
	Isuction=7;
	Idischarge=8;

	/*  Start compression chambers at +9
	 *  If more indices are needed, (for injection ports or other reasons) 
	 *  the 9 can be increased, or the new chambers can be added 
	 *  after the compression chambers indices 
	 */

	if (Ic1[0]!=NULL) // If Ic1 is still undefined
	{
		for (i=0;i<N;i++)
		{
			Ic1[i]=9+i;
			Ic2[i]=9+i+N;
		}
		NCV=Ic2[i-1]+1;//+1 since need to count 0-th index
	}
    
    //~ // Define the indices anyway
	//~ if (!Iinjection)
	//~ {
    //~ Iinjection=(int *)calloc(1,sizeof(int));
	//~ if (geo->flowModels.injection != NO_FLOW)
	//~ {	
		//~ Iinjection[0]=NCV;
		//~ NCV+=1;
	//~ }
	//~ else
		//~ Iinjection[0]=-1;
	//~ }
}

void sortAnglesCW(double *t1, double *t2)
{
	double temp;
	//Sort angles so that t2>t1 in a clockwise sense
	// idea from http://stackoverflow.com/questions/242404/sort-four-points-in-clockwise-order
	// more description: http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm
	
	/* If the signed area of the triangle formed between the points on a unit circle
	and the originis positive, the angles are sorted counterclockwise. Otherwise, the angles
	are sorted in a counter-clockwise manner.  Here we want the angles to be sorted CW, so
	if area is negative, swap angles
	
	Area obtained from the cross product of a vector from origin 
	to 1 and a vector to point 2, so use right hand rule to get 
	sign of cross product with unit length
	*/

	if (cos(*t1)*sin(*t2)-cos(*t2)*sin(*t1)>0)
	{	
		//Swap angles
		temp=*t1;
		*t1=*t2;
		*t2=temp;
	}
	while (*t1 > *t2)
	{
		// Make t2 bigger than t1
		*t2+=2*PI;
	}
}
void sortAnglesCCW(double *t1, double *t2)
{
	double temp;
	//Sort angles so that t2>t1 in a counter-clockwise sense
	// idea from http://stackoverflow.com/questions/242404/sort-four-points-in-clockwise-order
	// more description: http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm
	
	/* If the signed area of the triangle formed between the points on a unit circle with angles t1 and t2
	and the origin is positive, the angles are sorted counterclockwise. Otherwise, the angles
	are sorted in a counter-clockwise manner.  Here we want the angles to be sorted CCW, so
	if area is negative, swap angles
	
	Area obtained from the cross product of a vector from origin 
	to 1 and a vector to point 2, so use right hand rule to get 
	sign of cross product with unit length
	*/

	if (cos(*t1)*sin(*t2)-cos(*t2)*sin(*t1)<0)
	{	
		//Swap angles
		temp=*t1;
		*t1=*t2;
		*t2=temp;
	}
	while (*t1 > *t2)
	{
		// Make t2 bigger than t1
		*t2+=2*PI;
	}
}

double phi_s_sa(struct geoVals * geo,double theta)
{
	return (PI-(geo->phi.phi_fi0)+(geo->phi.phi_fo0))/((geo->phi.phi_fie)-(geo->phi.phi_oo0)-PI)*sin(theta)+(geo->phi.phi_fie)-PI;
}

double theta_d(struct geoVals *geo)
{
	return geo->phi.phi_fie-geo->phi.phi_fos-2*PI*nC_Max(geo)-PI;
}

double phi_d_dd(struct geoVals * geo,double theta)
{
	int iter;
	double phi_os,phi_o0,phi_ie,phi_i0,change,eps,f;
	double x1,x2,x3,y1,y2,phi,alpha;
	
	phi_os=geo->phi.phi_fos;
	phi_o0=geo->phi.phi_fo0;

	phi_ie=geo->phi.phi_fie;
	phi_i0=geo->phi.phi_fi0;
	alpha=PI-phi_i0+phi_o0;

	
	//Use secant method to calculate the involute angle at break
	eps=1e-8;
	change=999;
	iter=1;
	while ((iter<=3 || change>eps) && iter<100)
	{
		if (iter==1){x1=geo->phi.phi_fis; phi=x1;}
		if (iter==2){x2=geo->phi.phi_fis+0.1; phi=x2;}
		if (iter>2) {phi=x2;}

			f=1+cos(phi-phi_os)-(phi_os-phi_o0)*sin(phi-phi_os)+alpha*sin(phi-phi_ie+theta);
			

		if (iter==1){y1=f;}
		if (iter==2){y2=f;}
		if (iter>2)
		{
			y2=f;
			x3=x2-y2/(y2-y1)*(x2-x1);
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
		// If the value is still less than the starting angle
		// after 20 iterations
		if (iter>20 && x3<geo->phi.phi_fis)
		{
			return geo->phi.phi_fis;
		}
	}

	if (x3>geo->phi.phi_fis)
		return x3;
	else
	{
		return geo->phi.phi_fis;
	}
	//return phi_os+PI;
}
double polyArea(double *x,double*y,int N)
{
	// Takes in a non-closed set of curves which form a polygon and
	// determines the area between the curves

	// Technically this formula gives a negative area if the 
	// points are oriented clockwise, but the absolute value is
	// taken to eliminate the problem
	
	double sum=0;
	int i=0;
	for (i=0;i<N-2;i++)
	{
		sum+=x[i]*y[i+1]-x[i+1]*y[i];
	}
	//Close the two polygons
	sum+=x[N-1]*y[0]-x[0]*y[N-1];
	return fabs(sum/2.0);
}

int nC(struct geoVals *geo,double theta)
{
	return (int)((geo->phi.phi_fie-theta-geo->phi.phi_fos-PI)/(2.0*PI));
}
int nC_Max(struct geoVals *geo)
{
	return (int)((geo->phi.phi_fie-geo->phi.phi_fos-PI)/(2.0*PI));
}

// ****************************************************
//                  Volumes and Centroids 
// ****************************************************

double Vs_Wang(struct geoVals *geo ,double theta)
{
	double h,a,r,phi_e,alpha_i,alpha_o;
	
	h=geo->hs;
	a=geo->rb;
	r=geo->ro;
	phi_e=geo->phi.phi_fie;
	alpha_i=geo->phi.phi_fi0;
	alpha_o=geo->phi.phi_fo0;
	
	return h/2.0*a*r*(2.0*theta*phi_e-theta*theta-theta*(alpha_i+alpha_o+PI)+2.0*(1.0-cos(theta))-2.0*(phi_e-PI)*sin(theta)-PI/4*sin(2*theta));
}

double fxA_FixedInvolute(double rb, double phi,double phi0)
{
	// Anti-derivative term into which is substituted for integrations between the origin and an involute of the fixed scroll
    return rb*rb*rb/3.0*(4.0*(pow(phi-phi0,2)-2.0)*sin(phi)+(phi0-phi)*(pow(phi-phi0,2)-8.0)*cos(phi));
}
double fyA_FixedInvolute(double rb, double phi, double phi0)
{
	// Anti-derivative term into which is substituted for integrations between the origin and an involute of the fixed scroll
    return rb*rb*rb/3.0*((phi0-phi)*(pow(phi-phi0,2)-8.0)*sin(phi)-4.0*(pow(phi-phi0,2)-2.0)*cos(phi));
}
void Vs1_calcs(double theta,struct geoVals *geo,double *Vs,double *dVs,double *cx,double *cy)
{
	double B,h,ro,rb,phi_e,phi_o0,phi_i0,phi_ie,b,D,B_prime, VO,dVO,cx_O,cy_O,VIa,dVIa,cx_Ia,cy_Ia,VIb,dVIb,cx_Ib,cy_Ib,VIc,dVIc,cx_Ic,cy_Ic,cx_I,cy_I;
	h=geo->hs;
	ro=geo->ro;
	rb=geo->rb;
	phi_e=geo->phi.phi_fie;
	phi_ie=geo->phi.phi_fie;
	phi_o0=geo->phi.phi_oo0;
	phi_i0=geo->phi.phi_oi0;

	b=(-phi_o0+phi_e-PI);
    D=ro/rb*((phi_i0-phi_e)*sin(theta)-cos(theta)+1)/(phi_e-phi_i0);
    B=1.0/2.0*(sqrt(b*b-4.0*D)-b);
    B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_e-phi_i0)*sqrt(b*b-4*D));
    
    VO=h*rb*rb/6.0*(pow(phi_e-phi_i0,3)-pow(phi_e-theta-phi_i0,3));
    dVO=h*rb*rb/2.0*(pow(phi_e-theta-phi_i0,2));
    cx_O=h/VO*(fxA_FixedInvolute(rb,phi_ie,phi_i0)-fxA_FixedInvolute(rb,phi_ie-theta,phi_i0));
    cy_O=h/VO*(fyA_FixedInvolute(rb,phi_ie,phi_i0)-fyA_FixedInvolute(rb,phi_ie-theta,phi_i0));
    
    VIa=h*rb*rb/6.0*(pow(phi_e-PI+B-phi_o0,3)-pow(phi_e-PI-theta-phi_o0,3));
    dVIa=h*rb*rb/2.0*(pow(phi_e-PI+B-phi_o0,2)*B_prime+pow(phi_e-PI-theta-phi_o0,2));
    cx_Ia=h/VIa*(fxA_FixedInvolute(rb,phi_ie-PI+B,phi_o0)-fxA_FixedInvolute(rb,phi_ie-PI-theta,phi_o0));
    cy_Ia=h/VIa*(fyA_FixedInvolute(rb,phi_ie-PI+B,phi_o0)-fyA_FixedInvolute(rb,phi_ie-PI-theta,phi_o0));
    
    VIb=h*rb*ro/2.0*((B-phi_o0+phi_e-PI)*sin(B+theta)+cos(B+theta));
    dVIb=h*rb*ro*(B_prime+1)/2.0*((phi_e-PI+B-phi_o0)*cos(B+theta)-sin(B+theta));
    cx_Ib=1.0/3.0*(-rb*(B-phi_o0+phi_e-PI)*sin(B+phi_e)-rb*cos(B+phi_e)-ro*sin(theta-phi_e));
    cy_Ib=1.0/3.0*(-rb*sin(B+phi_e)+rb*(B-phi_o0+phi_e-PI)*cos(B+phi_e)-ro*cos(theta-phi_e));
    
    VIc=h*rb*ro/2.0;
    dVIc=0;
    cx_Ic=1.0/3.0*(rb*(-theta-phi_o0+phi_e-PI)*sin(theta-phi_e)-ro*sin(theta-phi_e)-rb*cos(theta-phi_e));
    cy_Ic=1.0/3.0*(rb*sin(theta-phi_e)+rb*(-theta-phi_o0+phi_e-PI)*cos(theta-phi_e)-ro*cos(theta-phi_e));
    
    cx_I=-(cx_Ia*VIa+cx_Ib*VIb-cx_Ic*VIc)/(VIa+VIb-VIc)+ro*cos(phi_ie-PI/2.0-theta);
    cy_I=-(cy_Ia*VIa+cy_Ib*VIb-cy_Ic*VIc)/(VIa+VIb-VIc)+ro*sin(phi_ie-PI/2.0-theta);
    
    *Vs=VO-(VIa+VIb-VIc)+1e-9;
    *dVs=dVO-(dVIa+dVIb-dVIc);
    *cx=(cx_O*VO-cx_I*(VIa+VIb-VIc))/ *Vs;
    *cy=(cy_O*VO-cy_I*(VIa+VIb-VIc))/ *Vs;
}

void Vc1_calcs(double theta,struct geoVals *geo,int alpha,double *V,double *dV,double *cx,double *cy)
{
	double h,ro,rb,phi_ie,phi_o0,phi_i0,psi;
	h=geo->hs;
	ro=geo->ro;
	rb=geo->rb;
	phi_ie=geo->phi.phi_fie;
	phi_o0=geo->phi.phi_oo0;
	phi_i0=geo->phi.phi_oi0;

	*V=-PI*h*rb*ro*(2*theta+4*alpha*PI-2*phi_ie-PI+phi_i0+phi_o0);
	*dV=-2.0*PI*h*rb*ro;
	psi=rb/3.0*(3.0*theta*theta+6.0*phi_o0*theta+3.0*phi_o0*phi_o0+PI*PI-15.0+(theta+phi_o0)*(12.0*PI*alpha-6.0*phi_ie)+3.0*phi_ie*phi_ie+12.0*PI*alpha*(PI*alpha-phi_ie))/(2.0*theta+phi_o0-2.0*phi_ie+phi_i0+4.0*PI*alpha-PI);
	*cx=-2.0*rb*cos(theta-phi_ie)-psi*sin(theta-phi_ie);
	*cy=+2.0*rb*sin(theta-phi_ie)-psi*cos(theta-phi_ie);
}
void Vd1_calcs(double theta,struct geoVals *geo,double *V,double *dV,double *cx,double *cy)
{
	double h,ro,rb,phi_ie,phi_o0,phi_i0,phi_os,VO,dVO,cx_O,cy_O,VIa,dVIa,cx_Ia,cy_Ia,VIb,dVIb,cx_Ib,cy_Ib,VIc,dVIc,cx_Ic,cy_Ic,VId,dVId,cx_Id,cy_Id,VI,dVI,cx_I,cy_I;
	h=geo->hs;
	ro=geo->ro;
	rb=geo->rb;
	phi_ie=geo->phi.phi_fie;
	phi_o0=geo->phi.phi_oo0;
	phi_i0=geo->phi.phi_oi0;
	phi_os=geo->phi.phi_fos;

	VO=h*rb*rb/6.0*(pow(phi_ie-theta-2.0*PI*nC(geo,theta)-phi_i0,3)-pow(phi_os+PI-phi_i0,3));
	dVO=-h*rb*rb/2.0*(pow(phi_ie-theta-2.0*PI*nC(geo,theta)-phi_i0,2));
	cx_O=h/VO*(fxA_FixedInvolute(rb,phi_ie-theta-2.0*PI*nC(geo,theta),phi_i0)-fxA_FixedInvolute(rb,phi_os+PI,phi_i0));
	cy_O=h/VO*(fyA_FixedInvolute(rb,phi_ie-theta-2.0*PI*nC(geo,theta),phi_i0)-fyA_FixedInvolute(rb,phi_os+PI,phi_i0));
    
	VIa=h*rb*rb/6.0*(pow(phi_ie-theta-2.0*PI*nC(geo,theta)-PI-phi_o0,3)-pow(phi_os-phi_o0,3));
	dVIa=-h*rb*rb/2.0*(pow(phi_ie-theta-2.0*PI*nC(geo,theta)-PI-phi_o0,2));
	cx_Ia=h/VIa*(fxA_FixedInvolute(rb,phi_ie-theta-2.0*PI*nC(geo,theta)-PI,phi_o0)-fxA_FixedInvolute(rb,phi_os,phi_o0));
	cy_Ia=h/VIa*(fyA_FixedInvolute(rb,phi_ie-theta-2.0*PI*nC(geo,theta)-PI,phi_o0)-fyA_FixedInvolute(rb,phi_os,phi_o0));

	VIb=h*rb*ro/2.0*((phi_os-phi_o0)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie));
	dVIb=h*rb*ro/2.0*((phi_os-phi_o0)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie));
	cx_Ib=1.0/3.0*(-ro*sin(theta-phi_ie)+rb*(phi_os-phi_o0)*sin(phi_os)+rb*cos(phi_os));
	cy_Ib=1.0/3.0*(-ro*cos(theta-phi_ie)-rb*(phi_os-phi_o0)*cos(phi_os)+rb*sin(phi_os));

	VIc=h*rb*ro/2.0;
	dVIc=0;
	cx_Ic=1.0/3.0*((rb*(-theta+phi_ie-phi_o0-2*PI*nC(geo,theta)-PI)-ro)*sin(theta-phi_ie)-rb*cos(theta-phi_ie));
	cy_Ic=1.0/3.0*((rb*(-theta+phi_ie-phi_o0-2*PI*nC(geo,theta)-PI)-ro)*cos(theta-phi_ie)+rb*sin(theta-phi_ie));

	VId= h*rb*ro/2.0*((phi_os-phi_i0+PI)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie)+1);
	dVId=h*rb*ro/2.0*((phi_os-phi_i0+PI)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie));
	cx_Id=(rb*(2*phi_os-phi_o0-phi_i0+PI)*sin(phi_os)-2*(ro*sin(theta-phi_ie)-rb*cos(phi_os)))/3.0;
	cy_Id=(-2*(ro*cos(theta-phi_ie)-rb*sin(phi_os))-rb*(2*phi_os-phi_o0-phi_i0+PI)*cos(phi_os))/3.0;

	VI=VIa+VIb+VIc+VId;
	dVI=dVIa+dVIb+dVIc+dVId;
	cx_I=-(cx_Ia*VIa+cx_Ib*VIb+cx_Ic*VIc+cx_Id*VId)/VI+ro*cos(phi_ie-PI/2.0-theta);
	cy_I=-(cy_Ia*VIa+cy_Ib*VIb+cy_Ic*VIc+cy_Id*VId)/VI+ro*sin(phi_ie-PI/2.0-theta);

	*V=VO-VI;
	*dV=dVO-dVI;
	*cx=(cx_O*VO-cx_I*VI)/ *V;
	*cy=(cy_O*VO-cy_I*VI)/ *V;
}
void Vdd_calcs(double theta,struct geoVals *geo,double *V,double *dV,double *cx,double *cy)
{
	double hs,ro,rb,phi_ie,phi_o0,phi_i0,phi_is,phi_os,VOa,dVOa,VOb,dVOb,VOc,dVOc,
		VIa,dVIa,VIb,dVIb,xa1,ya1,ra1,ta1_2,ta1_1,xa2,ya2,ra2,ta2_1,ta2_2,m_line,
		b_line,t1_line,t2_line,xoos,yoos,x1l,y1l,x2l,ml,om;
	hs=geo->hs;
	ro=geo->ro;
	rb=geo->rb;
	phi_ie=geo->phi.phi_fie;
	phi_o0=geo->phi.phi_oo0;
	phi_i0=geo->phi.phi_oi0;
	phi_is=geo->phi.phi_fis;
	phi_os=geo->phi.phi_fos;
	xa1=geo->disc.xa_arc1;
	ya1=geo->disc.ya_arc1;
	ra1=geo->disc.ra_arc1;
	ta1_1=geo->disc.t1_arc1; //TODO: fix angle definitions
	ta1_2=geo->disc.t2_arc1;
	xa2=geo->disc.xa_arc2;
	ya2=geo->disc.ya_arc2;
	ra2=geo->disc.ra_arc2;
	ta2_2=geo->disc.t2_arc2;
	ta2_1=geo->disc.t1_arc2;
	m_line=geo->disc.m_line;
	b_line=geo->disc.b_line;
	t1_line=geo->disc.t1_line;
	t2_line=geo->disc.t2_line;
	om=phi_ie-PI/2.0-theta;

	coords_inv(phi_os,  theta, geo,"oo",&xoos,&yoos);
    
    //#################### Oa portion ####################        
    VOa=hs*((-(ra1*(cos(ta1_2)*(ya1-yoos)-sin(ta1_2)*(xa1-xoos)-ra1*ta1_2))/2)-(-(ra1*(cos(ta1_1)*(ya1-yoos)-sin(ta1_1)*(xa1-xoos)-ra1*ta1_1))/2));
    dVOa=-hs*ra1*ro/2.0*((sin(om)*sin(ta1_2)+cos(om)*cos(ta1_2))-(sin(om)*sin(ta1_1)+cos(om)*cos(ta1_1)));
    
    //#################### Ob portion ####################
    x1l=t1_line; 
    y1l=m_line*t1_line+b_line;
    VOb=hs/2.0*((ro*xoos-ro*x1l)*sin(om)-(ro*cos(om)-2.0*x1l)*yoos+y1l*(ro*cos(om)-2.0*xoos));
    dVOb=ro*hs/2.0*(ro-yoos*sin(om)-xoos*cos(om)-y1l*sin(om)-x1l*cos(om));

	//#################### Oc portion ####################

	VOc=rb*hs/6*(
           3*ro*(phi_os-phi_i0+PI)*sin(theta+phi_os-phi_ie)
          +3*ro*cos(theta+phi_os-phi_ie)
          +3*ro*(phi_is-phi_i0)*sin(theta+phi_is-phi_ie)
          +3*ro*cos(theta+phi_is-phi_ie)
          +3*rb*((phi_is-phi_i0)*(phi_os-phi_o0)+1)*sin(phi_os-phi_is)
          -3*rb*(phi_os-phi_o0-phi_is+phi_i0)*cos(phi_os-phi_is)
          +rb*(pow(phi_os+PI-phi_i0,3)-pow(phi_is-phi_i0,3))+3*ro);
    dVOc=rb*hs*ro/2*(
          (phi_os-phi_i0+PI)*cos(theta+phi_os-phi_ie)
          -sin(theta+phi_os-phi_ie)
          +(phi_is-phi_i0)*cos(theta+phi_is-phi_ie)
          -sin(theta+phi_is-phi_ie));

    //#################### Ia portion ####################
    VIa=hs*ra2/2.0*(xa2*(sin(ta2_2)-sin(ta2_1))
                 -ya2*(cos(ta2_2)-cos(ta2_1))
                 -rb*(sin(ta2_2-phi_os)-sin(ta2_1-phi_os))
        -rb*(phi_os-phi_o0)*(cos(ta2_2-phi_os)-cos(ta2_1-phi_os))
                 +ra2*(ta2_2-ta2_1)  );
    dVIa=0.0;
    
    //#################### Ib portion #####################
    x1l=t1_line;
    x2l=t2_line;
    y1l=m_line*t1_line+b_line;
    ml=m_line;
    VIb=-hs*(x2l-x1l)/2.0*(rb*ml*(cos(phi_os)+(phi_os-phi_o0)*sin(phi_os))+b_line-rb*(sin(phi_os)-(phi_os-phi_o0)*cos(phi_os)));
    dVIb=0;
    
    *cx=ro*cos(om)/2.0;
    *cy=ro*sin(om)/2.0;
    *V =2.0*( VOa +VOb +VOc -VIa -VIb);
    *dV=2.0*(dVOa+dVOb+dVOc-dVIa-dVIb);
}

double dVs(struct geoVals *geo,double theta)
{
	double h,a,r,phi_e,alpha_i,alpha_o;

	h=geo->hs;
	a=geo->rb;
	r=geo->ro;
	phi_e=geo->phi.phi_fie;
	alpha_i=geo->phi.phi_fi0;
	alpha_o=geo->phi.phi_fo0;

	return h/2.0*a*r*(2.0*phi_e-2.0*theta-(alpha_i+alpha_o+PI)+2.0*sin(theta)-2.0*(phi_e-PI)*cos(theta)-PI/2.0*cos(2.0*theta));
}

double Vc(struct geoVals *geo,double theta,int alpha)
{
	// Convenience function return the volume directly
	double V,dV,cx,cy;
	Vc1_calcs(theta,geo,alpha,&V,&dV,&cx,&cy);
	return V;
}
double Vdisp(struct geoVals *geo)
{
	return -2.0*PI*(geo->hs)*(geo->rb)*(geo->ro)*(3.0*PI-2.0*geo->phi.phi_fie+(geo->phi.phi_fi0+geo->phi.phi_fo0));
}

double Vratio(struct geoVals *geo)
{
	return (3.0*PI-2.0*geo->phi.phi_fie+(geo->phi.phi_fi0+geo->phi.phi_fo0))/(-2.0*geo->phi.phi_fos-3.0*PI+geo->phi.phi_fi0+geo->phi.phi_fo0);
}

// This is the new method, but it doesn't work very well, so sticking 
// with the old method for now

	//double Vsa_integrate(struct geoVals *geo,double theta)
	//{
	//	double phi1, phi2,b,D,B,B_prime,rb,ro,phi_0;
	//	double phi_ie,phi_o0,phi_i0,phi_is,phi_os,V_Isa;
	//	phi_ie=geo->phi.phi_fie;
	//	phi_o0=geo->phi.phi_oo0;
	//	phi_i0=geo->phi.phi_oi0;
	//	phi_is=geo->phi.phi_fis;
	//	phi_os=geo->phi.phi_fos;
	//	ro=geo->ro;
	//	rb=geo->rb;
	//	phi_0=geo->phi.phi_fo0;
	//
	//	b=(-phi_o0+phi_ie-PI);
	//    D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0);
	//    B=1.0/2.0*(sqrt(b*b-4.0*D)-b);
	//    B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_ie-phi_i0)*sqrt(b*b-4*D));
	//	V_Isa=(geo->hs)*rb*rb/6.0*(pow(geo->phi.phi_foe-phi_0,3)-pow(phi_ie-PI+B-phi_0,3));
	//	//Volume of SA chamber is equal to area of empty walls minus the outer volume of the scroll set
	//	return (geo->hs)*PI*pow(geo->wall.r,2)-2*V_Isa;
	//}
	//double dVsa_integrate(struct geoVals *geo,double theta)
	//{
	//	double b,D,B,B_prime,rb,ro,phi_0;
	//	double phi_ie,phi_o0,phi_i0,phi_is,phi_os;
	//	phi_ie=geo->phi.phi_fie;
	//	phi_o0=geo->phi.phi_oo0;
	//	phi_i0=geo->phi.phi_oi0;
	//	phi_is=geo->phi.phi_fis;
	//	phi_os=geo->phi.phi_fos;
	//	ro=geo->ro;
	//	rb=geo->rb;
	//	phi_0=geo->phi.phi_fo0;
	//
	//	b=(-phi_o0+phi_ie-PI);
	//    D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0);
	//    B=1.0/2.0*(sqrt(b*b-4.0*D)-b);
	//    B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_ie-phi_i0)*sqrt(b*b-4*D));
	//
	//	return (geo->hs)*rb*rb*pow(phi_ie-PI+B-phi_0,2)*B_prime;
	//}

//And the old method for Vsa
double Vsa_integrate(struct geoVals *geo,double theta)
{
	double phi1, phi2,phi_0;
	double xos,yos,rb,A2,A1;

	rb=geo->rb;
	phi_0=geo->phi.phi_fo0;

	// For the fixed scroll part
	coords_inv(phi_s_sa(geo,theta),theta,geo,"oo",&xos,&yos);
	phi2=geo->phi.phi_foe;
	phi1=phi_s_sa(geo,theta);
	A2=rb/2.0*(rb*pow(phi2-phi_0,3)/3.0+cos(phi2)*(yos+(phi2-phi_0)*xos)-sin(phi2)*(xos-(phi2-phi_0)*yos));
	A1=rb/2.0*(rb*pow(phi1-phi_0,3)/3.0+cos(phi1)*(yos+(phi1-phi_0)*xos)-sin(phi1)*(xos-(phi1-phi_0)*yos));

	//Volume of SA chamber is equal to area of empty walls minus the outer volume of the scroll set
	return (geo->hs)*( PI*pow(geo->wall.r,2)-2*(A2-A1) );
}
double dVsa_integrate(struct geoVals *geo,double theta)
{
	return (Vsa_integrate(geo,theta+1e-5)-Vsa_integrate(geo,theta))/1e-5;
}


//~ //***************************************************************
//~ //                        LEAKAGE FUNCTIONS  
//~ //***************************************************************

//~ void addRadialLeak(struct flowVecVals *flowVec,struct geoVals *geo, int CV1, int CV2, double phi_max,double phi_min, int radialFlowModel, int *count)
//~ {
	//~ flowVec->CV1[*count]=CV1;
	//~ flowVec->CV2[*count]=CV2;
	//~ flowVec->A[*count]=(geo->delta_radial)*(geo->rb)*(1.0/2.0*(phi_max*phi_max-phi_min*phi_min)-(geo->phi.phi_fi0)*(phi_max-phi_min));
	//~ flowVec->flowModel[*count]=radialFlowModel;
	//~ *count=*count+1;
//~ }
//~ void addFlankLeak(struct flowVecVals *flowVec,struct geoVals *geo, int CV1, int CV2, int flankFlowModel, int *count)
//~ {
	//~ flowVec->CV1[*count]=CV1;
	//~ flowVec->CV2[*count]=CV2;
	//~ flowVec->A[*count]=(geo->delta_flank)*(geo->hs);
	//~ flowVec->flowModel[*count]=flankFlowModel;
	//~ *count=*count+1;
//~ }

void addPrimaryFlow(struct geoVals *geo, double theta, char *pathKey, double *A)
{
	double x_fie,y_fie,x_oob,y_oob,x_fis,y_fis,x_oos,y_oos;

	if (!strcmp(pathKey,"s-sa"))
	{
		coords_inv(geo->phi.phi_fie,theta,geo,"fi",&x_fie,&y_fie);
		coords_inv(phi_s_sa(geo,theta),theta,geo,"oo",&x_oob,&y_oob);
		*A=(geo->hs)*sqrt((x_fie-x_oob)*(x_fie-x_oob)+(y_fie-y_oob)*(y_fie-y_oob));
	}
	else if (!strcmp(pathKey,"d-dd"))
	{
		coords_inv(phi_d_dd(geo,theta),theta,geo,"fi",&x_fis,&y_fis);
		coords_inv(geo->phi.phi_oos,theta,geo,"oo",&x_oos,&y_oos);
		*A=(geo->hs)*sqrt((x_fis-x_oos)*(x_fis-x_oos)+(y_fis-y_oos)*(y_fis-y_oos));
		//~ if ((*A)>1000)
		//~ {
			//~ printf("Error in d-dd area %g>1000\n",*A);
			//~ printf("phi_d_dd(geo,theta): %g phi_d_dd(geo,theta-0.001): %g\n",phi_d_dd(geo,theta),phi_d_dd(geo,theta-0.001));
		//~ }
	}
    else
    {
        *A=-1;
    }
    
    return;
	//~ if (!strcmp(path,"discharge"))
	//~ {
		//~ /* Interpolate the pre-calculated discharge port blockage using Lagrange quadratic interpolation */
		//~ flowVec->A[*count]=geo->disc.Cd*interpVec((double *)geo->disc.thetaAdisc,(double *)geo->disc.Adisc,theta,NTHETA_ADISC);
	//~ }    
}
//~ void addPrimaryFlow(struct flowVecVals *flowVec,struct geoVals *geo, int CV1, int CV2, double theta,char *path,int FlowModel, int *count)
//~ {
	//~ double x_fie,y_fie,x_oob,y_oob,x_fis,y_fis,x_oos,y_oos,theta_hat,theta_0,d,r,A,alpha;

	//~ flowVec->CV1[*count]=CV1;
	//~ flowVec->CV2[*count]=CV2;

	//~ if (!strcmp(path,"s-sa"))
	//~ {
		//~ coords_inv(geo->phi.phi_fie,theta,geo,"fi",&x_fie,&y_fie);
		//~ coords_inv(phi_s_sa(geo,theta),theta,geo,"oo",&x_oob,&y_oob);
		//~ flowVec->A[*count]=(geo->hs)*sqrt((x_fie-x_oob)*(x_fie-x_oob)+(y_fie-y_oob)*(y_fie-y_oob));
	//~ }

	//~ if (!strcmp(path,"discharge"))
	//~ {
		//~ /* Interpolate the pre-calculated discharge port blockage using Lagrange quadratic interpolation */
		//~ flowVec->A[*count]=geo->disc.Cd*interpVec((double *)geo->disc.thetaAdisc,(double *)geo->disc.Adisc,theta,NTHETA_ADISC);
	//~ }

	//~ if (!strcmp(path,"d-dd"))
	//~ {
		//~ coords_inv(phi_d_dd(geo,theta),theta,geo,"fi",&x_fis,&y_fis);
		//~ coords_inv(geo->phi.phi_oos,theta,geo,"oo",&x_oos,&y_oos);
		//~ flowVec->A[*count]=(geo->hs)*sqrt((x_fis-x_oos)*(x_fis-x_oos)+(y_fis-y_oos)*(y_fis-y_oos));
		//~ if ((flowVec->A[*count])>1000)
		//~ {
			//~ printf("Error in d-dd area %g>1000\n",flowVec->A[*count]);
			//~ printf("phi_d_dd(geo,theta): %g phi_d_dd(geo,theta-0.001): %g\n",phi_d_dd(geo,theta),phi_d_dd(geo,theta-0.001));
		//~ }
	//~ }

	//~ if (!strcmp(path,"injection"))
	//~ {
		//~ // Calculate the opening of the injection port
		//~ // theta_hat is a shifted crank angle theta_hat==0 is 
		//~ // when the port is completely blocked and contact is
		//~ // right at the center of the port
		//~ r=0.001;   // Temporary value
		//~ theta_0=0; // The full-closed angle
		//~ theta_hat=theta-theta_0;
		//~ if (theta_hat<0)
			//~ theta_hat+=2*PI;
		//~ else if (theta_hat>2*PI)
			//~ theta_hat-=2*PI;
		//~ d=2*geo->rb*(1-cos(theta_hat));
		//~ if (d<r)
		//~ {
			//~ alpha=2*acos((r-d)/r);
			//~ A=r*r*(alpha-sin(alpha))/2.0;
		//~ }
		//~ else if (r<=d && d<=2*r)
		//~ {
			//~ alpha=2*acos((d-r)/r);
			//~ A=PI*r*r-r*r*(alpha-sin(alpha))/2.0;
		//~ }
		//~ else
		//~ {
			//~ A=PI*r*r;
		//~ }
		//~ flowVec->A[*count]=A;
	//~ }

	//~ if (!strcmp(path,"suction"))
	//~ {
		//~ flowVec->A[*count]=geo->suct.A_sa_suction;
	//~ }

	//~ flowVec->flowModel[*count]=FlowModel;
	//~ *count=*count+1;
//~ }
//~ int buildFlowVec(struct flowVecVals *flowVec, struct geoVals *geo,double theta,int useDDD,int LeftDischarge)
//~ {
	//~ int i, count=0, radialFlowModel,flankFlowModel,s_saFlowModel,suctionFlowModel,dischargeFlowModel,
		//~ d_ddFlowModel,alpha,Nc,injectionFlowModel,success;
	//~ double phi_ie;

	//~ LoadCVIndices(geo);
	//~ phi_ie=(geo->phi.phi_fie);
	//~ Nc=nC(geo,theta);

	//~ // Flush out any old values in flowVec
	//~ success=flushFlowVec(flowVec);
	//~ if (success!=true){printf("Failed flushing flowVec\n");}
	
	//~ //Load flow model values
	//~ radialFlowModel=geo->flowModels.radial;
	//~ flankFlowModel=geo->flowModels.flank;
	//~ s_saFlowModel=geo->flowModels.s_sa;
	//~ suctionFlowModel=geo->flowModels.suction;
	//~ d_ddFlowModel=geo->flowModels.d_dd;
	//~ dischargeFlowModel=geo->flowModels.discharge;
	//~ injectionFlowModel=geo->flowModels.injection;

	//~ if ((phi_ie-theta > phi_s_sa(geo,theta)) && Nc==0)
	//~ {
		//~ printf("Volume ratio too small and cannot solve");
	//~ }
	
	//~ //Fill radial values
	//~ if (Nc>0)
	//~ {
		//~ //Handle the suction area angles
		//~ if (phi_ie-theta > phi_s_sa(geo,theta))
		//~ {
			//~ addRadialLeak(flowVec,geo, Is1, Isa, phi_ie,phi_ie-theta,radialFlowModel, &count); 
			//~ addRadialLeak(flowVec,geo, Is2, Isa, phi_ie,phi_ie-theta,radialFlowModel, &count); 

			//~ addRadialLeak(flowVec,geo, Ic1[0], Isa, phi_ie-theta, phi_s_sa(geo,theta),radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Ic2[0], Isa, phi_ie-theta, phi_s_sa(geo,theta),radialFlowModel,&count); 

			//~ addRadialLeak(flowVec,geo, Ic1[0], Is2, phi_s_sa(geo,theta),phi_ie-theta-PI,radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Ic2[0], Is1, phi_s_sa(geo,theta),phi_ie-theta-PI,radialFlowModel,&count); 
		//~ }
		//~ else
		//~ {
			//~ addRadialLeak(flowVec,geo, Is1, Isa, phi_ie,phi_s_sa(geo,theta),radialFlowModel, &count); 
			//~ addRadialLeak(flowVec,geo, Is2, Isa, phi_ie,phi_s_sa(geo,theta),radialFlowModel, &count); 

			//~ addRadialLeak(flowVec,geo, Is2, Is1, phi_s_sa(geo,theta),phi_ie-theta,radialFlowModel,&count);

			//~ addRadialLeak(flowVec,geo, Ic1[0], Is2, phi_ie-theta,phi_ie-theta-PI,radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Ic2[0], Is1, phi_ie-theta,phi_ie-theta-PI,radialFlowModel,&count); 
		//~ }
		//~ //Compression Chamber leakages
		//~ for (alpha=1;alpha<=Nc;alpha++)
		//~ {
			                            
			//~ addRadialLeak(flowVec,geo, Ic1[alpha-1], Ic2[alpha-1], phi_ie-theta-2.0*PI*(alpha-1.0)-PI,phi_ie-theta-2.0*PI*(alpha),radialFlowModel,&count); 

			//~ //                               ^^^^ alpha-1 since C uses 0-based indexing
			//~ if (alpha>1)
			//~ {
				//~ addRadialLeak(flowVec,geo, Ic1[alpha-1], Ic2[alpha-2], phi_ie-theta-2.0*PI*(alpha-1.0),phi_ie-theta-2.0*PI*(alpha-1.0)-PI,radialFlowModel,&count); 
				//~ addRadialLeak(flowVec,geo, Ic2[alpha-1], Ic1[alpha-2], phi_ie-theta-2.0*PI*(alpha-1.0),phi_ie-theta-2.0*PI*(alpha-1.0)-PI,radialFlowModel,&count); 
			//~ }
		//~ }
		//~ //Discharge chamber leakages
		//~ // Caveat: does not take discharge region (across arcs or lines) into account
		//~ if (useDDD)
		//~ {
			//~ addRadialLeak(flowVec,geo, Iddd, Ic2[Nc-1], phi_ie-theta-2.0*PI*Nc,phi_ie-theta-2.0*PI*Nc-PI,radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Iddd, Ic1[Nc-1], phi_ie-theta-2.0*PI*Nc,phi_ie-theta-2.0*PI*Nc-PI,radialFlowModel,&count); 
		//~ }
		//~ else
		//~ {
			//~ addRadialLeak(flowVec,geo, Id1, Ic2[Nc-1], phi_ie-theta-2.0*PI*Nc,phi_ie-theta-2.0*PI*Nc-PI,radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Id2, Ic1[Nc-1], phi_ie-theta-2.0*PI*Nc,phi_ie-theta-2.0*PI*Nc-PI,radialFlowModel,&count); 
			//~ if (phi_ie-theta-2.0*PI*Nc-PI>geo->phi.phi_fis)
			//~ {
			//~ addRadialLeak(flowVec,geo, Id2, Id1, phi_ie-theta-2.0*PI*Nc-PI,geo->phi.phi_fis,radialFlowModel,&count); 
			//~ }
		//~ }
	//~ }
	//~ else
	//~ {
		//~ addRadialLeak(flowVec,geo, Is1, Isa, phi_ie, phi_s_sa(geo,theta),radialFlowModel,&count); 
		//~ addRadialLeak(flowVec,geo, Is2, Isa, phi_ie, phi_s_sa(geo,theta),radialFlowModel,&count); 
		//~ addRadialLeak(flowVec,geo, Is1, Is2, phi_s_sa(geo,theta), phi_ie-theta,radialFlowModel,&count); 
		//~ if (useDDD)
		//~ {
			//~ addRadialLeak(flowVec,geo, Iddd, Is2, phi_ie-theta, phi_ie-theta-PI,radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Iddd, Is1, phi_ie-theta, phi_ie-theta-PI,radialFlowModel,&count); 
		//~ }
		//~ else
		//~ {
			//~ addRadialLeak(flowVec,geo, Id1, Is2, phi_ie-theta, phi_ie-theta-PI,radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Id2, Is1, phi_ie-theta, phi_ie-theta-PI,radialFlowModel,&count); 
			//~ addRadialLeak(flowVec,geo, Id1, Id2, phi_ie-theta-PI, phi_ie-theta-2.0*PI,radialFlowModel,&count); 
		//~ }
		
	//~ }
	
	//~ //Fill flank values
	//~ if (Nc>0)
	//~ {
		//~ addFlankLeak(flowVec,geo,Is1,Ic1[0],flankFlowModel,&count);
		//~ addFlankLeak(flowVec,geo,Is2,Ic2[0],flankFlowModel,&count);

		//~ if (Nc>1)
		//~ {
			//~ for (alpha=2;alpha<=Nc;alpha++)
			//~ {
				//~ addFlankLeak(flowVec,geo,Ic1[alpha-1],Ic1[alpha-2],flankFlowModel,&count);
				//~ addFlankLeak(flowVec,geo,Ic2[alpha-1],Ic2[alpha-2],flankFlowModel,&count);
			//~ }
		//~ }
		//~ if (useDDD)
		//~ {
			//~ addFlankLeak(flowVec,geo,Iddd,Ic1[Nc-1],flankFlowModel,&count);
			//~ addFlankLeak(flowVec,geo,Iddd,Ic2[Nc-1],flankFlowModel,&count);
		//~ }
		//~ else
		//~ {
			//~ addFlankLeak(flowVec,geo,Id1,Ic1[Nc-1],flankFlowModel,&count);
			//~ addFlankLeak(flowVec,geo,Id2,Ic2[Nc-1],flankFlowModel,&count);
		//~ }
	//~ }
	//~ else
	//~ {
		//~ if (useDDD)
		//~ {
			//~ addFlankLeak(flowVec,geo,Is1,Iddd,flankFlowModel,&count);
			//~ addFlankLeak(flowVec,geo,Is2,Iddd,flankFlowModel,&count);
		//~ }
		//~ else
		//~ {
			//~ addFlankLeak(flowVec,geo,Is1,Id1,flankFlowModel,&count);
			//~ addFlankLeak(flowVec,geo,Is2,Id2,flankFlowModel,&count);
		//~ }
	//~ }

	//~ //Fill primary values
	//~ addPrimaryFlow(flowVec,geo,Isa,Isuction,theta,"suction",suctionFlowModel,&count);
	//~ addPrimaryFlow(flowVec,geo,Is1,Isa,theta,"s-sa",s_saFlowModel,&count);
	//~ addPrimaryFlow(flowVec,geo,Is2,Isa,theta,"s-sa",s_saFlowModel,&count);
	//~ if (useDDD==1)
	//~ {
		//~ addPrimaryFlow(flowVec,geo,Iddd,Idischarge,theta,"discharge",dischargeFlowModel,&count);
	//~ }
	//~ else
	//~ {
		//~ addPrimaryFlow(flowVec,geo,Idd,Id1,theta,"d-dd",d_ddFlowModel,&count);
		//~ addPrimaryFlow(flowVec,geo,Idd,Id2,theta,"d-dd",d_ddFlowModel,&count);
		//~ addPrimaryFlow(flowVec,geo,Idd,Idischarge,theta,"discharge",dischargeFlowModel,&count);
	//~ }
	//~ // Fill injection values
	//~ if (injectionFlowModel != NO_FLOW)
	//~ {
		//~ addPrimaryFlow(flowVec,geo,Iinjection[0],Ic1[0],theta,"injection",injectionFlowModel,&count);
		//~ addPrimaryFlow(flowVec,geo,Iinjection[0],Ic2[0],theta,"injection",injectionFlowModel,&count);
	//~ }

	//~ /* Check whether you are at the discharge angle, and if so, 
	//~ change all the CV indices to the appropriate names */
	//~ if (LeftDischarge==true)
	//~ {
		//~ for (i=0;i<count;i++)
		//~ {
			//~ if (flowVec->CV1[i]==Iddd)
				//~ flowVec->CV1[i]=Idd;
			//~ if (flowVec->CV2[i]==Iddd)
				//~ flowVec->CV2[i]=Idd;
			//~ if (flowVec->CVup[i]==Iddd)
				//~ flowVec->CVup[i]=Idd;

			//~ if (flowVec->CV1[i]==Ic1[nC_Max(geo)-1])
				//~ flowVec->CV1[i]=Id1;
			//~ if (flowVec->CV2[i]==Ic1[nC_Max(geo)-1])
				//~ flowVec->CV2[i]=Id1;
			//~ if (flowVec->CVup[i]==Ic1[nC_Max(geo)-1])
				//~ flowVec->CVup[i]=Id1;

			//~ if (flowVec->CV1[i]==Ic2[nC_Max(geo)-1])
				//~ flowVec->CV1[i]=Id2;
			//~ if (flowVec->CV2[i]==Ic2[nC_Max(geo)-1])
				//~ flowVec->CV2[i]=Id2;
			//~ if (flowVec->CVup[i]==Ic2[nC_Max(geo)-1])
				//~ flowVec->CVup[i]=Id2;
		//~ }
	
	//~ }
	//~ //Update the length of the flowVec
	//~ flowVec->N=count;
	//~ return true;
//~ }
//~ int flushFlowVec(struct flowVecVals *flowVec)
//~ {
	//~ int i;
	//~ for (i=0;i<NflowVec_MAX;i++)
	//~ {
		//~ flowVec->CV1[i]=0;
		//~ flowVec->CV2[i]=0;
		//~ flowVec->CVup[i]=0;
		//~ flowVec->A[i]=0;
		//~ flowVec->flowModel[i]=0;
		//~ flowVec->mdot[i]=0;
		//~ flowVec->mdot_L[i]=0;
		//~ flowVec->h_up[i]=0;
		//~ flowVec->h_down[i]=0;
		//~ flowVec->T_up[i]=0;
		//~ flowVec->T_down[i]=0;
		//~ flowVec->p_up[i]=0;
		//~ flowVec->p_down[i]=0;
		//~ flowVec->xL[i]=0;
		//~ flowVec->Ed[i]=0;
		//~ flowVec->Re[i]=0;
		//~ flowVec->Ma[i]=0;
	//~ }
	//~ return true;
//~ }

void setDiscGeo(struct geoVals *geo,char *Type)
{
	/*
	This function sets the discharge geometry parameters depending 
	on the type of discharge geometry desired.  The following string
	values for the variable Type are defined:

	1) "ArcLineArc", which takes the radius of arc 1 and arc 2 from the 
		disc input file
	2) "ArcLineArc-PMP", which takes the radius of arc 2 from the
		disc input file and determines the radius of arc 1 by first
		assuming that the phi_os=phi_is+pi and then solving for the 
		PMP arc radius
	3) "2Arc", which takes the radius of arc 2 from the disc input file
		and then calculates the radius of arc 1
	4) "2Arc-PMP", which does not require any information from the disc
		input file
	*/
	double a,b,c,x_is,y_is,x_os,y_os,nx_is,ny_is,nx_os,ny_os,dx,dy,r1,r2,r2_max,xarc1,xarc2,yarc1,yarc2,
		alpha,beta,d,L,xint,yint,xe,ye;
	
	// Common code for both ALA and 2Arc solutions
	coords_inv(geo->phi.phi_fis,0,geo,"fi",&x_is,&y_is);
	coords_inv(geo->phi.phi_fos,0,geo,"fo",&x_os,&y_os);
	coords_norm(geo->phi.phi_fis,"fi",&nx_is,&ny_is);
	coords_norm(geo->phi.phi_fos,"fo",&nx_os,&ny_os);
	dx=x_is-x_os;
	dy=y_is-y_os;

	// Maxiumum possible value of r2 for both perfect meshing and not perfect meshing
	a=cos(geo->phi.phi_fos-geo->phi.phi_fis)+1.0;
	b=geo->ro*a-dx*(sin(geo->phi.phi_fos)-sin(geo->phi.phi_fis))+dy*(cos(geo->phi.phi_fos)-cos(geo->phi.phi_fis));
	c=1.0/2.0*(2.0*dx*sin(geo->phi.phi_fis)*geo->ro-2.0*dy*cos(geo->phi.phi_fis)*geo->ro-dy*dy-dx*dx);
	if (geo->phi.phi_fos>geo->phi.phi_fis+PI)
		r2_max=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
	else //difference is equal to PI, it can't be less (collision otherwise)
		r2_max=-c/b;

	// First determine r2 if perfect meshing is being used for 2Arc solution
	if  (!strcmp(Type,"2Arc-PMP"))
		r2=r2_max;
	else // r2 value is passed in for either ArcLineArc w/ or w/o PMP, or 2Arc w/o PMP
		r2=geo->disc.ra_arc2;

	// If ALA, or ALA PMP, 2Arc PMP check that r2 is less than max possible
	if ( (!strcmp(Type,"ArcLineArc-PMP") || !strcmp(Type,"ArcLineArc") || !strcmp(Type,"2Arc-PMP")) 
		&& (r2>r2_max+1e-12))
	{
		printf("r2 value calculated [%g] is greater than max r2 possible [%g]\n",r2,r2_max);
		exit(EXIT_FAILURE);
	}

	// Then determine r1
	if (!strcmp(Type,"2Arc") || !strcmp(Type,"2Arc-PMP"))
	{
		// Covers both PMP and non-PMP cases
		r1=((1.0/2*dy*dy+1.0/2*dx*dx+r2*dx*sin(geo->phi.phi_fos)-r2*dy*cos(geo->phi.phi_fos))
        /(r2*cos(geo->phi.phi_fos-geo->phi.phi_fis)+dx*sin(geo->phi.phi_fis)-dy*cos(geo->phi.phi_fis)+r2));
	}
	else if (!strcmp(Type,"ArcLineArc-PMP"))
		r1=r2+geo->ro;
	else
		r1=geo->disc.ra_arc1;

	// Coordinates of the centers of the arcs
	xarc2 =  x_os+nx_os*r2;
	yarc2 =  y_os+ny_os*r2;
		// Inner starting normal has negative sign since 
		// you want the outward pointing unit normal vector
	xarc1 =  x_is-nx_is*r1;
	yarc1 =  y_is-ny_is*r1;

	if (!strcmp(Type,"ArcLineArc") || !strcmp(Type,"ArcLineArc-PMP"))
	{
        geo->disc.xa_arc2=xarc2;
        geo->disc.ya_arc2=yarc2;
        geo->disc.ra_arc2=r2;
        geo->disc.t2_arc2=atan2(y_os-yarc2,x_os-xarc2);
    
        geo->disc.xa_arc1=xarc1;
        geo->disc.ya_arc1=yarc1;
        geo->disc.ra_arc1=r1;
        geo->disc.t2_arc1=atan2(y_is-yarc1,x_is-xarc1);
                
        alpha=atan2(yarc2-yarc1,xarc2-xarc1);
        d=sqrt(pow(yarc2-yarc1,2)+pow(xarc2-xarc1,2));
		beta=acos((r1+r2)/d);
        L=sqrt(d*d-(r1+r2)*(r1+r2));
        geo->disc.t1_arc1=alpha+beta;

		xint=xarc1+r1*cos(alpha+beta)+L*sin(alpha+beta);
		yint=yarc1+r1*sin(alpha+beta)-L*cos(alpha+beta);
        geo->disc.t1_arc2=atan2(yint-yarc2,xint-xarc2);

		/* Shift the angles to ensure that going from t1 to t2 proceeds
		 in a counter-clockwise fashion, and takes less than 1 revolution */
        while (geo->disc.t2_arc2<geo->disc.t1_arc2)
            geo->disc.t2_arc2=geo->disc.t2_arc2+2.0*PI;
        while (geo->disc.t2_arc1<geo->disc.t1_arc1)
            geo->disc.t2_arc1=geo->disc.t2_arc1+2.0*PI;

        /*
        line given by y=m*t+b between the intersection points of the two circles
        */
        geo->disc.m_line=-1/tan(alpha+beta);
        geo->disc.t1_line=xarc1+r1*cos(geo->disc.t1_arc1);
        geo->disc.t2_line=xarc2+r2*cos(geo->disc.t1_arc2);
        geo->disc.b_line=yarc1+r1*sin(alpha+beta)-geo->disc.m_line*geo->disc.t1_line;
	}
        
    if (!strcmp(Type,"2Arc") || !strcmp(Type,"2Arc-PMP"))
	{
        geo->disc.xa_arc2=xarc2;
        geo->disc.ya_arc2=yarc2;
        geo->disc.ra_arc2=r2;
        geo->disc.t1_arc2=atan2(yarc1-yarc2,xarc1-xarc2);
        geo->disc.t2_arc2=atan2(y_os-yarc2,x_os-xarc2);
		/* Shift the angles to ensure that going from t1 to t2 proceeds
		 in a counter-clockwise fashion, and takes less than 1 revolution */
        while (geo->disc.t2_arc2<geo->disc.t1_arc2)
            geo->disc.t2_arc2=geo->disc.t2_arc2+2.0*PI;
    
        geo->disc.xa_arc1=xarc1;
        geo->disc.ya_arc1=yarc1;
        geo->disc.ra_arc1=r1;
        geo->disc.t2_arc1=atan2(y_is-yarc1,x_is-xarc1);
        geo->disc.t1_arc1=atan2(yarc2-yarc1,xarc2-xarc1);
		/* Shift the angles to ensure that going from t1 to t2 proceeds
		 in a counter-clockwise fashion, and takes less than 1 revolution */
        while (geo->disc.t2_arc1<geo->disc.t1_arc1)
            geo->disc.t2_arc1=geo->disc.t2_arc1+2.0*PI;
        
		/*
        line given by y=m*t+b with one element at the intersection
        point (with b=0, m=y/t)
        */
        geo->disc.b_line=0.0;
        geo->disc.t1_line=xarc2+r2*cos(geo->disc.t1_arc2);
        geo->disc.t2_line=geo->disc.t1_line;
        geo->disc.m_line=(yarc2+r2*sin(geo->disc.t1_arc2))/geo->disc.t1_line;
	}
	
	// If a negative value for disc_R is used in discharge input
	// file, use that as a fraction of the largest possible port
	// which is equal to the radius of arc 1
	if (geo->disc.R<0.0)
	{
		geo->disc.R=-geo->disc.R*geo->disc.ra_arc1;
		geo->disc.x0=geo->disc.xa_arc1;
		geo->disc.y0=geo->disc.ya_arc1;
		printf("Discharge port radius is now: %g [m]\n",geo->disc.R);
	}

	// If a negative value for wall radius is passed in, determine 
	// the near minimal wall radius that can ensure no scroll-wall contact
	// At theta=pi, the scrolls fill the largest possible circle whose radius is 
	// the distance from the outer ending angles of both scrolls with a center 
	// half the way to the origin for the orbiting scroll wrap
	// The factor of 1.03 gives a bit of breathing room 
	if (geo->wall.r<0.0)
	{
		geo->wall.x0=geo->ro/2.0*cos(geo->phi.phi_fie-PI/2-PI);
		geo->wall.y0=geo->ro/2.0*sin(geo->phi.phi_fie-PI/2-PI);
		coords_inv(geo->phi.phi_fie,PI,geo,"fo",&xe,&ye);
		geo->wall.r=1.03*sqrt(pow(geo->wall.x0-xe,2)+pow(geo->wall.y0-ye,2));
		printf("Shell wall radius is now: %g [m]\n",geo->wall.r);
	}
}

//~ // **************************************************************
//~ //                   DISCHARGE PORT BLOCKAGE CODE
//~ // **************************************************************
//~ double A_disc(struct geoVals *geo,double theta)
//~ {
	//~ gpc_polygon disc_hole,diff,lobe;
	//~ //FILE *fp;
	//~ double x,y,*xdiff=NULL,*ydiff=NULL,*t=NULL,*phi_v=NULL,om,A=0;
	//~ int i,k,I,Narc1,Nline,Narc2;

	//~ om=geo->phi.phi_fie-PI/2.0-theta;
    
	//~ //...............................................
	//~ //             Discharge port polygon
	//~ //...............................................

	//~ disc_hole.num_contours=1;
	//~ disc_hole.hole=NULL;
	//~ disc_hole.contour=(gpc_vertex_list*)malloc(1*sizeof(gpc_vertex_list));
	//~ disc_hole.contour->num_vertices=1000;
	//~ disc_hole.contour->vertex=(gpc_vertex*)malloc(1000*sizeof(gpc_vertex));

	//~ t=linspace(0,2*PI,1000);
	//~ for (i=0;i<1000;i++)
	//~ {
		//~ disc_hole.contour->vertex[i].x=geo->disc.x0+geo->disc.R*cos(t[i]);
		//~ disc_hole.contour->vertex[i].y=geo->disc.y0+geo->disc.R*sin(t[i]);
	//~ }
	//~ free(t);

	//~ lobe.num_contours=1;
	//~ lobe.hole=0;
	//~ lobe.contour=(gpc_vertex_list*)malloc(1*sizeof(gpc_vertex_list));
	//~ lobe.contour->num_vertices=500;
	//~ lobe.contour->vertex=(gpc_vertex*)malloc(500*sizeof(gpc_vertex));

	//~ // ..............................................
	//~ // Depending on type of discharge, build the polygon
	//~ // ..............................................
	
	//~ if (!strcmp(geo->disc.Type,"Arc-Arc") || !strcmp(geo->disc.Type,"2Arc") || !strcmp(geo->disc.Type,"2Arc-PMP"))
	//~ {	Narc1=150; Nline=0; Narc2=150;}
	//~ if (!strcmp(geo->disc.Type,"Arc"))
	//~ {	Narc1=300; Nline=0; Narc2=0;}
	//~ if (!strcmp(geo->disc.Type,"Arc-Line-Arc") || !strcmp(geo->disc.Type,"ArcLineArc") || !strcmp(geo->disc.Type,"ArcLineArc-PMP"))
	//~ {	Narc1=100; Nline=100; Narc2=100;}

	//~ //...............................................
	//~ //Proceeding clockwise around the orbiting scroll
	//~ //...............................................

	//~ // First is the outer involute segment with 250 elements
	//~ phi_v=linspace(geo->phi.phi_oos+PI,geo->phi.phi_oos,100);
	//~ for (i=0;i<100;i++)
	//~ {
		//~ coords_inv(phi_v[i],theta,geo,"oo",&x,&y);
		//~ lobe.contour->vertex[i].x=x;
		//~ lobe.contour->vertex[i].y=y;
	//~ }
	//~ free(phi_v);

	//~ //Next the small arc
	//~ if (Narc2>0)
	//~ {
		//~ t=linspace(geo->disc.t2_arc2,geo->disc.t1_arc2,Narc1);
		//~ for (i=0;i<Narc1;i++)
		//~ {
			//~ lobe.contour->vertex[i+100].x=-geo->disc.xa_arc2-geo->disc.ra_arc2*cos(t[i])+geo->ro*cos(om);
			//~ lobe.contour->vertex[i+100].y=-geo->disc.ya_arc2-geo->disc.ra_arc2*sin(t[i])+geo->ro*sin(om);
		//~ }
		//~ free(t);
	//~ }

	//~ //Next the line segment
	//~ if (Nline>0)
	//~ {
		//~ I=100+Narc2;
		//~ t=linspace(geo->disc.t2_line,geo->disc.t1_line,Nline);
		//~ for (i=0;i<Nline;i++)
		//~ {	
			//~ lobe.contour->vertex[i+I].x=-t[i]+geo->ro*cos(om);
			//~ lobe.contour->vertex[i+I].y=-geo->disc.m_line*t[i]-geo->disc.b_line+geo->ro*sin(om);
		//~ }
		//~ free(t);
	//~ }

	//~ //Next the big arc
	//~ if (Narc1>0)
	//~ {
		//~ I=100+Narc2+Nline;
		//~ t=linspace(geo->disc.t1_arc1,geo->disc.t2_arc1,Narc1);
		//~ for (i=0;i<Narc1;i++)
		//~ {
			//~ lobe.contour->vertex[i+I].x=-geo->disc.xa_arc1-geo->disc.ra_arc1*cos(t[i])+geo->ro*cos(om);
			//~ lobe.contour->vertex[i+I].y=-geo->disc.ya_arc1-geo->disc.ra_arc1*sin(t[i])+geo->ro*sin(om);
		//~ }
		//~ free(t);
	//~ }

	//~ // Last is the outer involute segment with 100 elements
	//~ I=100+Narc2+Nline+Narc1;
	//~ phi_v=linspace(geo->phi.phi_ois,geo->phi.phi_ois+PI/20.0,100);
	//~ for (i=0;i<100;i++)
	//~ {
		//~ coords_inv(phi_v[i],theta,geo,"oi",&x,&y);
		//~ lobe.contour->vertex[i+400].x=x;
		//~ lobe.contour->vertex[i+400].y=y;
	//~ }
	//~ free(phi_v);
	
	//~ //fp=fopen("disc_hole.csv", "wb");
	//~ //gpc_write_polygon(fp,0,&disc_hole);
	//~ //fclose(fp);

	//~ //fp=fopen("lobe.csv", "wb");
	//~ //gpc_write_polygon(fp,0,&lobe);
	//~ //fclose(fp);

	//~ gpc_polygon_clip(GPC_DIFF,&disc_hole,&lobe,&diff);

	//~ for (k=0;k<diff.num_contours;k++)
	//~ {
		//~ xdiff=(double *)calloc(diff.contour[k].num_vertices,sizeof(double));
		//~ ydiff=(double *)calloc(diff.contour[k].num_vertices,sizeof(double));
		//~ for (i=0;i<diff.contour[k].num_vertices;i++)
		//~ {
			//~ xdiff[i]=diff.contour[k].vertex[i].x;
			//~ ydiff[i]=diff.contour[k].vertex[i].y;
		//~ }
		//~ A+=polyArea(xdiff,ydiff,diff.contour[k].num_vertices-1);
		//~ free(xdiff);
		//~ free(ydiff);
	//~ }	

	//~ //fp=fopen("disc_hole.csv", "wb");
	//~ //gpc_write_polygon(fp,0,&disc_hole);
	//~ //fclose(fp);

	//~ //fp=fopen("lobe.csv", "wb");
	//~ //gpc_write_polygon(fp,0,&lobe);
	//~ //fclose(fp);

	//~ //fp=fopen("diff.csv", "wb");
	//~ //gpc_write_polygon(fp,0,&diff);
	//~ //fclose(fp);

	//~ // Free dynamically allocated polygon data
	//~ gpc_free_polygon(&disc_hole);
	//~ gpc_free_polygon(&lobe);
	//~ gpc_free_polygon(&diff);

	//~ return A;
//~ }

//~ void cleanUpGeo(void)
//~ {
	//~ free(Ic1);
	//~ Ic1=NULL;
	//~ free(Ic2);
	//~ Ic2=NULL;
	//~ free(Iinjection);
	//~ Iinjection=NULL;
//~ }