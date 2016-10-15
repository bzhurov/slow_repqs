//Repressilator
//true reducted model, 6 equations

#include "auto_f2c.h"
#include "math.h"

int func (integer ndim, const doublereal *u, const integer *icp,
	const doublereal *par, integer ijac,
	doublereal *f, doublereal *dfdu, doublereal *dfdp)
{

	/* System generated locals */
	integer dfdu_dim1 = ndim, dfdp_dim1 = ndim;
	//variables
	double M1 = u[0];
	double M2 = u[1];
	double M3 = u[2];
	double MI = u[3];
	double P1 = u[4];
	double P2 = u[5];
	double P3 = u[6];
	double PI = u[7];
	double S = u[8];
	double R = u[9];
	double D1 = u[10];
	double D2 = u[11];
	double D3 = u[12];
	double DR = u[13];
	double G1 = u[14];
	double G2 = u[15];
	double GI = u[16];
	double G3 = u[17];
	double G3I = u[18];

	//parameters
	double rm = par[0];
	double omega = par[1];
	double se = par[2];
	double pR = par[3];
	double rp = 0.1;
	double dm = 0.0033;
	double dp = 0.0033;
	double delta = 0.0;
	double kd = 0.025;
	double kud = 0.5;
	double kr = 0.012;
	double kur = 0.9;
	double gt = 20;
	double rmI = 0.1;
	double rs = 0.0067;
	double ds = 1.67e-05;
	double eta = 0.0067;
	double ka = 0.001667;
	double kua = 0.166667;
	double kRd = 0.000833;
	double kRud = 0.016667;
	double kRr = 0.000833;
	double kRur = 0.166667;

	//Help parameters

	//System
	f[0] = G1*rm - M1*dm;
	f[1] = G2*rm - M2*dm;
	f[2] = G3*rm - M3*dm + omega*rm*(-G3I + gt);
	f[3] = GI*rm*rmI - MI*dm;
	f[4] = 2*D1*kud + M1*rp - 2*pow(P1, 2)*kd - P1*dp;
	f[5] = 2*D2*kud + M2*rp - 2*pow(P2, 2)*kd - P2*dp;
	f[6] = 2*D3*kud + M3*rp - 2*pow(P3, 2)*kd - P3*dp;
	f[7] = MI*rp - PI*dp;
	f[8] = PI*rs + R*kua - S*ds - S*ka*(-R + pR) - eta*(S - se);
	f[9] = 2*DR*kRud - 2*pow(R, 2)*kRd - R*kua + S*ka*(-R + pR);
	f[10] = -D1*G2*kr - D1*GI*kr - D1*delta*dp - D1*kud + pow(P1, 2)*kd + kur*(-G2 + gt) + kur*(-GI + gt);
	f[11] = -D2*G3*kr - D2*delta*dp - D2*kud + pow(P2, 2)*kd + kur*(-G3 + gt);
	f[12] = -D3*G1*kr - D3*delta*dp - D3*kud + pow(P3, 2)*kd + kur*(-G1 + gt);
	f[13] = -DR*G3I*kRr - DR*delta*dp - DR*kRud + pow(R, 2)*kRd + kRur*(-G3I + gt);
	f[14] = -D3*G1*kr + kur*(-G1 + gt);
	f[15] = -D1*G2*kr + kur*(-G2 + gt);
	f[16] = -D1*GI*kr + kur*(-GI + gt);
	f[17] = -D2*G3*kr + kur*(-G3 + gt);
	f[18] = -DR*G3I*kRr + kRur*(-G3I + gt);

	//Jacobian

	if (ijac == 0)
	{
		return 0;
	}
	ARRAY2D(dfdu,0,0) = -dm;
	ARRAY2D(dfdu,0,1) = 0;
	ARRAY2D(dfdu,0,2) = 0;
	ARRAY2D(dfdu,0,3) = 0;
	ARRAY2D(dfdu,0,4) = 0;
	ARRAY2D(dfdu,0,5) = 0;
	ARRAY2D(dfdu,0,6) = 0;
	ARRAY2D(dfdu,0,7) = 0;
	ARRAY2D(dfdu,0,8) = 0;
	ARRAY2D(dfdu,0,9) = 0;
	ARRAY2D(dfdu,0,10) = 0;
	ARRAY2D(dfdu,0,11) = 0;
	ARRAY2D(dfdu,0,12) = 0;
	ARRAY2D(dfdu,0,13) = 0;
	ARRAY2D(dfdu,0,14) = rm;
	ARRAY2D(dfdu,0,15) = 0;
	ARRAY2D(dfdu,0,16) = 0;
	ARRAY2D(dfdu,0,17) = 0;
	ARRAY2D(dfdu,0,18) = 0;
	ARRAY2D(dfdu,1,0) = 0;
	ARRAY2D(dfdu,1,1) = -dm;
	ARRAY2D(dfdu,1,2) = 0;
	ARRAY2D(dfdu,1,3) = 0;
	ARRAY2D(dfdu,1,4) = 0;
	ARRAY2D(dfdu,1,5) = 0;
	ARRAY2D(dfdu,1,6) = 0;
	ARRAY2D(dfdu,1,7) = 0;
	ARRAY2D(dfdu,1,8) = 0;
	ARRAY2D(dfdu,1,9) = 0;
	ARRAY2D(dfdu,1,10) = 0;
	ARRAY2D(dfdu,1,11) = 0;
	ARRAY2D(dfdu,1,12) = 0;
	ARRAY2D(dfdu,1,13) = 0;
	ARRAY2D(dfdu,1,14) = 0;
	ARRAY2D(dfdu,1,15) = rm;
	ARRAY2D(dfdu,1,16) = 0;
	ARRAY2D(dfdu,1,17) = 0;
	ARRAY2D(dfdu,1,18) = 0;
	ARRAY2D(dfdu,2,0) = 0;
	ARRAY2D(dfdu,2,1) = 0;
	ARRAY2D(dfdu,2,2) = -dm;
	ARRAY2D(dfdu,2,3) = 0;
	ARRAY2D(dfdu,2,4) = 0;
	ARRAY2D(dfdu,2,5) = 0;
	ARRAY2D(dfdu,2,6) = 0;
	ARRAY2D(dfdu,2,7) = 0;
	ARRAY2D(dfdu,2,8) = 0;
	ARRAY2D(dfdu,2,9) = 0;
	ARRAY2D(dfdu,2,10) = 0;
	ARRAY2D(dfdu,2,11) = 0;
	ARRAY2D(dfdu,2,12) = 0;
	ARRAY2D(dfdu,2,13) = 0;
	ARRAY2D(dfdu,2,14) = 0;
	ARRAY2D(dfdu,2,15) = 0;
	ARRAY2D(dfdu,2,16) = 0;
	ARRAY2D(dfdu,2,17) = rm;
	ARRAY2D(dfdu,2,18) = -omega*rm;
	ARRAY2D(dfdu,3,0) = 0;
	ARRAY2D(dfdu,3,1) = 0;
	ARRAY2D(dfdu,3,2) = 0;
	ARRAY2D(dfdu,3,3) = -dm;
	ARRAY2D(dfdu,3,4) = 0;
	ARRAY2D(dfdu,3,5) = 0;
	ARRAY2D(dfdu,3,6) = 0;
	ARRAY2D(dfdu,3,7) = 0;
	ARRAY2D(dfdu,3,8) = 0;
	ARRAY2D(dfdu,3,9) = 0;
	ARRAY2D(dfdu,3,10) = 0;
	ARRAY2D(dfdu,3,11) = 0;
	ARRAY2D(dfdu,3,12) = 0;
	ARRAY2D(dfdu,3,13) = 0;
	ARRAY2D(dfdu,3,14) = 0;
	ARRAY2D(dfdu,3,15) = 0;
	ARRAY2D(dfdu,3,16) = rm*rmI;
	ARRAY2D(dfdu,3,17) = 0;
	ARRAY2D(dfdu,3,18) = 0;
	ARRAY2D(dfdu,4,0) = rp;
	ARRAY2D(dfdu,4,1) = 0;
	ARRAY2D(dfdu,4,2) = 0;
	ARRAY2D(dfdu,4,3) = 0;
	ARRAY2D(dfdu,4,4) = -4*P1*kd - dp;
	ARRAY2D(dfdu,4,5) = 0;
	ARRAY2D(dfdu,4,6) = 0;
	ARRAY2D(dfdu,4,7) = 0;
	ARRAY2D(dfdu,4,8) = 0;
	ARRAY2D(dfdu,4,9) = 0;
	ARRAY2D(dfdu,4,10) = 2*kud;
	ARRAY2D(dfdu,4,11) = 0;
	ARRAY2D(dfdu,4,12) = 0;
	ARRAY2D(dfdu,4,13) = 0;
	ARRAY2D(dfdu,4,14) = 0;
	ARRAY2D(dfdu,4,15) = 0;
	ARRAY2D(dfdu,4,16) = 0;
	ARRAY2D(dfdu,4,17) = 0;
	ARRAY2D(dfdu,4,18) = 0;
	ARRAY2D(dfdu,5,0) = 0;
	ARRAY2D(dfdu,5,1) = rp;
	ARRAY2D(dfdu,5,2) = 0;
	ARRAY2D(dfdu,5,3) = 0;
	ARRAY2D(dfdu,5,4) = 0;
	ARRAY2D(dfdu,5,5) = -4*P2*kd - dp;
	ARRAY2D(dfdu,5,6) = 0;
	ARRAY2D(dfdu,5,7) = 0;
	ARRAY2D(dfdu,5,8) = 0;
	ARRAY2D(dfdu,5,9) = 0;
	ARRAY2D(dfdu,5,10) = 0;
	ARRAY2D(dfdu,5,11) = 2*kud;
	ARRAY2D(dfdu,5,12) = 0;
	ARRAY2D(dfdu,5,13) = 0;
	ARRAY2D(dfdu,5,14) = 0;
	ARRAY2D(dfdu,5,15) = 0;
	ARRAY2D(dfdu,5,16) = 0;
	ARRAY2D(dfdu,5,17) = 0;
	ARRAY2D(dfdu,5,18) = 0;
	ARRAY2D(dfdu,6,0) = 0;
	ARRAY2D(dfdu,6,1) = 0;
	ARRAY2D(dfdu,6,2) = rp;
	ARRAY2D(dfdu,6,3) = 0;
	ARRAY2D(dfdu,6,4) = 0;
	ARRAY2D(dfdu,6,5) = 0;
	ARRAY2D(dfdu,6,6) = -4*P3*kd - dp;
	ARRAY2D(dfdu,6,7) = 0;
	ARRAY2D(dfdu,6,8) = 0;
	ARRAY2D(dfdu,6,9) = 0;
	ARRAY2D(dfdu,6,10) = 0;
	ARRAY2D(dfdu,6,11) = 0;
	ARRAY2D(dfdu,6,12) = 2*kud;
	ARRAY2D(dfdu,6,13) = 0;
	ARRAY2D(dfdu,6,14) = 0;
	ARRAY2D(dfdu,6,15) = 0;
	ARRAY2D(dfdu,6,16) = 0;
	ARRAY2D(dfdu,6,17) = 0;
	ARRAY2D(dfdu,6,18) = 0;
	ARRAY2D(dfdu,7,0) = 0;
	ARRAY2D(dfdu,7,1) = 0;
	ARRAY2D(dfdu,7,2) = 0;
	ARRAY2D(dfdu,7,3) = rp;
	ARRAY2D(dfdu,7,4) = 0;
	ARRAY2D(dfdu,7,5) = 0;
	ARRAY2D(dfdu,7,6) = 0;
	ARRAY2D(dfdu,7,7) = -dp;
	ARRAY2D(dfdu,7,8) = 0;
	ARRAY2D(dfdu,7,9) = 0;
	ARRAY2D(dfdu,7,10) = 0;
	ARRAY2D(dfdu,7,11) = 0;
	ARRAY2D(dfdu,7,12) = 0;
	ARRAY2D(dfdu,7,13) = 0;
	ARRAY2D(dfdu,7,14) = 0;
	ARRAY2D(dfdu,7,15) = 0;
	ARRAY2D(dfdu,7,16) = 0;
	ARRAY2D(dfdu,7,17) = 0;
	ARRAY2D(dfdu,7,18) = 0;
	ARRAY2D(dfdu,8,0) = 0;
	ARRAY2D(dfdu,8,1) = 0;
	ARRAY2D(dfdu,8,2) = 0;
	ARRAY2D(dfdu,8,3) = 0;
	ARRAY2D(dfdu,8,4) = 0;
	ARRAY2D(dfdu,8,5) = 0;
	ARRAY2D(dfdu,8,6) = 0;
	ARRAY2D(dfdu,8,7) = rs;
	ARRAY2D(dfdu,8,8) = -ds - eta - ka*(-R + pR);
	ARRAY2D(dfdu,8,9) = S*ka + kua;
	ARRAY2D(dfdu,8,10) = 0;
	ARRAY2D(dfdu,8,11) = 0;
	ARRAY2D(dfdu,8,12) = 0;
	ARRAY2D(dfdu,8,13) = 0;
	ARRAY2D(dfdu,8,14) = 0;
	ARRAY2D(dfdu,8,15) = 0;
	ARRAY2D(dfdu,8,16) = 0;
	ARRAY2D(dfdu,8,17) = 0;
	ARRAY2D(dfdu,8,18) = 0;
	ARRAY2D(dfdu,9,0) = 0;
	ARRAY2D(dfdu,9,1) = 0;
	ARRAY2D(dfdu,9,2) = 0;
	ARRAY2D(dfdu,9,3) = 0;
	ARRAY2D(dfdu,9,4) = 0;
	ARRAY2D(dfdu,9,5) = 0;
	ARRAY2D(dfdu,9,6) = 0;
	ARRAY2D(dfdu,9,7) = 0;
	ARRAY2D(dfdu,9,8) = ka*(-R + pR);
	ARRAY2D(dfdu,9,9) = -4*R*kRd - S*ka - kua;
	ARRAY2D(dfdu,9,10) = 0;
	ARRAY2D(dfdu,9,11) = 0;
	ARRAY2D(dfdu,9,12) = 0;
	ARRAY2D(dfdu,9,13) = 2*kRud;
	ARRAY2D(dfdu,9,14) = 0;
	ARRAY2D(dfdu,9,15) = 0;
	ARRAY2D(dfdu,9,16) = 0;
	ARRAY2D(dfdu,9,17) = 0;
	ARRAY2D(dfdu,9,18) = 0;
	ARRAY2D(dfdu,10,0) = 0;
	ARRAY2D(dfdu,10,1) = 0;
	ARRAY2D(dfdu,10,2) = 0;
	ARRAY2D(dfdu,10,3) = 0;
	ARRAY2D(dfdu,10,4) = 2*P1*kd;
	ARRAY2D(dfdu,10,5) = 0;
	ARRAY2D(dfdu,10,6) = 0;
	ARRAY2D(dfdu,10,7) = 0;
	ARRAY2D(dfdu,10,8) = 0;
	ARRAY2D(dfdu,10,9) = 0;
	ARRAY2D(dfdu,10,10) = -G2*kr - GI*kr - delta*dp - kud;
	ARRAY2D(dfdu,10,11) = 0;
	ARRAY2D(dfdu,10,12) = 0;
	ARRAY2D(dfdu,10,13) = 0;
	ARRAY2D(dfdu,10,14) = 0;
	ARRAY2D(dfdu,10,15) = -D1*kr - kur;
	ARRAY2D(dfdu,10,16) = -D1*kr - kur;
	ARRAY2D(dfdu,10,17) = 0;
	ARRAY2D(dfdu,10,18) = 0;
	ARRAY2D(dfdu,11,0) = 0;
	ARRAY2D(dfdu,11,1) = 0;
	ARRAY2D(dfdu,11,2) = 0;
	ARRAY2D(dfdu,11,3) = 0;
	ARRAY2D(dfdu,11,4) = 0;
	ARRAY2D(dfdu,11,5) = 2*P2*kd;
	ARRAY2D(dfdu,11,6) = 0;
	ARRAY2D(dfdu,11,7) = 0;
	ARRAY2D(dfdu,11,8) = 0;
	ARRAY2D(dfdu,11,9) = 0;
	ARRAY2D(dfdu,11,10) = 0;
	ARRAY2D(dfdu,11,11) = -G3*kr - delta*dp - kud;
	ARRAY2D(dfdu,11,12) = 0;
	ARRAY2D(dfdu,11,13) = 0;
	ARRAY2D(dfdu,11,14) = 0;
	ARRAY2D(dfdu,11,15) = 0;
	ARRAY2D(dfdu,11,16) = 0;
	ARRAY2D(dfdu,11,17) = -D2*kr - kur;
	ARRAY2D(dfdu,11,18) = 0;
	ARRAY2D(dfdu,12,0) = 0;
	ARRAY2D(dfdu,12,1) = 0;
	ARRAY2D(dfdu,12,2) = 0;
	ARRAY2D(dfdu,12,3) = 0;
	ARRAY2D(dfdu,12,4) = 0;
	ARRAY2D(dfdu,12,5) = 0;
	ARRAY2D(dfdu,12,6) = 2*P3*kd;
	ARRAY2D(dfdu,12,7) = 0;
	ARRAY2D(dfdu,12,8) = 0;
	ARRAY2D(dfdu,12,9) = 0;
	ARRAY2D(dfdu,12,10) = 0;
	ARRAY2D(dfdu,12,11) = 0;
	ARRAY2D(dfdu,12,12) = -G1*kr - delta*dp - kud;
	ARRAY2D(dfdu,12,13) = 0;
	ARRAY2D(dfdu,12,14) = -D3*kr - kur;
	ARRAY2D(dfdu,12,15) = 0;
	ARRAY2D(dfdu,12,16) = 0;
	ARRAY2D(dfdu,12,17) = 0;
	ARRAY2D(dfdu,12,18) = 0;
	ARRAY2D(dfdu,13,0) = 0;
	ARRAY2D(dfdu,13,1) = 0;
	ARRAY2D(dfdu,13,2) = 0;
	ARRAY2D(dfdu,13,3) = 0;
	ARRAY2D(dfdu,13,4) = 0;
	ARRAY2D(dfdu,13,5) = 0;
	ARRAY2D(dfdu,13,6) = 0;
	ARRAY2D(dfdu,13,7) = 0;
	ARRAY2D(dfdu,13,8) = 0;
	ARRAY2D(dfdu,13,9) = 2*R*kRd;
	ARRAY2D(dfdu,13,10) = 0;
	ARRAY2D(dfdu,13,11) = 0;
	ARRAY2D(dfdu,13,12) = 0;
	ARRAY2D(dfdu,13,13) = -G3I*kRr - delta*dp - kRud;
	ARRAY2D(dfdu,13,14) = 0;
	ARRAY2D(dfdu,13,15) = 0;
	ARRAY2D(dfdu,13,16) = 0;
	ARRAY2D(dfdu,13,17) = 0;
	ARRAY2D(dfdu,13,18) = -DR*kRr - kRur;
	ARRAY2D(dfdu,14,0) = 0;
	ARRAY2D(dfdu,14,1) = 0;
	ARRAY2D(dfdu,14,2) = 0;
	ARRAY2D(dfdu,14,3) = 0;
	ARRAY2D(dfdu,14,4) = 0;
	ARRAY2D(dfdu,14,5) = 0;
	ARRAY2D(dfdu,14,6) = 0;
	ARRAY2D(dfdu,14,7) = 0;
	ARRAY2D(dfdu,14,8) = 0;
	ARRAY2D(dfdu,14,9) = 0;
	ARRAY2D(dfdu,14,10) = 0;
	ARRAY2D(dfdu,14,11) = 0;
	ARRAY2D(dfdu,14,12) = -G1*kr;
	ARRAY2D(dfdu,14,13) = 0;
	ARRAY2D(dfdu,14,14) = -D3*kr - kur;
	ARRAY2D(dfdu,14,15) = 0;
	ARRAY2D(dfdu,14,16) = 0;
	ARRAY2D(dfdu,14,17) = 0;
	ARRAY2D(dfdu,14,18) = 0;
	ARRAY2D(dfdu,15,0) = 0;
	ARRAY2D(dfdu,15,1) = 0;
	ARRAY2D(dfdu,15,2) = 0;
	ARRAY2D(dfdu,15,3) = 0;
	ARRAY2D(dfdu,15,4) = 0;
	ARRAY2D(dfdu,15,5) = 0;
	ARRAY2D(dfdu,15,6) = 0;
	ARRAY2D(dfdu,15,7) = 0;
	ARRAY2D(dfdu,15,8) = 0;
	ARRAY2D(dfdu,15,9) = 0;
	ARRAY2D(dfdu,15,10) = -G2*kr;
	ARRAY2D(dfdu,15,11) = 0;
	ARRAY2D(dfdu,15,12) = 0;
	ARRAY2D(dfdu,15,13) = 0;
	ARRAY2D(dfdu,15,14) = 0;
	ARRAY2D(dfdu,15,15) = -D1*kr - kur;
	ARRAY2D(dfdu,15,16) = 0;
	ARRAY2D(dfdu,15,17) = 0;
	ARRAY2D(dfdu,15,18) = 0;
	ARRAY2D(dfdu,16,0) = 0;
	ARRAY2D(dfdu,16,1) = 0;
	ARRAY2D(dfdu,16,2) = 0;
	ARRAY2D(dfdu,16,3) = 0;
	ARRAY2D(dfdu,16,4) = 0;
	ARRAY2D(dfdu,16,5) = 0;
	ARRAY2D(dfdu,16,6) = 0;
	ARRAY2D(dfdu,16,7) = 0;
	ARRAY2D(dfdu,16,8) = 0;
	ARRAY2D(dfdu,16,9) = 0;
	ARRAY2D(dfdu,16,10) = -GI*kr;
	ARRAY2D(dfdu,16,11) = 0;
	ARRAY2D(dfdu,16,12) = 0;
	ARRAY2D(dfdu,16,13) = 0;
	ARRAY2D(dfdu,16,14) = 0;
	ARRAY2D(dfdu,16,15) = 0;
	ARRAY2D(dfdu,16,16) = -D1*kr - kur;
	ARRAY2D(dfdu,16,17) = 0;
	ARRAY2D(dfdu,16,18) = 0;
	ARRAY2D(dfdu,17,0) = 0;
	ARRAY2D(dfdu,17,1) = 0;
	ARRAY2D(dfdu,17,2) = 0;
	ARRAY2D(dfdu,17,3) = 0;
	ARRAY2D(dfdu,17,4) = 0;
	ARRAY2D(dfdu,17,5) = 0;
	ARRAY2D(dfdu,17,6) = 0;
	ARRAY2D(dfdu,17,7) = 0;
	ARRAY2D(dfdu,17,8) = 0;
	ARRAY2D(dfdu,17,9) = 0;
	ARRAY2D(dfdu,17,10) = 0;
	ARRAY2D(dfdu,17,11) = -G3*kr;
	ARRAY2D(dfdu,17,12) = 0;
	ARRAY2D(dfdu,17,13) = 0;
	ARRAY2D(dfdu,17,14) = 0;
	ARRAY2D(dfdu,17,15) = 0;
	ARRAY2D(dfdu,17,16) = 0;
	ARRAY2D(dfdu,17,17) = -D2*kr - kur;
	ARRAY2D(dfdu,17,18) = 0;
	ARRAY2D(dfdu,18,0) = 0;
	ARRAY2D(dfdu,18,1) = 0;
	ARRAY2D(dfdu,18,2) = 0;
	ARRAY2D(dfdu,18,3) = 0;
	ARRAY2D(dfdu,18,4) = 0;
	ARRAY2D(dfdu,18,5) = 0;
	ARRAY2D(dfdu,18,6) = 0;
	ARRAY2D(dfdu,18,7) = 0;
	ARRAY2D(dfdu,18,8) = 0;
	ARRAY2D(dfdu,18,9) = 0;
	ARRAY2D(dfdu,18,10) = 0;
	ARRAY2D(dfdu,18,11) = 0;
	ARRAY2D(dfdu,18,12) = 0;
	ARRAY2D(dfdu,18,13) = -G3I*kRr;
	ARRAY2D(dfdu,18,14) = 0;
	ARRAY2D(dfdu,18,15) = 0;
	ARRAY2D(dfdu,18,16) = 0;
	ARRAY2D(dfdu,18,17) = 0;
	ARRAY2D(dfdu,18,18) = -DR*kRr - kRur;
	//Jacobian for parameters

	if (ijac == 1) 
	{
		return 0;
	}


	ARRAY2D(dfdp,0,0) = G1;
	ARRAY2D(dfdp,0,1) = 0;
	ARRAY2D(dfdp,0,2) = 0;
	ARRAY2D(dfdp,0,3) = 0;
	ARRAY2D(dfdp,1,0) = G2;
	ARRAY2D(dfdp,1,1) = 0;
	ARRAY2D(dfdp,1,2) = 0;
	ARRAY2D(dfdp,1,3) = 0;
	ARRAY2D(dfdp,2,0) = G3 - omega*(G3I - gt);
	ARRAY2D(dfdp,2,1) = rm*(-G3I + gt);
	ARRAY2D(dfdp,2,2) = 0;
	ARRAY2D(dfdp,2,3) = 0;
	ARRAY2D(dfdp,3,0) = GI*rmI;
	ARRAY2D(dfdp,3,1) = 0;
	ARRAY2D(dfdp,3,2) = 0;
	ARRAY2D(dfdp,3,3) = 0;
	ARRAY2D(dfdp,4,0) = 0;
	ARRAY2D(dfdp,4,1) = 0;
	ARRAY2D(dfdp,4,2) = 0;
	ARRAY2D(dfdp,4,3) = 0;
	ARRAY2D(dfdp,5,0) = 0;
	ARRAY2D(dfdp,5,1) = 0;
	ARRAY2D(dfdp,5,2) = 0;
	ARRAY2D(dfdp,5,3) = 0;
	ARRAY2D(dfdp,6,0) = 0;
	ARRAY2D(dfdp,6,1) = 0;
	ARRAY2D(dfdp,6,2) = 0;
	ARRAY2D(dfdp,6,3) = 0;
	ARRAY2D(dfdp,7,0) = 0;
	ARRAY2D(dfdp,7,1) = 0;
	ARRAY2D(dfdp,7,2) = 0;
	ARRAY2D(dfdp,7,3) = 0;
	ARRAY2D(dfdp,8,0) = 0;
	ARRAY2D(dfdp,8,1) = 0;
	ARRAY2D(dfdp,8,2) = eta;
	ARRAY2D(dfdp,8,3) = -S*ka;
	ARRAY2D(dfdp,9,0) = 0;
	ARRAY2D(dfdp,9,1) = 0;
	ARRAY2D(dfdp,9,2) = 0;
	ARRAY2D(dfdp,9,3) = S*ka;
	ARRAY2D(dfdp,10,0) = 0;
	ARRAY2D(dfdp,10,1) = 0;
	ARRAY2D(dfdp,10,2) = 0;
	ARRAY2D(dfdp,10,3) = 0;
	ARRAY2D(dfdp,11,0) = 0;
	ARRAY2D(dfdp,11,1) = 0;
	ARRAY2D(dfdp,11,2) = 0;
	ARRAY2D(dfdp,11,3) = 0;
	ARRAY2D(dfdp,12,0) = 0;
	ARRAY2D(dfdp,12,1) = 0;
	ARRAY2D(dfdp,12,2) = 0;
	ARRAY2D(dfdp,12,3) = 0;
	ARRAY2D(dfdp,13,0) = 0;
	ARRAY2D(dfdp,13,1) = 0;
	ARRAY2D(dfdp,13,2) = 0;
	ARRAY2D(dfdp,13,3) = 0;
	ARRAY2D(dfdp,14,0) = 0;
	ARRAY2D(dfdp,14,1) = 0;
	ARRAY2D(dfdp,14,2) = 0;
	ARRAY2D(dfdp,14,3) = 0;
	ARRAY2D(dfdp,15,0) = 0;
	ARRAY2D(dfdp,15,1) = 0;
	ARRAY2D(dfdp,15,2) = 0;
	ARRAY2D(dfdp,15,3) = 0;
	ARRAY2D(dfdp,16,0) = 0;
	ARRAY2D(dfdp,16,1) = 0;
	ARRAY2D(dfdp,16,2) = 0;
	ARRAY2D(dfdp,16,3) = 0;
	ARRAY2D(dfdp,17,0) = 0;
	ARRAY2D(dfdp,17,1) = 0;
	ARRAY2D(dfdp,17,2) = 0;
	ARRAY2D(dfdp,17,3) = 0;
	ARRAY2D(dfdp,18,0) = 0;
	ARRAY2D(dfdp,18,1) = 0;
	ARRAY2D(dfdp,18,2) = 0;
	ARRAY2D(dfdp,18,3) = 0;


	return 0;
}


int stpnt (integer ndim, doublereal t, doublereal *u, doublereal *par)
{

	//init params
	par[0] = 0.0;
	par[1] = 1;
	par[2] = 0;
	par[3] = 50;

	//init variables
	u[0] = 0;
	u[1] = 0;
	u[2] = 0;
	u[3] = 0;
	u[4] = 0;
	u[5] = 0;
	u[6] = 0;
	u[7] = 0;
	u[8] = 0;
	u[9] = 0;
	u[10] = 0;
	u[11] = 0;
	u[12] = 0;
	u[13] = 0;
	u[14] = 0;
	u[15] = 0;
	u[16] = 0;
	u[17] = 0;
	u[18] = 0;
	return 0;
}

int pvls (integer ndim, const doublereal *u, doublereal *par)
{ return 0;}


int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc)
{ return 0;}


int icnd (integer ndim, const doublereal *par, const integer *icp,
	integer nint, const doublereal *u, const doublereal *uold,
	const doublereal *udot, const doublereal *upold, integer ijac,
	doublereal *fi, doublereal *dint)
{ return 0;}


int fopt (integer ndim, const doublereal *u, const integer *icp,
	const doublereal *par, integer ijac,
	doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{ return 0; }
