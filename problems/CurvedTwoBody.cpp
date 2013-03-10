#include "CurvedTwoBody.h"

#include <math.h>
#include <stdexcept>
#include <iostream>
using namespace std;

#define CUBE(x) ((x)*(x)*(x))
#define SQUARE(x) ((x)*(x))
#define SQR(x) SQUARE(x)

#define Ith(v, i) (NV_Ith_S((v), (i)-1))
#define IJth(v, i, j) (DENSE_ELEM(v, (i)-1, (j)-1))

#define CAPTURE 101
#define X 0
#define Y 1
#define Z 2

int CurvedTwoBody::rhsFunc( 
		realtype t, N_Vector y, N_Vector ydot, void *myself)
{
	CurvedTwoBody* me = (CurvedTwoBody*) myself;
	
	static realtype m1, m2, R;
	
	m1 = me->_mass_1;
	m2 = me->_mass_2;
	R = me->_sphere_radius;
	#define G 1

	static realtype theta1, theta2, phi1, phi2, theta1prime, theta2prime, 
					phi1prime, phi2prime;

	theta1 = Ith(y,1);
	theta2 = Ith(y,2);
	phi1 = Ith(y,3);
	phi2 = Ith(y,4);
	theta1prime = Ith(y,5);
	theta2prime = Ith(y,6);
	phi1prime = Ith(y,7);
	phi2prime = Ith(y,8);
 
 //calculate unit vectors
 	static realtype r1Hat[3];
    r1Hat[0] = sin(theta1) * cos(phi1);
    r1Hat[1] = sin(theta1) * sin(phi1);
    r1Hat[2] = cos(theta1);
    
    static realtype r2Hat[3];
    r2Hat[0] = sin(theta2) * cos(phi2);
    r2Hat[1] = sin(theta2) * sin(phi2);
    r2Hat[2] = cos(theta2);
 
    static realtype theta1Hat[3];
    theta1Hat[0] = cos(theta1) * cos(phi1);
    theta1Hat[1] = cos(theta1) * sin(phi1);
    theta1Hat[2] = -sin(theta1);
    
    static realtype theta2Hat[3];
    theta2Hat[0] = cos(theta2) * cos(phi2);
    theta2Hat[1] = cos(theta2) * sin(phi2);
    theta2Hat[2] = -sin(theta2);
    
    static realtype phi1Hat[3];
    phi1Hat[0] = -sin(phi1);
    phi1Hat[1] = cos(phi1);
    phi1Hat[2] = 0;
    
    static realtype phi2Hat[3];
    phi2Hat[0] = -sin(phi2);
    phi2Hat[1] = cos(phi2);
    phi2Hat[2] = 0;

    // calculate forces
    static realtype r1Hat_r2Hat;
	r1Hat_r2Hat = r1Hat[0] * r2Hat[0] + r1Hat[1] * r2Hat[1] + r1Hat[2] * r2Hat[2];	//dot product of r1Hat and r2Hat
    static realtype delta;
	delta = acos(r1Hat_r2Hat);	// delta is angle between particles
	
	if (R*delta < me->_capture_radius) {
		me->setError(CAPTURE, "Bodies Collided.");
		return -1;
	}

	static realtype dShort, dLong, fMagnitude;
    dShort = R*delta;		//short distance between the two particles
    dLong = R*(2*M_PI - delta);	//long distance between the two particles
    fMagnitude = G*m1*m2 * (1/SQR(dShort) - 1/SQR(dLong));	//magnitude of the gravitational force with wrap around

    static realtype p1Vec[3]; //Gram - Schmidt process
    p1Vec[X] = R*r2Hat[X] - R*r1Hat_r2Hat * r1Hat[X];
    p1Vec[Y] = R*r2Hat[Y] - R*r1Hat_r2Hat * r1Hat[Y];
    p1Vec[Z] = R*r2Hat[Z] - R*r1Hat_r2Hat * r1Hat[Z];
	
	static realtype p2Vec[3];
    p2Vec[X] = R*r1Hat[X] - R*r1Hat_r2Hat * r2Hat[X];
    p2Vec[Y] = R*r1Hat[Y] - R*r1Hat_r2Hat * r2Hat[Y];
    p2Vec[Z] = R*r1Hat[Z] - R*r1Hat_r2Hat * r2Hat[Z];
	
	static realtype p1Mag, p2Mag;
    p1Mag = sqrt (SQR(p1Vec[X]) +SQR(p1Vec[Y])+ SQR(p1Vec[Z]));
	p2Mag = sqrt (SQR(p2Vec[X]) +SQR(p2Vec[Y])+ SQR(p2Vec[Z]));    
    
	static realtype f1Hat[3];
    f1Hat[X] = p1Vec[X] / p1Mag;
    f1Hat[Y] = p1Vec[Y]/p1Mag;
    f1Hat[Z] = p1Vec[Z]/p1Mag;
    
	static realtype f2Hat[3]; // force components in theta direction
    f2Hat[X] = p2Vec[X] / p2Mag;
    f2Hat[Y] = p2Vec[Y]/p2Mag;
    f2Hat[Z] = p2Vec[Z]/p2Mag;
    
	static realtype f1Vec[3];
    f1Vec[X] = f1Hat[X] * fMagnitude;
    f1Vec[Y] = f1Hat[Y] * fMagnitude;
    f1Vec[Z] = f1Hat[Z] * fMagnitude;

    static realtype f2Vec[3];
    f2Vec[X] = f2Hat[X] * fMagnitude;
    f2Vec[Y] = f2Hat[Y] * fMagnitude;
    f2Vec[Z] = f2Hat[Z] * fMagnitude;
	
	static realtype aVec1[3];
	aVec1[X] = f1Vec[X]/m1;
	aVec1[Y] = f1Vec[Y]/m1;
	aVec1[Z] = f1Vec[Z]/m1;
	
	static realtype aVec2[3];
	aVec2[X] = f2Vec[X]/m2;
	aVec2[Y] = f2Vec[Y]/m2;
	aVec2[Z] = f2Vec[Z]/m2;
	
	static realtype a1MechLin[2]; //linear mechanical accelerations projections
	a1MechLin[0] = aVec1[X]*theta1Hat[X] + aVec1[Y]*theta1Hat[Y] + aVec1[Z]*theta1Hat[Z]; //acceleration in thetaHat
	a1MechLin[1] = aVec1[X]*phi1Hat[X] + aVec1[Y]*phi1Hat[Y] + aVec1[Z]*phi1Hat[Z];	//acceleration in phiHat
	
	static realtype a2MechLin[2]; //linear mechanical accelerations projections
	a2MechLin[0] = aVec2[X]*theta2Hat[X] + aVec2[Y]*theta2Hat[Y] + aVec2[Z]*theta2Hat[Z]; //acceleration in thetaHat
	a2MechLin[1] = aVec2[X]*phi2Hat[X] + aVec2[Y]*phi2Hat[Y] + aVec2[Z]*phi2Hat[Z];	//acceleration in phiHat

    static realtype a1Mechanical[2]; //angular mechanical acceleration projections
	a1Mechanical[0] = a1MechLin[0]/R; //acceleration in thetaHat
	a1Mechanical[1] = a1MechLin[1]/(R*sin(theta1));	//acceleration in phiHat

    static realtype a2Mechanical[2]; // angular mechanical acceleration projections
	a2Mechanical[0] = a2MechLin[0]/R; //acceleration in thetaHat
	a2Mechanical[1] = a2MechLin[1]/(R*sin(theta2));	//acceleration in phiHat         
    
    static realtype a1Curvature[2]; // effective acceleration due to curvature
    a1Curvature[0] = cos(theta1)*sin(theta1)*SQR(phi1prime);
    a1Curvature[1] = -2*(1/tan(theta1))*theta1prime*phi1prime;
    
    static realtype a2Curvature[2]; // effective acceleration due to curvature
    a2Curvature[0] = cos(theta2)*sin(theta2)*SQR(phi2prime);
    a2Curvature[1] = -2*(1/tan(theta2))*theta2prime*phi2prime;
    
    static realtype aTheta1, aPhi1;
	aTheta1 = a1Mechanical[0] + a1Curvature[0];
    aPhi1 = a1Mechanical[1] + a1Curvature[1];

	static realtype aTheta2, aPhi2;
    aTheta2 = a2Mechanical[0] + a2Curvature[0];    
    aPhi2 = a2Mechanical[1] + a2Curvature[1];
	
	Ith(ydot,1) = Ith(y,5);
	Ith(ydot,2) = Ith(y,6);	
	Ith(ydot,3) = Ith(y,7);	
	Ith(ydot,4) = Ith(y,8);	
	Ith(ydot,5) = aTheta1;	
	Ith(ydot,6) = aTheta2;	
	Ith(ydot,7) = aPhi1;	
	Ith(ydot,8) = aPhi2;	
	return 0; 
}

int CurvedTwoBody::denseJac( long int N, DenseMat J, realtype t, 
	N_Vector y, N_Vector fy, void *jac_data, N_Vector tmp1, 
	N_Vector tmp2, N_Vector tmp3 )
{
	CurvedTwoBody * me = (CurvedTwoBody *) jac_data;
	
	static realtype c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, 
					c14, c15, c16, c17, c18, c19, c20, c21, c22;
	
	static realtype s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, 
					s14, s15, s16, s17, s18, s19, s20, s21, s22;
	
	static realtype ac1, ct1, ct2, cc1, cc2;
	
	static realtype p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18;
	static realtype b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11;
	static realtype d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16;
	static realtype e1, e2, e3, e4, e5, e6, e7;
	static realtype f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, // f13, 
					f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, 
					f25, f26, f27, f28, f29, f30, f31, f32, f33, f34, f35, 
					f36, f37, f38, f39, f40, f41, f42, f43;
	
	static realtype g1, g2, g3, g4, g5, g6;
	
	static realtype h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11;

	static realtype sqrt2;
	sqrt2 = sqrt(RCONST(2.0));
	
	static realtype m1, m2, R;
	m1 = me->_mass_1;
	m2 = me->_mass_2;
	R = me->_sphere_radius;
	#define G 1

	static realtype theta1, theta2, phi1, phi2, theta1prime, theta2prime, 
					phi1prime, phi2prime;

	theta1 = Ith(y,1);		//theta1
	theta2 = Ith(y,2);		//theta2
	phi1 = Ith(y,3);			//phi1
	phi2 = Ith(y,4);			//phi2
	theta1prime = Ith(y,5);	//theta1 prime
	theta2prime = Ith(y,6);	//theta2 prime
	phi1prime = Ith(y,7);	//phi1 prime
	phi2prime = Ith(y,8);	//phi2 prime

	c1 = cos(theta1);
	c2 = cos(theta2);
	c3 = cos(phi1);
	c4 = cos(phi2);
	c5 = cos(phi1 - phi2);
	c6 = cos(2*theta1);
	c7 = cos(2*(theta1 - theta2));
	c8 = cos(2*theta2);
	c9 = cos(2*(theta1 + theta2));
	c10 = cos(2*(phi1 - phi2));
	c11 = cos(2*(theta1 + phi1 - phi2));
	c12 = cos(2*theta1 - 2*theta2 + phi1 - phi2);
	c13 = cos(2*(theta1 - theta2 + phi1 - phi2));
	c14 = cos(2*(theta2 + phi1 - phi2));
	c15 = cos(2*(theta1 + theta2 + phi1 - phi2));
	c16 = cos(2*theta1 + 2*theta2 + phi1 - phi2);
	c17 = cos(2*(theta1 - phi1 + phi2));
	c18 = cos(2*theta1 - 2*theta2 - phi1 + phi2);
	c19 = cos(2*(theta1 - theta2 - phi1 + phi2));
	c20 = cos(2*(theta2 - phi1 + phi2));
	c21 = cos(2*(theta1 + theta2 - phi1 + phi2));
	c22 = cos(2*theta1 + 2*theta2 - phi1 + phi2);
	
	s1 = sin(theta1);
	s2 = sin(theta2);
	s3 = sin(phi1);
	s4 = sin(phi2);
	s5 = sin(2*theta1);
	s6 = sin(2*(theta1 - theta2));
	s7 = sin(2*(theta1 + theta2));
	s8 = sin(2*(theta1 + phi1 - phi2));
	s9 = sin(2*theta1 - 2*theta2 + phi1 - phi2);
	s10 = sin(2*(theta1 - theta2 + phi1 - phi2));
	s11 = sin(2*(theta1 + theta2 + phi1 - phi2));
	s12 = sin(2*theta1 + 2*theta2 + phi1 - phi2);
	s13 = sin(2*(theta1 - phi1 + phi2));
	s14 = sin(2*theta1 - 2*theta2 - phi1 + phi2);
	s15 = sin(2*(theta1 - theta2 - phi1 + phi2));
	s16 = sin(2*(theta1 + theta2 - phi1 + phi2));
	s17 = sin(2*theta1 + 2*theta2 - phi1 + phi2);
	s18 = sin(2*theta2);
	s19 = sin(2*(theta2 + phi1 - phi2));
	s20 = sin(2*(theta2 - phi1 + phi2));
	s21 = sin(phi1 - phi2);
	s22 = sin(2*(phi1 - phi2));
	
	ac1 = acos(c1*c2 + c5*s1*s2);
	
	ct1 = 1/(tan(theta1));
	ct2 = 1/(tan(theta2));
	
	cc1 = 1/(sin(theta1));
	cc2 = 1/(sin(theta2));
	
	p1 = SQR(R);
	p2 = SQR(ac1);
	p3 = SQR(s1);
	p4 = SQR(s2);
	p5 = SQR(s21);
	p6 = SQR(phi1prime);
	p7 = SQR(phi2prime);
	
	b1 = c2*s1 - c1*c5*s2;
	b2 = c1*c2 + c5*s1*s2;
	b3 = -ac1+M_PI;
	b4 = ac1-2*M_PI;
	b5 = -20 + 4*c10 - 2*c11 + 4*c12 + c13 - 2*c14 + c15 - 4*c16 - 2*c17 + 4*c18 + c19 - 2*c20 + c21 - 4*c22 + 4*c6 + 6*c7 + 4*c8 + 6*c9;
	b6 = c2*c5*s1 - c1*s2;
	b7 = -2*ct2;
	b8 = -c1*c2*c5 - s1*s2;
	b9 = -2*s10 - 2*s11 + 8*s12 + 4*s13 - 8*s14 - 2*s15 - 2*s16 + 8*s17 - 8*s5 - 12*s6 - 12*s7 + 4*s8 - 8*s9;
	b10 = -c3*p3*R*s2*s21 - R*s2*s4;
	b11 = -c4*R*s2 + p3*R*s2*s21*s3;
	
	p8 = SQR(b1);
	p9 = SQR(b2);
	p10 = SQR(b4);
	p16 = SQR(b6);
	p11 = SQR(p1);
	
	d1 = -b2*c1*c3*R + b1*c3*R*s1;
	d2 = -b2*c3*R*s1 + c4*R*s2;
	d3 = b2*c1*R*s3 - b1*R*s1*s3;
	d4 = b2*R*s1*s3 - R*s2*s4;
	d5 = -8*b1*b3*G*m2*M_PI + (p1*sqrt(-b5*p1)*p10*p2*p6*s5)/(4*sqrt2);
	d6 = c3*R*s1 - b2*c4*R*s2;
	d7 = 8*b3*b6*G*m1*M_PI + (p1*sqrt(-b5*p1)*p10*p2*p7*s18)/(4*sqrt2);
	d8 = R*s1*s3 - b2*R*s2*s4;
	d9 = c1*c3*R + b1*c4*R*s2;
	d10 = -b2*c2*c4*R - b6*c4*R*s2;
	d11 = c2*c4*R - b6*c3*R*s1;
	d12 = b6*R*s1*s3 - c2*R*s4;
	d13 = c1*R*s3 + b1*R*s2*s4;
	d14 = c3*p3*R*s2*s21 + b2*R*s1*s3;
	d15 = -c4*p4*R*s1*s21 + b2*R*s2*s4;
	d16 = -b2*c4*R*s2 - p4*R*s1*s21*s4;
	
	p12 = SQR(d2);
	p13 = SQR(d4);
	p14 = SQR(d6);
	p15 = SQR(d8);
	p17 = (ac1*p2);
	p18 = (b4*p10);
	
	e1 = p12 + p13 + p1*p3*p8;
	e2 = ac1*p10 + b4*p2;
	e3 = sqrt(-b5*p1);
	e4 = sqrt(2 - 2*p9);
	e5 = sqrt(1 - p9);
	e6 = p10*p2;
	e7 = p14 + p15 + p1*p16*p4;
	
	f1 = sqrt(e1*e1*e1)*p1;
	f2 = sqrt(e7*e7*e7)*p1;
	f3 = 4*sqrt2*e3;
	f4 = s10 + s11 - 2*s12 + 2*s13 - 2*s14 - s15 - s16 + 2*s17 - 2*s19 + 2*s20 + 4*s22 - 2*s8 + 2*s9;
	f5 = -s10 + s11 - 4*s12 - 4*s14 - s15 + s16 - 4*s17 + 4*s18 - 2*s19 - 2*s20 - 6*s6 + 6*s7 - 4*s9;
	f6 = b1*c1*p1 - d4*R*s3;
	f7 = c4*d6*p4*s21 - d6*s3 + d8*p4*s21*s4;
	f8 = e5*p10*p17;
	f9 = e5*p18*p2;
	f10 = e5*e6*p1;
	f11 = c3*d8*R + f7*R - b6*c2*p1*p4*s21;
	f12 = sqrt(e7*e7*e7)*f10*p17*p18;
	//f13 = e5*e6*f1*p17*p18;
	f14 = 4*b3*e1*p4*p5;
	f15 = -f14/f8 - f14/f9 - (2*e1*p4*p5)/(e5*e6) + (2*b3*c5*cc1*e1*s2)/e6;
	f16 = (e6*f4*p11*p7*s18)/f3 - 8*b3*c2*G*m1*M_PI*s1*s21 - (8*b6*G*m1*M_PI*s1*s2*s21)/e5 + (ac1*e3*p1*p10*p7*s1*s18*s2*s21)/(2*e4) + (b4*e3*p1*p2*p7*s1*s18*s2*s21)/(2*e4);
	f17 = 4*b1*d7*e7;
	f18 = 4*b1*d5*e1;
	f19 = 4*b6*d5*e1;
	f20 = 4*b6*d7*e7;
	f21 = 4*d5*e1*s1*s2*s21;
	f22 = d1*d2 + d3*d4;
	f23 = 2*e7*f16;
	f24 = 4*d7*e7*s1*s2*s21;
	f25 = p10*p17 + p18*p2;
	f26 = 8*b1*G*m2*M_PI*s1*s2*s21;
	f27 = 8*b3*c1*G*m2*M_PI*s2*s21;
	f28 = e6*f4*p11*p6*s5;
	f29 = ac1*e3*p1*p10*p6*s1*s2*s21*s5;
	f30 = b4*e3*p1*p2*p6*s1*s2*s21*s5;
	f31 = 2*e7*p3*p5;
	f32 = 4*b3*e7*p3*p5;
	f33 = 2*b3*c5*cc2*e7*s1;
	f34 = 2*b3*cc1*s2*s21;
	f35 = e5*p17*p18;
	f36 = d11*d2 + d12*d4 + b1*b8*p1*p3;
	f37 = b10*d2 + b11*d4 - b1*c1*p1*p3*s2*s21;
	f38 = d15*d6 + d16*d8 + b6*c2*p1*p4*s1*s21;
	f39 = d14*d2 + b2*c3*d4*R*s1 + f6*p3*s2*s21;
	f40 = ac1*b6*e3*p1*p10;
	f41 = e6*f5*p11;
	f42 = b4*b6*e3*p1*p2;
	f43 = b2*e5*p1*p17;
	
	g1 = e5*f28 + f26*f3 - e5*f27*f3;
	g2 = e5*f3*(f29 + f30) + 2*e4*g1;
	g3 = e3*e5*e6;
	g4 = e6*f1*f35;
	g5 = p17*p18;
	g6 = b5*SQR(p1);
	
	h1 = f8 + f9;
	h2 = 4*e6*f2*f8*f9;
	h3 = e5*e6*f2*f8*f9;
	h4 = f8*f9*g3;
	h5 = 4*e4*e5*e6*f1*f3*f8*f9;
	h6 = f1*h4;
	h7 = f2*h4;
	h8 = e1*f8*f9*g2;
	h9 = e6*f24*h1;
	h10 = f31*f8*f9;
	h11 = e6*f32*h1;

	IJth(J,1,1) = 0;
	IJth(J,1,2) = 0;
	IJth(J,1,3) = 0;
	IJth(J,1,4) = 0;
	IJth(J,1,5) = 1;
	IJth(J,1,6) = 0;
	IJth(J,1,7) = 0;
	IJth(J,1,8) = 0;

	IJth(J,2,1) = 0;
	IJth(J,2,2) = 0;
	IJth(J,2,3) = 0;
	IJth(J,2,4) = 0;
	IJth(J,2,5) = 0;
	IJth(J,2,6) = 1;
	IJth(J,2,7) = 0;
	IJth(J,2,8) = 0;

	IJth(J,3,1) = 0;
	IJth(J,3,2) = 0;
	IJth(J,3,3) = 0;
	IJth(J,3,4) = 0;
	IJth(J,3,5) = 0;
	IJth(J,3,6) = 0;
	IJth(J,3,7) = 1;
	IJth(J,3,8) = 0;

	IJth(J,4,1) = 0;
	IJth(J,4,2) = 0;
	IJth(J,4,3) = 0;
	IJth(J,4,4) = 0;
	IJth(J,4,5) = 0;
	IJth(J,4,6) = 0;
	IJth(J,4,7) = 0;
	IJth(J,4,8) = 1;

	IJth(J,5,1) = -(16*d5*h4*(f22 + b1*b2*p1*p3 + c1*p1*p8*s1) + e1*SQR(e6)*f8*f9*p6*(4*sqrt2*c6*e5*g6 + b9*e4*p11*s5) + 4*e6*(2*f18*g3*h1 + e1*f8*f9*(32*b2*b3*e3*e5*G*m2*M_PI - 32*e3*G*m2*p8*M_PI + sqrt2*b1*e2*g6*p6*s5)))/(32*e6*h6);
	IJth(J,5,2) =(-(e1*e5*f3*(f40 + f42)*f8*f9*p6*s5) + e4*(-16*b1*b6*e1*f3*f8*f9*G*m2*M_PI + e5*(e6*f19*f3*h1 - 2*f8*f9*(d5*f3*f36 + 8*b3*b8*e1*f3*G*m2*M_PI - e1*f41*p6*s5))))/h5;
	IJth(J,5,3) = (-(e4*e5*f3*(2*d5*f39*f8*f9 + e6*f21*h1)) + h8)/h5;
	IJth(J,5,4) = (e4*e5*f3*(-2*d5*f37*f8*f9 + e6*f21*h1) - h8)/h5;
	IJth(J,5,5) = 0;
	IJth(J,5,6) = 0;
	IJth(J,5,7) = (e3*s5*phi1prime)/(4*sqrt2*sqrt(e1));
	IJth(J,5,8) = 0;

	IJth(J,6,1) = -(8*e3*(16*b1*b6*e7*f8*f9*G*m1*M_PI + e5*(e6*f17*h1 + 2*f8*f9*(d13*d7*d8 + d6*d7*d9 - b6*b8*d7*p1*p4 + 8*b3*b8*e7*G*m1*M_PI))) + sqrt2*e7*f8*f9*(4*b1*e2*g6 + b9*e5*e6*p11)*p7*s18)/(32*h7);
	IJth(J,6,2) = (2*sqrt2*c8*SQR(e3)*e7*f10*f8*f9*p7 + sqrt2*e5*e7*f41*f8*f9*p7*s18 + e3*(2*e7*f8*f9*(32*G*m1*p16*M_PI - sqrt2*(f40 + f42)*p7*s18) + 4*e5*(e6*f20*h1 + 2*f8*f9*(-(d10*d6*d7) + d7*s2*(-(c2*p1*p16) + b6*d8*R*s4) + b2*(b6*d7*p1*p4 - 8*b3*e7*G*m1*M_PI + c2*d7*d8*R*s4)))))/(16*h7);
	IJth(J,6,3) = (-h9 + f8*f9*(f23 - 2*d7*f11*s1))/h2;
	IJth(J,6,4) = (-((f23 + 2*d7*f38)*f8*f9) + h9)/h2;
	IJth(J,6,5) = 0;
	IJth(J,6,6) = 0;
	IJth(J,6,7) = 0;
	IJth(J,6,8) = (e3*s18*phi2prime)/(4*sqrt2*sqrt(e7));

	IJth(J,7,1) = (2*cc1*(2*b1*e1*G*g5*m2*M_PI*s2*s21 + 2*b3*G*m2*M_PI*(ct1*e1*f35 + f22*f35 + 2*b1*e1*p10*p17 + 2*b1*p12*p18*p2 + 2*b1*p13*p18*p2 + b1*f43*p18*p3 + 2*b1*p1*p18*p2*p3*p8 + c1*f35*p1*p8*s1)*s2*s21 + cc1*g4*theta1prime*phi1prime))/g4;
	IJth(J,7,2) = (-4*cc1*G*m2*M_PI*(b6*e1*g5*s2 + b3*(c2*e1*f35 + 2*b6*e1*f25*s2 - f35*f36*s2))*s21)/g4;
	IJth(J,7,3) = (2*(-(e6*f15) + f34*f39)*G*m2*M_PI)/(e6*f1);
	IJth(J,7,4) = (2*(e6*f15 + f34*f37)*G*m2*M_PI)/(e6*f1);
	IJth(J,7,5) = -2*ct1*phi1prime;
	IJth(J,7,6) = 0;
	IJth(J,7,7) = -2*ct1*theta1prime;
	IJth(J,7,8) = 0;

	IJth(J,8,1) = (-4*cc2*G*m1*M_PI*(b1*e7*g5*s1 + b3*(-(c1*e7*f35) + (2*b1*e7*f25 + d13*d8*f35 + d6*d9*f35 - b6*b8*f35*p1*p4)*s1))*s21)/f12;
	IJth(J,8,2) = (2*cc2*(2*b3*G*m1*M_PI*s1*s21*(-(d10*d6*f35) - ct2*e7*f35 + 2*b6*e7*p10*p17 + 2*b6*p14*p18*p2 + 2*b6*p15*p18*p2 + b6*f43*p18*p4 + 2*b6*p1*p16*p18*p2*p4 - c2*f35*p1*p16*s2 + b2*c2*d8*f35*R*s4 + b6*d8*e5*g5*R*s2*s4) + e7*g5*(2*b6*G*m1*M_PI*s1*s21 + cc2*sqrt(e7)*f10*theta2prime*phi2prime)))/f12;
	IJth(J,8,3) = (-2*G*m1*M_PI*(h10 + e5*(-(f33*f8*f9) + h11 + 2*b3*cc2*f11*f8*f9*p3*s21)))/h3;
	IJth(J,8,4) = (2*G*m1*M_PI*(h10 + e5*(h11 - f8*f9*(f33 + 2*b3*cc2*f38*s1*s21))))/h3;
	IJth(J,8,5) = 0;
	IJth(J,8,6) = b7*phi2prime;
	IJth(J,8,7) = 0;
	IJth(J,8,8) = b7*theta2prime;
	
	return 0; 
}

bool CurvedTwoBody::handleOption(int arg)
{
	switch (arg)
	{
		case 'P': // initial separation (\rho_0) 
			setInitialSeparation(atof(optarg));
			return true;
			
		case 'R': //sphere size
			setSphereRadius(atof(optarg));
			return true;

		case 'V': // difference in initial velocities of particle one and two
			setDeltaV(atof(optarg));
			return true;

		case 'c': // capture radius
			setCaptureRadius(atof(optarg));
			return true;

		case 'f': // mass fraction
			setMassFraction(atof(optarg));
			return true;
	}
}

void CurvedTwoBody::setInitialSeparation(realtype r)
{
	_initial_separation = r;
	updateInitialVector();
}

void CurvedTwoBody::setDeltaV(realtype s)
{
	_delta_v = s;
	updateInitialVector();
}

void CurvedTwoBody::setMassFraction(realtype f)
{
	if (f < 0 || f > 1)
		throw invalid_argument("CurvedTwoBody::setMassFraction(): mass fraction must be between 0 and 1.");
	
	_mass_1 = f;
	_mass_2 = 1 - f;
}

/**
 * Update the initial velocity vector to reflect the initial conditions
 * specified by _delta_v and _initial_separation.
 */
void CurvedTwoBody::updateInitialVector() 
{
	vector<realtype> v(getVectorSize());
	
	// v = (q1, q2, j1, j2, q1', q2', j1', j2')
	v.at(0) = (_initial_separation*_mass_2)/(_sphere_radius*(_mass_1+_mass_2));	// theta1 = theta1 initial
	v.at(1) = (_initial_separation*_mass_1)/(_sphere_radius*(_mass_1+_mass_2));		// theta2 = theta2 initial
	v.at(2) = 0;
	v.at(3) = M_PI;						// phi2 = phi2 initial

	v.at(4) = 0;
	v.at(5) = 0;
	v.at(6) = (_delta_v*_mass_2)/((_mass_1+_mass_2)*_sphere_radius*sin(v.at(0)));			// phi1prime = phi1prime initial
	v.at(7) = (_delta_v*_mass_1)/((_mass_1+_mass_2)*_sphere_radius*sin(v.at(1)));			// phi2prime = phi2prime initial

	setInitialVector(v);
}

void CurvedTwoBody::setCaptureRadius(realtype r)
{
	_capture_radius = r;
}

void CurvedTwoBody::setSphereRadius(realtype r)
{
	_sphere_radius = r;
	updateInitialVector();
}

void CurvedTwoBody::recordDataInit()
{
	cout << "t\tj1\tq1\tj2\tq2" << endl; 
}

void CurvedTwoBody::recordData(realtype t, N_Vector y)
{
	cout << t << "\t"
			<< Ith(y,3) << "\t"
			<< Ith(y,1) << "\t"
			<< Ith(y,4) << "\t"
			<< Ith(y,2)
			<< endl;
}
	
