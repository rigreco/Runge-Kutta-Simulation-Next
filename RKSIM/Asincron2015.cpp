/*  Runge Kutta simulator for Induction Motors
    This project is based on 1999-2000 Thesis work of Greco Riccardo.

	Asynchronous v1.0
	rksim.exe

    Copyright (C) 2015  Riccardo Greco rigreco.grc@gmail.com.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>

//#define rkfr 1.0e-3		/* RKF max resolution h<rk4r */

void gotoxy(int x, int y)
{
	COORD coord;
	coord.X = x;
	coord.Y = y;
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord);
}


//Global var

double p, rr, ls, rs, lm, lr, sn, gf, cr, j, vsr, vsi, iss, irs, coppia, e;

//

void f(double x, double y[6], double eq[6]);

// RUNGE KUTTA FHELBERG PARAMETERS

double acca[5] = { 1. / 4., 3. / 8., 12. / 13., 1., 1.2 };
double erre1[5] = { 1. / 4., 0., 0., 0., 0. };
double erre2[5] = { 3. / 32., 9. / 32., 0., 0., 0. };
double erre3[5] = { 1932. / 2197., -7200. / 2197., 7296. / 2197., 0., 0. };
double erre4[5] = { 439. / 216., -8., 3680. / 513., -845. / 4104., 0. };
double erre5[5] = { -8. / 27., 2., -3544. / 2565., 1859. / 4104., -11. / 40. };

double ips5[5] = { 25. / 216., 0., 1408. / 2565., 2197. / 4104., -1. / 5. };
double ips6[6] = { 16. / 135., 0., 6656. / 12825., 28561. / 56430., -9. / 50., 2. / 55. };

// ****************


void main(void)
{
	
	int xp, yp;
	char buf[20];

	double x, h, r11, r21, r31, r41, r51, r61,
		r12, r22, r32, r42, r52, r62, r13, r14, r15, r16,
		r23, r24, r25, r26, r33, r34, r35, r36, r43, r44,
		r45, r46, r53, r54, r55, r56, r63, r64, r65, r66;
	
	/* inizialization */
	double eq[6] = { 0., 0., 0., 0., 0., 0. };
	double y[7] = { 0., 0., 0., 0., 0., 0. };
	double y1[7] = { 0., 0., 0., 0., 0., 0. };
	//double err = 0.0;
	double time = 0.0;

	/* IM INPUT PARAMETERS */
	_cputs("Tri-phase asynchronous motor simulation by RKF4(5) method - Greco Riccardo \n\n");
	
	_cputs("Induction Motor Parameters \n\n");
	_cputs("Enter to confir defaut values \n");
	
	_cputs("Stator Voltage e [380.0 V] ");
	e = (double)atof(gets_s(buf));
	if (e == 0.0) { e = 380.0; }
	
	_cputs("Stator resistence rs [1.078 Ohm] ");
	rs = (double)atof(gets_s(buf));
	if (rs == 0.0) {rs = 1.078;}
	// printf("Valore=  %lf \n", rs);
	
	_cputs("Rotor resistence rr [0.898 Ohm] ");
	rr = (double)atof(gets_s(buf));
	if (rr == 0.0) { rr = 0.898; }

	_cputs("Stator inductance (ls=lr) [0.0093 H] ");
	ls = (double)atof(gets_s(buf));
	if (ls == 0.0) { ls = 0.0093; }
	
	lr = ls;

	_cputs("Magnetization inductance lm [0.236 H] ");
	lm = (double)atof(gets_s(buf));
	if (lm == 0.0) { lm = 0.236; }

	_cputs("Poles p [2.0] ");
	p = (double)atof(gets_s(buf));
	if (p == 0.0) { p = 2.0; }

	_cputs("Nominal slip sn [0.0257] ");
	sn = (double)atof(gets_s(buf));
	if (sn == 0.0) { sn = 0.0257; }
	
	_cputs("Inertia moment j [0.1 Kg*m2] ");
	j = (double)atof(gets_s(buf));
	if (j == 0.0) { j = 0.1; }

	// Charge

	_cputs("Load torque cr [0.0 N*m] ");
	cr = (double)atof(gets_s(buf));
	if (cr == 0.0) { cr = 0.0; }

	_cputs("Time simulation [1.0 s] ");
	time = (double)atof(gets_s(buf));
	if (time == 0.0) { time = 1.0; }

	_cputs("Press any key to start plotting ...");
	_getch();

	system("cls");

	/************* Setting Console ************/
	
	SetConsoleTitle("Asynchronous by Greco Riccardo"); // Set text of the console so you can find the window
	
	HWND hwnd = FindWindow(NULL, "Asynchronous by Greco Riccardo"); // Get the HWND The Main
	MoveWindow(hwnd, 50, 50, 800, 600, FALSE);
	HDC hdc = GetDC(hwnd);

	SetGraphicsMode(hdc, GM_COMPATIBLE);
	
	HANDLE  hConsole;
	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	/*****************************************/

	// Inizialization

	x = 0.0;
	
	h = 1e-4; /*Fix step integration*/

	
	// *********** Axes *****
		
	for (xp = 0; xp < 800; xp++)
	{
		SetPixel(hdc, xp, 300, RGB(255, 255, 255));
	}
	for (yp = 0; yp < 600; yp++)
	{
		SetPixel(hdc, 0, yp, RGB(255, 255, 255));
	}
	
	
	while (x<time)
	{
		// Currents symmatrical components

		iss = sqrt(y[0] * y[0] + y[1] * y[1]); /*Symmetric Stator current iss=sqrt(isr^2+isi^2) module of Is*/
		irs = sqrt(y[2] * y[2] + y[3] * y[3]); /*Symmetric Rotor current irs=sqrt(irr^2+iri^2) module of Ir*/
		
		/* Stator voltage */

		vsr = (e*M_SQRT2)*cos(2.0*M_PI*50.0*x); /*Real component of statoric Vs*/
		vsi = (e*M_SQRT2)*sin(2.0*M_PI*50.0*x); /*Immaginary component of statoric Vs*/

		// Electromagnetic torque (motor torque)

		coppia = 3. / 2.*lm*p*(y[1] * y[2] - y[0] * y[3]); /*(3/2)*lm*p*(isi*irr-isr*iri)*/
		
		
		/**************** Out to screen ******************/

		SetConsoleTextAttribute(hConsole, 4);
		gotoxy(30, 12);
		printf("Iss=  %f", iss);

		SetConsoleTextAttribute(hConsole, 5);
		gotoxy(30, 13);
		printf("Irs=  %f", irs);
		
		SetConsoleTextAttribute(hConsole, 1);
		gotoxy(30, 17);
		printf("Torque=  %f", coppia);
		
		SetConsoleTextAttribute(hConsole, 2);
		gotoxy(30, 10);
		printf("Pos=  %f", y[5]);
		
		SetConsoleTextAttribute(hConsole, 15);
		gotoxy(30, 8);
		printf("W= %f", y[4]);
		
		gotoxy(30, 20);
		printf("t= %f", x);

		/**********************************************/
		
		
		/* Begin 
		do
		{
		/* End */
		
		// RUNGE KUTTA

		// 0ø order

		f(x, y, eq);

		r11 = eq[0];
		r12 = eq[1];
		r13 = eq[2];
		r14 = eq[3];
		r15 = eq[4];
		r16 = eq[5];

		// 1ø order

		y1[0] = y[0] + erre1[0] * r11*h;
		y1[1] = y[1] + erre1[0] * r12*h;
		y1[2] = y[2] + erre1[0] * r13*h;
		y1[3] = y[3] + erre1[0] * r14*h;
		y1[4] = y[4] + erre1[0] * r15*h;
		y1[5] = y[5] + erre1[0] * r16*h;

		f(x + acca[0] * h, y1, eq);

		r21 = eq[0];
		r22 = eq[1];
		r23 = eq[2];
		r24 = eq[3];
		r25 = eq[4];
		r26 = eq[5];

		// 2ø order

		y1[0] = y[0] + erre2[0] * r11*h + erre2[1] * r21*h;
		y1[1] = y[1] + erre2[0] * r12*h + erre2[1] * r22*h;
		y1[2] = y[2] + erre2[0] * r13*h + erre2[1] * r23*h;
		y1[3] = y[3] + erre2[0] * r14*h + erre2[1] * r24*h;
		y1[4] = y[4] + erre2[0] * r15*h + erre2[1] * r25*h;
		y1[5] = y[5] + erre2[0] * r16*h + erre2[1] * r26*h;

		f(x + acca[1] * h, y1, eq);

		r31 = eq[0];
		r32 = eq[1];
		r33 = eq[2];
		r34 = eq[3];
		r35 = eq[4];
		r36 = eq[5];

		// 3ø order

		y1[0] = y[0] + erre3[0] * r11*h + erre3[1] * r21*h + erre3[2] * r31*h;
		y1[1] = y[1] + erre3[0] * r12*h + erre3[1] * r22*h + erre3[2] * r32*h;
		y1[2] = y[2] + erre3[0] * r13*h + erre3[1] * r23*h + erre3[2] * r33*h;
		y1[3] = y[3] + erre3[0] * r14*h + erre3[1] * r24*h + erre3[2] * r34*h;
		y1[4] = y[4] + erre3[0] * r15*h + erre3[1] * r25*h + erre3[2] * r35*h;
		y1[5] = y[5] + erre3[0] * r16*h + erre3[1] * r26*h + erre3[2] * r36*h;

		f(x + acca[2] * h, y1, eq);

		r41 = eq[0];
		r42 = eq[1];
		r43 = eq[2];
		r44 = eq[3];
		r45 = eq[4];
		r46 = eq[5];

		// 4ø order

		y1[0] = y[0] + erre4[0] * r11*h + erre4[1] * r21*h + erre4[2] * r31*h + erre4[3] * r41*h;
		y1[1] = y[1] + erre4[0] * r12*h + erre4[1] * r22*h + erre4[2] * r32*h + erre4[3] * r42*h;
		y1[2] = y[2] + erre4[0] * r13*h + erre4[1] * r23*h + erre4[2] * r33*h + erre4[3] * r43*h;
		y1[3] = y[3] + erre4[0] * r14*h + erre4[1] * r24*h + erre4[2] * r34*h + erre4[3] * r44*h;
		y1[4] = y[4] + erre4[0] * r15*h + erre4[1] * r25*h + erre4[2] * r35*h + erre4[3] * r45*h;
		y1[5] = y[5] + erre4[0] * r16*h + erre4[1] * r26*h + erre4[2] * r36*h + erre4[3] * r46*h;

		f(x + acca[3] * h, y1, eq);

		r51 = eq[0];
		r52 = eq[1];
		r53 = eq[2];
		r54 = eq[3];
		r55 = eq[4];
		r56 = eq[5];

		// 5ø order

		y1[0] = y[0] + erre5[0] * r11*h + erre5[1] * r21*h + erre5[2] * r31*h + erre5[3] * r41*h + erre5[4] * r51*h;
		y1[1] = y[1] + erre5[0] * r12*h + erre5[1] * r22*h + erre5[2] * r32*h + erre5[3] * r42*h + erre5[4] * r52*h;
		y1[2] = y[2] + erre5[0] * r13*h + erre5[1] * r23*h + erre5[2] * r33*h + erre5[3] * r43*h + erre5[4] * r53*h;
		y1[3] = y[3] + erre5[0] * r14*h + erre5[1] * r24*h + erre5[2] * r34*h + erre5[3] * r44*h + erre5[4] * r54*h;
		y1[4] = y[4] + erre5[0] * r15*h + erre5[1] * r25*h + erre5[2] * r35*h + erre5[3] * r45*h + erre5[4] * r55*h;
		y1[5] = y[5] + erre5[0] * r16*h + erre5[1] * r26*h + erre5[2] * r36*h + erre5[3] * r46*h + erre5[4] * r56*h;

		f(x + acca[4] * h, y1, eq);

		r61 = eq[0];
		r62 = eq[1];
		r63 = eq[2];
		r64 = eq[3];
		r65 = eq[4];
		r66 = eq[5];
		
		/* Begin 
		err = (fabs(erre1[0] * r11 + erre3[1] * r31 + erre4[2] * r41 + erre5[3] * r51 + erre6[4] * r61)*h);
		gotoxy(10, 8);
		printf("Errore=  %f", err);
		h = h / 2;

		}
		while (err>rkfr);
		/* End */
		
		// weighted average of operators 5^ order

		y[0] = y[0] + (ips5[0] * r11 + ips5[2] * r31 + ips5[3] * r41 + ips5[4] * r51)*h;

		y[1] = y[1] + (ips5[0] * r12 + ips5[2] * r32 + ips5[3] * r42 + ips5[4] * r52)*h;

		y[2] = y[2] + (ips5[0] * r13 + ips5[2] * r33 + ips5[3] * r43 + ips5[4] * r53)*h;

		y[3] = y[3] + (ips5[0] * r14 + ips5[2] * r34 + ips5[3] * r44 + ips5[4] * r54)*h;

		y[4] = y[4] + (ips5[0] * r15 + ips5[2] * r35 + ips5[3] * r45 + ips5[4] * r55)*h;

		y[5] = y[5] + (ips5[0] * r16 + ips5[2] * r36 + ips5[3] * r46 + ips5[4] * r56)*h;


		// Plot
					
		SetPixel(hdc, x * 1000, (-coppia) + 300, RGB(0, 0, 255)); // Blu (Torque)
		
		SetPixel(hdc, x * 1000, (-y[4]) + 300, RGB(255, 255, 255)); // White (Speed)
	
		SetPixel(hdc, x * 1000, (-y[5]) + 300, RGB(0, 255, 0)); // Green (Position)
					
		SetPixel(hdc, x * 1000, (-iss) + 300, RGB(255, 0, 0)); // Red

		SetPixel(hdc, x * 1000, (-irs) + 300, RGB(255, 0, 255)); // Magenta
				
		// **********

		x = x + h;
	}

	_getch();
	DestroyWindow(hwnd);

	return;

}


void f(double x, double y[6], double eq[6])

{
	double alfas, alfar, alfat, ks, kr, delta;



	alfas = rs / (ls + lm);
	alfar = rr / (lr + lm);
	alfat = ls + lm;
	ks = lm / (ls + lm);
	kr = lm / (lr + lm);
	delta = 1. / (1. - ks*kr);


	// Integrated equation sets (Asynchronous Math Model)


	// Isr real component of Is (Stator Current)

	eq[0] = (-alfas*y[0] + alfar*ks*y[2] + p*y[4] * ks*y[3] + p*y[4] * ks*kr*y[1] + (1. / alfat)*vsr)*delta;

	// Isi imaginary component of Is (Stator Current)

	eq[1] = (-kr*p*y[4] * ks*y[0] - p*y[4] * ks*y[2] + alfar*ks*y[3] - alfas*y[1] + (1. / alfat)*vsi)*delta;

	// Irr real component Ir (Rotor Current)

	eq[2] = (alfas*kr*y[0] - alfar*y[2] - p*y[4] * y[3] - p*y[4] * kr*y[1] - (kr / alfat)*vsr)*delta;

	// Iri imaginary component Ir (Rotor Current)

	eq[3] = (-alfar*y[3] + p*y[4] * y[2] + p*y[4] * kr*y[0] + alfas*kr*y[1] - (kr / alfat)*vsi)*delta;

	// W Angular Velocity (Mechanical balance equation)

	eq[4] = (3. / 2.*p*lm*(y[1] * y[2] - y[0] * y[3]) - cr) / j;

	// Position

	eq[5] = y[4];

}
