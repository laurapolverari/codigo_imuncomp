#include "de.h"
#include <math.h>
#include <cstdio>
#include <ctime>
#include <cvode/cvode.h>            /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>

#ifndef SUN_COMM_NULL
#define SUN_COMM_NULL NULL
#endif

#define NEQ   12               /* number of equations  */
#define RTOL  SUN_RCONST(1.0e-4) /* scalar relative tolerance           */
#define ATOL  SUN_RCONST(1.0e-4)

#define ZERO SUN_RCONST(0.0)

using namespace std;


// condições iniciais

real V0 = 18602.2, Ap0 = 1e4, Apm0 = 0.0, Thn0 = 1e3,
	 The0 = 0.0, Tkn0 = 5.0e3, Tke0 = 0.0, B0 = 2.5e5,
	 Ps0 = 0.0, Pl0 = 0.0, Bm0 = 0.0, A0 = 0.0;

// parâmetros

real pi_v = 0.0030666, c_v1 = 1.12583762e-9, c_v2 = 6.0e-1, k_v1 = 9.99415e-05,
	 k_v2 = 0.00484721, alpha_ap = 2.59993413e-1, beta_ap = 16.3662, c_ap1 = 4.08257620e0,
	 c_ap2 = 4.85488215e5, delta_apm = 2.83713e+08, alpha_th = 2.17e-1, beta_th = 144.577,
	 pi_th = 27119.2, delta_th = 0.0697364, alpha_tk = 0.0106936, beta_tk = 45.2902,
	 pi_tk = 1836.21, delta_tk = 0.0659293, alpha_b = 3.58e2, pi_b1 = 8.98e-5,
	 pi_b2 = 1.27e-8, beta_ps = 6.0e-6, beta_pl = 5.0e-6, beta_bm = 1.0e-6,
	 delta_ps = 2.5e0, delta_pl = 3.5e-1, gamma_bm = 9.75e-4, pi_bm1 = 8.11167858e0,
	 pi_bm2 = 3.79652355e3, pi_ps = 4.00409160e4, pi_pl = 20018.7, delta_a = 11738.9;


static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
	
    sunrealtype V = NV_Ith_S(y, 0);
    sunrealtype Ap = NV_Ith_S(y, 1);
    sunrealtype Apm = NV_Ith_S(y, 2);
    sunrealtype Thn = NV_Ith_S(y, 3);
    sunrealtype The = NV_Ith_S(y, 4);
    sunrealtype Tkn = NV_Ith_S(y, 5);
    sunrealtype Tke = NV_Ith_S(y, 6);
    sunrealtype B = NV_Ith_S(y, 7);
    sunrealtype Ps = NV_Ith_S(y, 8);
    sunrealtype Pl = NV_Ith_S(y, 9);
    sunrealtype Bm = NV_Ith_S(y, 10);
    sunrealtype A = NV_Ith_S(y, 11);

    /*dVdt =*/ NV_Ith_S(ydot, 0) = pi_v * V - ((c_v1 * V) / (c_v2 + V)) - k_v1 * V * A - k_v2 * V * Tke; // vírus
    
    /*dApdt =*/ NV_Ith_S(ydot, 1) = alpha_ap * (Ap0 - Ap) - beta_ap * Ap * ((c_ap1 * V) / (c_ap2 + V));  // APCs imaturas
    
    /*dApmdt =*/ NV_Ith_S(ydot, 2) = beta_ap * Ap * ((c_ap1 * V) / (c_ap2 + V)) - delta_apm * Apm;  // APCs maduras

    /*dThndt =*/ NV_Ith_S(ydot, 3) = alpha_th * (Thn0 - Thn) - beta_th * Apm * Thn;  // células T CD4+ imaturas

    /*dThedt =*/ NV_Ith_S(ydot, 4) = beta_th * Apm * Thn + pi_th * Apm * The - delta_th * The; // células T CD4+ maduras

    /*dTkndt =*/ NV_Ith_S(ydot, 5) = alpha_tk * (Tkn0 - Tkn) - beta_tk * Apm * Tkn; // células T CD8+ imaturas

    /*dTkedt =*/ NV_Ith_S(ydot, 6) = beta_tk * Apm * Tkn + pi_tk * Apm * Tke - delta_tk * Tke; // células T CD8+ maduras

    /*dBdt =*/ NV_Ith_S(ydot, 7) = alpha_b * (B0 - B) + pi_b1 * V * B + pi_b2 * The * B - beta_ps * Apm * B - beta_pl * The * B - beta_bm * The * B; // células B

    /*dPsdt =*/ NV_Ith_S(ydot, 8) = beta_ps * Apm * B - delta_ps * Ps; // plasmócitos de vida curta

    /*dPldt =*/ NV_Ith_S(ydot, 9) = beta_pl * The * B - delta_pl * Pl + gamma_bm * Bm; // plasmócitos de vida longa

    /*dBmdt =*/ NV_Ith_S(ydot, 10) = beta_bm * The * B + pi_bm1 * Bm * (1 - (Bm / pi_bm2)) - gamma_bm * Bm; // células B de memória

    /*dAdt =*/ NV_Ith_S(ydot, 11) = pi_ps * Ps + pi_pl * Pl - delta_a * A; // anticorpos

    return 0;
}

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}

static real modelcsv(real *x, ofstream& csv_file)
{
	SUNContext sunctx;
	sunrealtype t0, tf;
	N_Vector y;
	N_Vector abstol;
	SUNMatrix A;
	SUNLinearSolver LS;
	void* cvode_mem;
	int retval;

	y         = NULL;
	abstol    = NULL;
	A         = NULL;
	LS        = NULL;
	cvode_mem = NULL;

	t0 = 0.0;
	tf = 42.0;
	//tf = 100.0;

	retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
	if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

	y = N_VNew_Serial(NEQ, sunctx);
	if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }

	NV_Ith_S(y, 0) = V0;
	NV_Ith_S(y, 1) = Ap0;
	NV_Ith_S(y, 2) = Apm0;
	NV_Ith_S(y, 3) = Thn0;
	NV_Ith_S(y, 4) = The0;
	NV_Ith_S(y, 5) = Tkn0;
	NV_Ith_S(y, 6) = Tke0;
	NV_Ith_S(y, 7) = B0;
	NV_Ith_S(y, 8) = Ps0;
	NV_Ith_S(y, 9) = Pl0;
	NV_Ith_S(y, 10) = Bm0;
	NV_Ith_S(y, 11) = A0;

	pi_ps = x[0];
    printf("x[0] = %lf\n", x[0]);

	abstol = N_VNew_Serial(NEQ, sunctx);
	if (check_retval((void*)abstol, "N_VNew_Serial", 0)) { return (1); }

	cvode_mem = CVodeCreate(CV_BDF, sunctx);
  	if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

	retval = CVodeInit(cvode_mem, f, t0, y);
  	if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

	retval = CVodeSStolerances(cvode_mem, RTOL, ATOL);
    if (check_retval(&retval, "CVodeSStolerances", 1)) { return (1); }

	A = SUNDenseMatrix(NEQ, NEQ, sunctx);
	if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

	LS = SUNLinSol_Dense(y, A, sunctx);
	if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

	retval = CVodeSetLinearSolver(cvode_mem, LS, A);
	if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

	csv_file << "Time,V,Ap,Apm,Thn,The,Tkn,Tke,B,Ps,Pl,Bm,A\n";  // Write the header

	csv_file << t0 << "," << NV_Ith_S(y, 0) << "," << NV_Ith_S(y, 1) << "," << NV_Ith_S(y, 2) << "," << NV_Ith_S(y, 3)
					<< "," << NV_Ith_S(y, 4) << "," << NV_Ith_S(y, 5) << "," << NV_Ith_S(y, 6) << "," << NV_Ith_S(y, 7)
					<< "," << NV_Ith_S(y, 8) << "," << NV_Ith_S(y, 9) << "," << NV_Ith_S(y, 10) << "," << NV_Ith_S(y, 11) << "\n";  // Write data to CSV file

	sunrealtype t = t0;
	sunrealtype dt = 0.1;
	while (t < tf)
	{
		sunrealtype tret;
		retval = CVode(cvode_mem, t + dt, y, &tret, CV_NORMAL);
		
		if (check_retval(&retval, "CVode", 1)) { break; }
		if (retval == CV_SUCCESS)
		{
			t = tret;
			//printf("At time %g: The = %g, Tke = %g, V = %g, A = %g\n", t, NV_Ith_S(y, 4), NV_Ith_S(y, 6), NV_Ith_S(y, 0), NV_Ith_S(y, 11));
			csv_file << t << "," << NV_Ith_S(y, 0) << "," << NV_Ith_S(y, 1) << "," << NV_Ith_S(y, 2) << "," << NV_Ith_S(y, 3)
					<< "," << NV_Ith_S(y, 4) << "," << NV_Ith_S(y, 5) << "," << NV_Ith_S(y, 6) << "," << NV_Ith_S(y, 7)
					<< "," << NV_Ith_S(y, 8) << "," << NV_Ith_S(y, 9) << "," << NV_Ith_S(y, 10) << "," << NV_Ith_S(y, 11) << "\n";  // Write data to CSV file

		}

  	}

	N_VDestroy(y);            
	N_VDestroy(abstol);       
	CVodeFree(&cvode_mem);    
	SUNLinSolFree(LS);        
	SUNMatDestroy(A);         
	SUNContext_Free(&sunctx);

	return 0;
}

int main()
{  

	real *x = new real[1];
    x[0] = pi_ps;
	string var = "pi_ps";

    real porc1 = 0.2;
    real auxMin1 = 1-porc1;
  	real auxMax1 = 1+porc1;

    real porc2 = 0.1;
    real auxMin2 = 1-porc2;
  	real auxMax2 = 1+porc2;

    vector<real> val = {x[0]*auxMin1, x[0]*auxMin2, x[0], x[0]*auxMax2, x[0]*auxMax1};
	
	ofstream csv_file("model_results_min1_"  + var + ".csv");
    ofstream csv_file1("model_results_min2_" + var + ".csv");
    ofstream csv_file2("model_results_base_" + var + ".csv");
    ofstream csv_file3("model_results_max2_" + var + ".csv");
    ofstream csv_file4("model_results_max1_" + var + ".csv");

	real func = modelcsv(&val[0], csv_file);
    func = modelcsv(&val[1], csv_file1);
    func = modelcsv(&val[2], csv_file2);
    func = modelcsv(&val[3], csv_file3);
    func = modelcsv(&val[4], csv_file4);

    csv_file.close();
	csv_file1.close();
    csv_file2.close();
	csv_file3.close();
    csv_file4.close();

	return 0;
}