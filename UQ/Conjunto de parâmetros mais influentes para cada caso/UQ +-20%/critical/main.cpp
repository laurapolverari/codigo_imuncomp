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
#include <random>
#include <chrono>


#ifndef SUN_COMM_NULL
#define SUN_COMM_NULL NULL
#endif

#define TAM 9

#define NEQ   12               /* number of equations  */
#define RTOL  SUN_RCONST(1.0e-4) /* scalar relative tolerance           */
#define ATOL  SUN_RCONST(1.0e-4)

#define ZERO SUN_RCONST(0.0)

using namespace std;


// condições iniciais

real V0 = 19225.3, Ap0 = 1791.3, Apm0 = 0.0, Thn0 = 1e3,
	 The0 = 0.0, Tkn0 = 5.0e3, Tke0 = 0.0, B0 = 2.5e5,
	 Ps0 = 0.0, Pl0 = 0.0, Bm0 = 0.0, A0 = 0.0;

// parâmetros

real pi_v = 0.00167661, c_v1 = 1.12583762e-9, c_v2 = 6.0e-1, k_v1 = 1.20885e-05,
	 k_v2 = 0.02929, alpha_ap = 0.389196, beta_ap = 16.4079, c_ap1 = 4.08257620e0,
	 c_ap2 = 4.85488215e5, delta_apm = 4.27275e+07, alpha_th = 2.17e-1, beta_th = 5.99593,
	 pi_th = 2654.38, delta_th = 0.0528155, alpha_tk = 0.300002, beta_tk = 0.333305,
	 pi_tk = 2777.8, delta_tk = 0.0858949, alpha_b = 3.58e2, pi_b1 = 8.98e-5,
	 pi_b2 = 1.27e-8, beta_ps = 6.0e-6, beta_pl = 5.0e-6, beta_bm = 1.0e-6,
	 delta_ps = 2.5e0, delta_pl = 3.5e-1, gamma_bm = 9.75e-4, pi_bm1 = 8.11167858e0,
	 pi_bm2 = 3.79652355e3, pi_ps = 4.00409160e4, pi_pl = 302795, delta_a = 84788.3;

int simnum;

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


static real modelcsv(real *x, ofstream& csv_file, ofstream& csv_file1)
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

	retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
	if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

	y = N_VNew_Serial(NEQ, sunctx);
	if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }

	NV_Ith_S(y, 0) = V0;
	NV_Ith_S(y, 1) = x[8];
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

	k_v2 = x[0];
    beta_ap = x[1];
    delta_apm = x[2];
    pi_tk = x[3];
    delta_tk = x[4];
	pi_th = x[5];
	delta_th = x[6];
	alpha_ap = x[7];

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

	
	
	//csv_file << t0 << "," << NV_Ith_S(y, 0) << "," << NV_Ith_S(y, 4) << "," << NV_Ith_S(y, 6) << "," << NV_Ith_S(y, 11) << "," << simnum << "\n";  // Write data to CSV file

	vector<real> the_data, tke_data, virus_data, ant_data, t_data;
	t_data.push_back(t0);
	the_data.push_back(NV_Ith_S(y, 4));
	tke_data.push_back(NV_Ith_S(y, 6));
	virus_data.push_back(NV_Ith_S(y, 0));
	ant_data.push_back(NV_Ith_S(y, 11));

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

			t_data.push_back(t);
			the_data.push_back(NV_Ith_S(y, 4));
			tke_data.push_back(NV_Ith_S(y, 6));
			virus_data.push_back(NV_Ith_S(y, 0));
			ant_data.push_back(NV_Ith_S(y, 11));
			//csv_file << t << "," << NV_Ith_S(y, 0) << "," << NV_Ith_S(y, 4) << "," << NV_Ith_S(y, 6) << "," << NV_Ith_S(y, 11) << "," << simnum << "\n";  // Write data to CSV file

		}

  	}

	N_VDestroy(y);            
	N_VDestroy(abstol);       
	CVodeFree(&cvode_mem);    
	SUNLinSolFree(LS);        
	SUNMatDestroy(A);         
	SUNContext_Free(&sunctx);

	if(retval < 0)
	{
		simnum--;

		csv_file1 << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
			 << "," << x[4] << "," << x[5] << "," << x[6] << "," << x[7]
			 << "," << x[8] << "," << 1 << "\n";
	}
	else
	{
		for(int i = 0; i < the_data.size(); i++)
		{
			csv_file << t_data[i] << "," << virus_data[i] << "," << the_data[i] << "," << tke_data[i] << "," << ant_data[i] << "," << simnum << "\n";  // Write data to CSV file
		}

		csv_file1 << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
			 << "," << x[4] << "," << x[5] << "," << x[6] << "," << x[7]
			 << "," << x[8] << "," << 0 << "\n";
	}

	return 0;
}

int main()
{

	real porcentagem = 0.2;
	real auxMin = 1-porcentagem;
  	real auxMax = 1+porcentagem;

	simnum = 0;
	real func;
	ofstream csv_file("model_results_uq.csv");
	ofstream csv_file1("param_values_uq.csv");
	csv_file << "Time,V,The,Tke,A,SimNum\n";  // Write the header
	csv_file1 << "k_v2,beta_ap,delta_apm,pi_tk,delta_tk,pi_th,delta_th,alpha_ap,Ap0,Err\n";

	const int nrolls=10000;  // number of experiments
	// real low = delta_th*auxMin;
	// real high = delta_th*auxMax;
	// real *x = new real[1];

	vector<real> low = {k_v2*auxMin, beta_ap*auxMin, delta_apm*auxMin, pi_tk*auxMin,
							 delta_tk*auxMin, pi_th*auxMin, delta_th*auxMin, alpha_ap*auxMin,
							 Ap0*auxMin};

	vector<real> high = {k_v2*auxMax, beta_ap*auxMax, delta_apm*auxMax, pi_tk*auxMax,
							 delta_tk*auxMax, pi_th*auxMax, delta_th*auxMax, alpha_ap*auxMax,
							 Ap0*auxMax};

	real *x = new real[9];
	//printf("\n%f", model(x));

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

  	uniform_real_distribution<real> distribution_0(low[0],high[0]);
	uniform_real_distribution<real> distribution_1(low[1],high[1]);
	uniform_real_distribution<real> distribution_2(low[2],high[2]);
	uniform_real_distribution<real> distribution_3(low[3],high[3]);
	uniform_real_distribution<real> distribution_4(low[4],high[4]);
	uniform_real_distribution<real> distribution_5(low[5],high[5]);
	uniform_real_distribution<real> distribution_6(low[6],high[6]);
	uniform_real_distribution<real> distribution_7(low[7],high[7]);
	uniform_real_distribution<real> distribution_8(low[8],high[8]);


	while(true)
	{
		simnum++;
		if(simnum > nrolls)
		{
			break;
		}

		x[0] = distribution_0(generator);
		x[1] = distribution_1(generator);
		x[2] = distribution_2(generator);
		x[3] = distribution_3(generator);
		x[4] = distribution_4(generator);
		x[5] = distribution_5(generator);
		x[6] = distribution_6(generator);
		x[7] = distribution_7(generator);
		x[8] = distribution_8(generator);

		cout << simnum << ": x[0] = " << x[0] << "; " << "x[1] = " << x[1] << "; " << "x[2] = " << x[2] << "; " 
			<< "x[3] = " << x[3] << "; " << "x[4] = " << x[4] << "; " << "x[5] = " << x[5] << "; " 
			<< "x[6] = " << x[6] << "; " << "x[7] = " << x[7] << "; " << "x[8] = " << x[8] << endl;
		func = modelcsv(x, csv_file, csv_file1);
	}

	csv_file.close();
	csv_file1.close();

	return 0;
}