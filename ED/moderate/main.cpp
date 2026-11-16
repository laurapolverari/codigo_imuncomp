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

#define TAM 9

#define NEQ   12               /* number of equations  */
#define RTOL  SUN_RCONST(1.0e-4) /* scalar relative tolerance           */
#define ATOL  SUN_RCONST(1.0e-4)

#define ZERO SUN_RCONST(0.0)

using namespace std;


// condições iniciais
real V0 = 17417.4, Ap0 = 1e4, Apm0 = 0.0, Thn0 = 1e3,
	 The0 = 0.0, Tkn0 = 5.0e3, Tke0 = 0.0, B0 = 2.5e5,
	 Ps0 = 0.0, Pl0 = 0.0, Bm0 = 0.0, A0 = 0.0;

// parâmetros
real pi_v = 0.164336, c_v1 = 1.12583762e-9, c_v2 = 6.0e-1, k_v1 = 0.000123457,
	 k_v2 = 6.83881219e-4, alpha_ap = 2.59993413e-1, beta_ap = 4.29054542e7, c_ap1 = 4.08257620e0,
	 c_ap2 = 4.85488215e5, delta_apm = 1.10015322e5, alpha_th = 2.17e-1, beta_th = 33.4552,
	 pi_th = 1.0e-2, delta_th = 2.92e0, alpha_tk = 0.00502921, beta_tk = 14.1451,
	 pi_tk = 1.70684321e-6, delta_tk = 1.14072505e4, alpha_b = 3.58e2, pi_b1 = 8.98e-5,
	 pi_b2 = 1.27e-8, beta_ps = 6.0e-6, beta_pl = 5.0e-6, beta_bm = 1.0e-6,
	 delta_ps = 2.5e0, delta_pl = 3.5e-1, gamma_bm = 9.75e-4, pi_bm1 = 8.11167858e0,
	 pi_bm2 = 3.79652355e3, pi_ps = 4.00409160e4, pi_pl = 43628.9, delta_a = 18007;


vector<real> cd4_x, cd4_y, cd8_x, cd8_y, virus_x, virus_y, ant_x, ant_y;
vector<real> cd4_y_v, cd8_y_v, virus_y_v, ant_y_v;
int cd4_i, cd4_f, cd8_i, cd8_f, virus_i, virus_f, ant_i, ant_f;
vector<real> The_d_log, Tke_d_log, V_d_log, A_d_log; // relacionados aos dados utilizados
real mu_The, mu_Tke, mu_V, mu_A;
real sigma_The, sigma_Tke, sigma_V, sigma_A;
vector<real> Z_The_d, Z_Tke_d, Z_V_d, Z_A_d;


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

static real modelcsv(real *x)
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
	tf = 35.0;

	retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
	if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

	y = N_VNew_Serial(NEQ, sunctx);
	if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }

	/*
    Params description:

    x[0]=> k_v2
    x[1]=> beta_ap
    x[2]=> delta_apm
    x[3]=> pi_tk
    x[4]=> delta_tk
    x[5]=> pi_th
    x[6]=> delta_th
    x[7]=> alpha_ap
    x[8]=> Ap0
    */

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

	
	std::ofstream csv_file("model_results.csv");
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

	csv_file.close();

	N_VDestroy(y);            
	N_VDestroy(abstol);       
	CVodeFree(&cvode_mem);    
	SUNLinSolFree(LS);        
	SUNMatDestroy(A);         
	SUNContext_Free(&sunctx);

	return 0;
}

static real twonorm(vector<real>& v)
{
	real n = 0.0;

	for(int i = 0; i < v.size(); i++)
	{
		n += pow(v[i], 2);
	}

	n = sqrt(n);

	return n;

}

static void subt(vector<real>& v1, vector<real>& v2, vector<real>& vf)
{
	for(int i = 0; i < v1.size(); i++)
	{
		vf.push_back(v1[i] - v2[i]);
	}

	return;
}

real mean(vector<real>& v)
{
	real mean = 0.0;

	for(int i = 0; i < v.size(); i++)
	{
		mean += v[i];
	}

	mean = mean / float(v.size());

	return mean;
}

real standard_deviation(vector<real>& v)
{
	real std = 0.0;

	for(int i = 0; i < v.size(); i++)
	{
		std += pow(v[i], 2);
	}

	std = std / float(v.size() - 1);

	std = sqrt(std);

	return std;
}

static real model(real *x)
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
	tf = 35.0;

	retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
	if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

	y = N_VNew_Serial(NEQ, sunctx);
	if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }

	/*
    Params description:

    x[0]=> k_v2
    x[1]=> beta_ap
    x[2]=> delta_apm
    x[3]=> pi_tk
    x[4]=> delta_tk
    x[5]=> pi_th
    x[6]=> delta_th
    x[7]=> alpha_ap
    x[8]=> Ap0
    */

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


	vector<real> the_data, tke_data, virus_data, ant_data;

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

			the_data.push_back(NV_Ith_S(y, 4));
			tke_data.push_back(NV_Ith_S(y, 6));
			virus_data.push_back(NV_Ith_S(y, 0));
			ant_data.push_back(NV_Ith_S(y, 11));
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
		return 1e3;
	}

	
	vector<real> The_e, Tke_e, V_e, A_e;
	//vector<real> The_e, Tke_e;

	The_e = {the_data.begin() + cd4_i, the_data.begin() + cd4_f + 1};
	Tke_e = {tke_data.begin() + cd8_i, tke_data.begin() + cd8_f + 1};
	V_e = {virus_data.begin() + virus_i, virus_data.begin() + virus_f + 1};
	A_e = {ant_data.begin() + ant_i, ant_data.begin() + ant_f + 1};

	vector<real>::iterator The_valmin = min_element(The_e.begin(), The_e.end());
	real The_vmin = *The_valmin;

	vector<real>::iterator Tke_valmin = min_element(Tke_e.begin(), Tke_e.end());
	real Tke_vmin = *Tke_valmin;

	vector<real>::iterator V_valmin = min_element(V_e.begin(), V_e.end());
	real V_vmin = *V_valmin;
	
	vector<real>::iterator A_valmin = min_element(A_e.begin(), A_e.end());
	real A_vmin = *A_valmin;


	if(The_vmin <= 0 || Tke_vmin <= 0 || V_vmin <= 0 || A_vmin <= 0)
	{
		return 1e3;
	}

	// Transformação logarítmica do resultado do modelo

	float eps = 1.0;
	
	vector<real> The_m_log, Tke_m_log, V_m_log, A_m_log; // relacionados aos resultados do modelo

	for(int i = 0; i < The_e.size(); i++)
	{
		The_m_log.push_back(log10(The_e[i] + eps));
	}

	for(int i = 0; i < Tke_e.size(); i++)
	{
		Tke_m_log.push_back(log10(Tke_e[i] + eps));
	}

	for(int i = 0; i < V_e.size(); i++)
	{
		V_m_log.push_back(log10(V_e[i] + eps));
	}

	for(int i = 0; i < A_e.size(); i++)
	{
		A_m_log.push_back(log10(A_e[i] + eps));
	}

	vector<real> The_m_log_subt, Tke_m_log_subt, V_m_log_subt, A_m_log_subt;

	for(int i = 0; i < The_m_log.size(); i++)
	{
		The_m_log_subt.push_back(The_m_log[i] - mu_The);
	}

	for(int i = 0; i < Tke_m_log.size(); i++)
	{
		Tke_m_log_subt.push_back(Tke_m_log[i] - mu_Tke);
	}

	for(int i = 0; i < V_m_log.size(); i++)
	{
		V_m_log_subt.push_back(V_m_log[i] - mu_V);
	}

	for(int i = 0; i < A_m_log.size(); i++)
	{
		A_m_log_subt.push_back(A_m_log[i] - mu_A);
	}

	// Z-score

	vector<real> Z_The_m, Z_Tke_m, Z_V_m, Z_A_m;

	for(int i = 0; i < The_m_log_subt.size(); i++)
	{
		Z_The_m.push_back(The_m_log_subt[i] / sigma_The);
	}

	for(int i = 0; i < Tke_m_log_subt.size(); i++)
	{
		Z_Tke_m.push_back(Tke_m_log_subt[i] / sigma_Tke);
	}

	for(int i = 0; i < V_m_log_subt.size(); i++)
	{
		Z_V_m.push_back(V_m_log_subt[i] / sigma_V);
	}

	for(int i = 0; i < A_m_log_subt.size(); i++)
	{
		Z_A_m.push_back(A_m_log_subt[i] / sigma_A);
	}

	// Resíduos

	vector<real> res_The, res_Tke, res_V, res_A;

	subt(Z_The_m, Z_The_d, res_The);
	subt(Z_Tke_m, Z_Tke_d, res_Tke);
	subt(Z_V_m, Z_V_d, res_V);
	subt(Z_A_m, Z_A_d, res_A);

	// Erro

	real error = 0.0;

	for(int i = 0; i < res_The.size(); i++)
	{
		error += pow(res_The[i], 2) + pow(res_Tke[i], 2) + pow(res_V[i], 2) + pow(res_A[i], 2);
	}

	error = error / float(res_The.size());

	return error;
}

static real lagrange(real x, vector<real>& x_v, vector<real>& y_v)
{
	  real x0, x1, y0, y1, L0, L1, P;

	  for(int i = 0; i < x_v.size(); i++)
	  {
		if(x <= x_v[i])
		{
			x0 = x_v[i - 1];
			x1 = x_v[i];
			y0 = y_v[i - 1];
			y1 = y_v[i];
			break;
		}
	  }

	  L0 = (x - x1) / (x0 - x1);
      L1 = (x - x0) / (x1 - x0);

	  P = L0 * y0 + L1 * y1;

	  return P;
}

int main()
{
	real t0, tf;
	t0 = 0.0;
	tf = 35.0;
	real step = 0.1;
	vector<real> t;

	for(real v = t0; v < tf; v += step)
	{
		t.push_back(v);
	}

	
	// read files
	
	fstream fin;

    fin.open("CD4_Moderate.csv", ios::in);
    string line, word, temp;
	int cont = 0;
	int numdelinhas = 0;

	getline(fin, line);

	while (fin >> temp)
	{

		getline(fin, line);

		stringstream s(temp);

		while (cont < 2)
        {
            getline(s, word, ',');

			if (cont == 0)
			{
				cd4_x.push_back(stor(word));
			}
			else
			{
				cd4_y.push_back(stor(word));
			}

            cont += 1;
        }

		getline(fin, line);
		getline(fin, line);
		cont = 0;

	}

	fin.close();

	

    fin.open("CD8_Moderate.csv", ios::in);
	cont = 0;

	getline(fin, line);

	while (fin >> temp)
	{

		getline(fin, line);

		stringstream s(temp);

		while (cont < 2)
        {
            getline(s, word, ',');

			if (cont == 0)
			{
				cd8_x.push_back(stor(word));
			}
			else
			{
				cd8_y.push_back(stor(word));
			}

            cont += 1;
        }

		getline(fin, line);
		getline(fin, line);
		cont = 0;

	}

	fin.close();

	

    fin.open("Virus_Moderate.csv", ios::in);
	cont = 0;

	getline(fin, line);

	while (fin >> temp)
	{

		getline(fin, line);

		stringstream s(temp);

		while (cont < 2)
        {
            getline(s, word, ',');

			if (cont == 0)
			{
				virus_x.push_back(stor(word));
			}
			else
			{
				virus_y.push_back(pow(10, stor(word)));
			}

            cont += 1;
        }

		getline(fin, line);
		getline(fin, line);
		cont = 0;

	}

	fin.close();

	

    fin.open("SNIgG_Moderate.csv", ios::in);
	cont = 0;

	getline(fin, line);

	while (fin >> temp)
	{

		getline(fin, line);

		stringstream s(temp);

		while (cont < 2)
        {
            getline(s, word, ',');

			if (cont == 0)
			{
				ant_x.push_back(stor(word));
			}
			else
			{
				ant_y.push_back(pow(10, stor(word)));
			}

            cont += 1;
        }

		getline(fin, line);
		getline(fin, line);
		cont = 0;

	}

	fin.close();


	

	// cd4
	int tam;
	real a;

	a = round(cd4_x[0] * 10) / 10;
	cd4_i = int((round((a - t0) / step)));
	tam = cd4_x.size() - 1;
	if (cd4_x[tam] > tf)
	{
		cd4_f = t.size() - 1;
	}
	else
	{
		a = int(cd4_x[tam] * 10) / 10.0;
		cd4_f = int((round((a - t0) / step)));
	}

	cd4_y_v = vector<real>((cd4_f - cd4_i + 1), 0.0);
	vector<real> t_cd4;
	real val_x;

	t_cd4 = {t.begin() + cd4_i, t.begin() + cd4_f + 1};

	for(int i = 0; i < cd4_y_v.size(); i++)
	{
		val_x = t_cd4[i];
		cd4_y_v[i] = lagrange(val_x, cd4_x, cd4_y);
	}


	// cd8

	a = round(cd8_x[0] * 10) / 10;
	cd8_i = int((round((a - t0) / step)));
	tam = cd8_x.size() - 1;
	if (cd8_x[tam] > tf)
	{
		cd8_f = t.size() - 1;
	}
	else
	{
		a = int(cd8_x[tam] * 10) / 10.0;
		cd8_f = int((round((a - t0) / step)));
	}

	cd8_y_v = vector<real>((cd8_f - cd8_i + 1), 0.0);
	vector<real> t_cd8;

	t_cd8 = {t.begin() + cd8_i, t.begin() + cd8_f + 1};

	for(int i = 0; i < cd8_y_v.size(); i++)
	{
		val_x = t_cd8[i];
		cd8_y_v[i] = lagrange(val_x, cd8_x, cd8_y);
	}


	// virus

	a = round(virus_x[0] * 10) / 10;
	virus_i = int((round((a - t0) / step)));
	tam = virus_x.size() - 1;
	if (virus_x[tam] > tf)
	{
		virus_f = t.size() - 1;
	}
	else
	{
		a = int(virus_x[tam] * 10) / 10.0;
		virus_f = int((round((a - t0) / step)));
	}

	virus_y_v = vector<real>((virus_f - virus_i + 1), 0.0);
	vector<real> t_virus;

	t_virus = {t.begin() + virus_i, t.begin() + virus_f + 1};

	for(int i = 0; i < virus_y_v.size(); i++)
	{
		val_x = t_virus[i];
		virus_y_v[i] = lagrange(val_x, virus_x, virus_y);
	}


	// anticorpos

	a = round(ant_x[0] * 10) / 10;
	ant_i = int((round((a - t0) / step)));
	tam = ant_x.size() - 1;
	if (ant_x[tam] > tf)
	{
		ant_f = t.size() - 1;
	}
	else
	{
		a = int(ant_x[tam] * 10) / 10.0;
		ant_f = int((round((a - t0) / step)));
	}

	ant_y_v = vector<real>((ant_f - ant_i + 1), 0.0);
	vector<real> t_ant;

	t_ant = {t.begin() + ant_i, t.begin() + ant_f + 1};

	for(int i = 0; i < ant_y_v.size(); i++)
	{
		val_x = t_ant[i];
		ant_y_v[i] = lagrange(val_x, ant_x, ant_y);
	}

	// Transformação logarítmica dos dados
	float eps = 1.0;

	for(int i = 0; i < cd4_y_v.size(); i++)
	{
		The_d_log.push_back(log10(cd4_y_v[i] + eps));
	}

	for(int i = 0; i < cd8_y_v.size(); i++)
	{
		Tke_d_log.push_back(log10(cd8_y_v[i] + eps));
	}

	for(int i = 0; i < virus_y_v.size(); i++)
	{
		V_d_log.push_back(log10(virus_y_v[i] + eps));
	}

	for(int i = 0; i < ant_y_v.size(); i++)
	{
		A_d_log.push_back(log10(ant_y_v[i] + eps));
	}

	// Estatísticas

	vector<real> The_d_log_subt, Tke_d_log_subt, V_d_log_subt, A_d_log_subt;

	mu_The = mean(The_d_log);
	mu_Tke = mean(Tke_d_log);
	mu_V = mean(V_d_log);
	mu_A = mean(A_d_log);

	for(int i = 0; i < The_d_log.size(); i++)
	{
		The_d_log_subt.push_back(The_d_log[i] - mu_The);
	}

	for(int i = 0; i < Tke_d_log.size(); i++)
	{
		Tke_d_log_subt.push_back(Tke_d_log[i] - mu_Tke);
	}

	for(int i = 0; i < V_d_log.size(); i++)
	{
		V_d_log_subt.push_back(V_d_log[i] - mu_V);
	}

	for(int i = 0; i < A_d_log.size(); i++)
	{
		A_d_log_subt.push_back(A_d_log[i] - mu_A);
	}

	sigma_The = standard_deviation(The_d_log_subt);
	sigma_Tke = standard_deviation(Tke_d_log_subt);
	sigma_V = standard_deviation(V_d_log_subt);
	sigma_A = standard_deviation(A_d_log_subt);

	// Z-score

	for(int i = 0; i < The_d_log_subt.size(); i++)
	{
		Z_The_d.push_back(The_d_log_subt[i] / sigma_The);
	}

	for(int i = 0; i < Tke_d_log_subt.size(); i++)
	{
		Z_Tke_d.push_back(Tke_d_log_subt[i] / sigma_Tke);
	}

	for(int i = 0; i < V_d_log_subt.size(); i++)
	{
		Z_V_d.push_back(V_d_log_subt[i] / sigma_V);
	}

	for(int i = 0; i < A_d_log_subt.size(); i++)
	{
		Z_A_d.push_back(A_d_log_subt[i] / sigma_A);
	}




	
	int n_gens = 1000000;
	int pop_size = 50;
	int num_par = TAM;

	real porcentagem = 0.9;
	real auxMin = 1-porcentagem;
  	real auxMax = 1+porcentagem;

	 
	//  V0 = 9381.25, pi_v = 0.118744, beta_tk = 15.899,
	//  k_v1 = 0.000723997, k_v2 = 0.000954002e+1, beta_ap = 7749.01,
	//  delta_apm = 3.73984e+08, pi_tk = 2770.9, delta_tk = 0.123303,
	//  pi_pl = 28230.2, delta_a = 19314.8, alpha_tk = 0.0299076,
	//  beta_th = 19.2, pi_th = 1481.94e+1, delta_th = 0.0558038e+1,
	//  alpha_ap = 2.59993413e-1, Ap0 = 1e4;

	//  k_v2 = 0.000954002e+1, beta_ap = 7749.01, delta_apm = 3.73984e+08,
	//  pi_tk = 2770.9, delta_tk = 0.123303, pi_th = 1481.94e+1,
	//  delta_th = 0.0558038e+1, alpha_ap = 2.59993413e-1, Ap0 = 1e4;

	 k_v2 = 0.00825958, beta_ap = 776.954, delta_apm = 4.39069e+07,
	 pi_tk = 3079.3, delta_tk = 0.233863, pi_th = 4303.69,
	 delta_th = 0.156153, alpha_ap = 2.59993413e-1, Ap0 = 1e4;

	//  vector<real> l_bounds = {V0*auxMin, pi_v * auxMin, beta_tk*auxMin, k_v1*auxMin,
	// 						 k_v2*auxMin, beta_ap*auxMin, delta_apm*auxMin, pi_tk*auxMin,
	// 						 delta_tk*auxMin, pi_pl*auxMin, delta_a*auxMin, alpha_tk*auxMin,
	// 						 beta_th*auxMin, pi_th*auxMin, delta_th*auxMin, alpha_ap*auxMin,
	// 						 Ap0*auxMin};

	// vector<real> h_bounds = {V0*auxMax, pi_v * auxMax, beta_tk*auxMax, k_v1*auxMax,
	// 						 k_v2*auxMax, beta_ap*auxMax, delta_apm*auxMax, pi_tk*auxMax,
	// 						 delta_tk*auxMax, pi_pl*auxMax, delta_a*auxMax, alpha_tk*auxMax,
	// 						 beta_th*auxMax, pi_th*auxMax, delta_th*auxMax, alpha_ap*auxMax,
	// 						 Ap0*auxMax};

	vector<real> l_bounds = {k_v2*auxMin, beta_ap*auxMin, delta_apm*auxMin, pi_tk*auxMin,
							 delta_tk*auxMin, pi_th*auxMin, delta_th*auxMin, alpha_ap*auxMin,
							 Ap0*auxMin};

	vector<real> h_bounds = {k_v2*auxMax, beta_ap*auxMax, delta_apm*auxMax, pi_tk*auxMax,
							 delta_tk*auxMax, pi_th*auxMax, delta_th*auxMax, alpha_ap*auxMax,
							 Ap0*auxMax};
	
	printf("n_gens: %d, pop_size: %d, num_par: %d\n", n_gens, pop_size, num_par);
	//real *x = differential_evolution(false, model, &l_bounds[0], &h_bounds[0], n_gens, 0.7, 0.5, 1.0, pop_size, num_par, BEST1BIN, time(NULL), 0, 0.01);
	real *x = differential_evolution(false, model, &l_bounds[0], &h_bounds[0], n_gens, 0.7, 0.5, 1.0, pop_size, num_par, BEST1BIN, 1232, 0, 0.01);

	cout << "x: " << endl;
	for(int i = 0; i < num_par; i++)
	{
        	cout << " " << x[i];
    }
	cout << endl;

	printf("\n%f", model(x));

	real func = modelcsv(x);

	std::ofstream csv_file("param_values.csv");
	// csv_file << "V0, pi_v, beta_tk, k_v1, k_v2, beta_ap, delta_apm, pi_tk, delta_tk, pi_pl, delta_a, alpha_tk, beta_th, pi_th, delta_th, alpha_ap, Ap0\n";  // Write the header
	csv_file << "k_v2, beta_ap, delta_apm, pi_tk, delta_tk, pi_th, delta_th, alpha_ap, Ap0\n";  // Write the header


	// csv_file << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
	// 		 << "," << x[4] << "," << x[5] << "," << x[6] << "," << x[7]
	// 		 << "," << x[8] << "," << x[9] << "," << x[10] << "," << x[11]
	// 		 << ","	<< x[12] << "," << x[13] << "," << x[14] << "," << x[15]
	// 		 << "," << x[16] << "\n";  // Write data to CSV file

	csv_file << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
			 << "," << x[4] << "," << x[5] << "," << x[6] << "," << x[7]
			 << "," << x[8] << "\n";  // Write data to CSV file

	return 0;
}