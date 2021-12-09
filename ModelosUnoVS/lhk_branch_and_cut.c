#define TIEMPO_MAXIMO_EN_SEGUNDOS 1800
#define CORTAR_A_DISTANCIA 0
#define CPXSIZE_BITS_TEST_DISABLE

#define TIPO_VARIABLE CPX_BINARY

#include "C:/Program Files/IBM/ILOG/CPLEX_Studio201/cplex/include/ilcplex/cplex.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "leer_dimac.h"

static void
free_and_null(char** ptr);

static int resolver(CPXENVptr env, CPXLPptr lp);

void buscarSolucionInicial(int* coloreoValido);

int mejorSolucionValida = 0;

int resolver(CPXENVptr env, CPXLPptr lp) {
	int status = 0;
	int solstat = 0;

	double objval;
	double* x = NULL;

	/* Optimize the problem and obtain solution */

	status = CPXmipopt(env, lp);
	if (status) {
		fprintf(stderr, "Failed to optimize MIP.\n");
		goto TERMINATE;
	}

	solstat = CPXgetstat(env, lp);
	printf("Solution status %d.\n", solstat);

	status = CPXgetobjval(env, lp, &objval);
	if (status) {
		fprintf(stderr, "Failed to obtain objective value.\n");
		goto TERMINATE;
	}

	printf("Objective value %.10g\n", objval);

	int cur_numcols = CPXgetnumcols(env, lp);

	/* Allocate space for solution */

	x = (double*)malloc(cur_numcols * sizeof(double));

	if (x == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}

	status = CPXgetx(env, lp, x, 0, cur_numcols - 1);
	if (status) {
		fprintf(stderr, "Failed to obtain solution.\n");
		goto TERMINATE;
	}

	/* Write out the solution */
	FILE* f = fopen("results.txt", "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int j;
	for (j = 0; j < cur_numcols; j++) {
		if (fabs(x[j]) > 1e-10) {
			char* colname[1];
			char namestore[10];
			int surplus = 0;
			status = CPXgetcolname(env, lp, colname, namestore, 10, &surplus, j, j);
			if (status) {
				namestore[0] = 0;
				colname[0] = namestore;
			}
			printf("Column %3d (%10s):  Value = %17.10g\n", j, colname[0], x[j]);

			if (strlen(colname[0]) > 3) {
				if (strlen(colname[0]) == 6) {
					fprintf(f, "%c%c %c\n", colname[0][2], colname[0][3], colname[0][5]);
				}
				else {
					fprintf(f, "%c %c\n", colname[0][2], colname[0][4]);
				}
			}
		}
	}

	fclose(f);

TERMINATE:

	free_and_null((char**)&x);

	return (status);
}

/* This simple routine frees up the pointer *ptr, and sets *ptr
 to NULL */

static void free_and_null(char** ptr) {
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */

static int populatebyrow(CPXENVptr env, CPXLPptr lp, bool** matriz, int* pesos, int Nr_vert, int maximoColor) {
	int status = 0;

	int zero = 0;
	double* ind = NULL;
	char* tipos = NULL;
	int i, j, k;
	char** nombreColumnas = NULL;

	nombreColumnas = (char**)malloc(maximoColor * sizeof(char*));
	tipos = (char*)malloc(maximoColor * sizeof(char));
	ind = (double*)malloc(maximoColor * sizeof(double));
	if (tipos == NULL || ind == NULL || nombreColumnas == NULL) {
		status = CPXERR_NO_MEMORY;
		goto TERMINATE;
	}
	for (j = 0; j < maximoColor; j++) {
		tipos[j] = CPX_INTEGER;
		ind[j] = 1;
		nombreColumnas[j] = (char*)malloc(6 * sizeof(char));
		sprintf(nombreColumnas[j], "x_%d", j + 1);
	}

	status = CPXnewcols(env, lp, maximoColor, ind, NULL, NULL, tipos, nombreColumnas);

	for (j = 0; j < maximoColor; j++) {
		free_and_null((char**)&nombreColumnas[j]);
	}
	free_and_null((char**)&nombreColumnas);

	if (status)
		goto TERMINATE;

	nombreColumnas = (char**)malloc(maximoColor * Nr_vert * sizeof(char*));
	for (i = 0; i < Nr_vert; i++) {
		for (j = 0; j < maximoColor; j++) {
			nombreColumnas[i * maximoColor + j] = (char*)malloc(10 * sizeof(char));
			sprintf(nombreColumnas[i * maximoColor + j], "x_%d_%d", i + 1, j + 1);
		}
	}

	free_and_null((char**)&tipos);
	tipos = (char*)malloc(maximoColor * Nr_vert * sizeof(char));
	if (tipos == NULL) {
		status = CPXERR_NO_MEMORY;
		goto TERMINATE;
	}
	for (j = 0; j < maximoColor * Nr_vert; j++) {
		tipos[j] = TIPO_VARIABLE;
	}

	status = CPXnewcols(env, lp, maximoColor * Nr_vert, NULL, NULL, NULL, tipos, nombreColumnas);

	for (j = 0; j < maximoColor * Nr_vert; j++) {
		free_and_null((char**)&nombreColumnas[j]);
	}
	free_and_null((char**)&nombreColumnas);

	if (status)
		goto TERMINATE;

	int* posicionVariables = NULL;
	double* coeficienteVariables = NULL;
	double coeficienteTerminoIndependiente[1];
	coeficienteTerminoIndependiente[0] = 1;

	posicionVariables = (int*)malloc(Nr_vert * sizeof(int));
	coeficienteVariables = (double*)malloc(Nr_vert * sizeof(double));

	//Eliminacion Simetria por Xk>=Xk+1
	for (k = 0; k < maximoColor - 1; k++) {
		posicionVariables[0] = k;
		coeficienteVariables[0] = 1;
		posicionVariables[1] = k + 1;
		coeficienteVariables[1] = -1;
		status = CPXaddrows(env, lp, 0, 1, 2, NULL, "G", &zero, posicionVariables, coeficienteVariables, NULL, NULL);
		if (status)
			goto TERMINATE;

	}

	//Si el vertice i usa color k la variable de ese color debe pesar por lo menos el peso del vertice i
	for (k = 0; k < maximoColor; k++) {
		for (i = 0; i < Nr_vert; i++) {
			posicionVariables[0] = maximoColor + i * maximoColor + k;
			coeficienteVariables[0] = pesos[i];

			posicionVariables[1] = k;
			coeficienteVariables[1] = -1;
			status = CPXaddrows(env, lp, 0, 1, 2, NULL, "L", &zero, posicionVariables, coeficienteVariables, NULL, NULL);
			if (status)
				goto TERMINATE;
		}
	}

	//Todos los vertices deben estar coloreados
	for (i = 0; i < Nr_vert; i++) {
		for (k = 0; k < maximoColor; k++) {
			posicionVariables[k] = maximoColor + i * maximoColor + k;
			coeficienteVariables[k] = 1;
		}
		status = CPXaddrows(env, lp, 0, 1, maximoColor, coeficienteTerminoIndependiente, "E", &zero, posicionVariables, coeficienteVariables, NULL,
			NULL);
		if (status)
			goto TERMINATE;
	}

	//Incompatibilidades
	for (i = 0; i < Nr_vert; i++) {
		for (j = 0; j < Nr_vert; j++) {
			if (matriz[i][j]) {
				for (k = 0; k < maximoColor; k++) {
					posicionVariables[0] = maximoColor + i * maximoColor + k;
					posicionVariables[1] = maximoColor + j * maximoColor + k;

					coeficienteVariables[0] = coeficienteVariables[1] = 1;

					status = CPXaddrows(env, lp, 0, 1, 2, coeficienteTerminoIndependiente, "L", &zero, posicionVariables, coeficienteVariables, NULL, NULL);
					if (status)
						goto TERMINATE;
				}
			}
		}
	}

TERMINATE:

	free_and_null((char**)&posicionVariables);
	free_and_null((char**)&coeficienteVariables);

	free_and_null((char**)&ind);
	free_and_null((char**)&tipos);
	return (status);

} /* END populatebyrow */

int procesar_Tp(char* path, char* pathOutput) {
	printf("Inicio\n");

	read_graph_DIMACS_ascii(path);

	int status = 0;

	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;

	env = CPXopenCPLEX(&status);

	if (env == NULL) {
		char errmsg[CPXMESSAGEBUFSIZE];
		fprintf(stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "%s", errmsg);
		goto TERMINATE;
	}

	status = CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
	if (status != 0) {
		fprintf(stderr, "Failure to turn on screen indicator, error %d.\n", status);
		goto TERMINATE;
	}

	CPXsetintparam(env, CPXPARAM_MIP_Interval, 1000);

	CPXsetdblparam(env, CPX_PARAM_TILIM, TIEMPO_MAXIMO_EN_SEGUNDOS);
	CPXsetdblparam(env, CPX_PARAM_EPAGAP, CORTAR_A_DISTANCIA);

	lp = CPXcreateprob(env, &status, "tp_modelosuno");

	if (lp == NULL) {
		fprintf(stderr, "Failed to create LP.\n");
		goto TERMINATE;
	}

	///----
	int* coloreoValido = (int*)malloc(Nr_vert * sizeof(int));
	int maxColor = 7; // cantidad de lavados
	buscarSolucionInicial(coloreoValido);

	status = populatebyrow(env, lp, matriz, pesos, Nr_vert, maxColor);

	int cantidadVariables = maxColor + Nr_vert * maxColor;
	int* posicionVariables = (int*)malloc(cantidadVariables * sizeof(int));
	double* valorVariables = (double*)malloc(cantidadVariables * sizeof(double));

	int i, j;
	for (j = 0; j < maxColor; j++) {
		posicionVariables[j] = j;
		valorVariables[j] = 0;
	}

	for (i = 0; i < Nr_vert; i++) {
		for (j = 0; j < maxColor; j++) {
			posicionVariables[maxColor + i * maxColor + j] = maxColor + i * maxColor + j;
			if (coloreoValido[i] == j) {
				valorVariables[maxColor + i * maxColor + j] = 1;
				if (valorVariables[j] < pesos[i]) {
					valorVariables[j] = pesos[i];
				}
			}
			else {
				valorVariables[maxColor + i * maxColor + j] = 0;
			}
		}
	}

	int zero = 0;
	int one = 1;

	status = CPXaddmipstarts(env, lp, 1, cantidadVariables, &zero, posicionVariables, valorVariables, &one, NULL);
	if (status) {
		fprintf(stderr, "Error al cargar solucion inicial.\n");
		goto TERMINATE;
	}

	///----

	status = CPXwriteprob(env, lp, "tp.lp", NULL);
	if (status) {
		fprintf(stderr, "Failed to write LP to disk.\n");
		goto TERMINATE;
	}

	double inicio;
	status = CPXgettime(env, &inicio);
	if (status)
		goto TERMINATE;

	status = resolver(env, lp);
	if (status)
		goto TERMINATE;

	double fin;
	status = CPXgettime(env, &fin);
	if (status)
		goto TERMINATE;

	double objval;
	status = CPXgetobjval(env, lp, &objval);
	if (status)
		goto TERMINATE;

	double bestobjval;
	status = CPXgetbestobjval(env, lp, &bestobjval);
	if (status)
		goto TERMINATE;

	int solstat = CPXgetstat(env, lp);

	int nodecount = CPXgetnodecnt(env, lp);
	printf("\n");
	if (solstat != 101) {
		printf("Parece que no se llego al optimo, status: %d\n", solstat);
	}
	printf("Optimo: %f\n", objval);
	printf("Cota: %f\n", bestobjval);
	printf("Nodos: %d\n", nodecount);
	printf("Tiempo: %f\n", fin - inicio);

	/*=======================================================================*/

TERMINATE:

	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}

	/* Free the CPLEX environment, if necessary */

	if (env != NULL) {
		status = CPXcloseCPLEX(&env);

		/* Note that CPXcloseCPLEX produces no output, so the only
		 way to see the cause of the error is to use
		 CPXgeterrorstring.  For other CPLEX routines, the errors
		 will be seen if the CPXPARAM_ScreenOutput parameter is set to
		 CPX_ON */

		if (status) {
			char errmsg[CPXMESSAGEBUFSIZE];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	return (status);

} /* END main */

void buscarSolucionInicial(int* coloreoValido) {
	// carga del vector a mano segun el output de entrega anterior
	coloreoValido[1 - 1] = 5;
	coloreoValido[2 - 1] = 6;
	coloreoValido[3 - 1] = 3;
	coloreoValido[4 - 1] = 3;
	coloreoValido[5 - 1] = 2;
	coloreoValido[6 - 1] = 7;
	coloreoValido[7 - 1] = 6;
	coloreoValido[8 - 1] = 3;
	coloreoValido[9 - 1] = 3;
	coloreoValido[10 - 1] = 4;
	coloreoValido[11 - 1] = 5;
	coloreoValido[12 - 1] = 2;
	coloreoValido[13 - 1] = 4;
	coloreoValido[14 - 1] = 3;
	coloreoValido[15 - 1] = 3;
	coloreoValido[16 - 1] = 2;
	coloreoValido[17 - 1] = 2;
	coloreoValido[18 - 1] = 1;
	coloreoValido[19 - 1] = 3;
	coloreoValido[20 - 1] = 3;
	coloreoValido[21 - 1] = 2;
	coloreoValido[22 - 1] = 5;
	coloreoValido[23 - 1] = 6;
	coloreoValido[24 - 1] = 5;
	coloreoValido[25 - 1] = 1;
	coloreoValido[26 - 1] = 3;
	coloreoValido[27 - 1] = 1;
	coloreoValido[28 - 1] = 1;
	coloreoValido[29 - 1] = 1;
	coloreoValido[30 - 1] = 1;
	coloreoValido[31 - 1] = 1;
	coloreoValido[32 - 1] = 1;
	coloreoValido[33 - 1] = 1;
	coloreoValido[34 - 1] = 5;
	coloreoValido[35 - 1] = 4;
	coloreoValido[36 - 1] = 1;
	coloreoValido[37 - 1] = 3;
	coloreoValido[38 - 1] = 1;
	coloreoValido[39 - 1] = 1;
	coloreoValido[40 - 1] = 1;
	coloreoValido[41 - 1] = 1;
	coloreoValido[42 - 1] = 1;
	coloreoValido[43 - 1] = 1;
	coloreoValido[44 - 1] = 1;
	coloreoValido[45 - 1] = 5;
	coloreoValido[46 - 1] = 6;
	coloreoValido[47 - 1] = 2;
	coloreoValido[48 - 1] = 2;
	coloreoValido[49 - 1] = 6;
	coloreoValido[50 - 1] = 3;
	coloreoValido[51 - 1] = 3;
	coloreoValido[52 - 1] = 2;
	coloreoValido[53 - 1] = 2;
	coloreoValido[54 - 1] = 1;
	coloreoValido[55 - 1] = 3;
	coloreoValido[56 - 1] = 3;
	coloreoValido[57 - 1] = 2;
	coloreoValido[58 - 1] = 5;
	coloreoValido[59 - 1] = 2;
	coloreoValido[60 - 1] = 2;
	coloreoValido[61 - 1] = 3;
	coloreoValido[62 - 1] = 3;
	coloreoValido[63 - 1] = 2;
	coloreoValido[64 - 1] = 2;
	coloreoValido[65 - 1] = 1;
	coloreoValido[66 - 1] = 3;
	coloreoValido[67 - 1] = 3;
	coloreoValido[68 - 1] = 2;
	coloreoValido[69 - 1] = 2;
	coloreoValido[70 - 1] = 6;
	coloreoValido[71 - 1] = 5;
	coloreoValido[72 - 1] = 1;
	coloreoValido[73 - 1] = 3;
	coloreoValido[74 - 1] = 1;
	coloreoValido[75 - 1] = 1;
	coloreoValido[76 - 1] = 1;
	coloreoValido[77 - 1] = 1;
	coloreoValido[78 - 1] = 1;
	coloreoValido[79 - 1] = 1;
	coloreoValido[80 - 1] = 1;
	coloreoValido[81 - 1] = 5;
	coloreoValido[82 - 1] = 1;
	coloreoValido[83 - 1] = 1;
	coloreoValido[84 - 1] = 1;
	coloreoValido[85 - 1] = 1;
	coloreoValido[86 - 1] = 1;
	coloreoValido[87 - 1] = 1;
	coloreoValido[88 - 1] = 1;
	coloreoValido[89 - 1] = 1;
	coloreoValido[90 - 1] = 1;
	coloreoValido[91 - 1] = 1;
	coloreoValido[92 - 1] = 1;
	coloreoValido[93 - 1] = 6;
	coloreoValido[94 - 1] = 2;
	coloreoValido[95 - 1] = 4;
}
