
#include<stdio.h>
#include<math.h>
#define PI 3.14159265
#define MAX 50

typedef struct {
	double orientation;
	double thickness;
} laminate_data;

void section_break(); // introduces a section break - double \n

void input_mat_data(double *E, double *EE, double *v, double *vv, double *G); //input material data
void input_laminate(laminate_data laminate[], int *size); //input laminate orientation
double input_laminate_data(int i, int fn); //performs the printf and scanf_s fn for getting laminate data. 1 - orientation. 2 - thickness.
void LaminateDataOutput(laminate_data arr[], int size); //output values of a 1D array (of laminate orientation and thickness)
void Qcal(double Q[][3], double E1, double E2, double v12, double v21, double G12); //calculates the Q matrix
void output_MultiDArray(double[][3]); //output 2D or 3D array
void Tcal(double Tmat[][3], double ang); //calculates [T]
void QTcal(double QTmat[][3], double ang, double Q[][3]); //calculates [Qbar] matrix

double get_center(laminate_data arr[], int size); //calculate the geometrical center of the laminate
void get_interface_coord(double hi[], double center, laminate_data laminate[], int size); //calculates the ply interface coordinates h0, h1, h2 etc
void get_lam_matrix(double QTmat[][3][3], double M[][3], double hi[], int size, int power); //sends 3D [QT] matrix to be processed to get [A] [B] and [D] 
double cal_lam_mat_comp(double Q[][3][3], double h[], int i, int j, int size, int power); //calculates the [A] [B] and [D] matrix components

int check_laminate_sym(laminate_data laminate[], int size); //check if laminate is symmetrical. Return 1 if yes, 0 if no.

void input_failure_data(double* x1, double *x2, double *y1, double *y2, double *s1); //input material failure data
void input_thermal_data(double *a, double *aa, double *dtemp); //input material data - thermal related
void cal_alpha_mat(double alpha[], double T[][3], double a, double aa); //calculate alpha matrix for kth ply
void get_NT(double NT[], double QTmat[][3][3], double alpha[][3], double hi[], double dtemp, int size); //gets {NT}
double cal_NT(double QTmat[][3][3], double alpha[][3], double hi[], int colnum, int size); //cal {NT} components


void cat_orient_type(double orient_type[], laminate_data laminate[], double alpha_type[][3], double alpha [][3], int size, int *orient_size); //catagorise all oriented plies
int get_stress(double A[][3], double NT[], double NExt[], double orientation[], double QT[][3][3], double T[][3][3], double alpha[][3], double stress[][3], double delta_temp, int size); //get local and global stress for each orientation type. Returns local stress for each orientation type.
void check_failure(double stress[][3], double orientation[], double x, double xx, double y, double yy, double s, int size); //check for failure using local stress results found in get_stress. Max stress criterion and Tsai Wu criterion used.

int get_failure(double A[][3], double orientation[], double QT[][3][3], double T[][3][3], double stress_thermal[][3], double x, double xx, double y, double yy, double s, int size); //get Nx, Ny, Nxy failure values

void max_stress1(double stress[][3][3], double orientation[], double stress_thermal[][3], double x, double xx, double y, double yy, double s, int size); //obtain failure Nx, Ny, Nxy based on Max Stress Criterion
void get_lowest_N(double fail_arr[][3][6], double failure[], int rownum, int size); //get lowest magnitude {N} component

void TsaiWu1(double stress[][3][3], double orientation[], double stress_thermal [][3], double x, double xx, double y, double yy, double s, int size); //obtain failure Nx, Ny, Nxy based on Tsai Wu Failure Criterion

void TsaiWu_coefft(double x, double xx, double y, double yy, double s, double* F1, double* F2, double* F11, double* F22, double* F66, double* F12); //get Tsai Wu Coefft terms
double mat_mul(double A[][3], double B[], int col); //cal [3x3][3x1] matrix multipliation
void mat_inverse(double invr_mat[][3], double mat[][3], int* checknull); //inverse 3 x 3 matrix

void end_prog(); //to end prog

int main(void) {

	//Info on program
	printf("Composite matrix calculator. Plane stress assemptions being held. Assumes homogeneuos material. Max %d plies. \n", MAX);
	printf("Stress analysis function only works if laminate is symmetrical in layering. Assumes {M} = 0.\n");
	printf("Creator: Laurence Lu. Ver: 27.05.2021\n");
	section_break();


	//input material data
	double E1, E2, v12, v21, G12;
	input_mat_data(&E1, &E2, &v12, &v21, &G12);
	section_break();
	if (E1 == 0) {
		end_prog();
		return 0;
	}
	if (v12 * v21 == 1) {
		printf("Error! v12 x v21 = 1. \n");
		end_prog();
		return 0;
	}
	//check material data values
	//printf("Values of E1, E2, v12, v21 and G12 is %.4f, %.4f, %.4f, %.4f, %.4f respectively\n", E1, E2, v12, v21, G12);


	//input laminate data
	int size; //total number of plies
	laminate_data laminate[MAX] = { 0 }; //unable to do variable size array e.g array[size]
	input_laminate(laminate, &size);
	//safety feature to prevent out of bounds arr
	if (size > MAX) {
		end_prog();
		return 0;
	}
	printf("\n");

	//check laminate data values
	LaminateDataOutput(laminate, size);
	section_break();


	//calculate [Q] matrix
	double Q[3][3] = { 0 };
	Qcal(Q, E1, E2, v12, v21, G12);
	printf("[Q] matrix is (GPa): \n");
	output_MultiDArray(Q);
	section_break();


	//calculate [T] matrix
	double T[MAX][3][3] = { 0 };
	/* 3D array works like an array of 2D array
	   T[<T matrix corr to which ply>][3][3]
	   Same logic for QT matrix later */
	int count;
	for (count = 0;count < size;count++) {
		Tcal(T[count], laminate[count].orientation); //double Tmat[][3] is 2D array therefore, passing one element of an array of 2D arr --> T[count]. Same logic for [QT] later. Processing one 2D arr at a time.
	}

	//check [T] matrix
	count = 0;
	for (count = 0; count < size;count++) {
		printf("[T] for %.2f deg ply: \n", laminate[count].orientation);
		output_MultiDArray(T[count]);
		printf("\n");
	}
	printf("\n");


	//calculate [Qbar] matrix. Denote as [QT]
	double QT[MAX][3][3] = { 0 };
	count = 0;
	for (count = 0;count < size;count++) {
		QTcal(QT[count], laminate[count].orientation, Q);
	}

	//check [QT] matrix
	count = 0;
	for (count = 0;count < size;count++) {
		printf("[Qbar] for %.2f deg ply (GPa): \n", laminate[count].orientation);
		output_MultiDArray(QT[count]);
		printf("\n");
	}
	printf("\n");



	//Calculate geometrical center of laminate
	double center;
	center = get_center(laminate, size);
	
	//check center
	//printf("Center calculated is %.2f \n", center);

	//Cal ply interface coordinates
	double h[MAX+1] = { 0 };
	get_interface_coord(h, center, laminate, size);

	//check ply interface coordinates
	/* count = 0;
	printf("Ply interface coordinates: \n");
	for (count = 0;count <= size;count++) {
		printf("%.2f ", h[count]);
	}
	printf("\n"); */

	//Calculate [A]
	double A[3][3] = { 0 };
	get_lam_matrix(QT, A, h, size, 1);

	printf("[A] is (GPa mm): \n");
	output_MultiDArray(A);
	section_break();

	//Calculate [B]
	double B[3][3] = { 0 };
	get_lam_matrix(QT, B, h, size, 2);

	printf("[B] is (GPa mm^2): \n");
	output_MultiDArray(B);
	section_break();

	//Calculte [D]
	double D[3][3] = { 0 };
	get_lam_matrix(QT, D, h, size, 3);

	printf("[D] is (GPa mm^3): \n");
	output_MultiDArray(D);
	section_break();

	//check if laminate is symmetrical
	int sym_check;
	sym_check = check_laminate_sym(laminate, size);
	if (sym_check == 0) {
		//proceed to end program if non-symmetrical laminate
		end_prog();
		return 0;
	}
	
	//check to proceed for stress analysis
	char checkcont1;
	printf("Perform stress analysis? Y - yes. N - no. --- ");
	scanf_s(" %c", &checkcont1, 1);

	if (checkcont1 == 'Y' || checkcont1 == 'y') {

		double NT[3] = { 0 }; //placed here so var still valid after 'if' statement (likewise for delta temp & {alpha}). NT[0] = NxT, NT[1] = NyT, NT[2] = NxyT
		double delta_temp; //delta temp = room temp - curing temp
		delta_temp = 0;
		double alpha[MAX][3] = { 0 }; //AlphaX, AlphaY and AlphaXY of kth ply. Each ply stored in k-th row. AX in col 0. AY in col 1. AXY in col 2.

		char checkcont2;
		printf("Include thermal residual stresses? Y - yes. N - no. --- ");
		scanf_s(" %c", &checkcont2, 1);
		printf("\n");
		if (checkcont2 == 'Y' || checkcont2 == 'y') {
			double a1, a2; //a is alpha

			input_thermal_data(&a1, &a2, &delta_temp); //input thermal data

			//Cal {Alpha} for k-th ply 
			count = 0;
			for (count = 0; count < size;count++) {
				cal_alpha_mat(alpha[count], T[count], a1, a2); //performing one 1D arr at a time

				//check k-th ply Alpha
				printf("{Alpha} values for %.2f deg ply (x10^-6/C):\n", laminate[count].orientation);
				printf("%.4f %.4f %.4f\n", alpha[count][0], alpha[count][1], alpha[count][2]);
				printf("\n");
			}
			
			//Cal {NT} - Stress Resultants due to thermal
			get_NT(NT, QT, alpha, h, delta_temp, size);
		}
		//check {NT}
		printf("Value of Nx(T), Ny(T), Nxy(T) (MPa mm):\n");
		printf("%.5f %.5f %.5f", NT[0], NT[1], NT[2]);
		section_break();

		//check to see what function to perform - finding safety {N} range or Cal stresses in laminate
		int checkcont3;
		printf("Select type of stress analysis.\n");
		printf("1 - stresses in laminate. 2 - {N} safety range of laminate, assumes uniaxially loaded laminate, no moments. --- ");
		scanf_s("%d", &checkcont3);
		printf("\n");

		//catagorise all same orientated plies into a single array --> save computing power
		double orient_type[MAX / 2] = { 0 }; // [MAX / 2] as symmetric laminate, therefore, realistically, at most MAX/2 different types of unique orientation. Gets orientation types.
		double alpha_type[MAX / 2][3] = { 0 }; //used for cal mechanical strain in get_stress later
		int orient_size; //number of different types of orientations
		cat_orient_type(orient_type, laminate, alpha_type, alpha, size, &orient_size);

		double QT_fail_type[MAX / 2][3][3] = { 0 };
		double T_fail_type[MAX / 2][3][3] = { 0 };
		count = 0;
		for (count = 0;count < orient_size;count++) {
			QTcal(QT_fail_type[count], orient_type[count], Q);
			Tcal(T_fail_type[count], orient_type[count]);
		}
		/*check T_fail and QT_fail matrix along with types of orientation
			count = 0;
			for (count = 0;count < orient_size;count++) {
				printf("Orientation: %.2f\n", orient_type[count]);
				printf("[QT]:\n");
				output_MultiDArray(QT_fail_type[count]);
				printf("\n");
				printf("[T]:\n");
				output_MultiDArray(T_fail_type[count]);
				printf("{alpha}:\n");
				printf("%.4f %.4f %.4f\n", alpha_type[count][0], alpha_type[count][1], alpha_type[count][2]);
				printf("\n");
			}
			printf("Delta temp: %.1f\n", delta_temp);
			section_break(); */
		
		switch (checkcont3) {
		case 1: { //brackets w/in switch statement incld for to allow definition of variables within the switch statement
			//stress in laminate
			double stress_local[MAX / 2][3] = { 0 }; // stress_local[k][0] = 1-dir   stress_local[k][1] = 2-dir   stress_local[k][2] = 12-dir
			//input {N}
			double NExt[3] = { 0 };
			printf("NOTE: To calculate thermal residual stress, enter 0 for Nx, Ny and Nxy.\n");
			printf("Enter Nx (MPa mm): ");
			scanf_s("%lf", &NExt[0]);
			printf("Enter Ny (MPa mm): ");
			scanf_s("%lf", &NExt[1]);
			printf("Enter Nxy (MPa mm): ");
			scanf_s("%lf", &NExt[2]);
			printf("\n");

			int check;
			check = get_stress(A, NT, NExt, orient_type, QT_fail_type, T_fail_type, alpha_type, stress_local, delta_temp, orient_size);
			if (check != 0) {
				end_prog();
				return 0;
			}
			//check for failure
			char checkcont4;
			printf("Check for failure? Y- yes. N - no. --- ");
			scanf_s(" %c", &checkcont4, 1);
			printf("\n");
			if (checkcont4 == 'Y' || checkcont4 == 'y') {
				//input failure data of laminate
				double x, xx, y, yy, s; // x = X (+ve value), xx = X' (-ve value), y = Y (+ve value), yy = Y' (-ve value) and s = S (+ve value)
				input_failure_data(&x, &xx, &y, &yy, &s);
				check_failure(stress_local, orient_type, x, xx, y, yy, s, orient_size);
			}
			section_break();
			break;
		}
		case 2: { //brackets w/in switch statement incld for to allow definition of variables within the switch statement
			//{N} safety range

			//input failure data of laminate
			double x, xx, y, yy, s; // x = X (+ve value), xx = X' (-ve value), y = Y (+ve value), yy = Y' (-ve value) and s = S (+ve value)
			input_failure_data(&x, &xx, &y, &yy, &s);

			//get thermal stress residual
			double NExt[3] = { 0 };
			int check;
			double stress_local_thermal[MAX / 2][3] = { 0 }; //operates similar to stress_local in case 1
			if (checkcont2 == 'Y' || checkcont2 == 'y') { //if incld thermal resideual --> need to cal thermal residual stresses
				check = get_stress(A, NT, NExt, orient_type, QT_fail_type, T_fail_type, alpha_type, stress_local_thermal, delta_temp, orient_size);
				if (check != 0) {
					end_prog();
					return 0;
				}
			}

			check = get_failure(A, orient_type, QT_fail_type, T_fail_type, stress_local_thermal, x, xx, y, yy, s, orient_size);
			if (check != 0) {
				end_prog();
				return 0;
			}

			section_break();
			break;
		}
		}
		
	}

	end_prog();
	return 0;
}
//======================================================================================================================================================================================================================

//introduces a section break between sections
void section_break() {
	printf(" \n");
	printf(" \n");
}

//input material data
void input_mat_data(double *E, double *EE, double *v, double *vv, double *G) {
	//E = E1; EE = E2; v = v12; vv = v21; G = G12
	double E1, E2, v12, v21, G12;
	printf("Enter E1 value (GPa): ");
	scanf_s("%lf", &E1);

	printf("Enter E2 value (GPa): ");
	scanf_s("%lf", &E2);

	printf("Enter G12 value (GPa): ");
	scanf_s("%lf", &G12);

	*E = E1 * 1.0;
	*EE = E2 * 1.0;
	*G = G12 * 1.0;

	if (E1 == 0) {
		printf("Error! E1 cannot be zero!\n");
		return;
	}

	printf("Enter v12 value: ");
	scanf_s("%lf", &v12);
	*v = v12 * 1.0;

	v21 = (E2 * v12) / E1;
	*vv = v21;
}

//input laminate orientation and thickness
void input_laminate(laminate_data laminate[], int* size) {
	//input total number of plies / layers
	int s;
	printf("Enter total number of plies: ");
	scanf_s("%d", &s);
	*size = s;
	if (s > MAX) {
		printf("Exceed maximum allowed plies. Ending program!\n");
		printf("\n");
		return;
	}
	//input method of entering laminate data
	char ip_type;
	printf("Enter input method. S - Symmetrical laminate input, enter half layering. T - Total laminate input, enter full layering.\n");
	scanf_s(" %c", &ip_type, 1);
	printf("\n");
	int i;
	if (ip_type == 'S' || ip_type == 's') {
		int j = s - 1;
		printf("Enter data for first half of laminate in order.\n");

		for (i = 0;i < s / 2;i++) {
			//An exercise of filling in a symmetrical 1D array by just filling the first half and copying it to the second half
			laminate[i].orientation = input_laminate_data(i, 1);
			laminate[i].thickness = input_laminate_data(i, 2);

			laminate[j].orientation = laminate[i].orientation;
			laminate[j].thickness = laminate[i].thickness;
			j--;
		}
		if (s % 2 != 0) {
			//for odd number symmetricl laminates e.g [90/20/90]T
			printf("Ply %d orientation (deg): ", (s / 2) + 1);
			scanf_s("%lf", &laminate[(s / 2) ].orientation);
			printf("Ply %d thickness (mm): ", (s / 2) + 1);
			scanf_s("%lf", &laminate[(s / 2) ].thickness);
		}
	}
	else {
		//For filling in Total Laminate input
		printf("Enter ply data in order:\n");
		for (i = 0; i < s; i++) {
			laminate[i].orientation = input_laminate_data(i, 1);
			laminate[i].thickness = input_laminate_data(i, 2);
		}
	}
}

//performs the printf and scanf_s fn for getting laminate data. 1 - orientation. 2 - thickness.
double input_laminate_data(int i, int fn) {
	double rtn;
	switch (fn) {
	case 1: //ply orientation (deg)
		printf("Ply %d orientation (deg): ", i + 1);
		scanf_s("%lf", &rtn);
		return rtn;
		break;
	case 2: //ply thickness (mm)
		printf("Ply %d thickness (mm): ", i + 1);
		scanf_s("%lf", &rtn);
		return rtn;
		break;
	}
	return -1;
}

//output 1D arrary values (the laminate data of orientation and thickness)
void LaminateDataOutput(laminate_data arr[], int size) {
	printf("Total number of plies: %d\n", size);

	int i;
	printf("Ply orientation in deg (in order):\n");
	for (i = 0;i < size;i++) {
		printf("%.2f ", arr[i].orientation);
	}
	printf("\n");

	int j;
	printf("Ply thickness in mm (in order): \n");
	for (j = 0;j < size;j++) {
		printf("%.4f ", arr[j].thickness);
	}
}

//Calculates [Q] matrix
void Qcal(double Q[][3], double E1, double E2, double v12, double v21, double G12) {
	Q[0][0] = E1 / (1.0 - v12 * v21);
	Q[1][1] = E2 / (1.0 - v12 * v21);
	Q[2][2] = G12;

	Q[0][1] = (E2 * v12) / (1 - v12 * v21);
	Q[1][0] = Q[0][1];
}

//Output multi dimensional array --> if 2D, call the matrix name e.g Q; if 3D call matrix array e.g T[ply]
void output_MultiDArray(double arr[][3]) {
	int i, j;

	for (i = 0; i < 3;i++) {
		for (j = 0;j < 3;j++) {
			printf("%.4f ", arr[i][j]);
		}
		printf("\n");
	}
}

//calculates [T] matrix
void Tcal(double Tmat[][3], double ang) {
	double m, n;
	m = cos(ang * PI / 180);
	n = sin(ang * PI / 180);

	Tmat[0][0] = pow(m, 2);
	Tmat[0][1] = pow(n, 2);
	Tmat[0][2] = 2.0 * m * n;

	Tmat[1][0] = Tmat[0][1];
	Tmat[1][1] = Tmat[0][0];
	Tmat[1][2] = -2.0 * m * n;

	Tmat[2][0] = -n * m;
	Tmat[2][1] = n * m;
	Tmat[2][2] = pow(m, 2) - pow(n, 2);
}

//calculates [Qbar] matrix for nth ply
void QTcal(double QTmat[][3], double ang, double Q[][3]) {
	//cal m and n values
	double m, n;
	m = cos(ang * PI / 180);
	n = sin(ang * PI / 180);

	//cal actual [QT] component by component
	QTmat[0][0] = Q[0][0] * pow(m, 4) + Q[1][1] * pow(n, 4) + pow(n * m, 2) * (2 * Q[0][1] + 4 * Q[2][2]);
	QTmat[1][1] = Q[0][0] * pow(n, 4) + Q[1][1] * pow(m, 4) + pow(n * m, 2) * (2 * Q[0][1] + 4 * Q[2][2]);

	QTmat[0][1] = pow(m * n, 2) * (Q[0][0] + Q[1][1] - 4 * Q[2][2]) + (pow(m, 4) + pow(n, 4)) * Q[0][1];
	QTmat[1][0] = QTmat[0][1];

	QTmat[2][2] = pow(m * n, 2) * (Q[0][0] + Q[1][1] - 2 * Q[0][1]) + pow(m * m - n * n, 2) * Q[2][2];

	QTmat[0][2] = pow(m, 3) * n * Q[0][0] - m * pow(n, 3) * Q[1][1] + (m * pow(n, 3) - pow(m, 3) * n) * (Q[0][1] + 2 * Q[2][2]);
	QTmat[2][0] = QTmat[0][2];

	QTmat[1][2] = pow(n, 3) * m * Q[0][0] - n * pow(m, 3) * Q[1][1] + (n * pow(m, 3) - pow(n, 3) * m) * (Q[0][1] + 2 * Q[2][2]);
	QTmat[2][1] = QTmat[1][2];
}

//calculates the geometrical center of laminate
double get_center(laminate_data arr[], int size) {

	int i;
	double ctr;
	ctr = 0.0;
	// 1. Find total thickness
	for (i = 0;i < size;i++) {
		ctr = ctr + arr[i].thickness;
	}
	// 2. Divide total thickness by 2
	ctr = ctr / 2;

	return ctr;
}

//Calculate interface coordinate h0, h1, h2 etc
void get_interface_coord(double hi[], double center, laminate_data laminate[], int size) {

	/*
	h0 ----------------------- (-ctr)
	h1 ----------------------- (h0 + laminate[0].thickness)
	h2 ----------------------- (h1 + laminate[1].thickness)
	.
	.
	ctr-----------------------             |                                  Note: ply 1 = laminate[0]
	.                                      |                                        ply k = laminate[k-1]
	.                                      |
	h(k-2)--------------------            \ / +ve z-dir
	h(k-1)--------------------             v
	h(k)---------------------- (+ctr)
	
	*/

	hi[0] = -1.0 * center;

	int i;
	for (i = 1;i < size;i++) {
		hi[i] = hi[i - 1] + laminate[i - 1].thickness;
	}

	hi[size] = center;
}

//sends 3D QT matrix to be processed to get [A] [B] and [D]
void get_lam_matrix(double QTmat[][3][3], double M[][3], double hi[], int size, int power) {

	int i, j;
	//cal one component of the [A]/[B]/[D] matrix at a time
	for (i = 0;i < 3;i++) {
		for (j = 0; j < 3;j++) {
			M[i][j] = cal_lam_mat_comp(QTmat, hi, i, j, size, power);
		}
	}
}

//calculates individual [A] [B] and [D] components
double cal_lam_mat_comp(double Q[][3][3], double h[], int i, int j, int size, int power) {
	int counter;
	double M;
	//based on formula --> summation of products (see formula in ME4212 notes)
	M = 0;
	for (counter = 0; counter < size; counter++) {
		M += (1.0 / power) * (Q[counter][i][j]) * (pow(h[counter + 1], power) - pow(h[counter], power));
	}
	return M;
}

//check if laminate is symmetrical. Return 1 if yes, 0 if no
int check_laminate_sym(laminate_data laminate[], int size) {
	int i,j;
	j = size - 1;
	//An exercise of checking for symmetrical 1D arr
	for (i = 0;i < size;i++) {
		if (j <= i) {
			break; //improve speed of program
		}
		if (laminate[i].thickness != laminate[j].thickness || laminate[i].orientation != laminate[j].orientation) {
			return 0;
		}
		j--;
	}
	return 1;
}

//input material failure data
void input_failure_data(double* x1, double* x2, double* y1, double* y2, double* s1) {
	double x, xx, y, yy, s;

	printf("Enter X -- X-dir tensile stress failure (MPa): ");
	scanf_s("%lf", &x);

	printf("Enter X' -- X-dir compressive stress failure (MPa): ");
	scanf_s("%lf", &xx);
	if (xx > 0) {
		xx *= -1.0;
	}

	printf("Enter Y -- Y-dir tensile stress failure (MPa): ");
	scanf_s("%lf", &y);

	printf("Enter Y' -- Y-dir compressive stress failure (MPa): ");
	scanf_s("%lf", &yy);
	if (yy > 0) {
		yy *= -1.0;
	}

	printf("Enter S -- shear stress failure (MPa): ");
	scanf_s("%lf", &s);

	*x1 = x;
	*x2 = xx;
	*y1 = y;
	*y2 = yy;
	*s1 = s;

	printf("\n");
}


//input thermal data
void input_thermal_data(double *a, double *aa, double *dtemp) {
	double a1, a2, ctp, rtp;

	printf("Enter Alpha 1 value (x10^-6 /C): ");
	scanf_s("%lf", &a1);

	printf("Enter Alpha 2 value (x10^-6 /C): ");
	scanf_s("%lf", &a2);

	*a = a1;
	*aa = a2;

	printf("Enter Curing Temperature (deg C): ");
	scanf_s("%lf", &ctp);

	printf("Enter Room Temperature (deg C): ");
	scanf_s("%lf", &rtp);

	*dtemp = 1.0*(rtp - ctp);

	printf("\n");
}

//cal alpha matrix for k-th ply
void cal_alpha_mat(double alpha[], double T[][3], double a, double aa) {

	//see formula for why this is done

	//cal AlphaX col 0
	alpha[0] = T[0][0] * a + T[0][1] * aa;

	//cal AlphaY col 1
	alpha[1] = T[1][0] * a + T[1][1] * aa;

	//cal AlphaXY col 2 (note formula {a}=[T]{a1 a2} involves 0.5a(xy)
	alpha[2] = 2.0 * (T[2][0] * a * (-1.0) + T[2][1] * aa * (-1.0)); //(-1.0) as dealing with [T]^-1 which is basically for ang = X, sin(-x) = -sin(x) and cos(-x) = cos(x)
}

//get {NT} in MPa mm
void get_NT(double NT[], double QTmat[][3][3], double alpha[][3], double hi[], double dtemp, int size) {
	int counter;
	for (counter = 0;counter < 3;counter++) {
		NT[counter] = dtemp * cal_NT(QTmat, alpha, hi, counter, size) * pow(10, -6) * pow(10, 3); //similar to getting [A], [B], [D]. Cal one component of the matrix at a time; pow(10,-6) due to units of {alpha}; pow(10,3) due to GPa mm -> MPa mm conversion of [QT]
	}
}

//cal {NT} components
double cal_NT(double QTmat[][3][3], double alpha[][3], double hi[], int colnum, int size) {
	double m;
	int k;
	m = 0.0;
	for (k = 0;k < size;k++) {
		m += mat_mul(QTmat[k], alpha[k], colnum) * (hi[k + 1] - hi[k]); //for the summation part
	}
	return m;
}

//catagorise all oriented plies
void cat_orient_type(double orient_type[], laminate_data laminate[], double alpha_type[][3], double alpha[][3], int size, int *orient_size) {
	int i, j, k, check;
	j = 1;
	orient_type[0] = laminate[0].orientation;
	alpha_type[0][0] = alpha[0][0];
	alpha_type[0][1] = alpha[0][1];
	alpha_type[0][2] = alpha[0][2];
	for (i = 1;i < (size / 2) + 1;i++) {
		check = 0;
		for (k = 0; k < j; k++) {
			if (laminate[i].orientation == orient_type[k]) {
				check = 1;
				break;
			}
		}
		if (check == 0) {
			orient_type[j] = laminate[i].orientation;
			alpha_type[j][0] = alpha[i][0];
			alpha_type[j][1] = alpha[i][1];
			alpha_type[j][2] = alpha[i][2];
			j++;
		}
	}
	*orient_size = j;
}

//get local and global stress for each orientation type. Returns local stress for each orientation type
int get_stress(double A[][3], double NT[], double NExt[], double orientation[], double QT[][3][3], double T[][3][3], double alpha[][3], double stress[][3], double delta_temp, int size) {
	//get [A]^-1
	double A_inverse[3][3] = { 0 };
	int checknull;
	mat_inverse(A_inverse, A, &checknull);
	//safety feat to prevent division by 0 subsequently
	if (checknull == -1) {
		return -1;
	}
	

	//get ext {e}
	double e_ext[3] = { 0 }; //const throughout as laminate sym --> [B] = 0 --> {e}k = {eo} = {e}
	int i, k;
	for (i = 0;i < 3;i++) {
		e_ext[i] = mat_mul(A_inverse, NExt, i) * pow(10, -3); //input {N} is MPa mm
	}
	/*check ext{e}
	i = 0;
	for (i = 0;i < 3;i++) {
		printf("%.5f\n", e_ext[i]);
	}*/

	//get thermal {e}
	double e_thermal[3] = { 0 }; //const throughout as laminate sym --> [B] = 0 --> {e}k = {eo} = {e}
	i = 0;
	for (i = 0;i < 3;i++) {
		e_thermal[i] = mat_mul(A_inverse, NT, i) * pow(10, -3); //[NT] is in MPa mm
	}

	//get mechanical {e}
	double e_mech[MAX / 2][3] = { 0 };
	/*
		{e_mechancial} for k plies
			=
	ex		ey		y(xy)	-	k-th ply
	*/
	i = 0;
	k = 0;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			e_mech[k][i] = e_thermal[i] - delta_temp * alpha[k][i] * pow(10, -6);
		}
	}

	//get total {e} resulting in stress
	double e_total[MAX / 2][3] = { 0 };
	i = 0;
	k = 0;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			e_total[k][i] = e_mech[k][i] + e_ext[i];
			if (i == 2) {
				e_total[k][i] *= 0.5;
			}
		}
	}
	/*check total {e}
	i = 0;
	k = 0;
	for (k = 0;k < size;k++) {
		printf("For %.2f ply:\n", orientation[k]);
		for (i = 0;i < 3;i++) {
			printf("%.6f ", e_total[k][i]);
		}
		printf("\n");
	}*/

	//get stress wrt global coord
	double stress_global[MAX / 2][3] = { 0 };
	/*
	{stress} for k plies
			=
	x		y		xy	-	k-th ply
	*/
	i = 0;
	k = 0;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			stress_global[k][i] = mat_mul(QT[k], e_total[k], i) * pow(10, 3); //[QT] GPa --> MPa
		}
	}

	//get stress wrt local coord
	double stress_local[MAX / 2][3] = { 0 };
	i = 0;
	k = 0;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			stress_local[k][i] = mat_mul(T[k], stress_global[k], i);
			stress[k][i] = stress_local[k][i]; //send back to main function
		}
	}

	//output stress wrt global coord
	i = 0;
	k = 0;
	printf("STRESS WRT TO LAMINATE / GLOBAL COORD:\n");
	for (k = 0;k < size;k++) {
		printf("For %.2f ply:\n", orientation[k]);
		printf("In x-dir: %.2f MPa\n", stress_global[k][0]);
		printf("In y-dir: %.2f MPa\n", stress_global[k][1]);
		printf("In xy-dir: %.2f MPa\n", stress_global[k][2]);
		printf("\n");
	}
	printf("\n");
	//output stress wrt local coord
	i = 0;
	k = 0;
	printf("STRESS WRT TO PLY / LOCAL COORD:\n");
	for (k = 0;k < size;k++) {
		printf("For %.2f ply:\n", orientation[k]);
		printf("In 1-dir: %.2f MPa\n", stress_local[k][0]);
		printf("In 2-dir: %.2f MPa\n", stress_local[k][1]);
		printf("In 12-dir: %.2f MPa\n", stress_local[k][2]);
		printf("\n");
	}

	return 0;
}

//check for failure using local stress results found in get_stress. Max stress criterion and Tsai Wu criterion used.
void check_failure(double stress[][3], double orientation[], double x, double xx, double y, double yy, double s, int size) {
	//NOTE: stress[][3] refers to local stress of k-th ply
	/*
	{stress} for k plies
			=
	1		2		12	 - k-th ply
	*/

	//check using Max Stress Criterion
	int k;
	int check;
	
	printf("MAXIMUM STRESS CRITERION FAILURE CHECK:\n");
	for (k = 0;k < size;k++) {
		check = 0;
		printf("For %.2f deg plies: ", orientation[k]);
		if (stress[k][0] > x || stress[k][0] < xx) {
			printf("FAILED 1-dir. Stress: %.2f MPa exceed!  ", stress[k][0]);
			check = 1;
		}
		if (stress[k][1] > y || stress[k][1] < yy) {
			printf("FAILED 2-dir. Stress: %.2f MPa exceed!  ", stress[k][1]);
			check = 1;
		}
		if (fabs(stress[k][2]) > s) {
			printf("FAILED 12-dir. Stress: %.2f MPa exceed!  ", stress[k][2]);
			check = 1;
		}
		if (check == 0) {
			printf("PASSED");
		}
		printf("\n");
	}
	printf("\n");

	//check using Tsai Wu Failure Criterion
	//getting Tsai Wu coefft
	//xx and yy are auto converted to -ve values
	double F1, F2, F11, F22, F66, F12;
	TsaiWu_coefft(x, xx, y, yy, s, &F1, &F2, &F11, &F22, &F66, &F12);
	double TsaiWuValue;

	printf("TSAI WU FAILURE CRITERION CHECK:\n");
	printf("**Assumed F12 = -0.5sqrt(F11 x F22)**\n");
	k = 0;
	for (k = 0;k < size;k++) {
		printf("For %.2f deg plies: ", orientation[k]);
		TsaiWuValue = F1 * stress[k][0] + F2 * stress[k][1] + F11 * pow(stress[k][0], 2) + F22 * pow(stress[k][1], 2) + F66 * pow(stress[k][2], 2) + 2 * F12 * stress[k][0] * stress[k][1];
		if (TsaiWuValue < 1) {
			printf("PASSED");
		}
		else {
			printf("FAILED");
		}
		printf("\n");
	}
}



//get Nx, Ny, Nxy failure values (in MPa mm)
int get_failure(double A[][3], double orientation[], double QT[][3][3], double T[][3][3], double stress_thermal [][3], double x, double xx, double y, double yy, double s, int size) {
	//get [A]^-1
	double A_inverse[3][3] = { 0 };
	int checknull;
	mat_inverse(A_inverse, A, &checknull);
	//safety feat to prevent division by 0 subsequently
	if (checknull == -1) {
		return -1;
	}

	/*check A_inverse
	printf("[A]^-1 is: \n");
	int i, j;
	for (i = 0;i < 3;i++) {
		for (j = 0;j < 3;j++) {
			printf("%.6f ", A_inverse[i][j]);
		}
		printf("\n");
	}*/

	//get {e} (midplane strain). Note, since [B] = 0, symmetrical laminate, {e} of laminate = midplane strain

	/* Taking Nx, Ny and Nxy as 1, so will only be cal the coefft and manipulating the coefft to get the ans
	{e} and {stress} will be 3 x 3 matrix. Using {e} as an example ({stress} will follow an eqv method of storing the data):
			{e}
			 =
	ex		ey		.5y(xy)		-	Nx
	ex		ey		.5y(xy)		-	Ny
	ex		ey		.5y(xy)		-	Nxy
	*/

	double e[3][3] = { 0 };
	double N_temp[3] = { 0 };
	//cal one {N} component at a time
	int i, j, k;
	for (i = 0;i < 3;i++) {
		N_temp[0] = 0;
		N_temp[1] = 0;
		N_temp[2] = 0;

		N_temp[i] = 1;
		for (j = 0;j < 3;j++) {
			e[i][j] = mat_mul(A_inverse, N_temp, j); //no need x10^-3 as this will be cancelled off by the x10^3 of [QT] later
			if (j == 2) {
				e[i][j] *= 0.5;
			}
			//printf("%.6f ", e[i][j]);
		}
		//printf("\n");
	}
	
	/*
	{stress} = [orientation][{N}][{coord}]
	Therefore, for k-th ply
		{Stress}
			=
	X/1		Y/2		XY/12		-	Nx
	X/1		Y/2		XY/12		-	Ny
	X/1		Y/2		XY/12		-	Nxy
	*/
	
	//get {stress} wrt global coord
	double stress_global[MAX / 2][3][3] = { 0 }; //stress[<orientation type>][{N} = 3][{e} = 3]
	i = 0;
	j = 0;
	for (k = 0;k < size;k++) { //cal one ply, k-th ply, at a time
		for (i = 0;i < 3;i++) { //cal one {N} component at a time
			for (j = 0;j < 3;j++) {
				stress_global[k][i][j] = mat_mul(QT[k], e[i],j);
			}
		}
	}
	//get {stress} wrt local coord
	double stress_local[MAX / 2][3][3] = { 0 };
	i = 0;
	j = 0;
	k = 0;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			for (j = 0;j < 3;j++) {
				stress_local[k][i][j] = mat_mul(T[k], stress_global[k][i], j); //sending the Nx/Ny/Nxy of stress_global one {N} component at a time. i.e sending only the row arr of {stress} of k-th ply.
			}
		}
	}
	/*check {stress}
	i = 0;
	j = 0;
	k = 0;
	for (k = 0;k < size;k++) {
		printf("For %.2f ply, local stress in terms of Nx, Ny and Nxy (MPa):\n", orientation[k]);
		for (i = 0;i < 3;i++) {
			printf("%.6f Nx   %.6f Ny   %.6f Nxy\n", stress_local[k][0][i], stress_local[k][1][i], stress_local[k][2][i]);
		}
		printf("\n");
	}
	printf("\n");*/

	//check Max Stress Criterion
	max_stress1(stress_local, orientation, stress_thermal, x, xx, y, yy, s, size);

	//check Tsai-Wu Failure Criterion
	TsaiWu1(stress_local, orientation, stress_thermal, x, xx, y, yy, s, size);

	return 0;
}

//obtain failure Nx, Ny, Nxy based on Max Stress Criterion [MAX STRESS CRITERION]
void max_stress1(double stress[][3][3], double orientation[], double stress_thermal[][3], double x, double xx, double y, double yy, double s, int size) {
	/*
	{stress} = [orientation][{N}][{coord}]
	Therefore, for k-th ply
		{Stress}
			=
	1		2		12		-	Nx
	1		2		12		-	Ny
	1		2		12		-	Nxy
	*/

	double fail_arr[MAX/2][3][6] = { 0 };
	/*
	For k-th ply
		{Failure}
			=
	x	xx	y	yy	s	-s	-	Nx
	x	xx	y	yy	s	-s	-	Ny
	x	xx	y	yy	s	-s	-	Nxy
	RECALL: abs(t12) < s
	*/

	//need to adjust crit stresses to take into acct thermal residual stress. Remove thermal stresses from crit stress --> divide by coefft of {N}
	/*
		xx < stress(ext) + stress(thermal) < x
		xx - stress(thermal) < stress(ext) < x - stress(thermal)
		NOTE, stress(ext) is coefft of relevant {N} term
	*/
	int i, k;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			fail_arr[k][i][0] = (x - stress_thermal[k][0]) / stress[k][i][0];
			fail_arr[k][i][1] = (xx - stress_thermal[k][0]) / stress[k][i][0];
			fail_arr[k][i][2] = (y - stress_thermal[k][1]) / stress[k][i][1];
			fail_arr[k][i][3] = (yy - stress_thermal[k][1]) / stress[k][i][1];
			fail_arr[k][i][4] = (s - stress_thermal[k][2]) / stress[k][i][2];
			fail_arr[k][i][5] = fail_arr[k][i][4] * (-1.0);
		}
	}

	/*check fail_arr
	k = 0;
	for (k = 0;k < size;k++) {
		printf("For %.2f ply, safety range of Nx, Ny and Nxy (MPa mm):\n", orientation[k]);
		printf("%.2f %.2f (<- X and X')  %.2f %.2f (<-Y and Y')  %.2f (<-S) - Nx\n", fail_arr[k][0][0], fail_arr[k][0][1], fail_arr[k][0][2], fail_arr[k][0][3], fail_arr[k][0][4]);
		printf("%.2f %.2f (<- X and X')  %.2f %.2f (<-Y and Y')  %.2f (<-S) - Ny\n", fail_arr[k][1][0], fail_arr[k][1][1], fail_arr[k][1][2], fail_arr[k][1][3], fail_arr[k][1][4]);
		printf("%.2f %.2f (<- X and X')  %.2f %.2f (<-Y and Y')  %.2f (<-S) - Nxy\n", fail_arr[k][2][0], fail_arr[k][2][1], fail_arr[k][2][2], fail_arr[k][2][3], fail_arr[k][2][4]);
		printf("\n");
	} */
	
	//obtain lowest magnitude stress resultant (MPa mm) and the controlling laminate
	double fail_results[3][2] = { 0 }; //crit {N} components for overall laminate
	/*
		{fail_results}
			=
	+ve No.		-ve No.		-	Nx
	+ve No.		-ve No.		-	Ny
	+ve No.		-ve No.		-	Nxy
	*/

	i = 0;
	for (i = 0;i < 3;i++) { //doing one {N} component at a time
		get_lowest_N(fail_arr, fail_results[i], i, size);
	}

	//output results
	printf("**Assuming uniaxially loaded symmetrical laminate**\n");
	printf("MAXIMUM STRESS CRITERION - Stress Resultant estimated safety range (MPa mm):\n");
	printf("%.2f < Nx < %.2f\n", fail_results[0][1], fail_results[0][0]);
	printf("%.2f < Ny < %.2f\n", fail_results[1][1], fail_results[1][0]);
	printf("%.2f < Nxy < %.2f\n", fail_results[2][1], fail_results[2][0]);
	printf("\n");
	printf("\n");
}

//get lowest magnitude {N} component
void get_lowest_N(double fail_arr[][3][6], double failure[], int rownum, int size) {
	int k, i;
	//pre-initialise failure arr first using X and X' values to allow for comparison. One of the values will always be +ve and one will always be -ve
	if (fail_arr[0][rownum][0] > 0) {
		failure[0] = fail_arr[0][rownum][0];
		failure[1] = fail_arr[0][rownum][1];
	} 
	else {
		failure[0] = fail_arr[0][rownum][1];
		failure[1] = fail_arr[0][rownum][0];
	}

	for (k = 0;k < size;k++) {
		//Since for x xx, y yy, s -s, If one value is +ve, the other will be negative, will be comparing 2 at a time. See {failure} comments.
		for (i = 0;i < 3;i++) {
			if (fail_arr[k][rownum][0 + 2 * i] > 0) {
				//compare to see which is smaller mag
				//+ve
				if (fail_arr[k][rownum][0 + 2 * i] < failure[0]) {
					failure[0] = fail_arr[k][rownum][0 + 2 * i];
				}
				//-ve
				if (fail_arr[k][rownum][1 + 2 * i] > failure[1]) {
					failure[1] = fail_arr[k][rownum][1 + 2 * i];
				}
			}
			else { //i.e fail_arr[k][colnum][0 + 2 * i] < 0
				//compare to see which is smaller mag
				//+ve
				if (fail_arr[k][rownum][1 + 2 * i] < failure[0]) {
					failure[0] = fail_arr[k][rownum][1 + 2 * i];
				}
				//-ve
				if (fail_arr[k][rownum][0 + 2 * i] > failure[1]) {
					failure[1] = fail_arr[k][rownum][0 + 2 * i];
				}
			}
		}

	}
}

//obtain failure Nx, Ny, Nxy based on Tsai Wu Failure Criterion [TSAI WU FAILURE CRITERION]
void TsaiWu1(double stress[][3][3], double orientation[], double stress_thermal[][3], double x, double xx, double y, double yy, double s, int size) {
	//getting Tsai Wu coefft
	//xx and yy are auto converted to -ve values
	double F1, F2, F11, F22, F66, F12;
	TsaiWu_coefft(x, xx, y, yy, s, &F1, &F2, &F11, &F22, &F66, &F12);

	//obtaining the coefft for each {N} component
	double coefft[MAX / 2][3][3] = { 0 };
	/*
		{coefft} for k-th ply
			=
	^2		^1		^0		-	Nx
	^2		^1		^0		-	Ny
	^2		^1		^0		-	Nxy
	
	
	{stress} = [orientation][{N}][{coord}]
	Therefore, for k-th ply
		{Stress}
			=
	1		2		12		-	Nx
	1		2		12		-	Ny
	1		2		12		-	Nxy
	

	F1(stress1) + F2(stress2) + F11(stress1)^2 + F22(stress2)^2 + F66(stress12)^2 + 2F12(stress1 x stress2) - 1 < 0 --> for safety region
	*/

	int i, k;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			// NOTE: stress = stress(in terms of {N}) + thermal_stress Applicable to stress in 1-dir, 2-dir and 12-dir
			
			// ^2 coefft
			coefft[k][i][0] = F11 * pow(stress[k][i][0], 2) + F22 * pow(stress[k][i][1], 2) + F66 * pow(stress[k][i][2], 2) + 2 * F12 * stress[k][i][0] * stress[k][i][1];
			// ^1 coefft
			coefft[k][i][1] = F1 * stress[k][i][0] + F2 * stress[k][i][1] + F11 * (2 * stress[k][i][0] * stress_thermal[k][0]) + F22 * (2 * stress[k][i][1] * stress_thermal[k][1]) + F66 * (2 * stress[k][i][2] * stress_thermal[k][2]) + 2 * F12 * (stress[k][i][0] * stress_thermal[k][1] + stress[k][i][1] * stress_thermal[k][0]);
			// ^0 coefft
			coefft[k][i][2] = -1 + F1 * stress_thermal[k][0] + F2 * stress_thermal[k][1] + F11 * (pow(stress_thermal[k][0], 2)) + F22 * (pow(stress_thermal[k][1], 2)) + F66 * (pow(stress_thermal[k][2], 2)) + 2 * F12 * stress_thermal[k][0] * stress_thermal[k][1];

			//check to ensure quadratic eqn can be solved
			//if b^2-4ac < 0 --> prog end
			if ((pow(coefft[k][i][1], 2) - 4 * coefft[k][i][0] * coefft[k][i][2]) < 0) {
				printf("Error. Tsai Wu Quadratic equation results in imaginary solutions. Ending program!\n");
				return;
			}
		}
	}

	//get failure {N} then find governing {N} components
	double fail_arr[MAX / 2][3][2] = { 0 }; //crit {N} components for k-th ply
	/*
		{fail_arr} for k-th ply
			=
	+ve No.		-ve No.		-	Nx
	+ve No.		-ve No.		-	Ny
	+ve No.		-ve No.		-	Nxy
	*/
	k = 0;
	i = 0;
	double temp1, temp2;
	for (k = 0;k < size;k++) {
		for (i = 0;i < 3;i++) {
			temp1 = (-1.0 * coefft[k][i][1] + sqrt((pow(coefft[k][i][1], 2) - 4 * coefft[k][i][0] * coefft[k][i][2]))) / (2 * coefft[k][i][0]);
			temp2 = (-1.0 * coefft[k][i][1] - sqrt((pow(coefft[k][i][1], 2) - 4 * coefft[k][i][0] * coefft[k][i][2]))) / (2 * coefft[k][i][0]);
			if (temp1 > temp2) {
				fail_arr[k][i][0] = temp1;
				fail_arr[k][i][1] = temp2;
			}
			else {
				fail_arr[k][i][0] = temp2;
				fail_arr[k][i][1] = temp1;
			}
		}
	}

	double fail_results[3][2] = { 0 }; //crit {N} components for overall laminate
	/*
		{fail_results}
			=
	+ve No.		-ve No.		-	Nx
	+ve No.		-ve No.		-	Ny
	+ve No.		-ve No.		-	Nxy
	*/
	i = 0;
	k = 0;
	for (i = 0;i < 3;i++) { //doing one {N} component at a time
		//initialise first to allow comparison
		fail_results[i][0] = fail_arr[0][i][0];
		fail_results[i][1] = fail_arr[0][i][1];
		for (k = 1;k < size;k++) {
			//Comapre to see which has lower magnitude
			//+ve
			if (fabs(fail_results[i][0]) > fabs(fail_arr[k][i][0])) {
				fail_results[i][0] = fail_arr[k][i][0];
			}
			//-ve
			if (fabs(fail_results[i][1]) > fabs(fail_arr[k][i][1])) {
				fail_results[i][1] = fail_arr[k][i][1];
			}
		}

	}

	//output results
	printf("**Assuming uniaxially loaded symmetrical laminate**\n");
	printf("**Assumed F12 = -0.5sqrt(F11 x F22)**\n");
	printf("TSAI-WU FAILURE CRITERION - Stress Resultant estimated safety range (MPa mm):\n");
	printf("%.2f < Nx < %.2f\n", fail_results[0][1], fail_results[0][0]);
	printf("%.2f < Ny < %.2f\n", fail_results[1][1], fail_results[1][0]);
	printf("%.2f < Nxy < %.2f\n", fail_results[2][1], fail_results[2][0]);
	printf("\n");
}

// cal [3x3][3x1] matrix multipliation
double mat_mul(double A[][3], double B[], int col) {
	/*
		.		.									 B[0]
		.		.									 
	 x[col]  =  a[col][0] a[col][1] a[col][2]	x	 B[1]
		.		.
		.		.									 B[2]
	
	*/
	
	double x;
	x = A[col][0] * B[0] + A[col][1] * B[1] + A[col][2] * B[2];
	return x;
}

//inverse 3 x 3 matrix
void mat_inverse(double invr_mat[][3], double mat[][3], int *checknull) {
	double det;
	det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
	if (det == 0) {
		*checknull = -1;
		printf("\n");
		printf("No inverse! \n");
		printf("**Error in calculation**\n");
		printf("\n");
		return;
	}
	
	invr_mat[0][0] = (mat[1][1] * mat[2][2]) - (mat[1][2] * mat[2][1]);
	invr_mat[0][1] = (mat[0][2] * mat[2][1]) - (mat[0][1] * mat[2][2]);
	invr_mat[0][2] = (mat[0][1] * mat[1][2]) - (mat[0][2] * mat[1][1]);

	invr_mat[1][0] = (mat[1][2] * mat[2][0]) - (mat[1][0] * mat[2][2]);
	invr_mat[1][1] = (mat[0][0] * mat[2][2]) - (mat[0][2] * mat[2][0]);
	invr_mat[1][2] = (mat[1][0] * mat[0][2]) - (mat[0][0] * mat[1][2]);

	invr_mat[2][0] = (mat[1][0] * mat[2][1]) - (mat[2][0] * mat[1][1]);
	invr_mat[2][1] = (mat[2][0] * mat[0][1]) - (mat[0][0] * mat[2][1]);
	invr_mat[2][2] = (mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]);

	int i, j;
	for (i = 0;i < 3;i++) {
		for (j = 0;j < 3;j++) {
			invr_mat[i][j] *= (1.0 / det);
		}
	}
	*checknull = 0;
}

//get Tsai Wu Coefft terms
void TsaiWu_coefft(double x, double xx, double y, double yy, double s, double* F1, double* F2, double* F11, double* F22, double* F66, double* F12) {
	double a;
	a = fabs((1 / x) + (1 / xx));
	*F1 = a;

	a = fabs((1 / y) + (1 / yy));
	*F2 = a;
	
	a = (-1.0) / (x * xx);
	*F11 = a;
	
	a = (-1.0) / (y * yy);
	*F22 = a;
	
	a = 1 / pow(s, 2);
	*F66 = a;
	
	//Assumed:
	a = -0.5 * sqrt((*F11) * (*F22));
	*F12 = a;
}

//end prog funct
void end_prog() {
	//end program (prevent auto close of console window)
	//char ended;
	printf("Calculation completed! Press E and Enter to end and close program!\n");
	//scanf_s("%c", &ended,5); //IMPT, note the requirement to put char length (scanf_s req)
	char ended;
	ended = 'd';
	do {
		scanf_s(" %c", &ended, 1);
		if (ended == 'e') {
			ended = 'E';
		}
	} while (ended != 'E');
}