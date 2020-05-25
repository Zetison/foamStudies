#include <iostream>
#include <fstream>
#include "string"
#include "math.h"
using namespace std;
void printArray(std::ostream& out, double arr[][3], int rows, int cols, string tabs) 
{
	for (int i = 0; i < rows; i++)
	{
		out << tabs << "( ";
		for(int j = 0; j < cols; j++)
			out << arr[i][j] << " ";
		out << ")\n";
	}
}

void printHeader(std::ostream& out, string version, string format, string classType, string location, string object) 
{
	out << "/*--------------------------------*- C++ -*----------------------------------*\\";
	out << "\n| =========                 |                                                 |";
	out << "\n| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |";
	out << "\n|  \\\\    /   O peration     | Version:  v1912                                 |";
	out << "\n|   \\\\  /    A nd           | Website:  www.openfoam.com                      |";
	out << "\n|    \\\\/     M anipulation  |                                                 |\n";
	out << "\\*---------------------------------------------------------------------------*/";
	if (!version.empty())
	{
		out << "\nFoamFile\n\{\n";
		out << "\tversion     " << version << ";\n";
		out << "\tformat      " << format << ";\n";
		out << "\tclass       " << classType << ";\n";
		out << "\tlocation    \"" << location << "\";\n";
		out << "\tobject      " << object << ";\n}";
	}
	out << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
	out << "\n\n";
}

void printFooter(std::ostream& out) 
{
	out << "\n";
	out << "// ************************************************************************* //\n";
}

void insSubArray(double A[][3], int idx, double x, double y, double z){
	A[idx][0] = x;
	A[idx][1] = y;
	A[idx][2] = z;
}


int main(int argc, char **argv){
	double D         = stod(argv[1]);
	double nu        = stod(argv[2]);
	int M            = stoi(argv[3]);
	int nRe_0        = stoi(argv[4]);
	int np           = stoi(argv[5]);
	double omega_rot = stod(argv[6]);
	double u_inf     = stod(argv[7]);
	double delta_t   = stod(argv[8]);
	double end_time  = stod(argv[9]);

	double b = 64*D;
	double a = 2*b;
	double t = 4*D;
	double L = 2*M_PI*D;

	ofstream out1("./constant/parameters");
	out1.precision(15);	
	printHeader(out1, "", "ascii", "dictionary", "constant", "parameters");
	out1 << "nu " << nu << ";\n";
	out1 << "u_inf " << u_inf << ";\n";
	out1 << "end_time " << end_time << ";\n";
	out1 << "delta_t " << delta_t << ";\n";
	out1 << "omega_rot " << omega_rot << ";\n";
	out1 << "L " << L << ";\n";
	out1 << "D " << D << ";\n";
	printFooter(out1);
	out1.close();
	
	ofstream out2("./system/decomposeParDict");
	out2.precision(15);	
	printHeader(out2, "2.0", "ascii", "dictionary", "constant", "decomposeParDict");
	out2 << "numberOfSubdomains " << np << ";\n";
	out2 << "method simple;\n";
	out2 << "simpleCoeffs\n";
	out2 << "{\n";
	out2 << "\tn (" << np << " 1 1);\n";
	out2 << "}\n\n";
	printFooter(out2);
	out2.close();
		

	double Re = max(u_inf,omega_rot*D/2)*D/nu;  // Reynolds number
	int nRe = nRe_0*pow(2,(M>7)?(M-7):0);        // Roughly the number of elements
	                                             // withing D/sqrt(Re) outside the cylinder 
                                               // (should be at least 5)

	

	ofstream out("./system/blockMeshDict");
	out.precision(15);	
	printHeader(out, "2.0", "ascii", "dictionary", "system", "blockMeshDict");
	out << "vertices\n";
	out << "(\n";
	const int n = 14; 
	double points[2*n][3];
	int counter = 0;
	double theta;
	double r[2] = { D/2, t+D/2 };
	for (size_t i = 0; i < 2; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			theta = (j+1/2.0)*M_PI/2;
			insSubArray(points, counter++, r[i]*cos(theta), r[i]*sin(theta), 0);
		}
	}
	double x[3] = {-b/2, b/2, a-b/2 };
	double y[2] = {-b/2, b/2 };
	double R;
	for (size_t j = 0; j < 2; j++)
	{
		for (size_t i = 0; i < 3; i++)
		{
			R = sqrt(x[i]*x[i]+y[j]*y[j]);
			if (R > 1.1*r[1])
				insSubArray(points, counter++, x[i], y[j], 0);
		}
	}
	// Duplicate z points
	for (size_t i = 0; i < n; i++)
	{
		insSubArray(points, i+n, points[i][0], points[i][1], L);
	}
	
	printArray(out,points,2*n,3, "\t");
	out << ");\n";
	out << "\n";
	out << "blocks\n";
	out << "(\n";
	int N = 1 << M-1;
	const int nb = 9;
	double L_o = b/2-(D/2+t);
	double L_i = t;
	double h_o = b/N;
	double delta_s_i = D/sqrt(Re)/nRe;
	double beta = (L_i+L_o-delta_s_i+h_o)/(L_i+L_o);
	double N_o = log(h_o/(L_o*(1-beta)+h_o))/log(beta);
	double N_i = log((L_o*(1-beta)+h_o)/delta_s_i)/log(beta);
	double delta_e_i = delta_s_i*pow(beta,N_i-1);
	double delta_s_o = beta*delta_e_i;
	double delta_e_o = h_o/beta;
	double r_i = delta_e_i/delta_s_i;
	double r_o = delta_e_o/delta_s_o;
	int ref[5] = {(int)round(N_i),
	              N,
	              (int)round(N_o),
	              (int)round(N*(a-b)/b), 
	              1};
	int idx[n][4] = { {3, 7, 4, 0},      // 0
	                   {0, 4, 5, 1},     // 1
	                   {1, 5, 6, 2},     // 2
	                   {2, 6, 7 ,3},     // 3
	                   {7, 9, 12, 4},    // 4
	                   {4, 12, 11, 5},   // 5
	                   {5, 11, 8, 6},    // 6
	                   {6, 8, 9, 7},     // 7
	                   {9, 10, 13, 12}}; // 8
	for (size_t i = 0; i < nb; i++)
	{
		out << "\thex (";
		for (size_t k = 0; k < 2; k++)
		{
			for (size_t j = 0; j < 4; j++)
				out << idx[i][j]+n*k << " ";
		}
		out << ") (";
		if (i < 4)
			out << ref[0] << " " << ref[1] << " " << ref[4] << ") simpleGrading (" << r_i << " 1 1)\n"; 
		else if (i < 8)
			out << ref[2] << " " << ref[1] << " " << ref[4] << ") simpleGrading (" << r_o << " 1 1)\n"; 
		else
			out << ref[3] << " " << ref[1] << " " << ref[4] << ") simpleGrading (1 1 1)\n"; 
	}
	out << ");\n";
	out << "\n";
	out << "edges\n";
	out << "(\n";
	for (size_t k = 0; k < 2; k++)
	{
		for (size_t i = 0; i < 4; i++)
		{
			theta = (i+1)*M_PI/2;
			for (size_t j = 0; j < 2; j++)
			{
				out << "\tarc " << i+k*n+4*j << " ";
				if (i == 3)
					out << k*n+4*j;
				else 
					out << i+1+k*n+4*j;
				out << " (" << (D/2+j*t)*cos(theta) << " " << (D/2+j*t)*sin(theta) << " " << k*L << ")\n";
			}
		}
	}
	out << ");\n";
	out << "\n";
	out << "boundary\n";
	out << "(\n";
	out << "\tinlet\n";
	out << "\t{\n";
	out << "\t\ttype wall;\n";
	out << "\t\tfaces\n";
	out << "\t\t(\n";
	out << "\t\t\t( 11 8 22 25 )\n";
	out << "\t\t);\n";
	out << "\t}\n";
	
	out << "\toutlet\n";
	out << "\t{\n";
	out << "\t\ttype wall;\n";
	out << "\t\tfaces\n";
	out << "\t\t(\n";
	out << "\t\t\t( 10 13 27 24 )\n";
	out << "\t\t);\n";
	out << "\t}\n";
	
	out << "\twalls\n";
	out << "\t{\n";
	out << "\t\ttype wall;\n";
	out << "\t\tfaces\n";
	out << "\t\t(\n";
	out << "\t\t\t( 12 11 25 26 )\n";
	out << "\t\t\t( 13 12 26 27 )\n";
	out << "\t\t\t( 8 9 23 22 )\n";
	out << "\t\t\t( 9 10 24 23 )\n";
	out << "\t\t);\n";
	out << "\t}\n";
	
	out << "\tcylinder\n";
	out << "\t{\n";
	out << "\t\ttype wall;\n";
	out << "\t\tfaces\n";
	out << "\t\t(\n";
	out << "\t\t\t( 0 3 17 14 )\n";
	out << "\t\t\t( 1 0 14 15 )\n";
	out << "\t\t\t( 2 1 15 16 )\n";
	out << "\t\t\t( 3 2 16 17 )\n";
	out << "\t\t);\n";
	out << "\t}\n";
	
	out << "\tfrontAndBack\n";
	out << "\t{\n";
	out << "\t\ttype empty;\n";
	out << "\t\tfaces\n";
	out << "\t\t(\n";
	for (size_t k = 0; k < 2; k++)
	{
		for (size_t i = 0; i < nb; i++)
		{
			out << "\t\t\t( ";
			for (size_t j = 0; j < 4; j++)
			{
				if (k == 0)
				  out << idx[i][3-j] << " ";
				else
				  out << idx[i][j]+n*k << " ";
			}
			out << ")\n";
		}
	}
	out << "\t\t);\n";
	out << "\t}\n";
	out << ");\n";
	out << "\n";
	out << "mergePatchPairs\n";
	out << "(\n";
	out << ");\n";
	printFooter(out);
	out.close();
	return 0;
}
