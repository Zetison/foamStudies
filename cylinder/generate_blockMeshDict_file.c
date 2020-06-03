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
	double D   = SED_DIAM;
	int M      = SED_MESH;
	int N_Re_0 = SED_NRE_0;
	double bdD = SED_BDD;
	double Re  = SED_RE;
	double L   = SED_LENGTH;
	string str(argv[0]);
	bool createFile = str.compare("openfoam") == 0;

	double b = bdD*D;
	double a = 2*b;
	double t = 4*D;
	int N_Re = N_Re_0*pow(2,(M>7)?(M-7):0);        // Roughly the number of elements
	                                             // withing D/sqrt(Re) outside the cylinder 
                                               // (should be at least 5)

	int N = 1 << M-1;
	double delta_Nr = (a-b)/N;
	double f, df; // Find beta through Newton iterations
	double beta = 2.0;
	for (size_t i = 0; i < 100; i++)
	{
		f = (delta_Nr+(1-beta)*(b-D)/2)*(1-pow(beta,N_Re))-D/sqrt(Re)*(1-beta);
		df = (-(b-D)/2)*(1-pow(beta,N_Re)) - (delta_Nr+(1-beta)*(b-D)/2)*N_Re*pow(beta,N_Re-1) + D/sqrt(Re);
		beta -= f/df;
	}
	double delta_0 = delta_Nr + 0.5*(1-beta)*(b-D);
	int N_r = ceil(log(delta_Nr/delta_0)/log(beta));
	for (size_t i = 0; i < 100; i++)
	{
		f = D/sqrt(Re)*(1-pow(beta,N_r)) - (b-D)/2*(1-pow(beta,N_Re));
		df = -D/sqrt(Re)*N_r*pow(beta,N_r-1) + (b-D)/2*N_Re*pow(beta,N_Re-1);
		beta -= f/df;
	}
	if (!createFile)
	{
		cout << beta;
		return 0;
	}
	delta_0 = 0.5*(b-D)*(1-beta)/(1-pow(beta,N_r));
	double x_i = D/2;
	int N_t = 0;
	double delta = delta_0;
	while (x_i < D/2 + t)
	{
		x_i += delta;
		delta *= beta;
		N_t++;
	}
	t = x_i-D/2;
	if (createFile)
	{
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
			insSubArray(points, i+n, points[i][0], points[i][1], L);
		
		printArray(out,points,2*n,3, "\t");
		out << ");\n";
		out << "\n";
		out << "blocks\n";
		out << "(\n";
		const int nb = 9;
		int ref[5] = {N_t,
		              N,
		              N_r-N_t,
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
				for (size_t j = 0; j < 4; j++)
					out << idx[i][j]+n*k << " ";
			out << ") (";
			if (i < 4)
				out << ref[0] << " " << ref[1] << " " << ref[4] << ") simpleGrading (" << pow(beta,N_t-1) << " 1 1)\n"; 
			else if (i < 8)
				out << ref[2] << " " << ref[1] << " " << ref[4] << ") simpleGrading (" << pow(beta,N_r-N_t-1) << " 1 1)\n"; 
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
	}
	return 0;
}
