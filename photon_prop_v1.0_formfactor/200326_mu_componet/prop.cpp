#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

int main(){
	int N,M;
	long double C[1604][3];
	double S[1604][3];
	double P[1604][3];
	string fname[8];
	string c[4];

	fname[0] = "prop_xy_300";
	fname[1] = "prop_yz_300";

	c[0] = "0";
	c[1] = "1";
	c[2] = "2";
	c[3] = "3";

	for(int n=0; n<2;n++){
	ifstream fin(fname[n]+".dat");
		for(int row = 0; row<1604; ++row){
				for(int col=0; col<8; ++col){
					fin >> C[row][col];
				};
		};
		for(int j=0; j<4; j++){
			ofstream fout1("200326_mu_real"+fname[n]+"_0"+c[j]+".dat");
			for(int i = 0; i<1604 ; i=i+4){
				fout1 << i*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout2("200326_mu_real"+fname[n]+"_1"+c[j]+".dat");
			for(int i = 1; i<1604 ; i=i+4){
				fout2 << (i-1)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout3("200326_mu_real"+fname[n]+"_2"+c[j]+".dat");
			for(int i = 2; i<1604 ; i=i+4){
				fout3 << (i-2)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout4("200326_mu_real"+fname[n]+"_3"+c[j]+".dat");
			for(int i = 3; i<1604 ; i=i+4){
				fout4 << (i-3)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
		};
		for(int j=4; j<8; j++){
			ofstream fout5("200326_mu_imag"+fname[n]+"_0"+c[j-4]+".dat");
			for(int i = 0; i<1604 ; i=i+4){
				fout5 << i*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout6("200326_mu_imag"+fname[n]+"_1"+c[j-4]+".dat");
			for(int i = 1; i<1604 ; i=i+4){
				fout6 << (i-1)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout7("200326_mu_imag"+fname[n]+"_2"+c[j-4]+".dat");
			for(int i = 2; i<1604 ; i=i+4){
				fout7 << (i-2)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout8("200326_mu_imag"+fname[n]+"_3"+c[j-4]+".dat");
			for(int i = 3; i<1604 ; i=i+4){
				fout8 << (i-3)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
		};
	};
};

