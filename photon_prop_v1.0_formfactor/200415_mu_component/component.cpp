#include <iostream> 
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
using namespace std;

int main(){
// 25GeV/c < p_gamma < 40GeV/c
// 100MeV/c 刻み
	int N,M;
	long double C[1604][8];
	double S[1604][4];
	double P[1604][4];
	string fname[10];
	string c[4];

	fname[0] = "lep_xy_300";
	fname[1] = "lep_yz_300";
	fname[2] = "src_xy_300";
	fname[3] = "src_yz_300";
	fname[4] = "prop_xy_300";
	fname[5] = "prop_yz_300";
	fname[6] = "polten_xy_300";
	fname[7] = "polten_yz_300";
	fname[8] = "TT_xy_300";
	fname[9] = "TT_yz_300";

	c[0] = "0";
	c[1] = "1";
	c[2] = "2";
	c[3] = "3";

	for(int n=0; n<4;n++){
	ifstream fin(fname[n]+".dat");
		for(int row = 0; row<604; ++row){
				for(int col=0; col<4; ++col){
					fin >> C[row][col];
				};
		};
		for(int j=0; j<4; j++){
			ofstream fout1("200415_mu_"+fname[n]+"_0"+c[j]+".dat");
			for(int i = 0; i<604 ; i=i+4){
				fout1 << 25*pow(10,3.0)+i*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout2("200415_mu_"+fname[n]+"_1"+c[j]+".dat");
			for(int i = 1; i<604 ; i=i+4){
				fout2 << 25*pow(10,3.0)+(i-1)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout3("200415_mu_"+fname[n]+"_2"+c[j]+".dat");
			for(int i = 2; i<604 ; i=i+4){
				fout3 << 25*pow(10,3.0)+(i-2)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
			ofstream fout4("200415_mu_"+fname[n]+"_3"+c[j]+".dat");
			for(int i = 3; i<604 ; i=i+4){
				fout4 << 25*pow(10,3.0)+(i-3)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
			};
		};
	}
	for(int n=4; n<10;n++){
    ifstream fin(fname[n]+".dat");
        for(int row = 0; row<604; ++row){
                for(int col=0; col<8; ++col){
                    fin >> C[row][col];
                };  
        };  
        for(int j=0; j<4; j++){
            ofstream fout1("200415_mu_Re"+fname[n]+"_0"+c[j]+".dat");
            for(int i = 0; i<604 ; i=i+4){
                fout1 << 25*pow(10,3.0)+i*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
            ofstream fout2("200415_mu_Re"+fname[n]+"_1"+c[j]+".dat");
            for(int i = 1; i<604 ; i=i+4){
                fout2 << 25*pow(10,3.0)+(i-1)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
            ofstream fout3("200415_mu_Re"+fname[n]+"_2"+c[j]+".dat");
            for(int i = 2; i<604 ; i=i+4){
                fout3 << 25*pow(10,3.0)+(i-2)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
            ofstream fout4("200415_mu_Re"+fname[n]+"_3"+c[j]+".dat");
            for(int i = 3; i<604 ; i=i+4){
                fout4 << 25*pow(10,3.0)+(i-3)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
        };  
        for(int j=4; j<8; j++){
            ofstream fout5("200415_mu_Im"+fname[n]+"_0"+c[j-4]+".dat");
            for(int i = 0; i<604 ; i=i+4){
                fout5 << 25*pow(10,3.0)+i*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
            ofstream fout6("200415_mu_Im"+fname[n]+"_1"+c[j-4]+".dat");
            for(int i = 1; i<604 ; i=i+4){
                fout6 << 25*pow(10,3.0)+(i-1)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
            ofstream fout7("200415_mu_Im"+fname[n]+"_2"+c[j-4]+".dat");
            for(int i = 2; i<604 ; i=i+4){
                fout7 << 25*pow(10,3.0)+(i-2)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
            ofstream fout8("200415_mu_Im"+fname[n]+"_3"+c[j-4]+".dat");
            for(int i = 3; i<604 ; i=i+4){
                fout8 << 25*pow(10,3.0)+(i-3)*25 << "  " << setprecision(16) <<  C[i][j] << endl;
            };  
        };  
    };  
};

