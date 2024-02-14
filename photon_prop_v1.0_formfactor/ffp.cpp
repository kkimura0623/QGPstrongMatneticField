#include <iostream>
#include <fstream>
using namespace std;

int main(){
	int n=0;
	double pT[401],r[401],rate[401],n0r[401],n0i[401],n1r[401],n1i[401],n2r[401],n2i[401];
	ifstream fin0("epem_xy_300.dat");
	while(fin0>>pT[n]>>r[n]>>rate[n]>>n0r[n]>>n0i[n]>>n1r[n]>>n1i[n]>>n2r[n]>>n2i[n]){
		if(n==400) break;
		n++;
	};
	fin0.close();

	ofstream openfile0("200402_e_300_xy_ptrate.dat");
	for(int i=0; i<401; i++){
		openfile0 << pT[i] << "  " << rate[i] << endl;
	};
	openfile0.close();

	ofstream openfile1("200402_e_300_xy_ptn0r.dat");
	for(int i=0; i<401; i++){
		openfile1 << pT[i] << "  " << n0r[i] << endl;
	};
	openfile0.close();

	ofstream openfile2("200402_e_300_xy_ptn0i.dat");
	for(int i=0; i<401; i++){
		openfile2 << pT[i] << "  " << n0i[i] << endl;
	};
	openfile1.close();

	ofstream openfile3("200402_e_300_xy_ptn1r.dat");
	for(int i=0; i<401; i++){
		openfile3 << pT[i] << "  " << n1r[i] << endl;
	};
	openfile2.close();

	ofstream openfile4("200402_e_300_xy_ptn1i.dat");
	for(int i=0; i<401; i++){
		openfile4 << pT[i] << "  " << n1i[i] << endl;
	};
	openfile3.close();

	ofstream openfile5("200402_e_300_xy_ptn2r.dat");
	for(int i=0; i<401; i++){
		openfile5 << pT[i] << "  " << n2r[i] << endl;
	};
	openfile5.close();

	ofstream openfile6("200402_e_300_xy_ptn2i.dat");
	for(int i=0; i<401; i++){
		openfile6 << pT[i] << "  " << n2i[i] << endl;
	};
	openfile6.close();

	ofstream openfile7("200402_e_300_xy_rrate.dat");
	for(int i=0; i<401; i++){
		openfile7 << r[i] << "  " << rate[i] << endl;
	};
	openfile7.close();

	ofstream openfile8("200402_e_300_xy_rn0r.dat");
	for(int i=0; i<401; i++){
		openfile8 << r[i] << "  " << n0r[i] << endl;
	};
	openfile8.close();

	ofstream openfile9("200402_e_300_xy_rn0i.dat");
	for(int i=0; i<401; i++){
		openfile9 << r[i] << "  " << n0i[i] << endl;
	};
	openfile9.close();

	ofstream openfile10("200402_e_300_xy_rn1r.dat");
	for(int i=0; i<401; i++){
		openfile10 << r[i] << "  " << n1r[i] << endl;
	};
	openfile10.close();

	ofstream openfile11("200402_e_300_xy_rn1i.dat");
	for(int i=0; i<401; i++){
		openfile11 << r[i] << "  " << n1i[i] << endl;
	};
	openfile11.close();

	ofstream openfile12("200402_e_300_xy_rn2r.dat");
	for(int i=0; i<401; i++){
		openfile12 << r[i] << "  " << n2r[i] << endl;
	};
	openfile12.close();

	ofstream openfile13("200402_e_300_xy_rn2i.dat");
	for(int i=0; i<401; i++){
		openfile13 << r[i] << "  " << n2i[i] << endl;
	};
	openfile13.close();

	n = 0;
	ifstream fin1("epem_yz_300.dat");
	while(fin1>>pT[n]>>r[n]>>rate[n]>>n0r[n]>>n0i[n]>>n1r[n]>>n1i[n]>>n2r[n]>>n2i[n]){
		if(n==400) break;
		n++;
	};
	fin1.close();

	ofstream openfile14("200402_e_300_yz_ptrate.dat");
	for(int i=0; i<401; i++){
		openfile14 << pT[i] << "  " << rate[i] << endl;
		cout << pT[i] << endl;
	};
	openfile14.close();

	ofstream openfile15("200402_e_300_yz_ptn0r.dat");
	for(int i=0; i<401; i++){
		openfile15 << pT[i] << "  " << n0r[i] << endl;
	};
	openfile15.close();

	ofstream openfile16("200402_e_300_yz_ptn0i.dat");
	for(int i=0; i<401; i++){
		openfile16 << pT[i] << "  " << n0i[i] << endl;
	};
	openfile16.close();

	ofstream openfile17("200402_e_300_yz_ptn1r.dat");
	for(int i=0; i<401; i++){
		openfile17 << pT[i] << "  " << n1r[i] << endl;
	};
	openfile17.close();

	ofstream openfile18("200402_e_300_yz_ptn1i.dat");
	for(int i=0; i<401; i++){
		openfile18 << pT[i] << "  " << n1i[i] << endl;
	};
	openfile18.close();

	ofstream openfile19("200402_e_300_yz_ptn2r.dat");
	for(int i=0; i<401; i++){
		openfile19 << pT[i] << "  " << n2r[i] << endl;
	};
	openfile19.close();

	ofstream openfile20("200402_e_300_yz_ptn2i.dat");
	for(int i=0; i<401; i++){
		openfile20 << pT[i] << "  " << n2i[i] << endl;
	};
	openfile20.close();

	ofstream openfile21("200402_e_300_yz_rrate.dat");
	for(int i=0; i<401; i++){
		openfile21 << r[i] << "  " << rate[i] << endl;
	};
	openfile21.close();

	ofstream openfile22("200402_e_300_yz_rn0r.dat");
	for(int i=0; i<401; i++){
		openfile22 << r[i] << "  " << n0r[i] << endl;
	};
	openfile22.close();

	ofstream openfile23("200402_e_300_yz_rn0i.dat");
	for(int i=0; i<401; i++){
		openfile23 << r[i] << "  " << n0i[i] << endl;
	};
	openfile23.close();

	ofstream openfile24("200402_e_300_yz_rn1r.dat");
	for(int i=0; i<401; i++){
		openfile24 << r[i] << "  " << n1r[i] << endl;
	};
	openfile24.close();

	ofstream openfile25("200402_e_300_yz_rn1i.dat");
	for(int i=0; i<401; i++){
		openfile25 << r[i] << "  " << n1i[i] << endl;
	};
	openfile25.close();

	ofstream openfile26("200402_e_300_yz_rn2r.dat");
	for(int i=0; i<401; i++){
		openfile26 << r[i] << "  " << n2r[i] << endl;
	};
	openfile26.close();

	ofstream openfile27("200402_e_300_yz_rn2i.dat");
	for(int i=0; i<401; i++){
		openfile27 << r[i] << "  " << n2i[i] << endl;
	};
	openfile27.close();

};
