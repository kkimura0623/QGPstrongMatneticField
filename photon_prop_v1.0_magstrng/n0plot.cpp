#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(){
    ifstream fin("log");

    string line;

    while(getline(fin,line)){
        cout << "読み込んだ行: " << line << endl;
    }

    fin.close();

    return 0;

}