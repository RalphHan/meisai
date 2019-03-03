#include<bits/stdc++.h>
using namespace std;

int main(){
    int iter,num_of_block;
    ifstream fin("data\\process_korea123.dat");
    fin>>num_of_block>>iter;
    ofstream fout("data\\process_korea123.csv");
    fout<<"time";
    for(int i=0;i<num_of_block;i++)fout<<","<<i;
    fout<<endl;
    double Q;
    for(int t=1;t<=iter;t++){
        fout<<t;
        for(int i=0;i<num_of_block;i++){
        	fin>>Q;
            fout<<","<<Q;
        }
        fout<<endl;
    }
    fin.close();
    fout.close();
    return 0;
}
