#include<bits/stdc++.h>
using namespace std;
const int maxn=1000;
int pre,now,num_of_block;
double tot_population;
const double pi=acos(-1);
const double EARTH_RADIUS = 6378.137;
inline double rad(double d){
   return d * pi / 180.0;
}
struct Cell_Base{
    int mode;
    double x,y;
    string name;
    double DistTo(const Cell_Base&k1){
       double radLat1 = rad(y);
       double radLat2 = rad(k1.y);
       double a = radLat1 - radLat2;
       double b = rad(x) - rad(k1.x);
       double s = 2 * asin(sqrt(pow(sin(a/2),2) +
        cos(radLat1)*cos(radLat2)*pow(sin(b/2),2)));
       s*= EARTH_RADIUS;
       return s;
    }
    vector<double>D;
    vector<int>nei;
    double k,Dmin;
    double population;
}CB[maxn];
typedef double Cell;
Cell C[2][maxn];
inline double rd(){
    return rand()*1.0/RAND_MAX;
}
double mode[]={0,0.7,0.35,0.2};
double mode_init[]={0,0.1,0.05,0.01};
vector<int>choice;
void init(){
    //load_data
    ifstream fin("data\\loc.dat");
    fin>>num_of_block;
    int x;

    for(int i=0;i<num_of_block;i++){
        fin>>x;
//        if(x==1)
        choice.push_back(i);
        CB[i].k=mode[x];
        fin>>CB[i].x>>CB[i].y;
        C[0][i]=1e-3;
    }
    fin.close();
    fin.open("data\\korea_pop.dat");
    tot_population=0.0;
    for(int i=0;i<num_of_block;i++){
        fin>>CB[i].mode>>CB[i].name>>CB[i].population;
        tot_population+=CB[i].population;
    }
    random_shuffle(choice.begin(),choice.end());
    int up=min(num_of_block/20,(int)choice.size());
    for(int i=0;i<up;i++)C[0][choice[i]]=mode_init[CB[choice[i]].mode];
    fin.close();
    fin.open("data\\nei.dat");
    fin>>x;
    int num_of_nei;
    for(int i=0;i<num_of_block;i++){
        fin>>x>>num_of_nei;
        CB[i].Dmin=1e9;
        for(int j=0;j<num_of_nei;j++){
            fin>>x;
            CB[i].nei.push_back(x);
            double dist=CB[i].DistTo(CB[x]);
            CB[i].D.push_back(dist);
            CB[i].Dmin=min(CB[i].Dmin,dist);
        }
    }
    fin.close();
    pre=0,now=1;
}
void update(int id,Cell&c_new){
    Cell_Base&cb=CB[id];
    Cell&c_pre=C[pre][id];
    double n=0;
    double lambda=1.0;
    for(int i=0;i<cb.nei.size();i++){
        double Qnei=C[pre][cb.nei[i]];
        if(Qnei>max(c_pre,0.5))
            lambda*=(1+(Qnei-c_pre)/cb.D[i]*cb.Dmin);
    }
    c_new=(1-c_pre)*c_pre*cb.k*lambda+c_pre;
}

int main(){
    init();
    srand(time(NULL));
    int iter=40;
    ofstream fout("data\\process_korea123.dat");
    ofstream fout_tot("data\\process_korea123_tot.csv");
    fout<<num_of_block<<" "<<iter<<endl;
    fout_tot<<"choice";
    int up=min(num_of_block/20,(int)choice.size());
    for(int i=0;i<up;i++)fout_tot<<","<<CB[choice[i]].name;
    fout_tot<<endl<<endl;

    fout_tot<<"iter,process"<<endl;
    for(int t=1;t<=iter;t++){
        double Qtot=0;
        for(int i=0;i<num_of_block;i++){
            update(i,C[now][i]);
            fout<<C[now][i]<<" ";
            Qtot+=C[now][i]*CB[i].population;
        }
        fout<<endl;
        fout_tot<<t<<","<<Qtot/tot_population<<endl;
        swap(now,pre);
    }
    fout.close();
    fout_tot.close();
    return 0;
}
