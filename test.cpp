#include<bits/stdc++.h>
using namespace std;
const int maxn=107;
const double pi=acos(-1);
struct point{
    double x,y;
    point(double x=0.0,double y=0.0):x(x),y(y){}
    point operator -(const point&p)const{
        return point(x-p.x,y-p.y);
    }
    point operator *(const double& k){
        return point(x*k,y*k);
    }
    double l2(){
        return sqrt(x*x+y*y);
    }
}scatter[maxn];
double dis(const point &a,const point &b){
    return (a-b).l2();
}
inline double rd(){
    return 1.0*rand()/RAND_MAX;
}
int main(){
    freopen("test.txt","r",stdin);
    int T;
    cin>>T;
    srand(time(NULL));
    while(T--){
        int n;
        cin>>n;
        for(int i=0;i<n;i++){
            scanf("%lf%lf",&scatter[i].x,&scatter[i].y);
        }
        point O1=point(rd(),rd()),O2=point(rd(),rd());
        double step=1;
        int iter=0;
        while(step>1e-6){
            double maxx=-1;
            int loc=-1;
            for(int i=0;i<n;i++){
                double tmp=dis(O1,scatter[i])+dis(O2,scatter[i]);
                if(tmp>maxx){
                    maxx=tmp;
                    loc=i;
                }
            }
            double a=maxx/2;
            double c=dis(O1,O2)/2;
            double b=sqrt(a*a-c*c);
            double s=pi*a*b;
            iter++;
            if(iter%20==0){
                cout<<s<<endl;
            }
            double k1=b+a*a/b;
            double k2=-a/b/2;
            double O1A=dis(O1,scatter[loc]);
            double O2A=dis(O2,scatter[loc]);
            double u=scatter[loc].x,v=scatter[loc].y;
            double x1=k1*(O1.x-u)/O1A+k2*(O1.x-O2.x);
            double y1=k1*(O1.y-v)/O1A+k2*(O1.y-O2.y);
            double x2=k1*(O2.x-u)/O2A+k2*(O2.x-O1.x);
            double y2=k1*(O2.y-v)/O2A+k2*(O2.y-O1.y);
            double grad=sqrt(x1*x1+y1*y1+x2*x2+y2*y2);
            x1/=grad,y1/=grad,x2/=grad,y2/=grad;
            O1=O1-point(x1,y1)*step;
            O2=O2-point(x2,y2)*step;
            step*=0.98;
        }
    }
    return 0;
}
