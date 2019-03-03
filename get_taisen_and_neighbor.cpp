#include<bits/stdc++.h>
#include<unordered_map>
using namespace std;
#define fi first
#define se second
#define pb push_back
#define mp make_pair
#define rep(i,j,k) for (int i=(int)(j);i<=(int)(k);i++)
typedef double db;
const db eps=1e-8;
const db pi=acos(-1);
const double EARTH_RADIUS = 6378.137;
inline double rad(double d){
   return d * pi / 180.0;
}

inline int sign(db k) {
    if (k>eps) return 1;
    else if (k<-eps) return -1;
    return 0;
}
inline int cmp(db k1,db k2) {
    return sign(k1-k2);
}
inline int inmid(db k1,db k2,db k3) {
    return sign(k1-k3)*sign(k2-k3)<=0;   // k3 在 [k1,k2] 内
}
struct point {
    db x,y;
    point operator + (const point &k1) const {
        return (point) {
            k1.x+x,k1.y+y
        };
    }
    point operator - (const point &k1) const {
        return (point) {
            x-k1.x,y-k1.y
        };
    }
    point operator * (db k1) const {
        return (point) {
            x*k1,y*k1
        };
    }
    point operator / (db k1) const {
        return (point) {
            x/k1,y/k1
        };
    }
    int operator == (const point &k1) const {
        return cmp(x,k1.x)==0&&cmp(y,k1.y)==0;
    }
    // 逆时针旋转
    point turn(db k1) {
        return (point) {
            x*cos(k1)-y*sin(k1),x*sin(k1)+y*cos(k1)
        };
    }
    point turn90() {
        return (point) {
            -y,x
        };
    }
    bool operator < (const point& k1) const {
        int a=cmp(x,k1.x);
        if (a==-1) return 1;
        else if (a==1) return 0;
        else return cmp(y,k1.y)==-1;
    }
    db abs() {
        return sqrt(x*x+y*y);
    }
    db abs2() {
        return x*x+y*y;
    }
    db dis(const point& k1)const {
       double radLat1 = rad(y);
       double radLat2 = rad(k1.y);
       double a = radLat1 - radLat2;
       double b = rad(x) - rad(k1.x);
       double s = 2 * asin(sqrt(pow(sin(a/2),2) +
        cos(radLat1)*cos(radLat2)*pow(sin(b/2),2)));
       s*= EARTH_RADIUS;
       return s;
    }
    point unit() {
        db w=abs();
        return (point) {
            x/w,y/w
        };
    }
    void scan() {
        scanf("%lf%lf",&x,&y);
    }
    void print() {
        printf("%.11f %.11f\n",x,y);
    }
    db getw() {
        return atan2(y,x);
    }
    point getdel() {
        if (sign(x)==-1||(sign(x)==0&&sign(y)==-1)) return (*this)*(-1);
        else return (*this);
    }
    int getP() const {
        return sign(y)==1||(sign(y)==0&&sign(x)==-1);
    }
};
int inmid(const point& k1,const point& k2,const point& k3) {
    return inmid(k1.x,k2.x,k3.x)&&inmid(k1.y,k2.y,k3.y);   //点在矩形内
}
db cross(const point& k1,const point& k2) {
    return k1.x*k2.y-k1.y*k2.x;
}
db dot(const point& k1,const point& k2) {
    return k1.x*k2.x+k1.y*k2.y;
}
db rad(const point& k1,const point& k2) {
    return atan2(cross(k1,k2),dot(k1,k2));
}
// -pi -> pi
int compareangle (const point& k1,const point& k2) { //k1<k2返回true
    return k1.getP()<k2.getP()||(k1.getP()==k2.getP()&&sign(cross(k1,k2))>0);
}
point proj(const point& k1,const point& k2,const point& q) { // q 到直线 k1,k2 的投影
    point k=k2-k1;
    return k1+k*(dot(q-k1,k)/k.abs2());
}
point reflect(const point& k1,const point& k2,const point& q) {
    return proj(k1,k2,q)*2-q;
}
int clockwise(const point& k1,const point& k2,const point& k3) { // k1 k2 k3 逆时针 1 顺时针 -1 否则 0
    return sign(cross(k2-k1,k3-k1));
}
int checkLL(const point& k1,const point& k2,const point& k3,const point& k4) { // 求直线 (L) 线段 (S)k1,k2 和 k3,k4 的交点
    return cmp(cross(k3-k1,k4-k1),cross(k3-k2,k4-k2))!=0;
}
point getLL(const point&  k1,const point&  k2,const point&  k3,const point&  k4) {
    db w1=cross(k1-k3,k4-k3),w2=cross(k4-k3,k2-k3);
    return (k1*w2+k2*w1)/(w1+w2);
}
int intersect(db l1,db r1,db l2,db r2) { //区间交集
    if (l1>r1) swap(l1,r1);
    if (l2>r2) swap(l2,r2);
    return cmp(r1,l2)!=-1&&cmp(r2,l1)!=-1;
}
int checkSS(const point&  k1,const point&  k2,const point&  k3,const point&  k4) { //线段交点
    return intersect(k1.x,k2.x,k3.x,k4.x)&&intersect(k1.y,k2.y,k3.y,k4.y)&&
           sign(cross(k3-k1,k4-k1))*sign(cross(k3-k2,k4-k2))<=0&&
           sign(cross(k1-k3,k2-k3))*sign(cross(k1-k4,k2-k4))<=0;
}
db disSP(const point&  k1,const point&  k2,const point&  q) {
    point k3=proj(k1,k2,q);
    if (inmid(k1,k2,k3)) return q.dis(k3);
    else return min(q.dis(k1),q.dis(k2));
}
db disSS(const point&  k1,const point&  k2,const point&  k3,const point&  k4) {
    if (checkSS(k1,k2,k3,k4)) return 0;
    else return min(min(disSP(k1,k2,k3),disSP(k1,k2,k4)),min(disSP(k3,k4,k1),disSP(k3,k4,k2)));
}
int onS(const point&  k1,const point&  k2,const point&  q) {
    return inmid(k1,k2,q)&&sign(cross(k1-q,k2-k1))==0;   //点在线段上
}
struct circle {
    point o;
    db r;
    void scan() {
        o.scan();
        scanf("%lf",&r);
    }
    int inside(const point&  k)const {
        return cmp(r,o.dis(k));
    }
};
struct line {
    // p[0]->p[1]
    point p[2];
    line(const point&  k1,const point&  k2) {
        p[0]=k1;
        p[1]=k2;
    }
    point& operator [] (int k) {
        return p[k];
    }
    int include(const point&  k) {
        return sign(cross(p[1]-p[0],k-p[0]))>0;   //k在半平面内
    }
    point dir() {
        return p[1]-p[0];
    }
    line push() { // 向外 ( 左手边 ) 平移 eps
        const db eps = 1e-6;
        point delta=(p[1]-p[0]).turn90().unit()*eps;
        return {p[0]-delta,p[1]-delta};
    }
};
point getLL(line& k1,line& k2) {
    return getLL(k1[0],k1[1],k2[0],k2[1]);
}
int checkLL(line& k1,line& k2) {
    return checkLL(k1[0],k1[1],k2[0],k2[1]);
}
int parallel(line& k1,line& k2) {
    return sign(cross(k1.dir(),k2.dir()))==0;
}
int sameDir(line &k1,line &k2) {
    return parallel(k1,k2)&&sign(dot(k1.dir(),k2.dir()))==1;   //向量同方向
}
int operator < (line& k1,line& k2) { //按幅角比较
    if (sameDir(k1,k2)) return k2.include(k1[0]);
    return compareangle(k1.dir(),k2.dir());
}
int checkpos(line& k1,line& k2,line& k3) {
    return k3.include(getLL(k1,k2));
}
void getHL(vector<line> &L) { // 求半平面交 , 半平面是逆时针方向 , 输出（线）按照逆时针
    sort(L.begin(),L.end());
    deque<line> q;
    for (int i=0; i<(int)L.size(); i++) {
        if (i&&sameDir(L[i],L[i-1])) continue;
        while (q.size()>1&&!checkpos(q[q.size()-2],q[q.size()-1],L[i])) q.pop_back();
        while (q.size()>1&&!checkpos(q[1],q[0],L[i])) q.pop_front();
        q.push_back(L[i]);
    }
    while (q.size()>2&&!checkpos(q[q.size()-2],q[q.size()-1],q[0])) q.pop_back();
    while (q.size()>2&&!checkpos(q[1],q[0],q[q.size()-1])) q.pop_front();
    L.clear();
    for (int i=0; i<q.size(); i++) L.push_back(q[i]);
}
inline double rd(){
    return rand()*1.0/RAND_MAX;
}
struct point_hash{
    size_t operator()(const point &p)const{
        return (int)(p.x*1000000)^(int)(p.y*1000000);
    }
};
const int maxn=1000;
point po[maxn],pt[maxn];
vector<line> L;
unordered_map<point,int,point_hash>pt_mp;
int mode[maxn];
int main() {
    ofstream fout_nei("data\\nei.dat");
	ofstream fout_ts("data\\ts.dat");
	ofstream fout_loc("data\\loc.dat");
	srand(time(NULL));
	int n=0;
    ifstream fin_loc("data\\韩国城市信息（附人口）.csv");
    string lineStr;
    getline(fin_loc, lineStr);
    double x_max=-1e9,y_max=-1e9;
    double x_min=1e9,y_min=1e9;
    while (getline(fin_loc, lineStr)){
        istringstream ss(lineStr);
        //ss.ignore(100, '"');
        ss>>po[n].x;
        x_max=max(x_max,po[n].x);
        x_min=min(x_min,po[n].x);
        ss.get();
        ss>>po[n].y;
        y_max=max(y_max,po[n].y);
        y_min=min(y_min,po[n].y);
        ss.get();
        ss>>mode[n];
        n++;
    }
    fin_loc.close();
    fout_nei<<n<<endl;
	fout_ts<<n<<endl;
	fout_loc<<n<<endl;
    double x_gap=(x_max-x_min)/(n-1),y_gap=(y_max-y_min)/(n-1);
    x_max+=x_gap,x_min-=x_gap,y_max+=y_gap,y_min-=y_gap;
    for(int i=0;i<n;i++){
        pt[i].x=(po[i].x-x_min)/(x_max-x_min);
        pt[i].y=(po[i].y-y_min)/(y_max-y_min);
        fout_loc<<mode[i]<<" "<<pt[i].x<<" "<<pt[i].y<<endl;
        pt_mp[pt[i]]=i;
    }
    //double sum[]={0,0,0,0};int cnt[]={0,0,0,0};
    double sum=0;
    int cnt=0;
	for(int i=0;i<n;i++){
        L.clear();
        L.pb(line({0,0},{1,0}));
        L.pb(line({1,0},{1,1}));
        L.pb(line({1,1},{0,1}));
        L.pb(line({0,1},{0,0}));
		for(int j=0;j<n;j++){
			if(i!=j){
                point st=(pt[i]+pt[j])*0.5;
                point dir=(pt[j]-pt[i]).turn90();
                L.pb(line(st,st+dir));
			}
		}
        getHL(L);
        vector<int>nei_id;
        for(int j=0;j<L.size();j++){
            point rp=reflect(L[j].p[0],L[j].p[1],pt[i]);
            if(pt_mp.find(rp)!=pt_mp.end()){
                nei_id.pb(pt_mp[rp]);
            }
        }
        fout_nei<<i<<" "<<nei_id.size()<<endl;
        for(auto id : nei_id){
            fout_nei<<id<<" ";
            //if(mode[i]==mode[id]){
//                sum[mode[i]]+=po[i].dis(po[id]);
//                cnt[mode[i]]++;
            //}
            sum+=po[i].dis(po[id]);
            cnt++;
        }
        fout_nei<<endl;
        fout_ts<<i<<" "<<L.size()<<endl;
        L.pb(L[0]);
        for(int j=0;j<L.size()-1;j++){
            point xp=getLL(L[j],L[j+1]);
            fout_ts<<xp.x<<" "<<xp.y<<endl;
        }
        //cout<<nei_id.size()<<" "<<L.size()-1<<endl;
	}
//	cout<<"china east average distance:"<<sum[1]/cnt[1]<<"km"<<endl;
//	cout<<"china middle average distance:"<<sum[2]/cnt[2]<<"km"<<endl;
//	cout<<"china west average distance:"<<sum[3]/cnt[3]<<"km"<<endl;
    cout<<"korea average distance:"<<sum/cnt<<"km"<<endl;
	fout_loc.close();
	fout_nei.close();
	fout_ts.close();
    return 0;
}
