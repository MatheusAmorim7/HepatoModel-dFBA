#include "core/PhysiCell.h"
#include "modules/PhysiCell_settings.h"
#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

extern "C" void dpbtrf_(const char*,const int*,const int*,double*,const int*,int*);
extern "C" void dpbtrs_(const char*,const int*,const int*,const int*,const double*,const int*,double*,const int*,int*);
using V=std::vector<double>;
using R=std::array<double,4>;
void need(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
double total(const V& a){return std::accumulate(a.begin(),a.end(),0.);}
double maxabs(const V& a){double m=0;for(double x:a)m=std::max(m,std::abs(x));return m;}
bool valid(const V& a){for(double x:a)if(!std::isfinite(x)||x<0)return false;return true;}
double content(double p){need(std::isfinite(p)&&p>=0&&p<=150,"inlet_domain");double g=p*p*p+150*p;return .0013998848684210528*p+8*g/(g+23400);}
double saturation(double p){double g=p*p*p+150*p;return g/(g+23400);}
struct State {V b,e,p;R r;};
struct Diagnostics {double in=0,out=0,wall=0,water=0,root=0,body=0;};
struct Edge {int i,j;double k;};
struct Operator {
    V vb,ve,a,d,den,linear,hb;std::vector<int> node;std::vector<Edge> edges;
    V factor;double dt,m,qflow,wallG,D;int n=2500,kd=50,ld=51;
    Operator(const std::string& path,double seconds,double multiplicity,double diffusion=1000,double wall=.01,double flow=3.492):dt(seconds),m(multiplicity),qflow(flow),wallG(wall),D(diffusion){
        need(std::isfinite(dt)&&dt>0&&m>0&&std::isfinite(D)&&D>=0,"duration");
        std::ifstream f(path);need(bool(f),"input_open");int nb,ne;f>>nb>>ne;need(nb==48&&ne==2500,"input_shape");
        vb.resize(48);ve.resize(n);node.resize(n);V dummy(48);double x,y,z,qe,qb,p;
        for(int j=0;j<48;++j){f>>vb[j]>>qb>>p;need(vb[j]>0&&qflow*m*dt/60/vb[j]<=1,"CFL");}
        for(int i=0;i<n;++i){f>>ve[i]>>node[i]>>qe>>x>>y>>z;need(ve[i]>0&&node[i]>=-1&&node[i]<48,"volume_map");
            need(std::abs(x-(-490+20*(i%50)))<1e-10&&std::abs(y-(-490+20*(i/50)))<1e-10&&z==0,"grid_coordinates");}
        need(bool(f),"input_read");
        a.resize(n);d.assign(n,0);den.resize(n);linear.resize(48);hb.resize(48);V capacity(48,0);
        for(int i=0;i<n;++i)if(node[i]>=0)capacity[node[i]]+=ve[i];
        for(int j=0;j<48;++j){need(capacity[j]>0,"wall_support");linear[j]=vb[j]*.0013998848684210528;hb[j]=vb[j]*8;}
        for(int i=0;i<n;++i){a[i]=ve[i]*.002;
            if(node[i]>=0)d[i]=(dt/60)*m*wallG*ve[i]/capacity[node[i]];
            den[i]=a[i]+d[i];if(node[i]>=0)linear[node[i]]+=a[i]*d[i]/den[i];}
        for(int i=0;i<n;++i)for(int delta:{1,50})if(i+delta<n&&(delta==50||i%50!=49)){
            int j=i+delta;edges.push_back({i,j,D*(2*ve[i]*ve[j]/(ve[i]+ve[j]))/400});}
        need(edges.size()==4900,"edge_count");factor.assign(ld*n,0.);
        for(int i=0;i<n;++i)factor[ld*i]=ve[i];
        for(auto e:edges){double t=dt*e.k;factor[ld*e.i]+=t;factor[ld*e.j]+=t;factor[e.j-e.i+ld*e.i]-=t;}
        int info;char low='L';dpbtrf_(&low,&n,&kd,factor.data(),&ld,&info);need(info==0,"water_factorization");
    }
    State initial(const std::string& path){
        std::ifstream f(path);int nb,ne;f>>nb>>ne;State s;s.b.resize(48);s.p.resize(48);s.e.resize(n);double v,x,y,z;int nd;
        for(int i=0;i<48;++i)f>>v>>s.b[i]>>s.p[i];
        for(int i=0;i<n;++i)f>>v>>nd>>s.e[i]>>x>>y>>z;
        for(double& x:s.r)f>>x;need(bool(f),"initial_read");return s;
    }
    State propose(const State& old,bool opened,double inlet,Diagnostics& g){
        need(old.b.size()==48&&old.p.size()==48&&old.e.size()==size_t(n)&&valid(old.b)&&valid(old.e),"invalid_stock");
        for(double x:old.r)need(std::isfinite(x)&&x>=0,"invalid_reservoir");
        for(int j=0;j<48;++j)need(old.b[j]/vb[j]<=content(150),"blood_domain");
        for(int i=0;i<n;++i)need(old.e[i]/a[i]<=150+1e-10,"water_domain");
        double cin=content(inlet);State s=old;
        if(opened){need(s.r[1]>0&&std::abs(s.r[0]/s.r[1]-cin)<=1e-9,"feed_composition");
            double dv=qflow*m*dt/60;
            for(int j=0;j<48;++j){double enter=j%8==0?cin:old.b[j-1]/vb[j-1];double leave=old.b[j]/vb[j];s.b[j]+=dv*(enter-leave);
                if(j%8==0)g.in+=dv*enter;if(j%8==7)g.out+=dv*leave;}
            need(s.r[0]>=g.in&&s.r[1]>=6*dv,"feed_exhausted");s.r[0]-=g.in;s.r[1]-=6*dv;s.r[2]+=g.out;s.r[3]+=6*dv;
            need(s.r[1]>0&&std::abs(s.r[0]/s.r[1]-cin)<=1e-9,"returned_feed");}
        need(valid(s.b),"convection_stock");V before=s.b,aq(48,0),intercept(48,0);
        for(int i=0;i<n;++i)if(node[i]>=0){aq[node[i]]+=old.e[i];intercept[node[i]]+=a[i]*old.e[i]/den[i];}
        for(int j=0;j<48;++j){double target=s.b[j]+aq[j]-intercept[j],lo=0,hi=150,x=old.p[j];
            need(target>=0&&target<=150*linear[j]+hb[j]*saturation(150),"root_bracket");
            bool done=false;for(int it=0;it<100;++it){double rr=linear[j]*x+hb[j]*saturation(x)-target;
                if(std::abs(rr)<=1e-14){g.root=std::max(g.root,std::abs(rr));done=true;break;}
                if(rr<0)lo=x;else hi=x;x=(lo+hi)/2;}
            need(done,"root_iterations");s.p[j]=x;s.b[j]=vb[j]*content(x);}
        V debit(48,0);
        for(int i=0;i<n;++i)if(node[i]>=0){int j=node[i];double pe=(old.e[i]+d[i]*s.p[j])/den[i];s.e[i]=a[i]*pe;
            double transfer=d[i]*(s.p[j]-pe);debit[j]+=transfer;g.wall=std::max(g.wall,std::abs(s.e[i]-old.e[i]-transfer));}
        for(int j=0;j<48;++j)g.wall=std::max(g.wall,std::abs(s.b[j]-before[j]+debit[j]));
        // Solve for concentrations with accessible-volume capacities; never normalize.
        V wall=s.e,c=s.e;int nrhs=1,info;char low='L';dpbtrs_(&low,&n,&kd,&nrhs,factor.data(),&ld,c.data(),&n,&info);need(info==0,"water_solve");
        V net(n,0);for(auto e:edges){double t=dt*e.k*(c[e.i]-c[e.j]);net[e.i]-=t;net[e.j]+=t;}
        for(int i=0;i<n;++i){s.e[i]=ve[i]*c[i];g.water=std::max(g.water,std::abs(s.e[i]-wall[i]-net[i]));need(c[i]/.002<=150+1e-10,"solved_pressure");}
        g.body=total(s.b)+total(s.e)-total(old.b)-total(old.e)-g.in+g.out;
        need(valid(s.b)&&valid(s.e)&&g.wall<1e-11&&g.water<1e-11&&std::abs(g.body)<1e-10,"proposed_ledger");return s;
    }
};

void raw(std::ostream& f,const V& v){f.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));need(bool(f),"output_write");}
void commit(BioFVM::Microenvironment& env,const Operator& op,const State& s){for(int i=0;i<2500;++i)env.density_vector(i)[0]=s.e[i]/op.ve[i];}
double roundtrip(BioFVM::Microenvironment& env,const Operator& op,const State& s){double err=0;for(int i=0;i<2500;++i)err=std::max(err,std::abs(env.density_vector(i)[0]*op.ve[i]-s.e[i]));return err;}
int main(int argc,char** argv){try{
    need(argc==7,"arguments input output_prefix dt_s multiplicity valve mode");omp_set_num_threads(1);
    std::string input=argv[1],prefix=argv[2],mode=argv[6];double dt=std::stod(argv[3]),m=std::stod(argv[4]);bool opened=std::stoi(argv[5])!=0;
    Operator op(input,dt,m);State s=op.initial(input);auto& env=BioFVM::microenvironment;
    env.resize_densities(1);env.set_density(0,"oxygen","mM",0.,0.);
    env.resize_space(-500.,500.,-500.,500.,-10.,10.,20.,20.,20.);
    need(env.number_of_voxels()==2500,"native_grid");
    for(int i=0;i<2500;++i){auto& p=env.mesh.voxels[i].center;need(p[0]==-490+20*(i%50)&&p[1]==-490+20*(i/50)&&p[2]==0,"native_voxel_order");}
    commit(env,op,s);PhysiCell::PhysiCell_globals.current_time=0;
    if(mode=="faults"){
        int passed=0;
        for(int fault=0;fault<6;++fault){State bad=s;double inlet=65;
            if(fault==0)bad.r[0]*=.9;if(fault==1){bad.r[1]=1e-9;bad.r[0]=content(65)*1e-9;}
            if(fault==2){bad.r[0]=0;bad.r[1]=0;}if(fault==3)bad.e[0]=-1;
            if(fault==4)inlet=-1;if(fault==5)inlet=151;
            bool rejected=false;try{Diagnostics d;op.propose(bad,true,inlet,d);}catch(const std::exception&){rejected=true;}
            need(rejected&&PhysiCell::PhysiCell_globals.current_time==0&&roundtrip(env,op,s)<1e-15,"fault_nonmutation");++passed;}
        bool cfl=false;try{Operator bad(input,60,m);}catch(const std::exception&){cfl=true;}need(cfl,"CFL_fault");++passed;
        State eq=s;for(int j=0;j<48;++j){eq.p[j]=35;eq.b[j]=op.vb[j]*content(35);}for(int i=0;i<2500;++i)eq.e[i]=op.a[i]*35;
        Diagnostics g;State out=op.propose(eq,false,65,g);double error=0;for(int i=0;i<2500;++i)error=std::max(error,std::abs(out.e[i]-eq.e[i]));need(error<1e-12,"equilibrium");++passed;
        std::ofstream report(prefix+"_faults.json");report<<"{\"passed\":"<<passed<<",\"field_and_clock_unchanged\":true}\n";return 0;
    }
    need(mode=="run"||mode=="no_diffusion"||mode=="no_wall_flow","mode");
    Operator actual(input,dt,m,mode=="no_diffusion"?0:1000,mode=="no_wall_flow"?0:.01,mode=="no_wall_flow"?0:3.492);
    std::ofstream states(prefix+"_states.f64",std::ios::binary),ledger(prefix+"_ledger.csv");ledger<<std::setprecision(17)<<"time_s,feed_pmol,feed_nL,waste_pmol,waste_nL,blood_pmol,water_pmol,input_pmol,output_pmol,mass_error,fluid_error,wall_error,water_error,root_error,body_error,field_error,clock_min\n";
    std::array<double,11> times={{0,.03,.06,.12,.3,.6,1.5,6,15,30,60}};int next=0;
    auto snapshot=[&](){V e(2500);for(int i=0;i<2500;++i)e[i]=env.density_vector(i)[0]*actual.ve[i];raw(states,s.b);raw(states,e);raw(states,s.p);raw(states,V(s.r.begin(),s.r.end()));};
    snapshot();++next;int steps=std::lround(60/dt);need(std::abs(steps*dt-60)<1e-10,"duration_grid");
    double mass=total(s.b)+total(s.e)+s.r[0]+s.r[2],fluid=total(actual.vb)+total(actual.ve)+s.r[1]+s.r[3],cin=0,cout=0,drift=0;
    for(int step=1;step<=steps;++step){
        // The actual BioFVM density vector supplies the old aqueous state.
        for(int i=0;i<2500;++i)s.e[i]=env.density_vector(i)[0]*actual.ve[i];
        Diagnostics d;State proposed=actual.propose(s,opened,65,d);
        double merr=total(proposed.b)+total(proposed.e)+proposed.r[0]+proposed.r[2]-mass;
        double ferr=total(actual.vb)+total(actual.ve)+proposed.r[1]+proposed.r[3]-fluid;
        need(std::abs(merr)<1e-8&&std::abs(ferr)<1e-8,"global_ledger");
        commit(env,actual,proposed);s=std::move(proposed);PhysiCell::PhysiCell_globals.current_time=step*dt/60;
        double ferror=roundtrip(env,actual,s);need(ferror<1e-12,"native_field_roundtrip");drift=std::max(drift,ferror);cin+=d.in;cout+=d.out;
        ledger<<step*dt;for(double x:s.r)ledger<<','<<x;
        ledger<<','<<total(s.b)<<','<<total(s.e)<<','<<cin<<','<<cout<<','<<merr<<','<<ferr<<','<<d.wall<<','<<d.water<<','<<d.root<<','<<d.body<<','<<ferror<<','<<PhysiCell::PhysiCell_globals.current_time<<'\n';
        if(next<11&&std::abs(step*dt-times[next])<1e-10){snapshot();++next;}
    }
    need(next==11&&PhysiCell::PhysiCell_globals.current_time==1,"native_clock_samples");
    std::ofstream report(prefix+"_native.json");report<<std::setprecision(17)<<"{\"clock_min\":"<<PhysiCell::PhysiCell_globals.current_time<<",\"steps\":"<<steps<<",\"agents\":0,\"grid_voxels\":2500,\"field_roundtrip_error_pmol\":"<<drift<<",\"standard_diffusion_calls\":0,\"cell_scheduler_calls\":0,\"LP_calls\":0}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"REJECTED "<<e.what()<<std::endl;return 1;}}
