#include <memory>
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
struct Cache {ClpSimplex solver;std::vector<double> lo,hi;std::vector<int> members;};
struct MeasuredJoint:JointTemplate {
    double assembly=0,loading=0,solving=0,certificates=0;
    long long builds=0,hits=0,iterations=0,calls=0;
    std::map<int,std::unique_ptr<Cache>> cache;
    MeasuredJoint(LP& lp,const std::vector<ResourceField>& fields,int k,double dt):JointTemplate(lp,fields,k,dt){}
    JointAnswer run(LP& lp,int voxel,const std::vector<int>& members,const Matrix& lower,const Matrix& upper,const std::vector<double>& q,bool warm){
        auto t=Clock::now();std::vector<double> lo(nv,0),hi(nv,INF),rhs=rh;
        for(int c=0;c<k;++c)for(int j=0;j<n;++j){lo[c*(2*n+NF)+j]=lower[members[c]][j];hi[c*(2*n+NF)+j]=upper[members[c]][j];}
        for(int f=0;f<NF;++f){double reserve=q[f]>0?std::min(1e-12,q[f]*.5):0;rhs[stock_start+f]=q[f]-reserve;}
        assembly+=elapsed(t);t=Clock::now();Cache* entry=nullptr;std::unique_ptr<Cache> fresh;
        if(warm&&cache.count(voxel)){
            entry=cache.at(voxel).get();need(entry->members==members&&entry->lo==lo&&entry->hi==hi,"Cache_fixed_membership_and_bounds_required");
            for(int f=0;f<NF;++f)entry->solver.setRowUpper(stock_start+f,rhs[stock_start+f]);++hits;
        }else{
            fresh.reset(new Cache());entry=fresh.get();entry->lo=lo;entry->hi=hi;entry->members=members;
            entry->solver.setLogLevel(0);entry->solver.setPrimalTolerance(1e-10);entry->solver.setDualTolerance(1e-10);
            entry->solver.loadProblem(matrix,lo.data(),hi.data(),cost.data(),rl.data(),rhs.data());entry->solver.setOptimizationDirection(1);++builds;
            if(warm)cache.emplace(voxel,std::move(fresh));
        }
        loading+=elapsed(t);t=Clock::now();auto& p=entry->solver;p.dual();++calls;iterations+=p.numberIterations();solving+=elapsed(t);t=Clock::now();
        need(p.isProvenOptimal(),"Replay_LP_not_optimal_no_fallback");JointAnswer ans;ans.ok=true;ans.primary=p.objectiveValue();ans.v.assign(n,0);ans.cells.assign(k,std::vector<double>(n));auto x=p.primalColumnSolution();
        for(int c=0;c<k;++c)for(int j=0;j<n;++j){double v=x[c*(2*n+NF)+j];need(std::isfinite(v)&&v>=lo[c*(2*n+NF)+j]-2e-7&&v<=hi[c*(2*n+NF)+j]+2e-7,"Individual_bounds");ans.cells[c][j]=v;ans.v[j]+=v/k;}
        for(auto& values:ans.cells)for(auto& row:lp.sparse_stoichiometry){long double s=0;for(auto term:row)s+=term.second*values[term.first];need(std::abs(s)<2e-7,"Individual_Sv");}
        certificates+=elapsed(t);return ans;
    }
};

using CacheMap=std::map<int,std::unique_ptr<MeasuredJoint>>;
JointAnswer cached_select(MeasuredJoint& problem,LP& lp,int voxel,const std::vector<int>& members,const Matrix& lo,const Matrix& hi,const std::vector<double>& q){
    try{return problem.run(lp,voxel,members,lo,hi,q,true);}
    catch(const std::exception& e){
        need(std::string(e.what())=="Replay_LP_not_optimal_no_fallback","Unexpected_cache_failure");
        need(problem.cache.at(voxel)->solver.isProvenPrimalInfeasible(),"Unclassified_cache_failure");
        return JointAnswer();
    }
}
void cache_report(const CacheMap& problems,bool use_cache,long long calls,const std::string& path){
    long long solves=0,builds=0,hits=0,iterations=0;
    for(auto& item:problems){solves+=item.second->calls;builds+=item.second->builds;hits+=item.second->hits;iterations+=item.second->iterations;}
    need(solves==(use_cache?calls:0),"Cache_solve_count");
    std::ofstream f(path);f<<"{\"enabled\":"<<(use_cache?"true":"false")<<",\"joint_calls\":"<<calls<<",\"warm_calls\":"<<solves<<",\"builds\":"<<builds<<",\"hits\":"<<hits<<",\"reported_iterations\":"<<iterations<<"}\n";
}
