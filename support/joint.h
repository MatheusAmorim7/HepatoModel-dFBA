// Independent intracellular blocks share only the finite extracellular budget.
struct JointAnswer : Answer {std::vector<std::vector<double>> cells;};
struct JointTemplate {
    int n,k,nv,stock_start; CoinPackedMatrix matrix;
    std::vector<double> rl,rh,cost;
    JointTemplate(LP& lp,const std::vector<ResourceField>& fields,int K,double dt):n(lp.n),k(K),nv(K*(2*lp.n+20)),matrix(false,0,0){
        matrix.setDimensions(0,nv);cost.assign(nv,0);
        auto add=[&](CoinPackedVector& row,double lo,double hi){matrix.appendRow(row);rl.push_back(lo);rh.push_back(hi);};
        for(int c=0;c<k;++c){int off=c*(2*n+20);
            for(auto& terms:lp.sparse_stoichiometry){CoinPackedVector row;for(auto p:terms)row.insert(off+p.first,p.second);add(row,0,0);}
            for(int sign:{1,-1})for(int j=0;j<n;++j){CoinPackedVector row;row.insert(off+j,sign);row.insert(off+n+j,-1);add(row,-INF,0);}
            for(int f=0;f<20;++f){CoinPackedVector row;row.insert(off+fields[f].f.column,-1);row.insert(off+2*n+f,-1);add(row,-INF,0);}
            for(int j=0;j<n;++j)cost[off+n+j]=1./k;
        }
        stock_start=rh.size();
        for(int f=0;f<20;++f){CoinPackedVector row;for(int c=0;c<k;++c)row.insert(c*(2*n+20)+2*n+f,.12*dt/60);add(row,-INF,0);}
    }
    JointAnswer solve(LP& lp,const std::vector<int>& members,const std::vector<std::vector<double>>& lower,const std::vector<std::vector<double>>& upper,const std::vector<double>& q){
        std::vector<double> lo(nv,0),hi(nv,INF),rhs=rh;
        for(int c=0;c<k;++c)for(int j=0;j<n;++j){lo[c*(2*n+20)+j]=lower[members[c]][j];hi[c*(2*n+20)+j]=upper[members[c]][j];}
        for(int f=0;f<20;++f){double reserve=q[f]>0?std::min(1e-12,q[f]*.5):0;rhs[stock_start+f]=q[f]-reserve;}
        ClpSimplex p;p.setLogLevel(0);p.setPrimalTolerance(1e-10);p.setDualTolerance(1e-10);
        p.loadProblem(matrix,lo.data(),hi.data(),cost.data(),rl.data(),rhs.data());p.setOptimizationDirection(1);p.dual();
        JointAnswer ans;if(!p.isProvenOptimal()){need(p.isProvenPrimalInfeasible(),"Unclassified_joint_solver_failure");return ans;}
        ans.ok=true;ans.primary=p.objectiveValue();ans.v.assign(n,0);ans.cells.assign(k,std::vector<double>(n));auto x=p.primalColumnSolution();
        for(int c=0;c<k;++c)for(int j=0;j<n;++j){double v=x[c*(2*n+20)+j];need(std::isfinite(v)&&v>=lo[c*(2*n+20)+j]-2e-7&&v<=hi[c*(2*n+20)+j]+2e-7,"Joint_individual_bound");ans.cells[c][j]=v;ans.v[j]+=v/k;}
        for(auto& values:ans.cells)for(auto& row:lp.sparse_stoichiometry){long double s=0;for(auto p:row)s+=p.second*values[p.first];need(std::abs(s)<2e-7,"Joint_individual_Sv");}
        return ans;
    }
};
