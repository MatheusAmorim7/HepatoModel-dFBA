// Uniform groups share inventories, but retain separate stoichiometric problems.
struct GroupAnswer: Answer {
    double weighted_L1=0;
    std::vector<std::vector<double>> group_v;
};

class GroupLP {
    LP& reference;
    int n,G,nv,eidx;
    std::vector<int> counts;
    std::vector<double> weights,rowlo,rowhi;
    CoinPackedMatrix base;
    int col(int g,int j){return g*2*n+j;}
public:
    GroupLP(LP& lp,const std::vector<int>& sizes):reference(lp),n(lp.n),G(sizes.size()),
        nv(G*2*n+1),eidx(G*2*n),counts(sizes),base(false,0,0){
        int total=0;for(int count:counts){need(count>0,"Emptygroup");total+=count;weights.push_back(count/22.);}
        need(total==22&&G>=2&&G<=3,"Registeredgroups");base.setDimensions(0,nv);
        for(int g=0;g<G;++g){
            for(auto& terms:reference.sparse_stoichiometry){CoinPackedVector row;
                for(auto term:terms)row.insert(col(g,term.first),term.second);
                base.appendRow(row);rowlo.push_back(0);rowhi.push_back(0);}
            for(int sign:{1,-1})for(int j=0;j<n;++j){CoinPackedVector row;
                row.insert(col(g,j),sign);row.insert(col(g,n+j),-1);
                base.appendRow(row);rowlo.push_back(-INF);rowhi.push_back(0);}
        }
    }

    GroupAnswer select(const std::vector<Field>& fields,const std::vector<double>& stocks,
                       double remaining,double dt,double cap){
        CoinPackedMatrix matrix(base);auto rl=rowlo,rh=rowhi;
        const double a=SCALE*dt/60;int water=-1;
        for(size_t f=0;f<fields.size();++f){if(fields[f].name=="oxygen")continue;
            CoinPackedVector row;
            for(int g=0;g<G;++g)row.insert(col(g,fields[f].column),-a*weights[g]);
            bool coupled=fields[f].name=="alanine"||fields[f].name=="glutamine";
            if(coupled)row.insert(eidx,-1);
            if(fields[f].name=="water"){row.insert(eidx,1);water=f;}
            bool reachable=stocks[f]>0||(remaining>0&&coupled);
            double reserve=stocks[f]>0?std::min(RESERVE,stocks[f]*.5):(reachable?RESERVE:0);
            matrix.appendRow(row);rl.push_back(-INF);rh.push_back(stocks[f]-reserve);
        }
        need(water>=0,"Missinggroupwater");
        CoinPackedVector oxygen;
        for(int g=0;g<G;++g)oxygen.insert(col(g,reference.oi),-SCALE/60*weights[g]);
        matrix.appendRow(oxygen);rl.push_back(-INF);rh.push_back(cap);
        std::vector<double> lo(nv,0),hi(nv,INF),cost(nv,0);
        for(int g=0;g<G;++g)for(int j=0;j<n;++j){
            lo[col(g,j)]=reference.lower[j];hi[col(g,j)]=reference.upper[j];
            if(j==reference.oi)lo[col(g,j)]=std::max(lo[col(g,j)],-cap*60/SCALE);
            cost[col(g,n+j)]=weights[g];
        }
        double emax=std::min(remaining,std::min(stocks.back(),stocks[water]));
        need(emax>=0,"Negativegroupdeliverybudget");
        if(emax>0)emax-=std::min(RESERVE,.5*emax);hi[eidx]=emax;
        ClpSimplex problem;problem.setLogLevel(0);problem.setPrimalTolerance(1e-10);problem.setDualTolerance(1e-10);
        problem.loadProblem(matrix,lo.data(),hi.data(),cost.data(),rl.data(),rh.data());
        problem.setOptimizationDirection(1);problem.dual();GroupAnswer result;
        if(!problem.isProvenOptimal()){
            need(problem.isProvenPrimalInfeasible(),"Unclassifiedgroupprimaryfailure");return result;}
        result.primary=problem.objectiveValue();CoinPackedVector face;
        for(int g=0;g<G;++g)for(int j=0;j<n;++j)face.insert(col(g,n+j),weights[g]);
        problem.addRow(face.getNumElements(),face.getIndices(),face.getElements(),-INF,result.primary+1e-9);
        for(int j=0;j<nv;++j)problem.setObjectiveCoefficient(j,j==eidx?1:0);
        problem.primal();need(problem.isProvenOptimal(),"Groupsecondaryfailure");
        auto x=problem.primalColumnSolution();result.ok=true;result.extent=x[eidx];
        result.v.assign(n,0);result.group_v.assign(G,std::vector<double>(n));
        for(int g=0;g<G;++g){
            for(int j=0;j<n;++j){double v=x[col(g,j)];result.group_v[g][j]=v;
                need(std::isfinite(v)&&v>=lo[col(g,j)]-2e-7&&v<=hi[col(g,j)]+2e-7,"Groupboundcertificate");
                result.weighted_L1+=weights[g]*std::abs(v);}
            for(auto& row:reference.sparse_stoichiometry){long double balance=0;
                for(auto term:row)balance+=term.second*result.group_v[g][term.first];
                need(std::abs(balance)<2e-7,"GroupSvcertificate");}
        }
        // This fixed count-weighted map is the coupling definition, not a repair.
        for(int j=0;j<n;++j){long double total=0;
            for(int g=0;g<G;++g)total+=counts[g]*result.group_v[g][j];
            result.v[j]=total/22;}
        need(result.weighted_L1<=result.primary+2e-7,"GroupL1facecertificate");
        return result;
    }
};
