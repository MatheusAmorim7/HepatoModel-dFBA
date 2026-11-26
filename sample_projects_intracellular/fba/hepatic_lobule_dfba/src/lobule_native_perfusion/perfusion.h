// Same tissue graph and finite compartment as the sealed bath, with explicit inlet/outlet.
#include "../lobule_extension/bath.h"
struct PerfusionTransport:BathTransport {
    double F,Cin;
    PerfusionTransport(double dt,double flow,double inlet):BathTransport(1,true,2*dt),F(flow),Cin(inlet) {
        bb+=h*F;schur+=h*F;need(F>=0&&Cin>=0,"Nonnegative_source_parameters");
    }
    std::vector<double> advance(const std::vector<double>& q,double& incoming,double& outgoing,double& residual) {
        need(q.size()==2501,"Perfusion_shape");auto C=solve(std::vector<double>(q.begin(),q.begin()+n));
        double b=q[n]+h*F*Cin;for(double value:C)b+=h*g*value;b/=schur;
        for(int i=0;i<n;++i)C[i]+=z[i]*b;C.push_back(b);
        std::vector<double> next(2501),div(2501,0);
        for(auto edge:edges){double d=h*g*(C[edge.second]-C[edge.first]);div[edge.first]+=d;div[edge.second]-=d;}
        for(int i=0;i<n;++i){double d=h*g*(C[n]-C[i]);div[i]+=d;div[n]-=d;}
        incoming=h*F*Cin;outgoing=h*F*C[n];div[n]+=incoming-outgoing;residual=0;long double mass=0;
        for(int i=0;i<2501;++i){next[i]=C[i]*(i<n?.008:20);need(std::isfinite(next[i])&&next[i]>=0,"Negative_source_proposal_no_clip");
            residual=std::max(residual,std::abs(next[i]-q[i]-div[i]));mass+=(long double)next[i]-q[i];}
        need(residual<1e-10&&std::abs(double(mass)-incoming+outgoing)<1e-9,"Source_ledger");return next;
    }
};
