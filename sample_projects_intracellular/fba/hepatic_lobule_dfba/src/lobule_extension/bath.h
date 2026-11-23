extern "C" {
void dpbtrf_(char*,int*,int*,double*,int*,int*);
void dpbtrs_(char*,int*,int*,int*,double*,int*,double*,int*,int*);
}
class BathTransport {
public:
    int n=2500,kd=50,ld=51;double h,ratio,g=3.504,bb=0,schur=0;
    bool enabled;std::vector<double> band,z;std::vector<std::pair<int,int>> edges;
    BathTransport(double R,bool E,double dt):h(dt/2),ratio(R),enabled(E),band(51*2500,0){
        for(int i=0;i<n;++i)band[ld*i]=.008;
        if(enabled){
            auto edge=[&](int i,int j){edges.push_back({i,j});band[ld*i]+=h*g;band[ld*j]+=h*g;band[(j-i)+ld*i]-=h*g;};
            for(int y=0;y<50;++y)for(int x=0;x<50;++x){int i=x+50*y;if(x<49)edge(i,i+1);if(y<49)edge(i,i+50);}
            if(ratio>0)for(int i=0;i<n;++i)band[ld*i]+=h*g;
        }
        char low='L';int info=0;dpbtrf_(&low,&n,&kd,band.data(),&ld,&info);need(info==0,"Transport_band_factorization");
        if(ratio>0){need(enabled,"Bath_requires_contacts");z=solve(std::vector<double>(n,h*g));bb=20*ratio+h*n*g;schur=bb;
            for(double value:z)schur-=h*g*value;need(schur>0,"Positive_bath_Schur");}
    }
    std::vector<double> solve(std::vector<double> b){char low='L';int one=1,info=0;dpbtrs_(&low,&n,&kd,&one,band.data(),&ld,b.data(),&n,&info);need(info==0,"Transport_band_solve");return b;}
    std::vector<double> half(const std::vector<double>& q,double& residual,double& mass){
        need(q.size()==size_t(n+(ratio>0)),"Transport_shape");auto C=solve(std::vector<double>(q.begin(),q.begin()+n));
        if(ratio>0){double b=q[n];for(double value:C)b+=h*g*value;b/=schur;for(int i=0;i<n;++i)C[i]+=z[i]*b;C.push_back(b);}
        std::vector<double> next(C.size()),div(C.size(),0);
        for(auto edge:edges){double d=h*g*(C[edge.second]-C[edge.first]);div[edge.first]+=d;div[edge.second]-=d;}
        if(ratio>0)for(int i=0;i<n;++i){double d=h*g*(C[n]-C[i]);div[i]+=d;div[n]-=d;}
        residual=0;long double sum=0;
        for(size_t i=0;i<C.size();++i){next[i]=C[i]*(i<size_t(n)?.008:20*ratio);need(std::isfinite(next[i])&&next[i]>=0,"Negative_raw_transport_no_clipping");
            residual=std::max(residual,std::abs(next[i]-q[i]-div[i]));sum+=(long double)next[i]-q[i];}
        mass=double(sum);need(residual<1e-10&&std::abs(mass)<1e-9,"Transport_nodal_and_total_certificate");return next;
    }
};
