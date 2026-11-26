// Phase C: serial-compartment sinusoid transport for ONE field (backward Euler, banded tissue + Schur complement over K blood nodes).
// Requires `void need(bool,const char*)` to be declared before inclusion (as in the sealed chain).
#pragma once
#include <cmath>
#include <vector>
extern "C" {
void dpbtrf_(char*,int*,int*,double*,int*,int*);
void dpbtrs_(char*,int*,int*,int*,double*,int*,double*,int*,int*);
void dgetrf_(int*,int*,double*,int*,int*,int*);
void dgetrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*);
}
class SerialTransport {
public:
    static const int N=2500,NX=50,K=48,SEG=8,SECT=6,KD=50,LD=51;
    double dt,ge,gw,Fs,V,vox=.008;
    std::vector<int> node_of;std::vector<double> band,W,S;std::vector<int> piv;
    SerialTransport(double step,double edge_conductance,double wall_conductance,double sector_flow,double node_volume,const std::vector<int>& nodes)
        :dt(step),ge(edge_conductance),gw(wall_conductance),Fs(sector_flow),V(node_volume),node_of(nodes),band(LD*N,0),W(N*K,0),S(K*K,0),piv(K){
        need((int)node_of.size()==N,"Node_map_size");std::vector<int> m(K,0);for(int i=0;i<N;++i){need(node_of[i]>=0&&node_of[i]<K,"Node_range");++m[node_of[i]];}
        for(int k=0;k<K;++k)need(m[k]>0,"Empty_node");need(dt>0&&ge>=0&&gw>=0&&Fs>=0&&V>0,"Transport_parameters");
        for(int i=0;i<N;++i)band[LD*i]=vox/dt+gw;
        auto edge=[&](int i,int j){band[LD*i]+=ge;band[LD*j]+=ge;band[(j-i)+LD*i]-=ge;};
        for(int y=0;y<NX;++y)for(int x=0;x<NX;++x){int i=x+NX*y;if(x<NX-1)edge(i,i+1);if(y<NX-1)edge(i,i+NX);}
        char low='L';int n=N,kd=KD,ld=LD,info=0;dpbtrf_(&low,&n,&kd,band.data(),&ld,&info);need(info==0,"Tissue_band_factorization");
        for(int k=0;k<K;++k){std::vector<double> col(N,0);for(int i=0;i<N;++i)if(node_of[i]==k)col[i]=1;col=solve_tissue(col);for(int i=0;i<N;++i)W[i+N*k]=col[i];}
        for(int k=0;k<K;++k){S[k+K*k]=V/dt+Fs+gw*m[k];if(k%SEG!=0)S[k+K*(k-1)]=-Fs;}
        for(int i=0;i<N;++i)for(int l=0;l<K;++l)S[node_of[i]+K*l]-=gw*gw*W[i+N*l];
        int k2=K,lda=K;dgetrf_(&k2,&k2,S.data(),&lda,piv.data(),&info);need(info==0,"Schur_factorization");
    }
    std::vector<double> solve_tissue(std::vector<double> b) const{char low='L';int n=N,kd=KD,one=1,ld=LD,ldb=N,info=0;dpbtrs_(&low,&n,&kd,&one,const_cast<double*>(band.data()),&ld,b.data(),&ldb,&info);need(info==0,"Tissue_band_solve");return b;}
    // c_old: N concentrations (mM); b_old: K node concentrations (mM); Cin: inlet concentration (mM). Returns new amounts via c,b (mM).
    void advance(const std::vector<double>& c_old,const std::vector<double>& b_old,double Cin,std::vector<double>& c,std::vector<double>& b,double& incoming,double& outgoing,double& residual) const{
        need(c_old.size()==(size_t)N&&b_old.size()==(size_t)K&&Cin>=0,"Advance_shape");
        std::vector<double> rc(N);for(int i=0;i<N;++i)rc[i]=vox/dt*c_old[i];auto y=solve_tissue(rc);
        std::vector<double> r(K);for(int k=0;k<K;++k)r[k]=V/dt*b_old[k]+(k%SEG==0?Fs*Cin:0.);for(int i=0;i<N;++i)r[node_of[i]]+=gw*y[i];
        char tr='N';int k2=K,one=1,lda=K,ldb=K,info=0;auto f=S;dgetrs_(&tr,&k2,&one,f.data(),&lda,const_cast<int*>(piv.data()),r.data(),&ldb,&info);need(info==0,"Schur_solve");
        b=r;c=y;for(int l=0;l<K;++l)for(int i=0;i<N;++i)c[i]+=gw*W[i+N*l]*b[l];
        // independent certificate of the discrete equations (amounts, pmol) and global mass ledger
        std::vector<double> rt(N,0),rn(K,0);long double before=0,after=0;
        for(int i=0;i<N;++i){rt[i]=vox*(c[i]-c_old[i])-dt*gw*(b[node_of[i]]-c[i]);before+=vox*c_old[i];after+=vox*c[i];}
        auto add=[&](int i,int j){double d=dt*ge*(c[j]-c[i]);rt[i]-=d;rt[j]+=d;};
        for(int y2=0;y2<NX;++y2)for(int x=0;x<NX;++x){int i=x+NX*y2;if(x<NX-1)add(i,i+1);if(y2<NX-1)add(i,i+NX);}
        residual=0;for(double v:rt)residual=std::max(residual,std::abs(v));
        for(int k=0;k<K;++k){rn[k]=V*(b[k]-b_old[k])-dt*(Fs*((k%SEG==0?Cin:b[k-1])-b[k]));before+=V*b_old[k];after+=V*b[k];}
        for(int i=0;i<N;++i)rn[node_of[i]]-=dt*gw*(c[i]-b[node_of[i]]);
        for(double v:rn)residual=std::max(residual,std::abs(v));
        incoming=dt*Fs*Cin*SECT;outgoing=0;for(int k=SEG-1;k<K;k+=SEG)outgoing+=dt*Fs*b[k];
        double ledger=double(after-before)-incoming+outgoing;
        need(residual<1e-10&&std::abs(ledger)<1e-9,"Transport_equation_and_ledger_certificate");
        for(double v:c)need(std::isfinite(v)&&v>=0,"Negative_tissue_no_clip");for(double v:b)need(std::isfinite(v)&&v>=0,"Negative_node_no_clip");
    }
};
