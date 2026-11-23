// Isolated PhysiCell agents with joint uniform group Clp dFBA and FV gas.
#include "src/b0_dims.h"
#include "core/PhysiCell.h"
#include "modules/PhysiCell_standard_modules.h"
#include "addons/dFBA/src/dfba_intracellular.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <omp.h>

extern "C" {
void dgetrf_(int*,int*,double*,int*,int*,int*);
void dgetrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*);
}
using namespace PhysiCell;
using namespace BioFVM;
using PhysiCelldFBA::dFBAIntracellular;
const double SCALE=2.64,INF=1e100,RESERVE=1e-8;
void need(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
std::string bare(std::string id){return id.substr(0,2)=="R_" ? id.substr(2):id;}
struct Field{std::string name,rid;int column;double carbon,nitrogen;};

struct Gas {
    int n;double area,depth,D=12500,Cstar=.2,dt,dz,volume,g,gtop;
    std::vector<double> lu,z;std::vector<int> piv;
    Gas(int N,double A,double L,double h):n(N),area(A),depth(L),dt(h),dz(L/N),
        volume(A*dz*1e-6),g(D*A/dz*1e-6),gtop(2*D*A/dz*1e-6),lu(N*N,0),piv(N){
        auto at=[&](int i,int j)->double&{return lu[i+n*j];};
        for(int i=0;i<n;++i)at(i,i)=1;
        double x=dt*g/volume;
        for(int i=0;i<n-1;++i){at(i,i)+=x;at(i+1,i+1)+=x;at(i,i+1)-=x;at(i+1,i)-=x;}
        at(n-1,n-1)+=dt*gtop/volume;
        int info=0;dgetrf_(&n,&n,lu.data(),&n,piv.data(),&info);need(info==0,"FV factorization");
        z.assign(n,0);z[0]=dt/volume;z=solve(z);
    }
    std::vector<double> solve(std::vector<double> rhs){
        char trans='N';int one=1,info=0;
        dgetrs_(&trans,&n,&one,lu.data(),&n,piv.data(),rhs.data(),&n,&info);
        need(info==0,"FV solve");return rhs;
    }
    std::vector<double> baseline(std::vector<double> old){
        old.back()+=dt*gtop*Cstar/volume;return solve(old);
    }
    double half(){return dz/(2*D*area*1e-6);}
    double limit(double center,double resistance){
        need(center>=0,"Negative input oxygen");if(center==0)return 0;
        double lo=0,hi=std::min(center/resistance,6*SCALE/60);
        for(int k=0;k<90;++k){double mid=(lo+hi)/2,c=center-resistance*mid;
            double gap=mid-6*SCALE/60*c/(.01+c);if(gap>0)hi=mid;else lo=mid;}
        return (lo+hi)/2;
    }
    double endpoint_limit(const std::vector<double>& old){auto b=baseline(old);return limit(b[0],z[0]+half());}
};

struct Answer{bool ok=false;double primary=0,extent=0;std::vector<double> v;};
class LP {
public:
    dFBAModel model;int n,nv,oi;std::vector<double> lower,upper;std::vector<std::string> ids;
    CoinPackedMatrix base;std::vector<double> rowlo,rowhi;
    std::vector<std::vector<std::pair<int,double>>> sparse_stoichiometry;
    LP(const std::string& path):base(false,0,0){
        model.initModel(path.c_str());need(model.isInitialized(),"SBML load");
        n=model.reactions.size();nv=2*n+1;need(n==NR,"ExpectedNRreactions");base.setDimensions(0,nv);
        for(auto rx:model.reactions){ids.push_back(rx->getId());lower.push_back(rx->getLowerBound());upper.push_back(rx->getUpperBound());}
        oi=model.getReactionIndex("R_EX_HC00017_s");
        {const char* e=getenv("ATPM_EXPECTED");double ex=e?atof(e):5.0;need(model.getReactionLowerBound("R_CULTURE_ATPM")==ex && model.getReactionUpperBound("R_CULTURE_ATPM")==ex,"ATPM contract");}
        need(model.getReactionLowerBound("R_MEDIUM_ALAGLN_HYD")==0 && model.getReactionUpperBound("R_MEDIUM_ALAGLN_HYD")==0,"Doublecounted cleavage");
        for(auto met:model.metabolites){CoinPackedVector row;std::vector<std::pair<int,double>> terms_for_audit;
            for(int j=0;j<n;++j){auto terms=model.reactions[j]->getMetabolites();auto it=terms.find(met->getId());if(it!=terms.end()){row.insert(j,it->second);terms_for_audit.push_back({j,it->second});}}
            sparse_stoichiometry.push_back(terms_for_audit);
            base.appendRow(row);rowlo.push_back(0);rowhi.push_back(0);}
        for(int sign:{1,-1})for(int j=0;j<n;++j){CoinPackedVector row;row.insert(j,sign);row.insert(n+j,-1);
            base.appendRow(row);rowlo.push_back(-INF);rowhi.push_back(0);}
    }
    Answer select(const std::vector<Field>& fields,const std::vector<double>& stocks,double remaining,double dt,double cap){
        CoinPackedMatrix matrix(base);auto rl=rowlo,rh=rowhi;
        double a=SCALE*dt/60;int water=-1;
        for(size_t f=0;f<fields.size();++f){if(fields[f].name=="oxygen")continue;
            CoinPackedVector row;row.insert(fields[f].column,-a);
            bool coupled=fields[f].name=="alanine"||fields[f].name=="glutamine";
            if(coupled)row.insert(2*n,-1);
            if(fields[f].name=="water"){row.insert(2*n,1);water=f;}
            bool reachable=stocks[f]>0||(remaining>0&&coupled);
            double reserve=stocks[f]>0 ? std::min(RESERVE,stocks[f]*.5):(reachable?RESERVE:0);
            matrix.appendRow(row);rl.push_back(-INF);rh.push_back(stocks[f]-reserve);
        }
        need(water>=0,"Missingwater");auto lo=lower,hi=upper;
        lo[oi]=std::max(lo[oi],-cap*60/SCALE);
        lo.resize(nv,0);hi.resize(nv,INF);
        double emax=std::min(remaining,std::min(stocks.back(),stocks[water]));
        need(emax>=0,"Negative delivery budget");if(emax>0)emax-=std::min(RESERVE,.5*emax);hi[2*n]=emax;
        std::vector<double> cost(nv,0);for(int j=0;j<n;++j)cost[n+j]=1;
        ClpSimplex problem;problem.setLogLevel(0);problem.setPrimalTolerance(1e-10);problem.setDualTolerance(1e-10);
        problem.loadProblem(matrix,lo.data(),hi.data(),cost.data(),rl.data(),rh.data());problem.setOptimizationDirection(1);
        problem.dual();Answer result;
        if(!problem.isProvenOptimal()){need(problem.isProvenPrimalInfeasible(),"Unclassifiedprimarysolverfailure");return result;}
        result.primary=problem.objectiveValue();CoinPackedVector face;
        for(int j=0;j<n;++j)face.insert(n+j,1);
        problem.addRow(face.getNumElements(),face.getIndices(),face.getElements(),-INF,result.primary+1e-9);
        for(int j=0;j<nv;++j)problem.setObjectiveCoefficient(j,j==2*n?1:0);
        problem.primal();need(problem.isProvenOptimal(),"Secondarysolverfailure");
        auto x=problem.primalColumnSolution();result.ok=true;result.extent=x[2*n];result.v.assign(x,x+n);
        for(int j=0;j<n;++j)need(std::isfinite(x[j])&&x[j]>=lo[j]-2e-7&&x[j]<=hi[j]+2e-7,"Fluxboundcertificate");
        for(auto met:model.metabolites){long double sum=0;
            for(int j=0;j<n;++j)sum+=model.reactions[j]->getStoichCoefficient(met->getId())*x[j];
            need(std::abs(sum)<2e-7,"Sv certificate");}
        return result;
    }
};

#include "groups.h"

int main(int argc,char**argv){
    try{
        need(argc==2,"ExpectedfixtureXML");pugi::xml_document doc;need(doc.load_file(argv[1]),"Readfixture");
        auto cfg=doc.child("native_2d_fixture");need(cfg,"Missingfixture");
        const int force_agent_failure_step=cfg.attribute("force_agent_failure_step").as_int(-1);
        int n=cfg.attribute("n").as_int();double dt=cfg.attribute("dt_min").as_double(),T=cfg.attribute("T_min").as_double();
        need(((n==16&&(dt==2.5||dt==1.25||dt==.625))||(n==32&&dt==1.25))&&T==1440&&cfg.attribute("agents").as_int()==22,"Unregisteredgridtimepopulation");
        std::string out=cfg.attribute("output").value(),sbml=cfg.attribute("sbml").value();
        PhysiCell_settings.omp_num_threads=1;PhysiCell_settings.folder=out;omp_set_num_threads(1);SeedRandom(0);
        LP lp(sbml);std::vector<int> group_counts,agent_group;
        for(auto group:cfg.children("group")){int g=group_counts.size(),count=group.attribute("count").as_int();
            need(group.attribute("id").as_int()==g,"Grouporder");group_counts.push_back(count);
            for(int k=0;k<count;++k)agent_group.push_back(g);}
        need(agent_group.size()==22,"Groupcount");GroupLP grouped(lp,group_counts);
        std::vector<Field> fields;std::vector<double> stocks;std::map<std::string,int> fi;
        for(auto node:cfg.children("field")){Field f;f.name=node.attribute("name").value();f.rid="R_"+std::string(node.attribute("reaction").value());
            f.column=lp.model.getReactionIndex(f.rid);need(f.column>=0,"Missingfieldreaction");
            f.carbon=node.attribute("C").as_double();f.nitrogen=node.attribute("N").as_double();fi[f.name]=fields.size();fields.push_back(f);stocks.push_back(node.attribute("stock_pmol").as_double());
            auto terms=lp.model.reactions[f.column]->getMetabolites();need(terms.size()==1&&terms.begin()->second==-1,"Exchangedirection");}
        need(fields.size()==NF,"AllNFexchangesrequired");stocks.push_back(cfg.child("dipeptide").attribute("stock_pmol").as_double());
        const auto initial=stocks;double used=0,access=cfg.attribute("access").as_double();need(access==1,"Registeredaccess1");
        double area=cfg.attribute("area_um2").as_double(),depth=cfg.attribute("depth_um").as_double();
        Gas gas(n,area,depth,dt),midgas(n,area,depth,dt/2);std::vector<double> oxygen(n,.2);
        double width=std::sqrt(area),V=area*depth*1e-6;
        for(size_t f=0;f<=fields.size();++f){std::string name=f==fields.size()?"dipeptide":fields[f].name;
            if(f==0)microenvironment.set_density(0,name,"mM",0,0);else microenvironment.add_density(name,"mM",0,0);}
        microenvironment.resize_space(0.,depth,0.,width,0.,width,n,1,1);
        for(int i=0;i<n;++i)microenvironment.mesh.voxels[i].volume=area*depth/n;
        create_cell_container_for_microenvironment(microenvironment,30.);initialize_default_cell_definition();
        cell_defaults.name="hepatocyte";cell_defaults.type=0;cell_defaults.functions.update_phenotype=nullptr;
        cell_defaults.functions.custom_cell_rule=nullptr;cell_defaults.functions.volume_update_function=nullptr;
        std::vector<Cell*> agents;std::set<const void*> models,reaction_states;
        for(int k=0;k<22;++k){Cell* cell=create_cell();cell->type_name="hepatocyte";cell->is_movable=false;
            need(cell->assign_position(microenvironment.mesh.voxels[0].center),"Agentplacement");
            auto intra=new dFBAIntracellular();intra->experimental_culture=true;intra->use_metabolic_death=false;intra->sbml_filename=sbml;
            intra->sbml_model.initModel(sbml.c_str());cell->phenotype.intracellular=intra;agents.push_back(cell);models.insert(intra);
            need(intra->sbml_model.reactions.size()==NR,"Agentmodelsize");
            for(auto rx:intra->sbml_model.reactions)reaction_states.insert(rx);}
        need(models.size()==22&&reaction_states.size()==22*NR,"Sharedintracellularstate");
        std::ofstream mapping(out+"/agent_groups.csv");mapping<<"cell_id,group,count,physical_scale\n";
        for(size_t k=0;k<agents.size();++k)mapping<<agents[k]->ID<<','<<agent_group[k]<<','<<group_counts[agent_group[k]]<<",0.12\n";
        auto serialize_agent_fluxes=[&](const std::string& filename,int step){
            std::ofstream file(out+"/"+filename);file<<std::setprecision(17)<<"last_committed_step,cell_id,reaction,flux,intrinsic_lower,intrinsic_upper\n";
            for(auto cell:agents){auto intra=static_cast<dFBAIntracellular*>(cell->phenotype.intracellular);
                for(auto rx:intra->sbml_model.reactions)file<<step<<','<<cell->ID<<','<<bare(rx->getId())<<','<<rx->getFluxValue()<<','<<rx->getLowerBound()<<','<<rx->getUpperBound()<<'\n';}
        };
        serialize_agent_fluxes("agent_fluxes_initial.csv",-1);
        std::ofstream flux(out+"/fluxes.csv"),ledger(out+"/ledger.csv"),profiles(out+"/oxygen.csv"),objectives(out+"/objectives.csv");
        for(auto stream:{&flux,&ledger,&profiles,&objectives})*stream<<std::setprecision(17);
        std::ofstream agent_audit(out+"/agent_audit.csv");agent_audit<<std::setprecision(17);
        std::ofstream group_fluxes(out+"/group_fluxes.csv"),equivalence(out+"/equivalence.csv");
        group_fluxes<<std::setprecision(17)<<"step,group,count,reaction,flux\n";
        equivalence<<std::setprecision(17)<<"step,aggregate_primary,group_primary,aggregate_extent_pmol,group_extent_pmol,max_mean_flux_difference,weighted_group_L1,mean_L1\n";
        agent_audit<<"step,cell_id,max_selected_difference,max_Sv,max_bound_error,ATPM,oxygen,urea,serine\n";
        int atp_col=lp.model.getReactionIndex("R_CULTURE_ATPM"),urea_col=lp.model.getReactionIndex("R_EX_HC00089_s"),serine_col=lp.model.getReactionIndex("R_EX_HC00068_s");
        flux<<"step,reaction,flux\n";ledger<<"step,field,before_pmol,after_pmol,cell_exchange_pmol,delivery_pmol,residual_pmol\n";
        profiles<<"step,voxel,mM\n";objectives<<"step,primary_L1,selected_L1,extent_pmol,oxygen_cap_pmol_min,bottom_face_mM,gas_gain_pmol,oxygen_residual_pmol\n";
        for(int i=0;i<n;++i)profiles<<-1<<','<<i<<','<<oxygen[i]<<'\n';
        auto expose=[&](){for(int i=0;i<n;++i){for(size_t f=0;f<fields.size();++f)microenvironment(i)[f]=fields[f].name=="oxygen"?oxygen[i]:stocks[f]/V;
            microenvironment(i)[fields.size()]=stocks.back()/V;}};expose();
        int completed=0;std::string status="completed",detail="";double maxres=0;
        for(int step=0;step<int(T/dt);++step){
            auto b=gas.baseline(oxygen),bm=midgas.baseline(oxygen);
            double cap=std::min(gas.limit(oxygen[0],gas.half()),std::min(gas.limit(b[0],gas.z[0]+gas.half()),midgas.limit(bm[0],midgas.z[0]+midgas.half())));
            auto reference=lp.select(fields,stocks,access*initial.back()-used,dt,cap);
            auto result=grouped.select(fields,stocks,access*initial.back()-used,dt,cap);
            need(reference.ok==result.ok,"Uniformgroupfeasibilitymismatch");
            if(!result.ok){status="infeasible";detail="Rejectednextstep";break;}
            try{
                need(std::abs(reference.primary-result.primary)<2e-6,"Uniformgroupprimarymismatch");
                need(std::abs(reference.extent-result.extent)<2e-6,"Uniformgroupdeliverymismatch");
                double e=result.extent;need(e>=0&&e<=access*initial.back()-used,"Deliverybudget");
                need(e<=stocks.back()&&e<=stocks[fi.at("water")],"Deliveryphysicalstock");
                need(std::abs(stocks.back()-(initial.back()-used))<=std::max(1e-12,1e-12*std::abs(initial.back()-used)),"Deliveryhistory");
                std::vector<double> next=stocks,delta(stocks.size(),0),delivery(stocks.size(),0);
                delivery.back()=-e;delivery[fi.at("water")]=-e;delivery[fi.at("alanine")]=e;delivery[fi.at("glutamine")]=e;
                double cn[2]={0,0};
                for(size_t f=0;f<fields.size();++f){if(fields[f].name=="oxygen")continue;
                    delta[f]=SCALE*dt/60*result.v[fields[f].column];
                    cn[0]+=fields[f].carbon*result.v[fields[f].column];cn[1]+=fields[f].nitrogen*result.v[fields[f].column];}
                need(std::abs(cn[0])<2e-5&&std::abs(cn[1])<2e-5,"ExternalCN");
                for(size_t f=0;f<stocks.size();++f){next[f]=(long double)stocks[f]+delivery[f]+delta[f];need(next[f]>=0&&std::isfinite(next[f]),"Negative raw stock; no repair");}
                double sink=-SCALE/60*result.v[lp.oi];need(sink>=0,"NegativeO2sink");
                std::vector<double> end(n),middle(n);for(int i=0;i<n;++i){end[i]=b[i]-gas.z[i]*sink;middle[i]=bm[i]-midgas.z[i]*sink;need(end[i]>=0,"Negative raw oxygen");}
                double bottom=end[0]-gas.half()*sink,midbottom=middle[0]-midgas.half()*sink,startbottom=oxygen[0]-gas.half()*sink;
                for(double c:{startbottom,midbottom,bottom})need(c>=-1e-12&&sink<=6*SCALE/60*c/(.01+c)+1e-9,"LocalMM");
                double gain=dt*gas.gtop*(.2-end.back()),before=0,after=0;
                for(int i=0;i<n;++i){before+=oxygen[i]*gas.volume;after+=end[i]*gas.volume;}
                double ores=after-before-gain+dt*sink;need(std::abs(ores)<1e-7,"O2ledger");
                double meanL1=0,maxmean=0;
                for(int j=0;j<lp.n;++j){meanL1+=std::abs(result.v[j]);maxmean=std::max(maxmean,std::abs(result.v[j]-reference.v[j]));}
                double L1=result.weighted_L1;
                // Buffer evidence until every certificate accepts this proposal.
                std::ostringstream flux_rows,ledger_rows,profile_rows,objective_rows,agent_rows,group_rows,equivalence_rows;
                for(auto stream:{&flux_rows,&ledger_rows,&profile_rows,&objective_rows,&agent_rows,&group_rows,&equivalence_rows})*stream<<std::setprecision(17);
                equivalence_rows<<step<<','<<reference.primary<<','<<result.primary<<','<<reference.extent<<','<<result.extent<<','<<maxmean<<','<<L1<<','<<meanL1<<'\n';
                for(size_t g=0;g<group_counts.size();++g)for(int j=0;j<lp.n;++j)
                    group_rows<<step<<','<<g<<','<<group_counts[g]<<','<<bare(lp.ids[j])<<','<<result.group_v[g][j]<<'\n';
                double proposed_maxres=maxres;
                for(int j=0;j<lp.n;++j)flux_rows<<step<<','<<bare(lp.ids[j])<<','<<result.v[j]<<'\n';
                for(size_t f=0;f<stocks.size();++f){if(f<fields.size()&&fields[f].name=="oxygen")continue;
                    double res=next[f]-stocks[f]-delta[f]-delivery[f];proposed_maxres=std::max(proposed_maxres,std::abs(res));need(std::abs(res)<1e-7,"Stockledger");
                    ledger_rows<<step<<','<<(f==fields.size()?"dipeptide":fields[f].name)<<','<<stocks[f]<<','<<next[f]<<','<<delta[f]<<','<<delivery[f]<<','<<res<<'\n';}
                for(int i=0;i<n;++i)profile_rows<<step<<','<<i<<','<<end[i]<<'\n';
                objective_rows<<step<<','<<result.primary<<','<<L1<<','<<e<<','<<cap<<','<<bottom<<','<<gain<<','<<ores<<'\n';
                std::vector<std::vector<double>> previous_fluxes;std::vector<std::string> previous_status;
                for(auto cell:agents){auto intra=static_cast<dFBAIntracellular*>(cell->phenotype.intracellular);
                    std::vector<double> values;for(auto rx:intra->sbml_model.reactions)values.push_back(rx->getFluxValue());
                    previous_fluxes.push_back(values);previous_status.push_back(intra->last_solution_status);}
                try{
                int agent_index=0;
                for(auto cell:agents){auto intra=static_cast<dFBAIntracellular*>(cell->phenotype.intracellular);intra->last_solution_status="last_committed_optimal";
                    const auto& selected=result.group_v[agent_group[agent_index]];
                    for(int j=0;j<lp.n;++j){need(intra->sbml_model.reactions[j]->getId()==lp.ids[j],"Agentreactionorder");intra->sbml_model.reactions[j]->setFluxValue(selected[j]);}
                    std::vector<double> observed(lp.n);double difference=0,bound_error=0,sv=0;
                    for(int j=0;j<lp.n;++j){auto rx=intra->sbml_model.reactions[j];observed[j]=rx->getFluxValue();
                        need(std::isfinite(observed[j]),"Agentnonfiniteflux");difference=std::max(difference,std::abs(observed[j]-selected[j]));
                        bound_error=std::max(bound_error,std::max(rx->getLowerBound()-observed[j],observed[j]-rx->getUpperBound()));}
                    for(auto& row:lp.sparse_stoichiometry){long double balance=0;for(auto term:row)balance+=term.second*observed[term.first];sv=std::max(sv,double(std::abs(balance)));}
                    need(difference<=1e-12&&sv<2e-7&&bound_error<2e-7,"Agentreadbackcertificate");
                    agent_rows<<step<<','<<cell->ID<<','<<difference<<','<<sv<<','<<bound_error<<','<<observed[atp_col]<<','<<observed[lp.oi]<<','<<observed[urea_col]<<','<<observed[serine_col]<<'\n';
                    need(!(agent_index++==0&&step==force_agent_failure_step),"ForcedAgentCertificateFailure");}
                }catch(...){
                    // Restore uncommitted model state; do not repair the proposal.
                    for(size_t k=0;k<agents.size();++k){auto intra=static_cast<dFBAIntracellular*>(agents[k]->phenotype.intracellular);
                        for(int j=0;j<lp.n;++j)intra->sbml_model.reactions[j]->setFluxValue(previous_fluxes[k][j]);
                        intra->last_solution_status=previous_status[k];}
                    throw;
                }
                stocks=next;used+=e;oxygen=end;maxres=proposed_maxres;expose();++completed;
                flux<<flux_rows.str();ledger<<ledger_rows.str();profiles<<profile_rows.str();
                objectives<<objective_rows.str();agent_audit<<agent_rows.str();
                group_fluxes<<group_rows.str();equivalence<<equivalence_rows.str();
            }catch(const std::exception& error){status="numerically_inadmissible";detail=error.what();break;}
        }
        std::ofstream state(out+"/final_stocks.csv"),agentfile(out+"/agents.csv");state<<std::setprecision(17)<<"field,stock_pmol\n";
        for(size_t f=0;f<stocks.size();++f)if(f==fields.size()||fields[f].name!="oxygen")state<<(f==fields.size()?"dipeptide":fields[f].name)<<','<<stocks[f]<<'\n';
        agentfile<<"cell_id,type,voxel,weight,intracellular_status\n";
        for(auto cell:agents){auto intra=static_cast<dFBAIntracellular*>(cell->phenotype.intracellular);agentfile<<cell->ID<<",hepatocyte,"<<cell->get_current_voxel_index()<<",1,"<<intra->last_solution_status<<'\n';}
        serialize_agent_fluxes("agent_fluxes_final.csv",completed-1);
        std::ofstream summary(out+"/summary.json");summary<<std::setprecision(17)<<"{\"status\":\""<<status<<"\",\"detail\":\""<<detail
            <<"\",\"completed_steps\":"<<completed<<",\"completed_min\":"<<completed*dt<<",\"agents\":22,\"independent_models\":"<<models.size()
            <<",\"independent_reaction_states\":"<<reaction_states.size()<<",\"groups\":"<<group_counts.size()<<",\"joint_group_selection\":true,\"cumulative_extent_pmol\":"<<used<<",\"maximum_stock_residual_pmol\":"<<maxres<<",\"custom_FV_transport\":true,\"biological_validation\":false}\n";
        std::cout<<status<<" completed_min="<<completed*dt<<" agents="<<agents.size()<<std::endl;return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 2;}
}
