#define main archived_metabolic_entry_not_called
#include "support/legacy_lp.cpp"
#undef main
#define main archived_acellular_entry_not_called
#include "support/oxygen_operator.cpp"
#undef main
#include <sys/wait.h>
#include <unistd.h>
#include <memory>
#include <map>

struct ResourceField {Field f;int density;bool uptake;};
#include "support/joint.h"
using Matrix=std::vector<V>;
static const double seconds=.0075,amount_scale=.12*seconds/3600;
static std::string prefix,mode;
static bool child_phase=false,scheduler_phase=false;
static long long lp_calls=0,writer_calls=0,phenotype_calls=0,mechanics_calls=0,scheduler_entries=0;
static int commits=0;
static Matrix lower,upper,q0;
static std::vector<int> vox;
static std::vector<ResourceField> fields;
struct Host {int id,v;V xyz;double volume;std::string zone;};
static std::vector<Host> hosts;
static Operator* oxygen_operator=nullptr;
static State accepted_oxygen;
void reject(const std::string& reason);

struct CoordinatedHost:dFBAIntracellular {
    bool need_update() override{return false;}
    void update(Cell*,Phenotype&,double) override{reject("individual_metabolic_update");}
    void update() override{reject("individual_metabolic_solve");}
};
Matrix actual_fields(){Matrix q(2500,V(20));for(int v=0;v<2500;++v)for(int f=0;f<20;++f)q[v][f]=microenvironment(v)[f]*oxygen_operator->ve[v];return q;}
Matrix actual_flux(){Matrix a(3000,V(179));for(int c=0;c<3000;++c){auto m=static_cast<CoordinatedHost*>((*all_cells)[c]->phenotype.intracellular);for(int j=0;j<179;++j)a[c][j]=m->sbml_model.reactions[j]->getFluxValue();}return a;}
V flatten(const Matrix& a){V x;for(auto& row:a)x.insert(x.end(),row.begin(),row.end());return x;}
V observed(){V a=flatten(actual_fields()),v=flatten(actual_flux());a.insert(a.end(),v.begin(),v.end());
    for(auto c:*all_cells){a.push_back(c->ID);a.push_back(c->get_current_voxel_index());a.insert(a.end(),c->position.begin(),c->position.end());a.push_back(c->phenotype.volume.total);}
    a.push_back(PhysiCell_globals.current_time);a.push_back(commits);return a;}
void dump(const std::string& path,const V& x){std::ofstream f(path,std::ios::binary);raw(f,x);f.close();need(bool(f),"closed_evidence");}
void dump(const std::string& path,const Matrix& x){dump(path,flatten(x));}
void reject(const std::string& reason){
    if(child_phase){dump(prefix+"_rejected_child.f64",observed());std::ofstream f(prefix+"_rejection.json");
        f<<"{\"reason\":\""<<reason<<"\",\"LP_calls\":"<<lp_calls<<",\"commits\":"<<commits<<",\"scheduler_entries\":"<<scheduler_entries<<",\"writer_calls\":"<<writer_calls<<",\"phenotype_calls\":"<<phenotype_calls<<",\"mechanics_calls\":"<<mechanics_calls<<"}\n";f.close();std::_Exit(42);}
    throw std::runtime_error(reason);
}
void pheno(Cell* c,Phenotype&,double){++phenotype_calls;if(c->ID!=0)return;
    if(mode=="callback_field")microenvironment(0)[0]+=1e-4;
    if(mode=="callback_nan")microenvironment(0)[0]=std::numeric_limits<double>::quiet_NaN();
    if(mode=="callback_flux")static_cast<CoordinatedHost*>(c->phenotype.intracellular)->sbml_model.reactions[0]->setFluxValue(1e4);
    if(mode=="callback_clock")PhysiCell_globals.current_time=.5;
    if(mode=="callback_identity")c->position[0]+=.1;
}
void mech(Cell*,Phenotype&,double){++mechanics_calls;}
void zero(V& a){std::fill(a.begin(),a.end(),0.);}
extern "C" void __real__ZN6BioFVM11Basic_Agent29simulate_secretion_and_uptakeEPNS_16MicroenvironmentEd(Basic_Agent*,Microenvironment*,double);
extern "C" void __wrap__ZN6BioFVM11Basic_Agent29simulate_secretion_and_uptakeEPNS_16MicroenvironmentEd(Basic_Agent* a,Microenvironment* env,double dt){
    if(!child_phase||!scheduler_phase)reject("uptake_outside_staged_scheduler");
    for(auto r:{a->secretion_rates,a->uptake_rates,a->net_export_rates})for(double x:*r)if(x!=0)reject("nonzero_native_rate");
    int v=static_cast<Cell*>(a)->get_current_voxel_index();V before=env->density_vector(v);
    __real__ZN6BioFVM11Basic_Agent29simulate_secretion_and_uptakeEPNS_16MicroenvironmentEd(a,env,dt);
    if(before!=env->density_vector(v))reject("native_writer_drift");++writer_calls;
}
extern "C" void __wrap__ZN6BioFVM16Microenvironment24simulate_diffusion_decayEd(Microenvironment*,double){reject("standard_diffusion_forbidden");}
extern "C" void __wrap__Z32custom_microenvironment_functiond(double){reject("legacy_source_forbidden");}

struct Proposal {Matrix q,flux;State oxygen;};
void write_proposal(const std::string& path,const Proposal& p){std::ofstream f(path,std::ios::binary);
    raw(f,flatten(p.q));raw(f,flatten(p.flux));raw(f,p.oxygen.b);raw(f,p.oxygen.e);raw(f,p.oxygen.p);raw(f,V(p.oxygen.r.begin(),p.oxygen.r.end()));f.close();need(bool(f),"proposal_close");}
V read_values(std::istream& f,int n){V v(n);f.read(reinterpret_cast<char*>(v.data()),n*sizeof(double));need(bool(f),"input_read");return v;}
Proposal parse_proposal(std::istream& f){Proposal p;p.q.resize(2500);p.flux.resize(3000);
    for(auto& r:p.q)r=read_values(f,20);for(auto& r:p.flux)r=read_values(f,179);
    p.oxygen.b=read_values(f,48);p.oxygen.e=read_values(f,2500);p.oxygen.p=read_values(f,48);auto r=read_values(f,4);std::copy(r.begin(),r.end(),p.oxygen.r.begin());need(f.peek()==EOF,"trailing_proposal");return p;}
std::string file_bytes(const std::string& path){std::ifstream f(path,std::ios::binary);need(bool(f),"snapshot_open");return std::string((std::istreambuf_iterator<char>(f)),std::istreambuf_iterator<char>());}
void check_hosts(){
    need(all_cells->size()==3000,"host_count");for(int c=0;c<3000;++c){auto a=(*all_cells)[c];auto& h=hosts[c];
        if(a->ID!=h.id||a->get_current_voxel_index()!=h.v||a->position!=h.xyz||a->phenotype.volume.total!=h.volume||a->phenotype.death.dead||a->is_movable||a->is_out_of_domain)reject("host_identity_change");
        auto m=static_cast<CoordinatedHost*>(a->phenotype.intracellular);for(int j=0;j<179;++j)need(m->sbml_model.reactions[j]->getLowerBound()==lower[c][j]&&m->sbml_model.reactions[j]->getUpperBound()==upper[c][j],"host_bound_mutation");}
}
void certify(LP& lp,const Proposal& p,const State& transported){
    Matrix change(2500,V(20,0)),gross=change;long double C=0,N=0,od=0;
    for(int c=0;c<3000;++c){for(int j=0;j<179;++j){double x=p.flux[c][j];need(std::isfinite(x)&&x>=lower[c][j]-2e-7&&x<=upper[c][j]+2e-7,"host_flux_bounds");}
        for(auto& row:lp.sparse_stoichiometry){long double s=0;for(auto t:row)s+=t.second*p.flux[c][t.first];need(std::abs(s)<=2e-7,"host_Sv");}
        for(int f=0;f<20;++f){double d=amount_scale*p.flux[c][fields[f].f.column];change[vox[c]][f]+=d;gross[vox[c]][f]+=std::max(0.,-d);C+=d*fields[f].f.carbon;N+=d*fields[f].f.nitrogen;if(f==0)od+=d;}}
    need(std::abs(C)<=1e-9&&std::abs(N)<=1e-9,"element_ledger");
    for(int v=0;v<2500;++v)for(int f=0;f<20;++f){double available=f==0?transported.e[v]:q0[v][f],reserve=available>0?std::min(1e-12,available*.5):0;
        need(gross[v][f]<=available-reserve+1e-12,"gross_stock");need(std::isfinite(p.q[v][f])&&p.q[v][f]>=0,"negative_stock");need(std::abs(p.q[v][f]-available-change[v][f])<=1e-12,"field_update");}
    need(p.oxygen.b==transported.b&&p.oxygen.p==transported.p&&p.oxygen.r==transported.r,"blood_reservoir_changed");
    for(int v=0;v<2500;++v)need(p.oxygen.e[v]==p.q[v][0],"oxygen_field_consistency");
    double before=total(accepted_oxygen.b)+total(accepted_oxygen.e)+accepted_oxygen.r[0]+accepted_oxygen.r[2];
    double after=total(p.oxygen.b)+total(p.oxygen.e)+p.oxygen.r[0]+p.oxygen.r[2];need(std::abs(after-before-double(od))<1e-8,"oxygen_ledger");
    need(std::abs(p.oxygen.r[1]+p.oxygen.r[3]-accepted_oxygen.r[1]-accepted_oxygen.r[3])<1e-10,"fluid_ledger");
}
int main(int argc,char** argv){try{
    need(argc==8,"bundle prefix transport hosts model fields mode");prefix=argv[2];mode=argv[7];omp_set_num_threads(1);SeedRandom(0);
    Operator op(argv[3],seconds,1);oxygen_operator=&op;accepted_oxygen=op.initial(argv[3]);LP lp(argv[5]);
    pugi::xml_document doc;need(doc.load_file(argv[6]),"field_XML");auto cfg=doc.document_element().child("finite_resource_control");
    for(auto n:cfg.children("field")){ResourceField f;f.f.name=n.attribute("name").value();f.f.rid="R_"+std::string(n.attribute("reaction").value());f.f.column=lp.model.getReactionIndex(f.f.rid);f.f.carbon=n.attribute("C").as_double();f.f.nitrogen=n.attribute("N").as_double();f.density=fields.size();f.uptake=n.attribute("uptake").as_bool();need(f.f.column>=0,"field_map");fields.push_back(f);}need(fields.size()==20,"field_count");
    std::ifstream input(argv[1],std::ios::binary);q0.resize(2500);for(auto& row:q0)row=read_values(input,20);lower.resize(3000,V(179));upper=lower;
    for(int c=0;c<3000;++c)for(int j=0;j<179;++j){auto b=read_values(input,2);lower[c][j]=b[0];upper[c][j]=b[1];}need(input.peek()==EOF,"input_tail");
    auto& env=microenvironment;env.resize_densities(20);for(int f=0;f<20;++f)env.set_density(f,fields[f].f.name,"mM",0.,0.);env.resize_space(-500.,500.,-500.,500.,-10.,10.,20.,20.,20.);
    auto container=create_cell_container_for_microenvironment(env,30.);initialize_default_cell_definition();cell_defaults.name="coordinated_host";build_cell_definitions_maps();setup_signal_behavior_dictionaries();auto& def=cell_defaults;
    def.is_movable=false;def.functions.volume_update_function=nullptr;def.functions.update_phenotype=pheno;def.functions.custom_cell_rule=mech;def.functions.update_velocity=nullptr;def.functions.update_migration_bias=nullptr;def.functions.contact_function=nullptr;def.functions.set_orientation=nullptr;
    zero(def.phenotype.secretion.secretion_rates);zero(def.phenotype.secretion.uptake_rates);zero(def.phenotype.secretion.net_export_rates);zero(def.phenotype.death.rates);for(auto& row:def.phenotype.cycle.data.transition_rates)zero(row);
    zero(def.phenotype.cell_transformations.transformation_rates);zero(def.phenotype.cell_interactions.live_phagocytosis_rates);zero(def.phenotype.cell_interactions.attack_rates);zero(def.phenotype.cell_interactions.fusion_rates);
    def.phenotype.cell_interactions.apoptotic_phagocytosis_rate=0;def.phenotype.cell_interactions.necrotic_phagocytosis_rate=0;def.phenotype.cell_interactions.other_dead_phagocytosis_rate=0;
    PhysiCell_settings.rules_enabled=false;PhysiCell_settings.disable_automated_spring_adhesions=true;default_microenvironment_options.calculate_gradients=false;default_microenvironment_options.track_internalized_substrates_in_each_agent=false;
    std::ifstream host_input(argv[4]);std::set<const void*> models,reactions,metabolites;std::vector<std::vector<int>> members(2500);
    for(int c=0;c<3000;++c){Host h;h.xyz.resize(3);host_input>>h.id>>h.v>>h.xyz[0]>>h.xyz[1]>>h.xyz[2]>>h.volume>>h.zone;need(bool(host_input)&&h.id==c,"host_input");hosts.push_back(h);vox.push_back(h.v);members[h.v].push_back(c);
        auto a=create_cell(cell_defaults);a->set_total_volume(h.volume);need(a->assign_position(h.xyz),"position");auto m=new CoordinatedHost();m->experimental_culture=true;m->use_metabolic_death=false;m->sbml_filename=argv[5];m->sbml_model=lp.model;
        for(int j=0;j<179;++j){auto r=m->sbml_model.reactions[j];r->setLowerBound(lower[c][j]);r->setUpperBound(upper[c][j]);r->setFluxValue(0);reactions.insert(r);}for(auto r:m->sbml_model.metabolites)metabolites.insert(r);models.insert(m);a->phenotype.intracellular=m;}
    need(models.size()==3000&&reactions.size()==537000&&metabolites.size()==525000,"distinct_models");check_hosts();
    for(int v=0;v<2500;++v)for(int f=0;f<20;++f)env(v)[f]=q0[v][f]/op.ve[v];PhysiCell_globals.current_time=0;
    V parent_before=observed();if(mode=="run")dump(prefix+"_parent_before.f64",parent_before);
    std::ifstream status("/proc/self/status");std::string line;int threads=0;while(std::getline(status,line))if(line.substr(0,8)=="Threads:")threads=std::stoi(line.substr(8));need(threads==1,"fork_requires_single_thread");
    int ack[2],ready[2];need(pipe(ack)==0&&pipe(ready)==0,"pipes");pid_t child=fork();need(child>=0,"fork");
    if(child==0){child_phase=true;close(ack[1]);close(ready[0]);
        if(mode=="standard_diffusion")env.simulate_diffusion_decay(seconds/60);
        State source=accepted_oxygen;if(mode=="transport_mutation")source.e[0]=-1;
        Diagnostics d;State transported;try{transported=op.propose(source,true,65,d);}catch(const std::exception&){reject("invalid_transport_proposal");}
        Proposal p;p.oxygen=transported;p.q=q0;for(int v=0;v<2500;++v)p.q[v][0]=transported.e[v];p.flux.assign(3000,V(179));
        std::map<int,std::unique_ptr<JointTemplate>> templates;
        for(int v=0;v<2500;++v)if(!members[v].empty()){int k=members[v].size();if(!templates.count(k))templates[k].reset(new JointTemplate(lp,fields,k,seconds/60));auto a=templates[k]->solve(lp,members[v],lower,upper,p.q[v]);++lp_calls;if(!a.ok)reject("joint_LP_infeasible");for(int i=0;i<k;++i)p.flux[members[v][i]]=a.cells[i];}
        for(int c=0;c<3000;++c)for(int f=0;f<20;++f)p.q[vox[c]][f]+=amount_scale*p.flux[c][fields[f].f.column];for(int v=0;v<2500;++v)p.oxygen.e[v]=p.q[v][0];
        certify(lp,p,transported);if(mode=="late_reject")reject("injected_precommit_rejection");
        for(int c=0;c<3000;++c){auto m=static_cast<CoordinatedHost*>((*all_cells)[c]->phenotype.intracellular);for(int j=0;j<179;++j)m->sbml_model.reactions[j]->setFluxValue(p.flux[c][j]);}
        for(int v=0;v<2500;++v)for(int f=0;f<20;++f)env(v)[f]=p.q[v][f]/op.ve[v];if(mode=="nonzero_uptake")(*all_cells)[0]->phenotype.secretion.uptake_rates[0]=1;
        scheduler_phase=true;++scheduler_entries;container->update_all_cells(seconds/60,2.,.1,seconds/60);scheduler_phase=false;
        auto measured=actual_fields();for(int v=0;v<2500;++v)for(int f=0;f<20;++f){if(!std::isfinite(measured[v][f]))reject("scheduler_nonfinite_field");if(std::abs(measured[v][f]-p.q[v][f])>1e-12)reject("scheduler_field_write");}
        if(actual_flux()!=p.flux)reject("scheduler_flux_write");if(PhysiCell_globals.current_time!=0)reject("scheduler_clock_write");check_hosts();
        need(writer_calls==3000&&phenotype_calls==3000&&mechanics_calls==3000,"scheduler_counts");need(container->num_divisions_in_current_step==0&&container->num_deaths_in_current_step==0,"no_death_division");
        write_proposal(prefix+"_candidate.f64",p);char yes='R';need(write(ready[1],&yes,1)==1,"ready_write");need(read(ack[0],&yes,1)==1&&yes=='A',"commit_ACK");
        ++commits;PhysiCell_globals.current_time=seconds/60;if(mode=="commit_tamper"||mode=="candidate_replace")p.q[0][3]+=.001;if(mode=="candidate_replace")write_proposal(prefix+"_candidate.f64",p);write_proposal(prefix+"_committed.f64",p);dump(prefix+"_child_committed_world.f64",observed());
        std::ofstream report(prefix+"_child.json");report<<"{\"agents\":3000,\"individual_models\":3000,\"joint_LP_calls\":"<<lp_calls<<",\"scheduler_calls\":1,\"scheduler_entries\":"<<scheduler_entries<<",\"writer_calls\":"<<writer_calls<<",\"phenotype_calls\":"<<phenotype_calls<<",\"mechanics_calls\":"<<mechanics_calls<<",\"commits\":1,\"clock_min\":"<<std::setprecision(17)<<PhysiCell_globals.current_time<<"}\n";report.close();std::_Exit(0);
    }
    close(ack[0]);close(ready[1]);char flag=0;ssize_t count=read(ready[0],&flag,1);bool candidate=count==1&&flag=='R';
    const std::string certified_payload=candidate?file_bytes(prefix+"_candidate.f64"):std::string();
    if(candidate){std::istringstream captured(certified_payload,std::ios::in|std::ios::binary);auto p=parse_proposal(captured);Diagnostics d;auto transported=op.propose(accepted_oxygen,true,65,d);certify(lp,p,transported);
        std::ofstream saved(prefix+"_certified.f64",std::ios::binary);saved.write(certified_payload.data(),certified_payload.size());saved.close();need(bool(saved),"certified_capture_close");char yes='A';need(write(ack[1],&yes,1)==1,"ACK_write");}
    close(ack[1]);close(ready[0]);int exit_status;need(waitpid(child,&exit_status,0)==child,"wait_child");bool success=WIFEXITED(exit_status)&&WEXITSTATUS(exit_status)==0;
    need(observed()==parent_before&&PhysiCell_globals.current_time==0&&commits==0,"parent_canonical_world_changed");check_hosts();
    need(success==candidate,"worker_protocol_failure");bool approved=success;std::string final_rejection;
    if(success){
        const std::string final_payload=file_bytes(prefix+"_committed.f64");
        if(certified_payload!=final_payload){approved=false;final_rejection="certified_payload_mismatch";}
        if(approved){std::istringstream final_stream(final_payload,std::ios::in|std::ios::binary);auto p=parse_proposal(final_stream);Diagnostics d;certify(lp,p,op.propose(accepted_oxygen,true,65,d));
            std::ifstream wf(prefix+"_child_committed_world.f64",std::ios::binary);V world=read_values(wf,605002);need(wf.peek()==EOF,"world_tail");
            for(double x:world)if(!std::isfinite(x)){approved=false;final_rejection="nonfinite_committed_world";break;}
            V q=flatten(p.q),f=flatten(p.flux);
            if(approved){for(int i=0;i<50000;++i)if(std::abs(world[i]-q[i])>1e-12){approved=false;final_rejection="committed_world_fields";break;}}
            if(approved){for(int i=0;i<537000;++i)if(world[50000+i]!=f[i]){approved=false;final_rejection="committed_world_flux";break;}}
            if(approved){for(int i=587000;i<605000;++i)if(world[i]!=parent_before[i]){approved=false;final_rejection="committed_world_identity";break;}}
            if(approved&&(world[605000]!=seconds/60||world[605001]!=1)){approved=false;final_rejection="committed_world_clock";}
        }
    }
    std::ofstream report(prefix+"_coordinator.json");report<<"{\"accepted_snapshot\":"<<(approved?"true":"false")<<",\"canonical_parent_unchanged\":true,\"parent_world_promoted\":false,\"persistent_resume\":false,\"parent_agents\":3000,\"parent_threads\":"<<threads<<",\"child_exit\":"<<(WIFEXITED(exit_status)?WEXITSTATUS(exit_status):-1)<<",\"rejection_reason\":\""<<final_rejection<<"\"}\n";report.close();return approved?0:42;
}catch(const std::exception& e){std::cerr<<"NATIVE_WORKER_FAILED "<<e.what()<<std::endl;return 1;}}
