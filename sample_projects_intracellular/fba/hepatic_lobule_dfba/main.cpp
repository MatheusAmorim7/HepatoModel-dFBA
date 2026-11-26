#define main archived_types_entry_not_called
#include "src/zonation_2d_native_groups/native.cpp"
#undef main
#include "custom_modules/custom.h"
#include <array>
#include <chrono>
#include <cstdlib>
#include <sys/resource.h>
struct ResourceField {Field f;int density;bool uptake;};
#include "src/lobule_cell_capacity/joint.h"
#include "src/lobule_native_perfusion/perfusion.h"
#include "src/vascular.h"
using Matrix=std::vector<std::vector<double>>;
#include "src/cadence/cache.h"
#include "src/cadence/sink.h"
struct PhaseClock {
    using Clock=std::chrono::steady_clock;
    Clock::time_point start=Clock::now(),last=start;
    std::string phase="config";
    std::map<std::string,double> seconds;
    void mark(const std::string& next){auto now=Clock::now();seconds[phase]+=std::chrono::duration<double>(now-last).count();last=now;phase=next;}
    void finish(const std::string& path){mark("unmeasured_profile_write");std::ofstream f(path);f<<std::setprecision(17)<<"{\"profiled_wall_s\":"<<std::chrono::duration<double>(last-start).count()<<",\"phases\":{";bool first=true;for(auto p:seconds){if(!first)f<<',';first=false;f<<'"'<<p.first<<"\":"<<p.second;}f<<"}}\n";}
};
std::string output;
long long passive_checks=0,unexpected_updates=0,writers=0,phenotypes=0,mechanics=0;
double writer_drift=0;
double diff(const std::vector<double>& a,const std::vector<double>& b){need(a.size()==b.size(),"Shape");double d=0;for(size_t i=0;i<a.size();++i)d=std::max(d,std::abs(a[i]-b[i]));return d;}
void raw(std::ostream& s,const std::vector<double>& a){s.write(reinterpret_cast<const char*>(a.data()),a.size()*sizeof(double));need(bool(s),"Evidence_write");}
void raw(std::ostream& s,const Matrix& a){for(auto& row:a)raw(s,row);}
std::vector<double> values(const char* text,int n){std::istringstream s(text);std::vector<double>a(n);for(double& x:a)need(bool(s>>x)&&std::isfinite(x),"Input");double tail;need(!(s>>tail),"Extra_input");return a;}
void forbidden(const char* reason){std::ofstream f(output+"/unexpected.json");f<<"{\"reason\":\""<<reason<<"\"}\n";f.close();std::cerr<<"UNEXPECTED "<<reason<<std::endl;std::exit(43);}

// The global coordinator owns selection/debits. Each cell retains its own SBML objects.
struct CoordinatedHost:dFBAIntracellular {
    bool need_update() override {++passive_checks;return false;}
    void update(Cell*,Phenotype&,double) override {++unexpected_updates;forbidden("individual_update");}
    void update() override {++unexpected_updates;forbidden("individual_optimize");}
};
extern "C" void __real__ZN6BioFVM11Basic_Agent29simulate_secretion_and_uptakeEPNS_16MicroenvironmentEd(Basic_Agent*,Microenvironment*,double);
extern "C" void __wrap__ZN6BioFVM11Basic_Agent29simulate_secretion_and_uptakeEPNS_16MicroenvironmentEd(Basic_Agent* a,Microenvironment* m,double dt){
    ++writers;for(auto p:{a->secretion_rates,a->uptake_rates,a->net_export_rates})for(double x:*p)if(x!=0)forbidden("second_native_debit");
    int voxel=static_cast<Cell*>(a)->get_current_voxel_index();auto before=(*m)(voxel);
    __real__ZN6BioFVM11Basic_Agent29simulate_secretion_and_uptakeEPNS_16MicroenvironmentEd(a,m,dt);
    writer_drift=std::max(writer_drift,diff(before,(*m)(voxel))*.008);if(writer_drift>1e-12)forbidden("native_writer_drift");
}
extern "C" void __wrap__ZN6BioFVM16Microenvironment24simulate_diffusion_decayEd(Microenvironment*,double){forbidden("duplicate_diffusion");}
extern "C" void __wrap__Z32custom_microenvironment_functiond(double){forbidden("legacy_source");}
std::vector<int> pc(3000,0),mc(3000,0);
void phenotype_marker(Cell* c,Phenotype&,double){++phenotypes;++pc.at(c->ID);}
void mechanics_marker(Cell* c,Phenotype&,double){++mechanics;++mc.at(c->ID);}
void zero(std::vector<double>& v){std::fill(v.begin(),v.end(),0.);}

int main(int argc,char** argv){try{
    if(argc==4&&std::string(argv[1])=="sink_fault")return sink_fault_test(argv[2],argv[3]);
    PhaseClock timing;
    need(argc==7&&load_PhysiCell_config_file(argv[1]),"Isolated_XML");omp_set_num_threads(1);
    std::string format=argv[6];need(format=="raw"||format=="gzip","Output_format");bool zip=format=="gzip";std::string suffix=zip?".gz":"";
    std::string label=argv[2];output=argv[3];bool reject=std::string(argv[4])=="reject_after_one";
    std::string cache_mode=argv[5];need(cache_mode=="cold"||cache_mode=="warm","Cache_mode");bool use_cache=cache_mode=="warm";
    auto cfg=physicell_config_root.child("finite_resource_control");need(cfg.attribute("metabolic_scheduler_control").as_bool(),"Not_root_dropin");
    auto item=cfg.find_child_by_attribute("case","label",label.c_str());need(item,"Case");
    double dt=item.attribute("dt").as_double(),F=item.attribute("F").as_double();int steps=item.attribute("steps").as_int();
    need(dt>0&&dt<=1&&steps>=1&&steps<=2000,"Bounded_cases");double DTP=item.attribute("dt_phenotype").as_double(2.),DTM=item.attribute("dt_mechanics").as_double(.1);const int RP=std::lround(DTP/dt),RM=std::lround(DTM/dt);need(RP>=1&&RM>=1&&std::abs(RP*dt-DTP)<1e-9&&std::abs(RM*dt-DTM)<1e-9,"Cadence_integers");
    timing.mark("model_template");std::string model=cfg.attribute("model").value();LP lp(model);long long setup_initial_solve_calls=1;GLUCOSE_REWARD=item.attribute("glucose_reward").as_double(0);if(GLUCOSE_REWARD!=0)GLUCOSE_COL=lp.model.getReactionIndex("R_EX_HC00040_s");
    timing.mark("cell_setup");setup_microenvironment();auto container=create_cell_container_for_microenvironment(microenvironment,30.);create_cell_types();load_cells_from_pugixml();
    need(all_cells->size()==3000&&microenvironment.number_of_voxels()==2500,"Population_mesh");
    PhysiCell_settings.rules_enabled=false;PhysiCell_settings.disable_automated_spring_adhesions=true;default_microenvironment_options.calculate_gradients=false;
    std::vector<ResourceField> fields;
    for(auto node:cfg.children("field")){ResourceField f;f.f.name=node.attribute("name").value();f.f.rid="R_"+std::string(node.attribute("reaction").value());f.f.column=lp.model.getReactionIndex(f.f.rid);f.density=microenvironment.find_density_index(f.f.name);f.uptake=node.attribute("uptake").as_bool();f.f.carbon=node.attribute("C").as_double();f.f.nitrogen=node.attribute("N").as_double();need(f.f.column>=0&&f.density>=0,"Mapping");fields.push_back(f);}
    need(fields.size()==NF&&fields[0].f.column==lp.oi,"Fields");
    std::vector<int> vox(3000);std::vector<std::vector<int>> members(2500);std::vector<std::array<double,4>> identity;
    std::set<const void*> models,reactions,metabolites;
    Matrix lower(3000,lp.lower),upper(3000,lp.upper);
    std::ofstream hosts(output+"/hosts.csv");hosts<<"cell_id,voxel,x,y,z\n"<<std::setprecision(17);
    timing.mark("host_installation_and_checks");
    for(int c=0;c<3000;++c){auto cell=(*all_cells)[c];auto& p=cell->phenotype;need(cell->ID==c&&p.intracellular==nullptr&&!p.death.dead&&!cell->is_out_of_domain,"Host_identity");
        vox[c]=cell->get_current_voxel_index();members[vox[c]].push_back(c);identity.push_back({cell->position[0],cell->position[1],cell->position[2],p.volume.total});
        hosts<<c<<','<<vox[c]<<','<<cell->position[0]<<','<<cell->position[1]<<','<<cell->position[2]<<'\n';
        cell->is_movable=false;cell->functions.update_phenotype=phenotype_marker;cell->functions.custom_cell_rule=mechanics_marker;cell->functions.volume_update_function=nullptr;cell->functions.contact_function=nullptr;cell->functions.update_velocity=nullptr;cell->functions.update_migration_bias=nullptr;
        zero(p.secretion.secretion_rates);zero(p.secretion.uptake_rates);zero(p.secretion.net_export_rates);zero(p.death.rates);for(auto& row:p.cycle.data.transition_rates)zero(row);
        zero(p.cell_transformations.transformation_rates);zero(p.cell_interactions.live_phagocytosis_rates);zero(p.cell_interactions.attack_rates);zero(p.cell_interactions.fusion_rates);p.cell_interactions.apoptotic_phagocytosis_rate=0;p.cell_interactions.necrotic_phagocytosis_rate=0;p.cell_interactions.other_dead_phagocytosis_rate=0;
        auto m=new CoordinatedHost();m->experimental_culture=true;m->use_metabolic_death=false;m->sbml_filename=model;
        m->sbml_model=lp.model;++setup_initial_solve_calls;
        need(m->sbml_model.isInitialized()&&m->sbml_model.reactions.size()==NR&&m->sbml_model.metabolites.size()==NM,"Individual_model_shape");
        for(int j=0;j<NR;++j){auto a=m->sbml_model.reactions[j],b=lp.model.reactions[j];need(a!=b&&a->getId()==b->getId()&&a->getMetabolites()==b->getMetabolites()&&a->getLowerBound()==b->getLowerBound()&&a->getUpperBound()==b->getUpperBound(),"Deep_model_identity");reactions.insert(a);}
        for(auto a:m->sbml_model.metabolites)metabolites.insert(a);models.insert(m);p.intracellular=m;
        if(c%500==499)std::cout<<"initialized "<<c+1<<" individual hosts"<<std::endl;
    }hosts.close();need(models.size()==3000&&reactions.size()==3000*(size_t)NR&&metabolites.size()==3000*(size_t)NM,"No_shared_model_objects");
    timing.mark("case_setup");int count=0;for(auto row:cfg.child("capacity_maps").children("profile")){int c=row.attribute("cell_id").as_int();need(c==count++,"Profile_order");for(std::string rid:{"r0054","r0353","r0396","r0123","r0487"}){int j=lp.model.getReactionIndex("R_"+rid);double factor=row.attribute(rid.c_str()).as_double();need(factor>0,"Factor");if(rid=="r0054")lower[c][j]*=factor;else upper[c][j]*=factor;}double gub=row.attribute("glucose_ub").as_double(-1),gcap=row.attribute("glycogen_cap").as_double(1.0);if(gub>=0){int jg2=lp.model.getReactionIndex("R_EX_HC00040_s");lower[c][jg2]=0;upper[c][jg2]=gub;}if(gcap!=1.0){upper[c][lp.model.getReactionIndex("R_GLYCOGEN_SUPPLY")]*=gcap;}double tg=row.attribute("task_glucose").as_double(1.0);need(tg>0,"Task_factor");if(tg!=1.0){int jg=lp.model.getReactionIndex("R_EX_HC00040_s"),jh=lp.model.getReactionIndex("R_GLYCOGEN_SUPPLY");lower[c][jg]*=tg;upper[c][jg]*=tg;upper[c][jh]*=tg;}}
    need(count==3000,"Profiles");for(int c=0;c<3000;++c)for(auto f:fields)if(!f.uptake)lower[c][f.f.column]=std::max(lower[c][f.f.column],0.);
    std::ofstream bounds(output+"/bounds.f64",std::ios::binary);for(int c=0;c<3000;++c)for(int j=0;j<NR;++j)raw(bounds,std::vector<double>{lower[c][j],upper[c][j]});bounds.close();
    Matrix initial(2500),initialflux(3000);count=0;for(auto row:cfg.child("initial_stocks").children("voxel")){need(row.attribute("id").as_int()==count,"Stock_order");initial[count++]=values(row.child_value(),NF);}need(count==2500,"Stocks");
    count=0;for(auto row:cfg.child("frozen_inputs").children("cell")){need(row.attribute("id").as_int()==count,"Flux_order");initialflux[count++]=values(row.child_value(),NR);}need(count==3000,"Fluxes");
auto q=initial;auto vc=cfg.child("vascular");need(vc,"Vascular_config");const double ge=vc.attribute("ge").as_double(),gw=item.attribute("gw").as_double(vc.attribute("gw").as_double()),plasma=item.attribute("plasma_flow_nL_min").as_double(vc.attribute("plasma_flow_nL_min").as_double()),Vb=item.attribute("blood_volume_nL").as_double(vc.attribute("blood_volume_nL").as_double());auto nm_d=values(vc.child("node_map").child_value(),2500);std::vector<int> node_map(nm_d.begin(),nm_d.end());std::vector<double> Cin(NF,0.);for(auto in:vc.children("inlet")){bool found=false;for(int f=0;f<NF;++f)if(fields[f].f.name==in.attribute("field").value()){Cin[f]=in.attribute("mM").as_double();found=true;}need(found,"Inlet_field");}std::vector<std::vector<double>> vessel(NF,std::vector<double>(SerialTransport::K,0.)),vessel_prop;for(int f=0;f<NF;++f)vessel[f].assign(SerialTransport::K,Cin[f]);const double Vnode=Vb/SerialTransport::K;
    auto actual_fields=[&](){Matrix a(2500,std::vector<double>(NF));for(int v=0;v<2500;++v)for(int f=0;f<NF;++f)a[v][f]=microenvironment(v)[fields[f].density]*.008;return a;};
    auto actual_flux=[&](){Matrix a(3000,std::vector<double>(NR));for(int c=0;c<3000;++c){auto m=static_cast<CoordinatedHost*>((*all_cells)[c]->phenotype.intracellular);for(int j=0;j<NR;++j)a[c][j]=m->sbml_model.reactions[j]->getFluxValue();}return a;};
    for(int c=0;c<3000;++c){auto m=static_cast<CoordinatedHost*>((*all_cells)[c]->phenotype.intracellular);for(int j=0;j<NR;++j){m->sbml_model.reactions[j]->setLowerBound(lower[c][j]);m->sbml_model.reactions[j]->setUpperBound(upper[c][j]);m->sbml_model.reactions[j]->setFluxValue(initialflux[c][j]);}}
    for(int v=0;v<2500;++v)for(int f=0;f<NF;++f)microenvironment(v)[fields[f].density]=q[v][f]/.008;
    CacheMap templates;int occupied=0;for(auto& group:members)if(!group.empty()){++occupied;if(!templates.count(group.size()))templates.emplace(group.size(),std::unique_ptr<MeasuredJoint>(new MeasuredJoint(lp,fields,group.size(),dt)));}
    need(occupied==1884,"Occupied_voxels");SerialTransport transport(dt,ge,gw,plasma/SerialTransport::SECT,Vnode,node_map);const double beta=item.attribute("o2_buffer").as_double(vc.attribute("o2_buffer").as_double(0.));SerialTransport transport_o2(dt,ge,gw,plasma/SerialTransport::SECT*(1+beta),Vnode*(1+beta),node_map);const int LBP=item.attribute("local_pool").as_int(0);std::vector<int> mnode(SerialTransport::K,0);for(int v=0;v<2500;++v)++mnode[node_map[v]];double lbp_roundoff=0,lbp_max_draw=0;
    BinarySink states(output+"/states.f64"+suffix,zip),oxy(output+"/vessel.f64"+suffix,zip),fluxes(output+"/fluxes.f64"+suffix,zip);std::ofstream ledger(output+"/ledger.csv");
    ledger<<"step,field,incoming_pmol,outgoing_pmol,metabolic_delta_pmol,ledger_error_pmol,transport_residual_pmol,LP_calls,scheduler_calls,passive_checks,writer_calls,scheduler_field_drift_pmol\n"<<std::setprecision(17);
    Matrix terminal_available;std::ofstream cadence(output+"/cadence.csv");cadence<<"step,phenotype_calls,mechanics_calls\n";
    long long calls=0,accepted=0,source_calls=0,scheduler_calls=0,clock_ticks=0;double max_field_drift=0;PhysiCell_globals.current_time=0;
    for(int step=0;step<steps;++step){timing.mark("snapshot");double dt_time=step*dt;need(PhysiCell_globals.current_time==dt_time,"Clock");std::vector<double> inc(NF,0.),outg(NF,0.),res(NF,0.),delta(NF,0.),errf(NF,0.);
        auto accepted_fields=actual_fields();auto accepted_flux=actual_flux();auto accepted_vessel=vessel;
        timing.mark("transport");auto available=q;vessel_prop=vessel;{std::vector<double> c_old(2500),c_new,b_new;for(int f=0;f<NF;++f){for(int v=0;v<2500;++v)c_old[v]=q[v][f]/.008;(f==0?transport_o2:transport).advance(c_old,vessel[f],Cin[f],c_new,b_new,inc[f],outg[f],res[f]);for(int v=0;v<2500;++v)available[v][f]=c_new[v]*.008;vessel_prop[f]=b_new;}}++source_calls;Matrix avail_lp=available;if(LBP)for(int v=0;v<2500;++v)for(int f=0;f<NF;++f)avail_lp[v][f]+=(f==0?Vnode*(1+beta):Vnode)*vessel_prop[f][node_map[v]]/mnode[node_map[v]];
        timing.mark("joint_lp_including_assembly_and_certificates");Matrix selected(3000,std::vector<double>(NR,0));bool ok=true;int failed=-1;
        for(int v=0;v<2500;++v){if(members[v].empty())continue;auto& problem=*templates.at(members[v].size());auto answer=use_cache?cached_select(problem,lp,v,members[v],lower,upper,avail_lp[v]):problem.solve(lp,members[v],lower,upper,avail_lp[v]);++calls;if(!answer.ok){ok=false;failed=v;break;}for(size_t k=0;k<members[v].size();++k)selected[members[v][k]]=answer.cells[k];}
        timing.mark("stock_validation");if(!ok||(reject&&step==1)){
            auto current=actual_fields();auto currentflux=actual_flux();need(current==accepted_fields&&currentflux==accepted_flux&&vessel==accepted_vessel&&scheduler_calls==accepted&&clock_ticks==accepted&&PhysiCell_globals.current_time==dt_time,"Atomic_precommit_rejection");
            std::ofstream a(output+"/rejected_fields.f64",std::ios::binary),b(output+"/rejected_fluxes.f64",std::ios::binary),d(output+"/rejected_available.f64",std::ios::binary);raw(a,current);raw(b,currentflux);raw(d,available);std::ofstream ox(output+"/rejected_vessel.f64",std::ios::binary);raw(ox,vessel);
            states.finish();oxy.finish();fluxes.finish();
            std::ofstream f(output+"/rejection.json");f<<"{\"kind\":\""<<(ok?"injected_coordinator_rejection":"local_LP_infeasible")<<"\",\"LP_calls\":"<<calls<<",\"failed_voxel\":"<<failed<<",\"commits\":"<<accepted<<",\"scheduler_calls\":"<<scheduler_calls<<",\"clock_ticks\":"<<clock_ticks<<",\"source_proposals\":"<<source_calls<<",\"setup_initialSolve_calls\":"<<setup_initial_solve_calls<<"}\n";
            cache_report(templates,use_cache,calls,output+"/cache.json");
            std::cout<<"Rejected entire uncommitted transaction; no repair"<<std::endl;return 42;
        }
        if(steps==201&&step==steps-1)terminal_available=available;
        auto next=available;Matrix gross(2500,std::vector<double>(NF,0));
        for(int c=0;c<3000;++c)for(int f=0;f<NF;++f){double d=selected[c][fields[f].f.column]*.12*dt/60;next[vox[c]][f]+=d;gross[vox[c]][f]+=std::max(0.,-d);delta[f]+=d;}
        std::vector<std::vector<double>> node_deficit(NF,std::vector<double>(SerialTransport::K,0.));
        for(int v=0;v<2500;++v)for(int f=0;f<NF;++f){double reserve=avail_lp[v][f]>0?std::min(1e-12,avail_lp[v][f]*.5):0;need(gross[v][f]<=avail_lp[v][f]-reserve+1e-12,"Gross_stock_no_repair");if(next[v][f]<0){need(LBP==1,"Negative_tissue_without_local_pool");node_deficit[f][node_map[v]]+=-next[v][f];next[v][f]=0;}}
        for(int f=0;f<NF;++f)for(int k=0;k<SerialTransport::K;++k){double Veff=(f==0?Vnode*(1+beta):Vnode);double amount=Veff*vessel_prop[f][k]-node_deficit[f][k];lbp_max_draw=std::max(lbp_max_draw,node_deficit[f][k]);if(amount<0){need(amount>-1e-12,"Negative_node_after_local_pool");lbp_roundoff=std::max(lbp_roundoff,-amount);amount=0;}vessel_prop[f][k]=amount/Veff;}
        for(int f=0;f<NF;++f){long double before=0,after=0;for(int v=0;v<2500;++v){before+=q[v][f];after+=next[v][f];}for(int k=0;k<SerialTransport::K;++k){before+=(f==0?Vnode*(1+beta):Vnode)*vessel[f][k];after+=(f==0?Vnode*(1+beta):Vnode)*vessel_prop[f][k];}errf[f]=double(after-before)-inc[f]+outg[f]-delta[f];need(std::abs(errf[f])<1e-9,"Global_ledger_per_field");}

        timing.mark("commit");
        // All groups and stocks passed: now commit fluxes and fields exactly once.
        for(int c=0;c<3000;++c){auto m=static_cast<CoordinatedHost*>((*all_cells)[c]->phenotype.intracellular);for(int j=0;j<NR;++j)m->sbml_model.reactions[j]->setFluxValue(selected[c][j]);}
        for(int v=0;v<2500;++v)for(int f=0;f<NF;++f)microenvironment(v)[fields[f].density]=next[v][f]/.008;
        auto snapshot=actual_fields();timing.mark("scheduler");++scheduler_calls;container->update_all_cells(dt_time,DTP,DTM,dt);timing.mark("post_scheduler_validation");auto measured=actual_fields();
        need(phenotypes==3000LL*(1+step/RP)&&mechanics==3000LL*(1+step/RM),"Horizon_aware_counts");
        for(int v=0;v<2500;++v)max_field_drift=std::max(max_field_drift,diff(snapshot[v],measured[v]));need(max_field_drift<=1e-12,"Scheduler_field_invariance");
        auto applied=actual_flux();need(applied==selected,"Scheduler_flux_invariance");
        need(all_cells->size()==3000&&container->num_divisions_in_current_step==0&&container->num_deaths_in_current_step==0,"No_host_loss");
        for(int c=0;c<3000;++c){auto cell=(*all_cells)[c];auto a=identity[c];need(pc[c]==1+step/RP&&mc[c]==1+step/RM,"Individual_callback_counts");need(cell->ID==c&&cell->position[0]==a[0]&&cell->position[1]==a[1]&&cell->position[2]==a[2]&&cell->get_current_voxel_index()==vox[c]&&cell->phenotype.volume.total==a[3]&&!cell->phenotype.death.dead&&models.count(cell->phenotype.intracellular)==1,"Host_model_invariance");auto m=static_cast<CoordinatedHost*>(cell->phenotype.intracellular);for(int j=0;j<NR;++j)need(m->sbml_model.reactions[j]->getLowerBound()==lower[c][j]&&m->sbml_model.reactions[j]->getUpperBound()==upper[c][j],"Scheduler_bound_invariance");}
        timing.mark("binary_output");raw(states,q);raw(states,available);raw(states,next);raw(oxy,vessel_prop);raw(fluxes,applied);
        timing.mark("bookkeeping_and_logs");q.swap(next);vessel.swap(vessel_prop);++accepted;PhysiCell_globals.current_time=(step+1)*dt;++clock_ticks;
        for(int f=0;f<NF;++f)ledger<<step<<','<<fields[f].f.name<<','<<inc[f]<<','<<outg[f]<<','<<delta[f]<<','<<errf[f]<<','<<res[f]<<','<<calls<<','<<scheduler_calls<<','<<passive_checks<<','<<writers<<','<<max_field_drift<<'\n';ledger.flush();
        cadence<<step<<','<<phenotypes<<','<<mechanics<<'\n';cadence.flush();
        std::cout<<label<<" step="<<accepted<<" jointLP="<<calls<<std::endl;
    }
    timing.mark("terminal_cold_reference");
    if(steps==201){auto applied=actual_flux();Matrix ref(3000,std::vector<double>(NR));double err=0;int checked=0;
        for(int v=0;v<2500;++v){if(members[v].empty())continue;auto a=templates.at(members[v].size())->solve(lp,members[v],lower,upper,terminal_available[v]);need(a.ok,"Terminal_cold_feasible");++checked;double obj=0;
            for(size_t k=0;k<members[v].size();++k){ref[members[v][k]]=a.cells[k];for(double x:applied[members[v][k]])obj+=std::abs(x)/members[v].size();}err=std::max(err,std::abs(obj-a.primary));need(std::abs(obj-a.primary)<2e-6,"Terminal_objective");}
        BinarySink f(output+"/terminal_fluxes.f64.gz",true);raw(f,ref);f.finish();std::ofstream r(output+"/terminal_reference.json");r<<std::setprecision(17)<<"{\"LP_calls\":"<<checked<<",\"objective_difference\":"<<err<<"}\n";}
    timing.mark("final_summary");need(calls==1884*steps&&passive_checks==3000*(steps-1)&&writers==3000*steps&&unexpected_updates==0&&phenotypes==3000LL*(1+(steps-1)/RP)&&mechanics==3000LL*(1+(steps-1)/RM),"Scheduler_counts");
    struct rusage usage;need(getrusage(RUSAGE_SELF,&usage)==0,"RSS");
    timing.mark("sink_close");states.finish();oxy.finish();fluxes.finish();timing.mark("completion_summary");{std::ofstream lb(output+"/lbp.json");lb<<std::setprecision(17)<<"{\"local_pool\":"<<LBP<<",\"max_node_draw_pmol\":"<<lbp_max_draw<<",\"max_roundoff_reset_pmol\":"<<lbp_roundoff<<"}\n";}
    std::ofstream summary(output+"/summary.json");summary<<std::setprecision(17)<<"{\"status\":\"completed\",\"accepted_steps\":"<<accepted<<",\"joint_LP_calls\":"<<calls<<",\"setup_initialSolve_calls\":"<<setup_initial_solve_calls<<",\"agents\":3000,\"individual_models\":3000,\"unique_reactions\":"<<3000*NR<<",\"unique_metabolites\":"<<3000*NM<<",\"scheduler_calls\":"<<scheduler_calls<<",\"clock_ticks\":"<<clock_ticks<<",\"source_proposals\":"<<source_calls<<",\"passive_checks\":"<<passive_checks<<",\"unexpected_updates\":"<<unexpected_updates<<",\"writer_calls\":"<<writers<<",\"scheduler_field_drift_pmol\":"<<max_field_drift<<",\"writer_drift_pmol\":"<<writer_drift<<",\"peak_RSS_kB\":"<<usage.ru_maxrss<<"}\n";
    summary.close();cache_report(templates,use_cache,calls,output+"/cache.json");timing.mark("final_evidence_flush");ledger.flush();timing.mark("cache_destruction");templates.clear();timing.finish(output+"/profile.json");
    return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<std::endl;return 1;}}
