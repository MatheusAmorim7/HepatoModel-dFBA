#ifndef PHYSICELL_SHARED_POOL_EVENT_CONTROL_H
#define PHYSICELL_SHARED_POOL_EVENT_CONTROL_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <sys/statvfs.h>
#include <coin/CoinPackedMatrix.hpp>
#include <coin/CoinPackedVector.hpp>
#include "../addons/dFBA/src/dfba_Reaction.h"

namespace SharedPoolEvents
{
using namespace PhysiCell;
using namespace BioFVM;
using PhysiCelldFBA::dFBAIntracellular;
const double INF = 1e100;
const double EPS = 1e-10;

inline void require(bool ok, const std::string& message)
{ if (!ok) throw std::runtime_error(message); }

inline std::string unprefix(const std::string& value)
{ return value.substr(0, 2) == "R_" ? value.substr(2) : value; }

inline std::string csv_set(const std::set<std::string>& values)
{
    std::ostringstream out;
    out << '"' << '[';
    bool first = true;
    for (const auto& value : values)
    {
        require(value.find_first_of("\"\\\n,") == std::string::npos, "Unsupported field name");
        if (!first) out << ',';
        out << "\"\"" << value << "\"\"";
        first = false;
    }
    out << ']' << '"';
    return out.str();
}

struct Port { int field, column; };
struct Group
{
    std::string name;
    std::vector<Cell*> cells;
    dFBAIntracellular* intra;
    double scale;
    int offset;
    std::vector<Port> ports;
};

class Control
{
    std::vector<Group> groups;
    std::vector<std::string> fields;
    std::vector<double> initial;
    std::vector<int> pool_rows;
    ClpSimplex lp;
    double pmol_per_mM;
    int oxygen;
    std::string folder;

    std::vector<double> stocks() const
    {
        std::vector<double> result(fields.size());
        for (size_t i=0; i<fields.size(); ++i) result[i]=microenvironment(0)[i]*pmol_per_mM;
        return result;
    }

    void snapshot(const std::string& name, double time)
    {
        PhysiCell_globals.current_time = time;
        save_PhysiCell_to_MultiCellDS_v2(folder+"/"+name, microenvironment, time);
    }

    bool solve()
    {
        lp.dual();
        if (lp.isProvenOptimal()) return true;
        require(lp.isProvenPrimalInfeasible(), "Unclassified Clp status: "+std::to_string(lp.status()));
        return false;
    }

public:
    explicit Control(pugi::xml_node root)
    {
        require(microenvironment.number_of_voxels()==1, "One voxel required");
        pmol_per_mM=microenvironment.mesh.voxels[0].volume*1e-6;
        require(std::abs(pmol_per_mM-20)<1e-10, "Expected 20 nL");
        require(std::abs(PhysiCell_settings.max_time-720)<1e-10, "Expected 720 minutes");
        folder=PhysiCell_settings.folder;
        auto user=root.child("user_parameters");
        require(!user.child("culture_finite_o2_enabled").text().as_bool(), "Disable the legacy oxygen source");
        require(user.child("culture_o2_reservoir_pmol").text().as_double()==32 &&
            user.child("culture_o2_feed_pmol_h").text().as_double()==3 &&
            user.child("culture_o2_target_pmol").text().as_double()==4, "Expected reference oxygen controller");
        for (auto variable : root.child("microenvironment_setup").children("variable"))
        {
            int f=microenvironment.find_density_index(variable.attribute("name").value());
            require(f>=0, "Unknown field");
            require(std::string(variable.attribute("units").value())=="mM", "Expected mM fields");
            require(microenvironment.diffusion_coefficients[f]==0 && microenvironment.decay_rates[f]==0,
                "No diffusion/decay permitted in this control");
            require(!variable.child("Dirichlet_boundary_condition").attribute("enabled").as_bool(),
                "Dirichlet not permitted");
        }
        fields=microenvironment.density_names;
        initial=stocks();
        oxygen=microenvironment.find_density_index("oxygen");
        require(oxygen>=0 && std::abs(initial[oxygen]-4)<1e-10, "Expected 4 pmol oxygen");
        std::map<std::string,std::vector<Cell*> > by_name;
        for (Cell* cell : *all_cells)
        {
            require(!cell->phenotype.death.dead, "Live initial cells required");
            by_name[cell->type_name].push_back(cell);
        }
        require(by_name.size()==3, "Expected three homogeneous groups");
        int columns=0;
        std::vector<double> lower,upper,objective;
        for (const std::string name : {"basal","local","tumor"})
        {
            Group group;
            group.name=name;
            group.cells=by_name.at(name);
            require(!group.cells.empty(), "Empty group");
            group.intra=dynamic_cast<dFBAIntracellular*>(group.cells[0]->phenotype.intracellular);
            require(group.intra && group.intra->experimental_culture, "Culture dFBA required");
            group.offset=columns;
            auto settings=root.child("cell_definitions").find_child_by_attribute("cell_definition","name",name.c_str())
                .child("phenotype").child("intracellular");
            auto culture=settings.child("settings").child("experimental_culture");
            double flux_scale=culture.child("flux_scale_mmol_gDW_h").text().as_double();
            double mass=culture.child("dry_mass_g").text().as_double();
            require(flux_scale>0 && mass>0, "Positive physical conversion required");
            group.scale=mass*flux_scale*1e9*group.cells.size();
            double expected=(name=="basal" ? 2.52 : name=="local" ? .12 : 1.);
            require(std::abs(group.scale-expected)<1e-10, "Aggregate physical scale changed");
            auto& reactions=group.intra->sbml_model.reactions;
            std::map<std::string,int> index;
            for (size_t j=0;j<reactions.size();++j)
            {
                index[reactions[j]->getId()]=j;
                lower.push_back(reactions[j]->getLowerBound());
                upper.push_back(reactions[j]->getUpperBound());
                objective.push_back(0.);
            }
            for (auto exchange : settings.child("transport_model").children("exchange"))
            {
                int f=microenvironment.find_density_index(exchange.attribute("substrate").value());
                require(f>=0, "Unknown port field");
                int j=index.at(exchange.child("fba_flux").text().as_string());
                double cap=exchange.child("allow_uptake").text().as_bool()
                    ? exchange.child("Vmax").text().as_double()/flux_scale : 0.;
                // Use XML doubles, not the legacy float-valued kinetic struct.
                lower[group.offset+j]=-cap;
                require(lower[group.offset+j]<=upper[group.offset+j], "Inconsistent port bounds");
                group.ports.push_back({f,group.offset+j});
            }
            for (Cell* cell : group.cells)
            {
                auto intra=dynamic_cast<dFBAIntracellular*>(cell->phenotype.intracellular);
                require(intra && intra->sbml_filename==group.intra->sbml_filename,
                    "Heterogeneous model cannot be grouped");
                require(intra->culture_dry_mass_g==group.intra->culture_dry_mass_g &&
                    intra->culture_flux_scale==group.intra->culture_flux_scale, "Heterogeneous conversion");
                require(cell->get_current_voxel_index()==0, "Different voxel cannot be grouped");
                auto& rr=intra->sbml_model.reactions;
                require(rr.size()==reactions.size(), "Heterogeneous reactions");
                for (size_t j=0;j<rr.size();++j)
                    require(rr[j]->getId()==reactions[j]->getId() && rr[j]->getLowerBound()==reactions[j]->getLowerBound()
                        && rr[j]->getUpperBound()==reactions[j]->getUpperBound(), "Heterogeneous bounds");
            }
            columns+=reactions.size();
            groups.push_back(group);
        }
        for (const auto& group : groups)
            for (size_t j=0;j<group.intra->sbml_model.reactions.size();++j)
            { lower.push_back(0.); upper.push_back(INF); objective.push_back(group.scale); }
        CoinPackedMatrix matrix(false,0,0);
        matrix.setDimensions(0,2*columns);
        std::vector<double> row_lower,row_upper;
        auto row=[&](const CoinPackedVector& values,double lo,double hi)
        { matrix.appendRow(values); row_lower.push_back(lo); row_upper.push_back(hi); };
        for (const auto& group : groups)
        {
            auto& model=group.intra->sbml_model;
            std::map<std::string,CoinPackedVector> balances;
            for (size_t j=0;j<model.reactions.size();++j)
                for (const auto& term : model.reactions[j]->getMetabolites())
                    balances[term.first].insert(group.offset+j,term.second);
            for (const auto& balance : balances) row(balance.second,0,0);
        }
        for (int j=0;j<columns;++j)
        {
            CoinPackedVector positive,negative;
            positive.insert(j,1); positive.insert(columns+j,-1);
            negative.insert(j,-1); negative.insert(columns+j,-1);
            row(positive,-INF,0); row(negative,-INF,0);
        }
        for (size_t f=0;f<fields.size();++f)
        {
            CoinPackedVector consumption;
            for (const auto& group : groups)
                for (const auto& port : group.ports)
                    if (port.field==int(f)) consumption.insert(port.column,-group.scale);
            pool_rows.push_back(row_lower.size());
            row(consumption,-INF,INF);
        }
        lp.loadProblem(matrix,lower.data(),upper.data(),objective.data(),row_lower.data(),row_upper.data());
        lp.setOptimizationDirection(1);
        lp.setPrimalTolerance(1e-9);
        lp.setDualTolerance(1e-9);
        lp.setLogLevel(0);
        std::ofstream definition(folder+"/native_definition.csv");
        definition << std::setprecision(17) << "cell_type,count,aggregate_scale,reaction,lower,upper\n";
        for (const auto& group : groups)
            for (size_t j=0;j<group.intra->sbml_model.reactions.size();++j)
                definition << group.name << ',' << group.cells.size() << ',' << group.scale << ','
                    << unprefix(group.intra->sbml_model.reactions[j]->getId()) << ','
                    << lower[group.offset+j] << ',' << upper[group.offset+j] << '\n';
    }

    void run(double maximum)
    {
        require(std::isfinite(maximum) && maximum>EPS && maximum<=2, "Interval must be in (0,2]");
        set_save_biofvm_mesh_as_matlab(true);
        set_save_biofvm_data_as_matlab(true);
        set_save_biofvm_cell_data(true);
        set_save_biofvm_cell_data_as_custom_matlab(true);
        snapshot("initial",0);
        std::ofstream segments(folder+"/segments.csv"), fluxes(folder+"/fluxes.csv");
        segments << std::setprecision(17);
        fluxes << std::setprecision(17) << "solution_id,time_min,status,cell_type,reaction,flux_normalized\n";
        segments << "start_min,end_min,duration_min,solution_id,status,feed_pmol_h,feed_regime,reservoir_before_pmol,reservoir_after_pmol,blocked,hits";
        for (auto f:fields) segments << ",before_" << f;
        for (auto f:fields) segments << ',' << f;
        segments << '\n';
        std::set<std::string> blocked;
        for (size_t f=0;f<fields.size();++f) if (initial[f]<=0) blocked.insert(fields[f]);
        double t=0,remaining=32,covered=0,first_failure=-1;
        bool empty=false,target=true,terminal=false;
        int solution_id=0,zero_steps=0;
        // Only this loop owns the medium. No cell, secretion, diffusion or old feed loop is run.
        while (t<720-EPS)
        {
            if (solution_id%100==0)
            {
                struct statvfs disk;
                require(statvfs("/mnt/c",&disk)==0, "Cannot check Windows free space");
                require(double(disk.f_bavail)*disk.f_frsize>10.*1024*1024*1024, "Windows disk below 10 GiB floor");
            }
            auto before=stocks();
            for (size_t f=0;f<fields.size();++f)
            {
                if (before[f]>1e-9) blocked.erase(fields[f]);
                lp.setRowUpper(pool_rows[f],blocked.count(fields[f]) ? (int(f)==oxygen && !empty ? 3. : 0.) : INF);
            }
            std::vector<double> rates(fields.size(),0);
            std::string status;
            if (t<2-EPS) status="startup";
            else if (terminal) status="infeasible";
            else if (solve()) status="optimal";
            else
            {
                double original=lp.rowUpper()[pool_rows[oxygen]];
                lp.setRowUpper(pool_rows[oxygen],INF);
                bool future=(!empty && solve());
                lp.setRowUpper(pool_rows[oxygen],original);
                require(!future, "Future oxygen recovery requires a separate explicit policy");
                terminal=true;
                first_failure=t;
                status="infeasible";
            }
            for (const auto& group : groups)
            {
                auto& rr=group.intra->sbml_model.reactions;
                for (size_t j=0;j<rr.size();++j)
                {
                    double value=status=="optimal" ? lp.primalColumnSolution()[group.offset+j] : 0.;
                    fluxes << solution_id << ',' << t << ',' << status << ',' << group.name << ',' << unprefix(rr[j]->getId()) << ',' << value << '\n';
                    for (Cell* cell : group.cells)
                    {
                        auto intra=static_cast<dFBAIntracellular*>(cell->phenotype.intracellular);
                        intra->sbml_model.reactions[j]->setFluxValue(value);
                        intra->last_solution_status=status;
                    }
                }
                if (status=="optimal") for (const auto& port : group.ports)
                    rates[port.field]+=group.scale*lp.primalColumnSolution()[port.column];
            }
            double feed=empty ? 0. : (target ? std::min(3.,std::max(0.,-rates[oxygen])) : 3.);
            double requested=std::min(maximum,720-t);
            if (t<2-EPS) requested=std::min(requested,2-t);
            std::vector<std::pair<double,std::string> > candidates={{requested,"clock"}};
            for (size_t f=0;f<fields.size();++f)
            {
                double slope=(rates[f]+(int(f)==oxygen ? feed : 0.))/60;
                if (!blocked.count(fields[f]) && slope<-1e-12)
                    candidates.push_back({std::max(0.,before[f])/(-slope),"zero:"+fields[f]});
            }
            if (feed>0) candidates.push_back({std::max(0.,remaining)*60/feed,"reservoir"});
            double oxygen_slope=(rates[oxygen]+feed)/60;
            if (!target && oxygen_slope>1e-12)
                candidates.push_back({std::max(0.,4-before[oxygen])/oxygen_slope,"target"});
            double duration=requested;
            for (const auto& candidate : candidates) duration=std::min(duration,candidate.first);
            std::set<std::string> hits;
            for (const auto& candidate : candidates) if (std::abs(candidate.first-duration)<=EPS) hits.insert(candidate.second);
            double reservoir_before=remaining;
            for (size_t f=0;f<fields.size();++f)
                microenvironment(0)[f]+=(rates[f]+(int(f)==oxygen ? feed : 0.))*duration/(60*pmol_per_mM);
            remaining-=feed*duration/60;
            auto after=stocks();
            segments << t << ',' << t+duration << ',' << duration << ',' << solution_id << ',' << status << ',' << feed << ','
                << (empty ? "off" : target ? "target" : "below_target") << ',' << reservoir_before << ',' << remaining << ','
                << csv_set(blocked) << ',' << csv_set(hits);
            for (double value:before) segments << ',' << value;
            for (double value:after) segments << ',' << value;
            segments << '\n';
            if (status=="optimal") covered+=duration;
            for (const auto& hit : hits)
            {
                if (hit.substr(0,5)=="zero:") blocked.insert(hit.substr(5));
                else if (hit=="reservoir") empty=true;
                else if (hit=="target") target=true;
            }
            if (after[oxygen]<4-1e-9) target=false;
            if (duration<=EPS) require(++zero_steps<12,"Zero-time event loop"); else zero_steps=0;
            t+=duration;
            ++solution_id;
            for (double value:after) require(std::isfinite(value) && value>=-1e-8,"Invalid native stock");
            require(remaining>=-1e-8 && after[oxygen]<=4+1e-8,"Invalid oxygen budget");
        }
        require(segments.good() && fluxes.good(),"Output write failed");
        segments.close(); fluxes.close();
        snapshot("final",t);
        std::ofstream summary(folder+"/native_summary.json");
        summary << std::setprecision(17) << "{\n\"end_min\":" << t << ",\n\"maximum_interval_min\":" << maximum
            << ",\n\"first_infeasible_min\":" << (first_failure<0 ? "null" : precise(first_failure))
            << ",\n\"feasible_minutes\":" << covered << ",\n\"feasible_percent\":" << covered/720*100
            << ",\n\"reservoir_remaining_pmol\":" << remaining << ",\n\"cells_frozen\":true,\n\"biological_validation\":false,\n\"all_groups_zero_on_joint_failure\":true,\n\"native_agents\":" << all_cells->size() << "\n}\n";
        require(summary.good(),"Summary write failed");
        std::cout << "NATIVE_EVENT_CONTROL_FINISHED " << t << " " << covered/720*100 << std::endl;
    }

    static std::string precise(double value)
    { std::ostringstream s; s<<std::setprecision(17)<<value; return s.str(); }
};
}
#endif
