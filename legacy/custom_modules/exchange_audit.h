#ifndef PHYSICELL_EXCHANGE_AUDIT_H
#define PHYSICELL_EXCHANGE_AUDIT_H
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <cmath>

// Opt-in observation only. No concentration, flux, parameter or phenotype edits.
inline void record_exchange_inventory(const char* phase)
{
    static const bool enabled = std::getenv("PHYSICELL_EXCHANGE_AUDIT") != nullptr;
    if (!enabled) return;
    static std::ofstream file(PhysiCell::PhysiCell_settings.folder + "/exchange_inventory.csv");
    static bool header = false;
    if (!header) {
        file << "time_min,phase,substrate,total_pmol,min_mM,max_mM,negative_voxels,nonfinite_voxels\n";
        header = true;
    }
    auto& env = BioFVM::microenvironment;
    for (int s=0; s<env.number_of_densities(); ++s) {
        double total=0., lo=std::numeric_limits<double>::infinity(), hi=-lo;
        int negative=0, nonfinite=0;
        for (int v=0; v<env.number_of_voxels(); ++v) {
            const double c=env.density_vector(v)[s];
            total += c * env.mesh.voxels[v].volume * 1e-6;
            lo=std::min(lo,c); hi=std::max(hi,c);
            if (c < -1e-12) ++negative;
            if (!std::isfinite(c)) ++nonfinite;
        }
        file << std::setprecision(17) << PhysiCell::PhysiCell_globals.current_time << ','
             << phase << ',' << env.density_names[s] << ',' << total << ',' << lo << ','
             << hi << ',' << negative << ',' << nonfinite << '\n';
    }
}
#endif
