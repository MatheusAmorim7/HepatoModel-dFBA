#ifndef DFBA_EXCHANGE_RATES_H
#define DFBA_EXCHANGE_RATES_H

// mmol/gDW/h * gDW /60 * 1e15um^3/L -> mM*um^3/min.
// This conversion must not impose a phenotype or change a solved flux.
inline double solved_exchange_biofvm_rate(double flux_mmol_gDW_h, double dry_mass_g)
{
    return flux_mmol_gDW_h * dry_mass_g / 60.0 * 1e15;
}

#endif
