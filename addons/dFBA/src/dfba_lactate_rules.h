#ifndef PHYSICELL_DFBA_LACTATE_RULES_H
#define PHYSICELL_DFBA_LACTATE_RULES_H

namespace PhysiCelldFBA
{
constexpr double lactate_oxygen_switch_mM = 0.06;

inline double lactate_auxiliary_uptake_oxygen_gate(double oxygen_mM)
{
    // Equal onset/saturation thresholds define a step, not an interpolation.
    return oxygen_mM > lactate_oxygen_switch_mM ? 1.0 : 0.0;
}
}

#endif
