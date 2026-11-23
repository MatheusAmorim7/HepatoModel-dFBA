#ifndef FINITE_POOL_FEED_H
#define FINITE_POOL_FEED_H
#include <algorithm>
#include <cmath>
#include <stdexcept>

class FinitePoolFeed {
    double left_, rate_, target_;
public:
    FinitePoolFeed(double amount, double rate_per_hour, double target)
        : left_(amount), rate_(rate_per_hour), target_(target) {
        if(!std::isfinite(amount)||!std::isfinite(rate_)||!std::isfinite(target_)||
           amount<0.||rate_<0.||target_<0.) throw std::invalid_argument("Invalid finite feed");
    }
    double remaining() const {return left_;}
    double advance(double stock, double dt_min) {
        if(!std::isfinite(stock)||!std::isfinite(dt_min)||stock<0.||dt_min<=0.)
            throw std::invalid_argument("Invalid feed stock or timestep");
        const double amount=std::min(left_,std::min(rate_*dt_min/60.,std::max(0.,target_-stock)));
        left_-=amount;
        return amount;
    }
};
#endif
