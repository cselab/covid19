#ifndef INTERVENTION_H
#define INTERVENTION_H

#ifdef CENTERED
template <class T>
T intervention(T R0 /* r0 before intervention */,
               double t /* time */,
               T kbeta /* reduction factor */,
               T tact /* intervention time */,
               T dtact /* intervention duration */)
{
    T r0;
    T t0 = tact - 0.5 * dtact;
    if (t < t0) {
        r0 = R0;
    } else if (t < tact + 0.5 * dtact) {
        r0 = (1. - (1 - kbeta) / dtact * (t - t0)) * R0;
    } else {
        r0 = kbeta * R0;
    }
    return r0;
}
#else
template <class T>
T intervention(T R0 /* r0 before intervention */,
               double t /* time */,
               T kbeta /* reduction factor */,
               T tact /* intervention time */,
               T dtact /* intervention duration */)
{
    T r0;
    T t0 = tact;
    if (t < t0) {
        r0 = R0;
    } else if (t < tact + dtact) {
        r0 = (1. - (1 - kbeta) / dtact * (t - t0)) * R0;
    } else {
        r0 = kbeta * R0;
    }
    return r0;
}
#endif

#ifdef CENTERED
template <class T>
T intervention_smooth(T R0 /* r0 before intervention */,
               double t /* time */,
               T kbeta /* reduction factor */,
               T tact /* intervention time */,
               T dtact /* intervention duration */)
{
    using std::log;
    T c = -2.0*log(0.025/0.975)/dtact;
    
    using std::exp;
    return R0-R0*(1-kbeta)/(1+exp(-c*(t-tact)));
}
#else
template <class T>
T intervention_smooth(T R0 /* r0 before intervention */,
               double t /* time */,
               T kbeta /* reduction factor */,
               T tact /* intervention time */,
               T dtact /* intervention duration */)
{
    using std::log;
    T c = -2.0*log(0.025/0.975)/dtact;
    
    using std::exp;
    return R0-R0*(1-kbeta)/(1+exp(-c*(t-tact-0.5*dtact)));
}
#endif

template <class T>
T intervention_step(T R0 /* r0 before intervention */,
               double t /* time */,
               T kbeta /* reduction factor */,
               T tact /* intervention time */)
{
    T r0;
    if (t < tact) {
        r0 = R0;
    } else { 
        r0 = R0*kbeta;
    }
    return r0;
}

   
template <class T>
T intervention_exp(T R0 /* r0 before intervention */,
               double t /* time */,
               T k /* reduction factor */,
               T tact /* intervention time */)
{
    T r0;
    if (t < tact) {
        r0 = R0;
    } else { 
        using std::exp;
        double dt = t-tact;     // careful: diff will not work
        r0 = R0*exp(-k*dt);
    }
    return r0;
}



#endif
