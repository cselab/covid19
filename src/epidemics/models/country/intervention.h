#ifndef INTERVENTION_H
#define INTERVENTION_H


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


template <class T>
T intervention_step(T R0 /* r0 before intervention */,
               double t /* time */,
               T kbeta /* reduction factor */,
               T tact /* intervention time */
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
               T tact /* intervention time */
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
