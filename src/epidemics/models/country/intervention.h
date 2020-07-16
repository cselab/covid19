#ifndef INTERVENTION_H
#define INTERVENTION_H


template <class T>
T intervention(T R0 /* r0 vefore intervention */, 
               double t /* time */,
               T k /* reduction factor */,
               T tact /* intervention time */,
               T dtact /* intervention duration */)
{
    T r0;
    double t0 = tact - 0.5 * dtact;
    if (t < t0) {
        r0 = R0;
    } else if (t < tact + 0.5 * dtact) {
        r0 = (1. - (t - t0) / dtact * (1. - k)) * R0;
    } else {
        r0 = k * R0;
    }
    return r0;
}

#endif
