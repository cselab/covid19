#include <backward.hpp>

namespace epidemics {

backward::SignalHandling sh;

void printStackTrace(FILE *f, int depth) {
    backward::StackTrace st;
    st.load_here(depth);
    st.skip_n_firsts(2);
    backward::Printer p;
    p.object = true;
    p.color_mode = backward::ColorMode::automatic;
    p.address = true;
    p.print(st, f);
}

}  // namespace epidemics
