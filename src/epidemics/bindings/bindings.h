#include <pybind11/pybind11.h>

#include <epidemics/utils/signal.h>

class SignalRAII {
public:
    SignalRAII() {
        check_signals_func = []() {
            // https://stackoverflow.com/questions/14707049/allowing-ctrl-c-to-interrupt-a-python-c-extension
            if (PyErr_CheckSignals() != 0)
                throw std::runtime_error("Signal received. Breaking.");
        };
    }
    ~SignalRAII() {
        check_signals_func = nullptr;
    }
};
