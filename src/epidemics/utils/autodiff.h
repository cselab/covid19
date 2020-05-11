#pragma once

#include <cmath>

namespace epidemics {

template <typename T>
auto sqr(T x) {
    return x * x;
}

// In the case we have AutoDiff<AutoDiff<...>, ...> this may be useful.
template <typename T>
auto inverse(T x) {
    return 1 / x;
}

/** First partial derivatives class. */
template <typename T, int N>
struct AutoDiff {
    AutoDiff() : AutoDiff(T{}) { }

    /* implicit */ AutoDiff(T value) {
        val() = value;
        for (int i = 0; i < N; ++i)
            d(i) = T{};
    }
    template <typename First, typename Second, typename ...Args>
    AutoDiff(First first, Second second, Args ...args) :
        v_{first, second, args...}
    { }

    explicit operator T() const {
        return v_[0];
    }

    /// Return the non-derivative value.
    const T &val() const { return v_[0]; }
    T &val() { return v_[0]; }

    /// Return the first derivative wrt the ith variable.
    const T &d(int i) const { return v_[1 + i]; }
    T &d(int i) { return v_[1 + i]; }

    AutoDiff operator-() const {
        AutoDiff result;
        result.val() = -val();
        for (int i = 0; i < N; ++i)
            result.d(i) = -d(i);
        return result;
    }

    friend AutoDiff operator+(const AutoDiff &a, const AutoDiff &b) {
        AutoDiff result = a;
        result.val() += b.val();
        for (int i = 0; i < N; ++i)
            result.d(i) += b.d(i);
        return result;
    }
    friend AutoDiff operator+(const AutoDiff &a, T b) {
        AutoDiff result = a;
        result.val() += b;
        return result;
    }
    friend AutoDiff operator-(const AutoDiff &a, const AutoDiff &b) {
        AutoDiff result = a;
        result.val() -= b.val();
        for (int i = 0; i < N; ++i)
            result.d(i) -= b.d(i);
        return result;
    }
    friend AutoDiff operator*(const AutoDiff &a, const AutoDiff &b) {
        AutoDiff result;
        result.val() = a.val() * b.val();
        for (int i = 0; i < N; ++i)
            result.d(i) = a.val() * b.d(i) + a.d(i) * b.val();
        return result;
    }
    friend AutoDiff operator/(const AutoDiff &a, const AutoDiff &b) {
        T inv_b0 = inverse(b.v_[0]);

        AutoDiff result;
        result.val() = a.val() * inv_b0;
        for (int i = 0; i < N; ++i)
            result.d(i) = (a.d(i) * b.val() - a.val() * b.d()) * sqr(inv_b0);
        return result;
    }

    AutoDiff &operator+=(const AutoDiff &b) {
        val() += b.val();
        for (int i = 0; i < N; ++i)
            d(i) += b.d(i);
        return *this;
    }
    AutoDiff &operator-=(const AutoDiff &b) {
        val() -= b.val();
        for (int i = 0; i < N; ++i)
            d(i) -= b.d(i);
        return *this;
    }
    AutoDiff &operator*=(const AutoDiff &b) {
        for (int i = 0; i < N; ++i)
            d(i) = val() * b.d(i) + d(i) * b.val();
        val() *= b.val();
        return *this;
    }
    AutoDiff &operator/=(const AutoDiff &b) {
        /* ... */
        return *this = (*this / b);
    }


    /* Too lazy to write other operators now... */
    friend bool operator<(const AutoDiff &a, const AutoDiff &b) {
        return a.val() < b.val();
    }
    friend bool operator<(const AutoDiff &a, const T &b) {
        return a.val() < b;
    }
    friend bool operator>(const AutoDiff &a, const AutoDiff &b) {
        return a.val() > b.val();
    }
    friend bool operator>(const AutoDiff &a, const T &b) {
        return a.val() > b;
    }

    friend AutoDiff abs(const AutoDiff &x) {
        return x >= 0 ? x : -x;
    }

private:
    T v_[N + 1];
};


template <typename T, typename ...Args>
AutoDiff<T, sizeof...(Args)> make_ad(T value, Args ...derivatives) {
    return AutoDiff<T, sizeof...(Args)>{value, derivatives...};
}

}  // namespace epidemics
