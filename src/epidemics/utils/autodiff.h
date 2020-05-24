#pragma once

#include <array>
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

template <typename T, size_t N_>
struct ADStaticStorage
{
    static constexpr size_t N() { return N_; }

    std::array<T, 1 + N_> v_;
};

/** First partial derivatives class. */
template <typename Derived, typename T, typename Storage>
struct AutoDiffBase : protected Storage
{
    using ValueType = T;

    explicit operator T() const {
        return this->v_[0];
    }

    /// Return the non-derivative value.
    const T &val() const { return this->v_[0]; }
    T &val() { return this->v_[0]; }

    /// Return the first derivative wrt the ith variable.
    const T &d(int i) const { return this->v_[1 + i]; }
    T &d(int i) { return this->v_[1 + i]; }

    // NOTE: If we had a span<> here, we could just have `span<T> d();`.
    const T *dBegin() const noexcept { return this->v_.data() + 1; }
    const T *dEnd() const noexcept { return this->v_.data() + this->N() + 1; }
    T *dBegin() noexcept { return this->v_.data() + 1; }
    T *dEnd() noexcept { return this->v_.data() + this->N() + 1; }

    Derived operator-() const {
        Derived result;
        result.val() = -val();
        for (size_t i = 0; i < this->N(); ++i)
            result.d(i) = -d(i);
        return result;
    }

    friend Derived operator+(const Derived &a, const Derived &b) {
        assert(a.N() == b.N());
        Derived result = a;
        result.val() += b.val();
        for (size_t i = 0; i < result.N(); ++i)
            result.d(i) += b.d(i);
        return result;
    }
    friend Derived operator-(const Derived &a, const Derived &b) {
        assert(a.N() == b.N());
        Derived result = a;
        result.val() -= b.val();
        for (size_t i = 0; i < result.N(); ++i)
            result.d(i) -= b.d(i);
        return result;
    }
    friend Derived operator*(const Derived &a, const Derived &b) {
        assert(a.N() == b.N());
        Derived result;
        result.val() = a.val() * b.val();
        for (size_t i = 0; i < result.N(); ++i)
            result.d(i) = a.val() * b.d(i) + a.d(i) * b.val();
        return result;
    }
    friend Derived operator/(const Derived &a, const Derived &b) {
        T inv_b0 = inverse(b.v_[0]);

        Derived result;
        result.val() = a.val() * inv_b0;
        for (size_t i = 0; i < result.N(); ++i)
            result.d(i) = (a.d(i) - result.val() * b.d(i)) * inv_b0;
        return result;
    }

    Derived &operator+=(const Derived &b) {
        val() += b.val();
        for (size_t i = 0; i < this->N(); ++i)
            d(i) += b.d(i);
        return derived();
    }
    Derived &operator-=(const Derived &b) {
        val() -= b.val();
        for (size_t i = 0; i < this->N(); ++i)
            d(i) -= b.d(i);
        return derived();
    }
    Derived &operator*=(const Derived &b) {
        for (size_t i = 0; i < this->N(); ++i)
            d(i) = val() * b.d(i) + d(i) * b.val();
        val() *= b.val();
        return derived();
    }
    Derived &operator/=(const Derived &b) {
        /* TODO */
        *this = (*this / b);
        return derived();
    }

    Derived inv() const {
        Derived result;
        result.val() = inverse(val());
        auto inv2 = sqr(result.val());
        for (size_t i = 0; i < result.N(); ++i)
            result.d(i) = d(i) * -inv2;
        return result;
    }

    /* Too lazy to write other operators now... */
    friend bool operator<(const Derived &a, const Derived &b) {
        return a.val() < b.val();
    }
    friend bool operator<(const Derived &a, const T &b) {
        return a.val() < b;
    }
    friend bool operator>(const Derived &a, const Derived &b) {
        return a.val() > b.val();
    }
    friend bool operator>(const Derived &a, const T &b) {
        return a.val() > b;
    }

    friend Derived abs(const Derived &x) {
        return x >= 0 ? x : -x;
    }

private:
    Derived &derived() { return *static_cast<Derived *>(this); }
    const Derived &derived() const { return *static_cast<const Derived *>(this); }
};

template <typename T, size_t N_>
struct AutoDiff : AutoDiffBase<AutoDiff<T, N_>, T, ADStaticStorage<T, N_>>
{
    static constexpr size_t kNumVariables = N_;

    AutoDiff() : AutoDiff(T{}) { }

    /* implicit */ AutoDiff(T value) {
        this->val() = value;
        for (size_t i = 0; i < this->N(); ++i)
            this->d(i) = T{};
    }
    /* implicit */ AutoDiff(T value, std::array<T, N_> der) {
        this->val() = value;
        for (size_t i = 0; i < this->N(); ++i)
            this->d(i) = der[i];
    }
};

template <typename T, typename ...Args>
AutoDiff<T, sizeof...(Args)> make_ad(T value, Args ...derivatives) {
    constexpr size_t N = sizeof...(Args); 
    return AutoDiff<T, N>{
        value, std::array<T, N>{static_cast<T>(derivatives)...}};
}

}  // namespace epidemics
