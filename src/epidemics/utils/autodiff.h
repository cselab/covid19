#pragma once

#include <array>
#include <cmath>
#include <vector>

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

protected:
    std::array<T, 1 + N_> v_;
};

template <typename T>
struct ADDynamicStorage
{
    ADDynamicStorage(std::vector<T> v) : v_(std::move(v)) { }

    size_t N() const noexcept { return v_.size() - 1; }

protected:
    std::vector<T> v_;
};

/** First partial derivatives class. */
template <typename Derived, typename T, typename Storage>
struct AutoDiffBase : Storage
{
    using Storage::Storage;
    using ValueType = T;

    explicit operator T() const {
        return val();
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

    // Unary operators.
    Derived operator+() const {
        return derived();
    }
    Derived operator-() const {
        Derived result = derived();
        result.val() = -result.val();
        for (size_t i = 0; i < result.N(); ++i)
            result.d(i) = -result.d(i);
        return result;
    }

    // Binary operator +.
    friend Derived operator+(Derived a, T b) {
        a.val() += b;
        return a;
    }
    friend Derived operator+(T a, Derived b) {
        b.val() += a;
        return b;
    }
    friend Derived operator+(Derived a, const Derived &b) {
        assert(a.N() == b.N());
        a.val() += b.val();
        for (size_t i = 0; i < a.N(); ++i)
            a.d(i) += b.d(i);
        return a;
    }

    // Binary operator -.
    friend Derived operator-(Derived a, T b) {
        a.val() -= b;
        return a;
    }
    friend Derived operator-(T a, Derived b) {
        b.val() = a - b.val();
        for (size_t i = 0; i < b.N(); ++i)
            b.d(i) = -b.d(i);
        return b;
    }
    friend Derived operator-(Derived a, const Derived &b) {
        assert(a.N() == b.N());
        a.val() -= b.val();
        for (size_t i = 0; i < a.N(); ++i)
            a.d(i) -= b.d(i);
        return a;
    }

    // Binary operator *.
    friend Derived operator*(Derived a, T b) {
        a.val() *= b;
        for (size_t i = 0; i < a.N(); ++i)
            a.d(i) *= b;
        return a;
    }
    friend Derived operator*(T a, Derived b) {
        b.val() *= a;
        for (size_t i = 0; i < b.N(); ++i)
            b.d(i) *= a;
        return b;
    }
    friend Derived operator*(Derived a, const Derived &b) {
        assert(a.N() == b.N());
        T aval = a.val();
        a.val() *= b.val();
        for (size_t i = 0; i < a.N(); ++i)
            a.d(i) = aval * b.d(i) + a.d(i) * b.val();
        return a;
    }

    // Binary operator /.
    friend Derived operator/(Derived a, T b) {
        return std::move(a) * inverse(b);
    }
    friend Derived operator/(T a, Derived b) {
        T inv_b0 = inverse(b.val());

        b.val() = a * inv_b0;
        auto factor = -b.val() * inv_b0;
        for (size_t i = 0; i < b.N(); ++i)
            b.d(i) *= factor;
        return b;
    }
    friend Derived operator/(Derived a, const Derived &b) {
        T inv_b0 = inverse(b.val());

        assert(a.N() == b.N());
        a.val() *= inv_b0;
        for (size_t i = 0; i < a.N(); ++i)
            a.d(i) = (a.d(i) - a.val() * b.d(i)) * inv_b0;
        return a;
    }

    // Assignment-modify operators.
    Derived &operator+=(T b) {
        val() += b;
        return derived();
    }
    Derived &operator+=(const Derived &b) {
        val() += b.val();
        for (size_t i = 0; i < this->N(); ++i)
            d(i) += b.d(i);
        return derived();
    }
    Derived &operator-=(T b) {
        val() -= b;
        return derived();
    }
    Derived &operator-=(const Derived &b) {
        val() -= b.val();
        for (size_t i = 0; i < this->N(); ++i)
            d(i) -= b.d(i);
        return derived();
    }
    Derived &operator*=(T b) {
        val() *= b;
        for (size_t i = 0; i < this->N(); ++i)
            d(i) *= b;
        return derived();
    }
    Derived &operator*=(const Derived &b) {
        for (size_t i = 0; i < this->N(); ++i)
            d(i) = val() * b.d(i) + d(i) * b.val();
        val() *= b.val();
        return derived();
    }
    Derived &operator/=(T b) {
        return (*this *= inverse(b));
    }
    Derived &operator/=(const Derived &b) {
        /* TODO */
        *this = (*this / b);
        return derived();
    }

    Derived inv() const {
        Derived result = derived();
        result.val() = inverse(result.val());
        auto minv2 = -sqr(result.val());
        for (size_t i = 0; i < result.N(); ++i)
            result.d(i) *= minv2;
        return result;
    }

    /* Too lazy to write other operators now... */
#define EPIDEMICS_AD_CMD_OPERATOR(OP) \
    friend bool operator OP(const Derived &a, const Derived &b) { \
        return a.val() OP b.val(); \
    } \
    friend bool operator OP(const Derived &a, T b) { \
        return a.val() OP b; \
    } \
    friend bool operator OP(T a, const Derived &b) { \
        return a OP b.val(); \
    }
    EPIDEMICS_AD_CMD_OPERATOR(<);
    EPIDEMICS_AD_CMD_OPERATOR(<=);
    EPIDEMICS_AD_CMD_OPERATOR(>);
    EPIDEMICS_AD_CMD_OPERATOR(>=);
#undef EPIDEMICS_AD_CMP_OPERATOR

    // Overloads of global functions.
    friend Derived abs(const Derived &x) {
        return x >= 0 ? x : -x;
    }

protected:
    struct size_tag { };

private:
    Derived &derived() { return *static_cast<Derived *>(this); }
    const Derived &derived() const { return *static_cast<const Derived *>(this); }
};

/// Autodifferentiation algebra with a compile-time number of derivatives.
template <typename T, size_t N_>
struct AutoDiff : AutoDiffBase<AutoDiff<T, N_>, T, ADStaticStorage<T, N_>>
{
private:
    using Base = AutoDiffBase<AutoDiff<T, N_>, T, ADStaticStorage<T, N_>>;
    using size_tag = typename Base::size_tag;

public:
    static constexpr size_t kNumVariables = N_;

    AutoDiff() : AutoDiff(T{}) { }
    AutoDiff(size_tag, size_t N) : AutoDiff(T{}) {
        (void)N;
        assert(N == N_);
    }

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

/// Autodifferentiation algebra with a runtime number of derivatives.
template <typename T>
struct DynamicAutoDiff : AutoDiffBase<DynamicAutoDiff<T>, T, ADDynamicStorage<T>>
{
private:
    using Base = AutoDiffBase<DynamicAutoDiff<T>, T, ADDynamicStorage<T>>;
    using size_tag = typename Base::size_tag;

public:
    DynamicAutoDiff(size_tag, size_t N) : Base(std::vector<T>(N)) { }
    DynamicAutoDiff(std::vector<T> v) : Base(std::move(v)) { }
};

template <typename T, typename ...Args>
AutoDiff<T, sizeof...(Args)> make_ad(T value, Args ...derivatives) {
    constexpr size_t N = sizeof...(Args);
    return AutoDiff<T, N>{
        value, std::array<T, N>{static_cast<T>(derivatives)...}};
}

template <typename T, typename ...Args>
DynamicAutoDiff<T> make_dynamic_ad(T value, Args ...derivatives) {
    return DynamicAutoDiff<T>(std::vector<T>{value, static_cast<T>(derivatives)...});
}

}  // namespace epidemics
