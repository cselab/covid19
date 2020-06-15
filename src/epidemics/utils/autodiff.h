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
    static constexpr size_t N() noexcept { return N_; }

protected:
    static constexpr void checkSize() noexcept { }
    static constexpr void checkSize(const ADStaticStorage &) noexcept { }

    T &get(size_t index) {
        assert(index <= N_);
        return v_[index];
    }
    const T &get(size_t index) const {
        assert(index <= N_);
        return v_[index];
    }
    T *data() noexcept { return v_.data(); }
    const T *data() const noexcept { return v_.data(); }

private:
    std::array<T, 1 + N_> v_;
};

template <typename T>
struct ADDynamicStorage
{
    ADDynamicStorage(std::vector<T> v) : v_(std::move(v)) { }
    ADDynamicStorage(T value, const std::vector<T> &d) : v_(d.size() + 1) {
        v_[0] = value;
        for (size_t i = 0; i < d.size(); ++i)
            v_[1 + i] = d[i];
    }

    size_t N() const noexcept { return v_.size() - 1; }
protected:
    void checkSize() const {
        assert(!v_.empty());
    }
    void checkSize(const ADDynamicStorage &other) const {
        (void)other;
        assert(v_.size() == other.v_.size());
    }

    T &get(size_t index) {
        assert(index < v_.size());
        return v_[index];
    }
    const T &get(size_t index) const {
        assert(index < v_.size());
        return v_[index];
    }
    T *data() noexcept { return v_.data(); }
    const T *data() const noexcept { return v_.data(); }

private:
    // TODO: A custom container with T* v_ and size_t N.
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
    const T &val() const { return this->get(0); }
    T &val() { return this->get(0); }

    /// Return the first derivative wrt the ith variable.
    const T &d(int i) const { return this->get(1 + i); }
    T &d(int i) { return this->get(1 + i); }

    // NOTE: If we had a span<> here, we could just have `span<T> d();`.
    const T *dBegin() const noexcept { return this->data() + 1; }
    const T *dEnd() const noexcept { return this->data() + this->N() + 1; }
    T *dBegin() noexcept { return this->data() + 1; }
    T *dEnd() noexcept { return this->data() + this->N() + 1; }

    // Unary operators.
    Derived operator+() const {
        return derived();
    }
    Derived operator-() const {
        this->checkSize();
        Derived result(size_tag{}, this->N());
        result.val() = -val();
        for (size_t i = 0; i < this->N(); ++i)
            result.d(i) = -d(i);
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
        a += b;
        return a;
    }
    Derived &operator+=(T b) {
        val() += b;
        return derived();
    }
    Derived &operator+=(const Derived &b) {
        this->checkSize(b);
        val() += b.val();
        for (size_t i = 0; i < this->N(); ++i)
            d(i) += b.d(i);
        return derived();
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
        a -= b;
        return a;
    }
    Derived &operator-=(T b) {
        val() -= b;
        return derived();
    }
    Derived &operator-=(const Derived &b) {
        this->checkSize(b);
        val() -= b.val();
        for (size_t i = 0; i < this->N(); ++i)
            d(i) -= b.d(i);
        return derived();
    }

    // Binary operator *.
    friend Derived operator*(Derived a, T b) {
        a *= b;
        return a;
    }
    friend Derived operator*(T a, Derived b) {
        b *= a;
        return b;
    }
    friend Derived operator*(Derived a, const Derived &b) {
        a *= b;
        return a;
    }
    Derived &operator*=(T b) {
        val() *= b;
        for (size_t i = 0; i < this->N(); ++i)
            d(i) *= b;
        return derived();
    }
    Derived &operator*=(const Derived &b) {
        this->checkSize(b);
        for (size_t i = 0; i < this->N(); ++i)
            d(i) = val() * b.d(i) + d(i) * b.val();
        val() *= b.val();
        return derived();
    }

    // Binary operator /.
    friend Derived operator/(Derived a, T b) {
        a /= b;
        return a;
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
        a /= b;
        return a;
    }
    Derived &operator/=(T b) {
        return *this *= inverse(b);
    }
    Derived &operator/=(const Derived &b) {
        this->checkSize(b);
        T inv_b0 = inverse(b.val());
        val() *= inv_b0;
        for (size_t i = 0; i < this->N(); ++i)
            d(i) = (d(i) - val() * b.d(i)) * inv_b0;
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

    friend Derived exp(Derived x) {
        using std::exp;
        x.val() = exp(x.val());
        for (size_t i = 0; i < x.N(); ++i)
            x.d(i) *= x.val();
        return x;
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

public:
    using size_tag = typename Base::size_tag;
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

public:
    using size_tag = typename Base::size_tag;

    DynamicAutoDiff(size_tag, size_t N) : Base(std::vector<T>(1 + N)) { }
    DynamicAutoDiff(std::vector<T> v) : Base(std::move(v)) { }
    DynamicAutoDiff(T value, const std::vector<T> &d) : Base(value, d) { }

    // Cannot construct a DynamicAutoDiff without num of derivatives.
    DynamicAutoDiff(T value) = delete;

    // NOTE: This is dangerous to use, we have it because of boost!
    DynamicAutoDiff() : Base(std::vector<T>{}) { }
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
