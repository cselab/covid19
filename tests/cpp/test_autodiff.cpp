#include <gtest/gtest.h>

#include <epidemics/utils/autodiff.h>

using namespace epidemics;

template <typename T, size_t N>
void assertEqAD(AutoDiff<T, N> a, AutoDiff<T, N> b) {
    ASSERT_EQ(a.val(), b.val());
    for (size_t i = 0; i < N; ++i)
        ASSERT_EQ(a.d(i), b.d(i));
}

TEST(AutoDiff, arithmetic) {
    AutoDiff<double, 2> x = make_ad(10.,  1., 0.);
    AutoDiff<double, 2> y = make_ad(100., 0., 1.);

    auto func1 = [](auto x, auto y) {
        return x * x + y * y * y;
    };
    assertEqAD(func1(x, 0), make_ad(100.,     20., 0.));
    assertEqAD(func1(0, y), make_ad(1000000., 0.,  30000.));
    assertEqAD(func1(x, y), make_ad(1000100., 20., 30000.));

    auto func2 = [](auto x, auto y) {
        return x * y;
    };
    assertEqAD(func2(x, y), make_ad(1000., 100., 10.));
}

TEST(AutoDiff, inv) {
    auto x = make_ad(8.0, 3.0);
    assertEqAD(x / x,   make_ad(1.0, 0.0));
    assertEqAD(x.inv(), 1 / x);
    assertEqAD(x.inv(), make_ad(1. / 8, -3. / 64));
}
