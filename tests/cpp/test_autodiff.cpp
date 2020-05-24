#include <gtest/gtest.h>

#include <epidemics/utils/autodiff.h>

using namespace epidemics;

template <typename AD>
static void assertEqAD(const AD &a, const AD &b) {
    ASSERT_EQ(a.N(), b.N());
    ASSERT_EQ(a.val(), b.val());
    for (size_t i = 0; i < a.N(); ++i)
        ASSERT_EQ(a.d(i), b.d(i));
}

template <typename Make>
static void test_arithmetic(Make make) {
    // Test unary +, -.
    assertEqAD(+make(100., 10., 20.), make(100., 10., 20.));
    assertEqAD(-make(100., 10., 20.), make(-100., -10., -20.));

    // Test binary +.
    assertEqAD(make(100., 10., 20.) + 1, make(101., 10., 20.));
    assertEqAD(1 + make(100., 10., 20.), make(101., 10., 20.));
    assertEqAD(make(100., 10., 20.) + make(10., 2., 3.), make(110., 12., 23.));

    // Test binary -.
    assertEqAD(make(100., 10., 20.) - 1, make(99., 10., 20.));
    assertEqAD(1 - make(100., 10., 20.), make(-99., -10., -20.));
    assertEqAD(make(100., 10., 20.) - make(10., 2., 3.), make(90., 8., 17.));

    // Test *.
    assertEqAD(make(100., 10., 20.) * 3, make(300., 30., 60.));
    assertEqAD(3 * make(100., 10., 20.), make(300., 30., 60.));
    assertEqAD(make(3., 4., 5.) * make(11., 12., 13.), make(33., 3*12+4*11, 3*13+5*11));

    // Test /.
    assertEqAD(make(100., 10., 20.) / 2, make(50., 5., 10.));
    assertEqAD(2 / make(8., 3., 5.), make(1. / 4, -3. / 32, -5. / 32));
    assertEqAD(make(8., 3., 5.).inv(), make(1. / 8, -3. / 64, -5. / 64));
    assertEqAD(make(3., 4., 5.) / make(8., 3., 5.), make(3 / 8., (4*8-3*3)/64., (5*8-5*3)/64.));
}

template <typename Make>
static void test_comparison(Make make) {
    ASSERT_TRUE(5 > make(4., 1.));
    ASSERT_FALSE(5 > make(5., 1.));
    ASSERT_FALSE(5 > make(6., 1.));

    ASSERT_TRUE(5 >= make(4., 1.));
    ASSERT_TRUE(5 >= make(5., 1.));
    ASSERT_FALSE(5 >= make(6., 1.));

    ASSERT_FALSE(5 < make(4., 1.));
    ASSERT_FALSE(5 < make(5., 1.));
    ASSERT_TRUE(5 < make(6., 1.));

    ASSERT_FALSE(5 <= make(4., 1.));
    ASSERT_TRUE(5 <= make(5., 1.));
    ASSERT_TRUE(5 <= make(6., 1.));


    ASSERT_FALSE(make(4., 1.) > 5);
    ASSERT_FALSE(make(5., 1.) > 5);
    ASSERT_TRUE(make(6., 1.) > 5);

    ASSERT_FALSE(make(4., 1.) >= 5);
    ASSERT_TRUE(make(5., 1.) >= 5);
    ASSERT_TRUE(make(6., 1.) >= 5);

    ASSERT_TRUE(make(4., 1.) < 5);
    ASSERT_FALSE(make(5., 1.) < 5);
    ASSERT_FALSE(make(6., 1.) < 5);

    ASSERT_TRUE(make(4., 1.) <= 5);
    ASSERT_TRUE(make(5., 1.) <= 5);
    ASSERT_FALSE(make(6., 1.) <= 5);


    ASSERT_FALSE(make(4., 1.) > make(5., 2.));
    ASSERT_FALSE(make(5., 1.) > make(5., 2.));
    ASSERT_TRUE(make(6., 1.) > make(5., 2.));

    ASSERT_FALSE(make(4., 1.) >= make(5., 2.));
    ASSERT_TRUE(make(5., 1.) >= make(5., 2.));
    ASSERT_TRUE(make(6., 1.) >= make(5., 2.));

    ASSERT_TRUE(make(4., 1.) < make(5., 2.));
    ASSERT_FALSE(make(5., 1.) < make(5., 2.));
    ASSERT_FALSE(make(6., 1.) < make(5., 2.));

    ASSERT_TRUE(make(4., 1.) <= make(5., 2.));
    ASSERT_TRUE(make(5., 1.) <= make(5., 2.));
    ASSERT_FALSE(make(6., 1.) <= make(5., 2.));
}

template <typename Make>
static void test_expression(Make make) {
    auto x = make(10.,  1., 0.);
    auto y = make(100., 0., 1.);

    auto func1 = [](const auto &x, const auto &y) {
        return x * x + y * y * y;
    };
    assertEqAD(func1(x, 0), make(100.,     20., 0.));
    assertEqAD(func1(0, y), make(1000000., 0.,  30000.));
    assertEqAD(func1(x, y), make(1000100., 20., 30000.));

    auto func2 = [](const auto &x, const auto &y) {
        return x * y;
    };
    assertEqAD(func2(x, y), make(1000., 100., 10.));
}

TEST(AutoDiff, arithmetic) {
    test_arithmetic(make_ad<double, double, double>);
    test_arithmetic(make_dynamic_ad<double, double, double>);
}

TEST(AutoDiff, comparison) {
    test_comparison(make_ad<double, double>);
    test_comparison(make_dynamic_ad<double, double>);
}

TEST(AutoDiff, expression) {
    test_expression(make_ad<double, double, double>);
    test_expression(make_dynamic_ad<double, double, double>);
}
