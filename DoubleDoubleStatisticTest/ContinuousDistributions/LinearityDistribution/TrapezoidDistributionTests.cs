using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.LinearityDistribution {
    [TestClass()]
    public class TrapezoidDistributionTests {
        readonly TrapezoidDistribution dist_a0b1c4d5 = new(a: 0, b: 1, c: 4, d: 5);
        readonly TrapezoidDistribution dist_a1b2c3d6 = new(a: 1, b: 2, c: 3, d: 6);
        readonly TrapezoidDistribution dist_a2b3c6d10 = new(a: 2, b: 3, c: 6, d: 10);

        TrapezoidDistribution[] Dists => [
            dist_a0b1c4d5 ,
            dist_a1b2c3d6 ,
            dist_a2b3c6d10,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"A={dist.A}");
                Console.WriteLine($"B={dist.B}");
                Console.WriteLine($"C={dist.C}");
                Console.WriteLine($"D={dist.D}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode - 1e-4), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode + 1e-4), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 7; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 7; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 7; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (TrapezoidDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.01, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (TrapezoidDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.NegativeInfinity)) && dist.PDF(ddouble.NegativeInfinity) >= 0d, "pdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MinValue)) && dist.PDF(ddouble.MinValue) >= 0d, "pdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MinValue / 2)) && dist.PDF(ddouble.MinValue / 2) >= 0d, "pdf(-lval / 2)");

                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.PositiveInfinity)) && dist.PDF(ddouble.PositiveInfinity) >= 0d, "pdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MaxValue)) && dist.PDF(ddouble.MaxValue) >= 0d, "pdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MaxValue / 2)) && dist.PDF(ddouble.MaxValue / 2) >= 0d, "pdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.PDF(ddouble.NaN)), "pdf(NaN)");

                Assert.IsTrue(dist.CDF(ddouble.NegativeInfinity, Interval.Lower) == 0d, "cdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue, Interval.Lower)) && dist.CDF(ddouble.MinValue, Interval.Lower) >= 0d, "cdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue / 2, Interval.Lower)) && dist.CDF(ddouble.MinValue / 2, Interval.Lower) >= 0d, "cdf(-lval / 2)");

                Assert.IsTrue(dist.CDF(ddouble.PositiveInfinity, Interval.Lower) == 1d, "cdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue, Interval.Lower)) && dist.CDF(ddouble.MaxValue, Interval.Lower) <= 1d, "cdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue / 2, Interval.Lower)) && dist.CDF(ddouble.MaxValue / 2, Interval.Lower) <= 1d, "cdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.CDF(ddouble.NaN, Interval.Lower)), "cdf(NaN)");

                Assert.IsTrue(dist.CDF(ddouble.NegativeInfinity, Interval.Upper) == 1d, "ccdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue, Interval.Upper)) && dist.CDF(ddouble.MinValue, Interval.Upper) <= 1d, "ccdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue / 2, Interval.Upper)) && dist.CDF(ddouble.MinValue / 2, Interval.Upper) <= 1d, "ccdf(-lval / 2)");

                Assert.IsTrue(dist.CDF(ddouble.PositiveInfinity, Interval.Upper) == 0d, "cdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue, Interval.Upper)) && dist.CDF(ddouble.MaxValue, Interval.Upper) >= 0d, "ccdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue / 2, Interval.Upper)) && dist.CDF(ddouble.MaxValue / 2, Interval.Upper) >= 0d, "ccdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.CDF(ddouble.NaN, Interval.Upper)), "ccdf(NaN)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(0d, Interval.Lower)), "quantile(0)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(1d, Interval.Lower)), "quantile(1)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(0d, Interval.Upper)), "cquantile(0)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(1d, Interval.Upper)), "cquantile(1)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.Epsilon, Interval.Lower)), "quantile(0+eps)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)), "quantile(1-eps)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.Epsilon, Interval.Upper)), "cquantile(0+eps)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)), "cquantile(1-eps)");
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist_a0b1c4d5 = [
                0.000000000000000000e+00,
                6.250000000000000000e-02,
                1.250000000000000000e-01,
                1.875000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                1.875000000000000555e-01,
                1.250000000000000000e-01,
                6.250000000000006939e-02,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a1b2c3d6 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                8.333333333333334259e-02,
                1.666666666666666852e-01,
                2.500000000000000000e-01,
                3.333333333333333703e-01,
                3.333333333333333703e-01,
                3.333333333333333703e-01,
                3.333333333333333703e-01,
                3.333333333333333703e-01,
                3.055555555555555802e-01,
                2.777777777777777901e-01,
                2.500000000000000000e-01,
                2.222222222222222654e-01,
                1.944444444444444753e-01,
                1.666666666666667129e-01,
                1.388888888888888951e-01,
                1.111111111111110911e-01,
                8.333333333333335646e-02,
                5.555555555555554553e-02,
                2.777777777777780746e-02,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a2b3c6d10 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                4.545454545454545581e-02,
                9.090909090909091161e-02,
                1.363636363636363535e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.818181818181818232e-01,
                1.704545454545454697e-01,
                1.590909090909090884e-01,
                1.477272727272727348e-01,
                1.363636363636363535e-01,
                1.250000000000000000e-01,
                1.136363636363636465e-01,
                1.022727272727272790e-01,
                9.090909090909091161e-02,
                7.954545454545454419e-02,
                6.818181818181817677e-02,
                5.681818181818182323e-02,
                4.545454545454545581e-02,
                3.409090909090908839e-02,
                2.272727272727272790e-02,
                1.136363636363636395e-02,
                0.000000000000000000e+00,
            ];

            foreach ((TrapezoidDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b1c4d5 , expected_dist_a0b1c4d5 ),
                (dist_a1b2c3d6 , expected_dist_a1b2c3d6 ),
                (dist_a2b3c6d10, expected_dist_a2b3c6d10),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.25, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist_a0b1c4d5 = [
                0.000000000000000000e+00,
                7.812500000000001735e-03,
                3.125000000000000694e-02,
                7.031249999999998612e-02,
                1.250000000000000000e-01,
                1.874999999999999722e-01,
                2.499999999999999722e-01,
                3.124999999999999445e-01,
                3.750000000000000555e-01,
                4.374999999999999445e-01,
                5.000000000000000000e-01,
                5.625000000000000000e-01,
                6.250000000000000000e-01,
                6.875000000000000000e-01,
                7.499999999999998890e-01,
                8.125000000000000000e-01,
                8.750000000000000000e-01,
                9.296875000000000000e-01,
                9.687500000000000000e-01,
                9.921875000000000000e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a1b2c3d6 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.041666666666666956e-02,
                4.166666666666667823e-02,
                9.375000000000000000e-02,
                1.666666666666666852e-01,
                2.500000000000000000e-01,
                3.333333333333333148e-01,
                4.166666666666666297e-01,
                5.000000000000001110e-01,
                5.798611111111109384e-01,
                6.527777777777776791e-01,
                7.187500000000000000e-01,
                7.777777777777776791e-01,
                8.298611111111111605e-01,
                8.750000000000000000e-01,
                9.131944444444444198e-01,
                9.444444444444444198e-01,
                9.687500000000000000e-01,
                9.861111111111111605e-01,
                9.965277777777777901e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a2b3c6d10 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                5.681818181818181976e-03,
                2.272727272727272790e-02,
                5.113636363636363952e-02,
                9.090909090909091161e-02,
                1.363636363636363535e-01,
                1.818181818181818232e-01,
                2.272727272727272652e-01,
                2.727272727272727071e-01,
                3.181818181818181768e-01,
                3.636363636363636465e-01,
                4.090909090909091161e-01,
                4.545454545454545303e-01,
                5.000000000000000000e-01,
                5.454545454545454142e-01,
                5.909090909090909394e-01,
                6.363636363636363535e-01,
                6.803977272727272929e-01,
                7.215909090909091717e-01,
                7.599431818181818787e-01,
                7.954545454545454142e-01,
                8.281250000000000000e-01,
                8.579545454545454142e-01,
                8.849431818181818787e-01,
                9.090909090909090606e-01,
                9.303977272727272929e-01,
                9.488636363636363535e-01,
                9.644886363636363535e-01,
                9.772727272727272929e-01,
                9.872159090909090606e-01,
                9.943181818181817677e-01,
                9.985795454545454142e-01,
                1.000000000000000000e+00,
            ];

            foreach ((TrapezoidDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b1c4d5 , expected_dist_a0b1c4d5 ),
                (dist_a1b2c3d6 , expected_dist_a1b2c3d6 ),
                (dist_a2b3c6d10, expected_dist_a2b3c6d10),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.25, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}