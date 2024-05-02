using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class NoncentralSnedecorFDistributionTests {
        readonly NoncentralSnedecorFDistribution dist_n_1_m_1_lambda_1 = new(n: 1, m: 1, lambda: 1);
        readonly NoncentralSnedecorFDistribution dist_n_2_m_1_lambda_2 = new(n: 2, m: 1, lambda: 2);
        readonly NoncentralSnedecorFDistribution dist_n_1_m_2_lambda_3 = new(n: 1, m: 2, lambda: 3);
        readonly NoncentralSnedecorFDistribution dist_n_2_m_2_lambda_4 = new(n: 2, m: 2, lambda: 4);
        readonly NoncentralSnedecorFDistribution dist_n_3_m_4_lambda_1 = new(n: 3, m: 4, lambda: 1);
        readonly NoncentralSnedecorFDistribution dist_n_4_m_2_lambda_2 = new(n: 4, m: 2, lambda: 2);
        readonly NoncentralSnedecorFDistribution dist_n_2_m_4_lambda_3 = new(n: 2, m: 4, lambda: 3);
        readonly NoncentralSnedecorFDistribution dist_n_4_m_4_lambda_4 = new(n: 4, m: 4, lambda: 4);
        readonly NoncentralSnedecorFDistribution dist_n_6_m_8_lambda_1 = new(n: 6, m: 8, lambda: 1);
        readonly NoncentralSnedecorFDistribution dist_n_8_m_6_lambda_2 = new(n: 8, m: 6, lambda: 2);
        readonly NoncentralSnedecorFDistribution dist_n_8_m_8_lambda_3 = new(n: 8, m: 8, lambda: 3);
        readonly NoncentralSnedecorFDistribution dist_n_10_m_10_lambda_4 = new(n: 10, m: 10, lambda: 4);

        NoncentralSnedecorFDistribution[] Dists => [
            dist_n_1_m_1_lambda_1,
            dist_n_2_m_1_lambda_2,
            dist_n_1_m_2_lambda_3,
            dist_n_2_m_2_lambda_4,
            dist_n_3_m_4_lambda_1,
            dist_n_4_m_2_lambda_2,
            dist_n_2_m_4_lambda_3,
            dist_n_4_m_4_lambda_4,
            dist_n_6_m_8_lambda_1,
            dist_n_8_m_6_lambda_2,
            dist_n_8_m_8_lambda_3,
            dist_n_10_m_10_lambda_4
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"D1={dist.N}");
                Console.WriteLine($"D2={dist.M}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                //Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            Assert.Inconclusive();

            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536 * 2);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
        public void IrregularValueTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            ddouble[] expected_dist_n_1_m_1_lambda_1 = [
            ];
            ddouble[] expected_dist_n_2_m_1_lambda_2 = [
            ];
            ddouble[] expected_dist_n_1_m_2_lambda_3 = [
            ];
            ddouble[] expected_dist_n_2_m_2_lambda_4 = [
            ];
            ddouble[] expected_dist_n_3_m_4_lambda_1 = [
            ];
            ddouble[] expected_dist_n_4_m_2_lambda_2 = [
            ];
            ddouble[] expected_dist_n_2_m_4_lambda_3 = [
            ];
            ddouble[] expected_dist_n_4_m_4_lambda_4 = [
            ];
            ddouble[] expected_dist_n_6_m_8_lambda_1 = [
            ];
            ddouble[] expected_dist_n_8_m_6_lambda_2 = [
            ];
            ddouble[] expected_dist_n_8_m_8_lambda_3 = [
            ];
            ddouble[] expected_dist_n_10_m_10_lambda_4 = [
            ];

            foreach ((NoncentralSnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_n_1_m_1_lambda_1,   expected_dist_n_1_m_1_lambda_1),
                 (dist_n_2_m_1_lambda_2,   expected_dist_n_2_m_1_lambda_2),
                 (dist_n_1_m_2_lambda_3,   expected_dist_n_1_m_2_lambda_3),
                 (dist_n_2_m_2_lambda_4,   expected_dist_n_2_m_2_lambda_4),
                 (dist_n_3_m_4_lambda_1,   expected_dist_n_3_m_4_lambda_1),
                 (dist_n_4_m_2_lambda_2,   expected_dist_n_4_m_2_lambda_2),
                 (dist_n_2_m_4_lambda_3,   expected_dist_n_2_m_4_lambda_3),
                 (dist_n_4_m_4_lambda_4,   expected_dist_n_4_m_4_lambda_4),
                 (dist_n_6_m_8_lambda_1,   expected_dist_n_6_m_8_lambda_1),
                 (dist_n_8_m_6_lambda_2,   expected_dist_n_8_m_6_lambda_2),
                 (dist_n_8_m_8_lambda_3,   expected_dist_n_8_m_8_lambda_3),
                 (dist_n_10_m_10_lambda_4, expected_dist_n_10_m_10_lambda_4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (ddouble.IsInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else if (expected > 0) {
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
            ddouble[] expected_dist_n_1_m_1_lambda_1 = [
            ];
            ddouble[] expected_dist_n_2_m_1_lambda_2 = [
            ];
            ddouble[] expected_dist_n_1_m_2_lambda_3 = [
            ];
            ddouble[] expected_dist_n_2_m_2_lambda_4 = [
            ];
            ddouble[] expected_dist_n_3_m_4_lambda_1 = [
            ];
            ddouble[] expected_dist_n_4_m_2_lambda_2 = [
            ];
            ddouble[] expected_dist_n_2_m_4_lambda_3 = [
            ];
            ddouble[] expected_dist_n_4_m_4_lambda_4 = [
            ];
            ddouble[] expected_dist_n_6_m_8_lambda_1 = [
            ];
            ddouble[] expected_dist_n_8_m_6_lambda_2 = [
            ];
            ddouble[] expected_dist_n_8_m_8_lambda_3 = [
            ];
            ddouble[] expected_dist_n_10_m_10_lambda_4 = [
            ];

            foreach ((NoncentralSnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_n_1_m_1_lambda_1,   expected_dist_n_1_m_1_lambda_1),
                 (dist_n_2_m_1_lambda_2,   expected_dist_n_2_m_1_lambda_2),
                 (dist_n_1_m_2_lambda_3,   expected_dist_n_1_m_2_lambda_3),
                 (dist_n_2_m_2_lambda_4,   expected_dist_n_2_m_2_lambda_4),
                 (dist_n_3_m_4_lambda_1,   expected_dist_n_3_m_4_lambda_1),
                 (dist_n_4_m_2_lambda_2,   expected_dist_n_4_m_2_lambda_2),
                 (dist_n_2_m_4_lambda_3,   expected_dist_n_2_m_4_lambda_3),
                 (dist_n_4_m_4_lambda_4,   expected_dist_n_4_m_4_lambda_4),
                 (dist_n_6_m_8_lambda_1,   expected_dist_n_6_m_8_lambda_1),
                 (dist_n_8_m_6_lambda_2,   expected_dist_n_8_m_6_lambda_2),
                 (dist_n_8_m_8_lambda_3,   expected_dist_n_8_m_8_lambda_3),
                 (dist_n_10_m_10_lambda_4, expected_dist_n_10_m_10_lambda_4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
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