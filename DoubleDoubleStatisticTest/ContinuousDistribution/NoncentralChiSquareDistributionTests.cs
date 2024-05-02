using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class NoncentralChiSquareDistributionTests {
        readonly NoncentralChiSquareDistribution dist_nu1lambda1 = new(nu: 1, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu2lambda1 = new(nu: 2, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu3lambda1 = new(nu: 3, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu4lambda1 = new(nu: 4, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu8lambda1 = new(nu: 8, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu16lambda1 = new(nu: 16, lambda: 1);
        readonly NoncentralChiSquareDistribution dist_nu32lambda1 = new(nu: 32, lambda: 1);

        readonly NoncentralChiSquareDistribution dist_nu1lambda2 = new(nu: 1, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu2lambda2 = new(nu: 2, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu3lambda2 = new(nu: 3, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu4lambda2 = new(nu: 4, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu8lambda2 = new(nu: 8, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu16lambda2 = new(nu: 16, lambda: 2);
        readonly NoncentralChiSquareDistribution dist_nu32lambda2 = new(nu: 32, lambda: 2);

        readonly NoncentralChiSquareDistribution dist_nu1lambda4 = new(nu: 1, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu2lambda4 = new(nu: 2, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu3lambda4 = new(nu: 3, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu4lambda4 = new(nu: 4, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu8lambda4 = new(nu: 8, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu16lambda4 = new(nu: 16, lambda: 4);
        readonly NoncentralChiSquareDistribution dist_nu32lambda4 = new(nu: 32, lambda: 4);

        readonly NoncentralChiSquareDistribution dist_nu1lambda8 = new(nu: 1, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu2lambda8 = new(nu: 2, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu3lambda8 = new(nu: 3, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu4lambda8 = new(nu: 4, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu8lambda8 = new(nu: 8, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu16lambda8 = new(nu: 16, lambda: 8);
        readonly NoncentralChiSquareDistribution dist_nu32lambda8 = new(nu: 32, lambda: 8);

        readonly NoncentralChiSquareDistribution dist_nu1lambda16 = new(nu: 1, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu2lambda16 = new(nu: 2, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu3lambda16 = new(nu: 3, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu4lambda16 = new(nu: 4, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu8lambda16 = new(nu: 8, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu16lambda16 = new(nu: 16, lambda: 16);
        readonly NoncentralChiSquareDistribution dist_nu32lambda16 = new(nu: 32, lambda: 16);

        readonly NoncentralChiSquareDistribution dist_nu1lambda32 = new(nu: 1, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu2lambda32 = new(nu: 2, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu3lambda32 = new(nu: 3, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu4lambda32 = new(nu: 4, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu8lambda32 = new(nu: 8, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu16lambda32 = new(nu: 16, lambda: 32);
        readonly NoncentralChiSquareDistribution dist_nu32lambda32 = new(nu: 32, lambda: 32);


        NoncentralChiSquareDistribution[] Dists => [
            dist_nu1lambda1,
            dist_nu2lambda1,
            dist_nu3lambda1,
            dist_nu4lambda1,
            dist_nu8lambda1,
            dist_nu16lambda1,
            dist_nu32lambda1,
            dist_nu1lambda2,
            dist_nu2lambda2,
            dist_nu3lambda2,
            dist_nu4lambda2,
            dist_nu8lambda2,
            dist_nu16lambda2,
            dist_nu32lambda2,
            dist_nu1lambda4,
            dist_nu2lambda4,
            dist_nu3lambda4,
            dist_nu4lambda4,
            dist_nu8lambda4,
            dist_nu16lambda4,
            dist_nu32lambda4,
            dist_nu1lambda8,
            dist_nu2lambda8,
            dist_nu3lambda8,
            dist_nu4lambda8,
            dist_nu8lambda8,
            dist_nu16lambda8,
            dist_nu32lambda8,
            dist_nu1lambda16,
            dist_nu2lambda16,
            dist_nu3lambda16,
            dist_nu4lambda16,
            dist_nu8lambda16,
            dist_nu16lambda16,
            dist_nu32lambda16,
            dist_nu1lambda32,
            dist_nu2lambda32,
            dist_nu3lambda32,
            dist_nu4lambda32,
            dist_nu8lambda32,
            dist_nu16lambda32,
            dist_nu32lambda32,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.Nu}");
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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

            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536 * 2);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);

                Console.WriteLine($"cdf(median)={cdf_median}");

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);
                ddouble ccdf_median = dist.CDF(dist.Median, Interval.Upper);

                Console.WriteLine($"ccdf(median)={ccdf_median}");

                Assert.IsTrue(ddouble.Abs(cdf_median + ccdf_median - 1) < 1e-28);

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-28);
                    }
                }

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-10"; p /= 10) {
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
            foreach (NoncentralChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-10"; p /= 10) {
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"quantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }


        [TestMethod()]
        public void IrregularValueTest() {
            foreach (NoncentralChiSquareDistribution dist in Dists) {
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
            ddouble[] expected_nu1lambda1 = [
            ];
            ddouble[] expected_nu2lambda4 = [
            ];
            ddouble[] expected_nu4lambda2 = [
            ];

            foreach ((NoncentralChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1lambda1, expected_nu1lambda1),
                (dist_nu2lambda4, expected_nu2lambda4),
                (dist_nu4lambda2, expected_nu4lambda2),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else if (ddouble.IsPositiveInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_nu1lambda1 = [
            ];
            ddouble[] expected_nu2lambda4 = [
            ];
            ddouble[] expected_nu4lambda2 = [
            ];

            foreach ((NoncentralChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1lambda1, expected_nu1lambda1),
                (dist_nu2lambda4, expected_nu2lambda4),
                (dist_nu4lambda2, expected_nu4lambda2),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }
    }
}