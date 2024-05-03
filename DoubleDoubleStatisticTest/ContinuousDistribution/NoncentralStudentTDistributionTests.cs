using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class NoncentralStudentTDistributionTests {
        readonly NoncentralStudentTDistribution dist_nu1mu0 = new(nu: 1, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu2mu0 = new(nu: 2, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu3mu0 = new(nu: 3, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu4mu0 = new(nu: 4, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu5mu0 = new(nu: 5, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu8mu0 = new(nu: 8, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu16mu0 = new(nu: 16, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu32mu0 = new(nu: 32, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu64mu0 = new(nu: 64, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu128mu0 = new(nu: 128, mu: 0);
        readonly NoncentralStudentTDistribution dist_nu129mu0 = new(nu: 129, mu: 0);

        readonly NoncentralStudentTDistribution dist_nu1mu1 = new(nu: 1, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu2mu1 = new(nu: 2, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu3mu1 = new(nu: 3, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu4mu1 = new(nu: 4, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu5mu1 = new(nu: 5, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu8mu1 = new(nu: 8, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu16mu1 = new(nu: 16, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu32mu1 = new(nu: 32, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu64mu1 = new(nu: 64, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu128mu1 = new(nu: 128, mu: 1);
        readonly NoncentralStudentTDistribution dist_nu129mu1 = new(nu: 129, mu: 1);

        readonly NoncentralStudentTDistribution dist_nu1mum1 = new(nu: 1, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu2mum1 = new(nu: 2, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu3mum1 = new(nu: 3, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu4mum1 = new(nu: 4, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu5mum1 = new(nu: 5, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu8mum1 = new(nu: 8, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu16mum1 = new(nu: 16, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu32mum1 = new(nu: 32, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu64mum1 = new(nu: 64, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu128mum1 = new(nu: 128, mu: -1);
        readonly NoncentralStudentTDistribution dist_nu129mum1 = new(nu: 129, mu: -1);

        readonly NoncentralStudentTDistribution dist_nu1mu2 = new(nu: 1, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu2mu2 = new(nu: 2, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu3mu2 = new(nu: 3, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu4mu2 = new(nu: 4, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu5mu2 = new(nu: 5, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu8mu2 = new(nu: 8, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu16mu2 = new(nu: 16, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu32mu2 = new(nu: 32, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu64mu2 = new(nu: 64, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu128mu2 = new(nu: 128, mu: 2);
        readonly NoncentralStudentTDistribution dist_nu129mu2 = new(nu: 129, mu: 2);

        readonly NoncentralStudentTDistribution dist_nu1mu4 = new(nu: 1, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu2mu4 = new(nu: 2, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu3mu4 = new(nu: 3, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu4mu4 = new(nu: 4, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu5mu4 = new(nu: 5, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu8mu4 = new(nu: 8, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu16mu4 = new(nu: 16, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu32mu4 = new(nu: 32, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu64mu4 = new(nu: 64, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu128mu4 = new(nu: 128, mu: 4);
        readonly NoncentralStudentTDistribution dist_nu129mu4 = new(nu: 129, mu: 4);

        readonly NoncentralStudentTDistribution dist_nu1mu16 = new(nu: 1, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu2mu16 = new(nu: 2, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu3mu16 = new(nu: 3, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu4mu16 = new(nu: 4, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu5mu16 = new(nu: 5, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu8mu16 = new(nu: 8, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu16mu16 = new(nu: 16, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu32mu16 = new(nu: 32, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu64mu16 = new(nu: 64, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu128mu16 = new(nu: 128, mu: 16);
        readonly NoncentralStudentTDistribution dist_nu129mu16 = new(nu: 129, mu: 16);


        NoncentralStudentTDistribution[] Dists => [
            dist_nu1mu0,
            dist_nu2mu0,
            dist_nu3mu0,
            dist_nu4mu0,
            dist_nu5mu0,
            dist_nu8mu0,
            dist_nu16mu0,
            dist_nu32mu0,
            dist_nu64mu0,
            dist_nu128mu0,
            dist_nu129mu0,

            dist_nu1mu1,
            dist_nu2mu1,
            dist_nu3mu1,
            dist_nu4mu1,
            dist_nu5mu1,
            dist_nu8mu1,
            dist_nu16mu1,
            dist_nu32mu1,
            dist_nu64mu1,
            dist_nu128mu1,
            dist_nu129mu1,

            dist_nu1mum1,
            dist_nu2mum1,
            dist_nu3mum1,
            dist_nu4mum1,
            dist_nu5mum1,
            dist_nu8mum1,
            dist_nu16mum1,
            dist_nu32mum1,
            dist_nu64mum1,
            dist_nu128mum1,
            dist_nu129mum1,

            dist_nu1mu2,
            dist_nu2mu2,
            dist_nu3mu2,
            dist_nu4mu2,
            dist_nu5mu2,
            dist_nu8mu2,
            dist_nu16mu2,
            dist_nu32mu2,
            dist_nu64mu2,
            dist_nu128mu2,
            dist_nu129mu2,

            dist_nu1mu4,
            dist_nu2mu4,
            dist_nu3mu4,
            dist_nu4mu4,
            dist_nu5mu4,
            dist_nu8mu4,
            dist_nu16mu4,
            dist_nu32mu4,
            dist_nu64mu4,
            dist_nu128mu4,
            dist_nu129mu4,

            dist_nu1mu16,
            dist_nu2mu16,
            dist_nu3mu16,
            dist_nu4mu16,
            dist_nu5mu16,
            dist_nu8mu16,
            dist_nu16mu16,
            dist_nu32mu16,
            dist_nu64mu16,
            dist_nu128mu16,
            dist_nu129mu16,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Nu={dist.Nu}");
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
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            Assert.Inconclusive();

            foreach (NoncentralStudentTDistribution dist in Dists) {
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
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-10, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 8192);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-24);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                //ignore
                if (ddouble.Abs(dist.Mu) >= 8d) {
                    continue;
                }

                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-24);
                    }
                }
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (NoncentralStudentTDistribution dist in Dists) {
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
            ddouble[] expected_nu1mu1 = [
            ];
            ddouble[] expected_nu2mum1 = [
            ];
            ddouble[] expected_nu4mu2 = [
            ];
            ddouble[] expected_nu5mu4 = [
            ];
            ddouble[] expected_nu8mu16 = [
            ];

            foreach ((NoncentralStudentTDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1mu1, expected_nu1mu1),
                (dist_nu2mum1, expected_nu2mum1),
                (dist_nu4mu2, expected_nu4mu2),
                (dist_nu5mu4, expected_nu5mu4),
                (dist_nu8mu16,  expected_nu8mu16)
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_nu1mu1 = [
            ];
            ddouble[] expected_nu2mum1 = [
            ];
            ddouble[] expected_nu4mu2 = [
            ];
            ddouble[] expected_nu5mu4 = [
            ];
            ddouble[] expected_nu8mu16 = [
            ];

            foreach ((NoncentralStudentTDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1mu1, expected_nu1mu1),
                (dist_nu2mum1, expected_nu2mum1),
                (dist_nu4mu2, expected_nu4mu2),
                (dist_nu5mu4, expected_nu5mu4),
                (dist_nu8mu16,  expected_nu8mu16)
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }
    }
}