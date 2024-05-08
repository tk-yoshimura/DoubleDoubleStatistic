using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class PowerDistributionTests {
        readonly PowerDistribution dist_k1a1 = new(k: 1, alpha: 1);
        readonly PowerDistribution dist_k1a2 = new(k: 1, alpha: 2);
        readonly PowerDistribution dist_k2a1 = new(k: 2, alpha: 1);
        readonly PowerDistribution dist_k2a2 = new(k: 2, alpha: 2);
        readonly PowerDistribution dist_k4a3 = new(k: 4, alpha: 3);

        PowerDistribution[] Dists => [
            dist_k1a1,
            dist_k1a2,
            dist_k2a1,
            dist_k2a2,
            dist_k4a3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (PowerDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"K={dist.K}");
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
            foreach (PowerDistribution dist in Dists) {
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
            foreach (PowerDistribution dist in Dists) {
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
            foreach (PowerDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (PowerDistribution dist in Dists) {
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
            foreach (PowerDistribution dist in Dists) {
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
            foreach (PowerDistribution dist in Dists) {
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
            foreach (PowerDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (PowerDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (PowerDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (PowerDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (PowerDistribution dist in Dists) {
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
            foreach (PowerDistribution dist in Dists) {
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

            foreach (PowerDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 5) * 0.1, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (PowerDistribution dist in Dists) {
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
            ddouble[] expected_dist_k1a1 = [
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ];
            ddouble[] expected_dist_k1a2 = [
                0.0,
                0.03125,
                0.0625,
                0.09375,
                0.125,
                0.15625,
                0.1875,
                0.21875,
                0.25,
                0.28125,
                0.3125,
                0.34375,
                0.375,
                0.40625,
                0.4375,
                0.46875,
                0.5,
                0.53125,
                0.5625,
                0.59375,
                0.625,
                0.65625,
                0.6875,
                0.71875,
                0.75,
                0.78125,
                0.8125,
                0.84375,
                0.875,
                0.90625,
                0.9375,
                0.96875,
                1.0,
                1.03125,
                1.0625,
                1.09375,
                1.125,
                1.15625,
                1.1875,
                1.21875,
                1.25,
                1.28125,
                1.3125,
                1.34375,
                1.375,
                1.40625,
                1.4375,
                1.46875,
                1.5,
                1.53125,
                1.5625,
                1.59375,
                1.625,
                1.65625,
                1.6875,
                1.71875,
                1.75,
                1.78125,
                1.8125,
                1.84375,
                1.875,
                1.90625,
                1.9375,
                1.96875,
                2.0,
            ];
            ddouble[] expected_dist_k2a1 = [
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ];
            ddouble[] expected_dist_k2a2 = [
                0.0,
                0.125,
                0.25,
                0.375,
                0.5,
                0.625,
                0.75,
                0.875,
                1.0,
                1.125,
                1.25,
                1.375,
                1.5,
                1.625,
                1.75,
                1.875,
                2.0,
                2.125,
                2.25,
                2.375,
                2.5,
                2.625,
                2.75,
                2.875,
                3.0,
                3.125,
                3.25,
                3.375,
                3.5,
                3.625,
                3.75,
                3.875,
                4.0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ];
            ddouble[] expected_dist_k4a3 = [
                0.0,
                0.046875,
                0.1875,
                0.421875,
                0.75,
                1.171875,
                1.6875,
                2.296875,
                3.0,
                3.796875,
                4.6875,
                5.671875,
                6.75,
                7.921875,
                9.1875,
                10.546875,
                12.0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ];

            foreach ((PowerDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1a1, expected_dist_k1a1),
                (dist_k1a2, expected_dist_k1a2),
                (dist_k2a1, expected_dist_k2a1),
                (dist_k2a2, expected_dist_k2a2),
                (dist_k4a3, expected_dist_k4a3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 64, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist_k1a1 = [
                0.0,
                0.015625,
                0.03125,
                0.046875,
                0.0625,
                0.078125,
                0.09375,
                0.109375,
                0.125,
                0.140625,
                0.15625,
                0.171875,
                0.1875,
                0.203125,
                0.21875,
                0.234375,
                0.25,
                0.265625,
                0.28125,
                0.296875,
                0.3125,
                0.328125,
                0.34375,
                0.359375,
                0.375,
                0.390625,
                0.40625,
                0.421875,
                0.4375,
                0.453125,
                0.46875,
                0.484375,
                0.5,
                0.515625,
                0.53125,
                0.546875,
                0.5625,
                0.578125,
                0.59375,
                0.609375,
                0.625,
                0.640625,
                0.65625,
                0.671875,
                0.6875,
                0.703125,
                0.71875,
                0.734375,
                0.75,
                0.765625,
                0.78125,
                0.796875,
                0.8125,
                0.828125,
                0.84375,
                0.859375,
                0.875,
                0.890625,
                0.90625,
                0.921875,
                0.9375,
                0.953125,
                0.96875,
                0.984375,
                1.0,
            ];
            ddouble[] expected_dist_k1a2 = [
                0.0,
                2.44140625e-4,
                9.765625e-4,
                0.002197265625,
                0.00390625,
                0.006103515625,
                0.0087890625,
                0.011962890625,
                0.015625,
                0.019775390625,
                0.0244140625,
                0.029541015625,
                0.03515625,
                0.041259765625,
                0.0478515625,
                0.054931640625,
                0.0625,
                0.070556640625,
                0.0791015625,
                0.088134765625,
                0.09765625,
                0.107666015625,
                0.1181640625,
                0.129150390625,
                0.140625,
                0.152587890625,
                0.1650390625,
                0.177978515625,
                0.19140625,
                0.205322265625,
                0.2197265625,
                0.234619140625,
                0.25,
                0.265869140625,
                0.2822265625,
                0.299072265625,
                0.31640625,
                0.334228515625,
                0.3525390625,
                0.371337890625,
                0.390625,
                0.410400390625,
                0.4306640625,
                0.451416015625,
                0.47265625,
                0.494384765625,
                0.5166015625,
                0.539306640625,
                0.5625,
                0.586181640625,
                0.6103515625,
                0.635009765625,
                0.66015625,
                0.685791015625,
                0.7119140625,
                0.738525390625,
                0.765625,
                0.793212890625,
                0.8212890625,
                0.849853515625,
                0.87890625,
                0.908447265625,
                0.9384765625,
                0.968994140625,
                1.0,
            ];
            ddouble[] expected_dist_k2a1 = [
                0.0,
                0.03125,
                0.0625,
                0.09375,
                0.125,
                0.15625,
                0.1875,
                0.21875,
                0.25,
                0.28125,
                0.3125,
                0.34375,
                0.375,
                0.40625,
                0.4375,
                0.46875,
                0.5,
                0.53125,
                0.5625,
                0.59375,
                0.625,
                0.65625,
                0.6875,
                0.71875,
                0.75,
                0.78125,
                0.8125,
                0.84375,
                0.875,
                0.90625,
                0.9375,
                0.96875,
                1.0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ];
            ddouble[] expected_dist_k2a2 = [
                0.0,
                9.765625e-4,
                0.00390625,
                0.0087890625,
                0.015625,
                0.0244140625,
                0.03515625,
                0.0478515625,
                0.0625,
                0.0791015625,
                0.09765625,
                0.1181640625,
                0.140625,
                0.1650390625,
                0.19140625,
                0.2197265625,
                0.25,
                0.2822265625,
                0.31640625,
                0.3525390625,
                0.390625,
                0.4306640625,
                0.47265625,
                0.5166015625,
                0.5625,
                0.6103515625,
                0.66015625,
                0.7119140625,
                0.765625,
                0.8212890625,
                0.87890625,
                0.9384765625,
                1.0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ];
            ddouble[] expected_dist_k4a3 = [
                0.0,
                2.44140625e-4,
                0.001953125,
                0.006591796875,
                0.015625,
                0.030517578125,
                0.052734375,
                0.083740234375,
                0.125,
                0.177978515625,
                0.244140625,
                0.324951171875,
                0.421875,
                0.536376953125,
                0.669921875,
                0.823974609375,
                1.0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ];

            foreach ((PowerDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1a1, expected_dist_k1a1),
                (dist_k1a2, expected_dist_k1a2),
                (dist_k2a1, expected_dist_k2a1),
                (dist_k2a2, expected_dist_k2a2),
                (dist_k4a3, expected_dist_k4a3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 64, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}