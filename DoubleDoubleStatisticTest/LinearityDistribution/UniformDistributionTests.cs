using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.LinearityDistribution {
    [TestClass()]
    public class UniformDistributionTests {
        readonly UniformDistribution dist_a0b1 = new(a: 0, b: 1);
        readonly UniformDistribution dist_a1b2 = new(a: 1, b: 2);
        readonly UniformDistribution dist_a2b4 = new(a: 2, b: 4);

        UniformDistribution[] Dists => [
            dist_a0b1,
            dist_a1b2,
            dist_a2b4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (UniformDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Min={dist.A}");
                Console.WriteLine($"Max={dist.B}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (UniformDistribution dist in Dists) {
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
        public void MedianTest() {
            foreach (UniformDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (UniformDistribution dist in Dists) {
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
            foreach (UniformDistribution dist in Dists) {
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
            foreach (UniformDistribution dist in Dists) {
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
            foreach (UniformDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (UniformDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (UniformDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (UniformDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (UniformDistribution dist in Dists) {
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
            foreach (UniformDistribution dist in Dists) {
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
        public void PDFExpectedTest() {
            ddouble[] expected_dist_a0b1 = [
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
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
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a1b2 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
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
            ddouble[] expected_dist_a2b4 = [
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
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
            ];

            foreach ((UniformDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b1, expected_dist_a0b1),
                (dist_a1b2, expected_dist_a1b2),
                (dist_a2b4, expected_dist_a2b4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.AreEqual(expected, actual, $"{dist} pdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist_a0b1 = [
                0.000000000000000000e+00,
                1.250000000000000000e-01,
                2.500000000000000000e-01,
                3.750000000000000000e-01,
                5.000000000000000000e-01,
                6.250000000000000000e-01,
                7.500000000000000000e-01,
                8.750000000000000000e-01,
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
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a1b2 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.250000000000000000e-01,
                2.500000000000000000e-01,
                3.750000000000000000e-01,
                5.000000000000000000e-01,
                6.250000000000000000e-01,
                7.500000000000000000e-01,
                8.750000000000000000e-01,
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
            ddouble[] expected_dist_a2b4 = [
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
                6.250000000000000000e-02,
                1.250000000000000000e-01,
                1.875000000000000000e-01,
                2.500000000000000000e-01,
                3.125000000000000000e-01,
                3.750000000000000000e-01,
                4.375000000000000000e-01,
                5.000000000000000000e-01,
                5.625000000000000000e-01,
                6.250000000000000000e-01,
                6.875000000000000000e-01,
                7.500000000000000000e-01,
                8.125000000000000000e-01,
                8.750000000000000000e-01,
                9.375000000000000000e-01,
                1.000000000000000000e+00,
            ];

            foreach ((UniformDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b1, expected_dist_a0b1),
                (dist_a1b2, expected_dist_a1b2),
                (dist_a2b4, expected_dist_a2b4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.AreEqual(expected, actual, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }
    }
}