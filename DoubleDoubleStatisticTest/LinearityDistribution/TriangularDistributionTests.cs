using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.LinearityDistribution {
    [TestClass()]
    public class TriangularDistributionTests {
        readonly TriangularDistribution dist_a0b4c1 = new(a: 0, b: 4, c: 1);
        readonly TriangularDistribution dist_a1b3c2 = new(a: 1, b: 3, c: 2);
        readonly TriangularDistribution dist_a2b6c3 = new(a: 2, b: 6, c: 3);

        TriangularDistribution[] Dists => [
            dist_a0b4c1,
            dist_a1b3c2,
            dist_a2b6c3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (TriangularDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"A={dist.A}");
                Console.WriteLine($"B={dist.B}");
                Console.WriteLine($"C={dist.C}");
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
            foreach (TriangularDistribution dist in Dists) {
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
            foreach (TriangularDistribution dist in Dists) {
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
            foreach (TriangularDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (TriangularDistribution dist in Dists) {
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
            foreach (TriangularDistribution dist in Dists) {
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
            foreach (TriangularDistribution dist in Dists) {
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
            foreach (TriangularDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (TriangularDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 7; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (TriangularDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 7; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (TriangularDistribution dist in Dists) {
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
            foreach (TriangularDistribution dist in Dists) {
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
            foreach (TriangularDistribution dist in Dists) {
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
            ddouble[] expected_dist_a0b4c1 = [
                0.000000000000000000e+00,
                6.250000000000000000e-02,
                1.250000000000000000e-01,
                1.875000000000000000e-01,
                2.500000000000000000e-01,
                3.125000000000000000e-01,
                3.750000000000000000e-01,
                4.375000000000000000e-01,
                5.000000000000000000e-01,
                4.791666666666666852e-01,
                4.583333333333333148e-01,
                4.375000000000000000e-01,
                4.166666666666666852e-01,
                3.958333333333333148e-01,
                3.750000000000000000e-01,
                3.541666666666666852e-01,
                3.333333333333333148e-01,
                3.125000000000000000e-01,
                2.916666666666666852e-01,
                2.708333333333333148e-01,
                2.500000000000000000e-01,
                2.291666666666666574e-01,
                2.083333333333333426e-01,
                1.875000000000000000e-01,
                1.666666666666666574e-01,
                1.458333333333333426e-01,
                1.250000000000000000e-01,
                1.041666666666666713e-01,
                8.333333333333332871e-02,
                6.250000000000000000e-02,
                4.166666666666666435e-02,
                2.083333333333333218e-02,
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
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a1b3c2 = [
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
                8.750000000000000000e-01,
                7.500000000000000000e-01,
                6.250000000000000000e-01,
                5.000000000000000000e-01,
                3.750000000000000000e-01,
                2.500000000000000000e-01,
                1.250000000000000000e-01,
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
            ddouble[] expected_dist_a2b6c3 = [
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
                4.791666666666666852e-01,
                4.583333333333333148e-01,
                4.375000000000000000e-01,
                4.166666666666666852e-01,
                3.958333333333333148e-01,
                3.750000000000000000e-01,
                3.541666666666666852e-01,
                3.333333333333333148e-01,
                3.125000000000000000e-01,
                2.916666666666666852e-01,
                2.708333333333333148e-01,
                2.500000000000000000e-01,
                2.291666666666666574e-01,
                2.083333333333333426e-01,
                1.875000000000000000e-01,
                1.666666666666666574e-01,
                1.458333333333333426e-01,
                1.250000000000000000e-01,
                1.041666666666666713e-01,
                8.333333333333332871e-02,
                6.250000000000000000e-02,
                4.166666666666666435e-02,
                2.083333333333333218e-02,
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

            foreach ((TriangularDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b4c1, expected_dist_a0b4c1),
                (dist_a1b3c2, expected_dist_a1b3c2),
                (dist_a2b6c3, expected_dist_a2b6c3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
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
            ddouble[] expected_dist_a0b4c1 = [
                0.000000000000000000e+00,
                3.906250000000000000e-03,
                1.562500000000000000e-02,
                3.515625000000000000e-02,
                6.250000000000000000e-02,
                9.765625000000000000e-02,
                1.406250000000000000e-01,
                1.914062500000000000e-01,
                2.500000000000000000e-01,
                3.111979166666666852e-01,
                3.697916666666666852e-01,
                4.257812500000000000e-01,
                4.791666666666666852e-01,
                5.299479166666666297e-01,
                5.781250000000000000e-01,
                6.236979166666666297e-01,
                6.666666666666666297e-01,
                7.070312500000000000e-01,
                7.447916666666666297e-01,
                7.799479166666666297e-01,
                8.125000000000000000e-01,
                8.424479166666666297e-01,
                8.697916666666666297e-01,
                8.945312500000000000e-01,
                9.166666666666666297e-01,
                9.361979166666666297e-01,
                9.531250000000000000e-01,
                9.674479166666666297e-01,
                9.791666666666666297e-01,
                9.882812500000000000e-01,
                9.947916666666666297e-01,
                9.986979166666666297e-01,
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
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a1b3c2 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                7.812500000000000000e-03,
                3.125000000000000000e-02,
                7.031250000000000000e-02,
                1.250000000000000000e-01,
                1.953125000000000000e-01,
                2.812500000000000000e-01,
                3.828125000000000000e-01,
                5.000000000000000000e-01,
                6.171875000000000000e-01,
                7.187500000000000000e-01,
                8.046875000000000000e-01,
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
            ddouble[] expected_dist_a2b6c3 = [
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
                3.906250000000000000e-03,
                1.562500000000000000e-02,
                3.515625000000000000e-02,
                6.250000000000000000e-02,
                9.765625000000000000e-02,
                1.406250000000000000e-01,
                1.914062500000000000e-01,
                2.500000000000000000e-01,
                3.111979166666666852e-01,
                3.697916666666666852e-01,
                4.257812500000000000e-01,
                4.791666666666666852e-01,
                5.299479166666666297e-01,
                5.781250000000000000e-01,
                6.236979166666666297e-01,
                6.666666666666666297e-01,
                7.070312500000000000e-01,
                7.447916666666666297e-01,
                7.799479166666666297e-01,
                8.125000000000000000e-01,
                8.424479166666666297e-01,
                8.697916666666666297e-01,
                8.945312500000000000e-01,
                9.166666666666666297e-01,
                9.361979166666666297e-01,
                9.531250000000000000e-01,
                9.674479166666666297e-01,
                9.791666666666666297e-01,
                9.882812500000000000e-01,
                9.947916666666666297e-01,
                9.986979166666666297e-01,
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

            foreach ((TriangularDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b4c1, expected_dist_a0b4c1),
                (dist_a1b3c2, expected_dist_a1b3c2),
                (dist_a2b6c3, expected_dist_a2b6c3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
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