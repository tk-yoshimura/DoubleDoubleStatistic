using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class GeometricDistributionTests {
        readonly GeometricDistribution dist_p25 = new(p: 0.25);
        readonly GeometricDistribution dist_p50 = new(p: 0.50);
        readonly GeometricDistribution dist_p875 = new(p: 0.875);
        readonly GeometricDistribution dist_p999 = new(p: 0.999);
        readonly GeometricDistribution dist_p001 = new(p: 0.001);

        GeometricDistribution[] Dists => [
            dist_p25,
            dist_p50,
            dist_p875,
            dist_p999,
            dist_p001,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (GeometricDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"P={dist.P}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            Assert.AreEqual(3.0, (double)dist_p25.Mean, 1e-10);
            Assert.AreEqual(1.0, (double)dist_p50.Mean, 1e-10);
            Assert.AreEqual(0.1428571428571428, (double)dist_p875.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(12.0, (double)dist_p25.Variance, 1e-10);
            Assert.AreEqual(2.0, (double)dist_p50.Variance, 1e-10);
            Assert.AreEqual(0.16326530612244897, (double)dist_p875.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(2.0207259421636903, (double)dist_p25.Skewness, 1e-10);
            Assert.AreEqual(2.1213203435596424, (double)dist_p50.Skewness, 1e-10);
            Assert.AreEqual(3.181980515339464, (double)dist_p875.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(6.083333333333333, (double)dist_p25.Kurtosis, 1e-10);
            Assert.AreEqual(6.5, (double)dist_p50.Kurtosis, 1e-10);
            Assert.AreEqual(12.125, (double)dist_p875.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(2.249340578475233, (double)dist_p25.Entropy, 1e-10);
            Assert.AreEqual(1.3862943611198906, (double)dist_p50.Entropy, 1e-10);
            Assert.AreEqual(0.4305944700073563, (double)dist_p875.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (GeometricDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int x = -1; x <= 5; x++) {
                    ddouble pdf = dist.PMF(x);

                    Console.WriteLine($"pmf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (GeometricDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 1000000).ToArray();

                for (int i = -1; i <= 50; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, (double)dist.PMF(i) * 0.2 + 1e-5, $"{dist},{i}");
                }

                Assert.AreEqual((double)dist.Mean, samples.Average(), (double)dist.Mean * 0.01, $"{dist},mean");
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (GeometricDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 100000).ToArray();

                GeometricDistribution? dist_fit = GeometricDistribution.Fit(samples);

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);

                Assert.AreEqual((double)dist.P, (double)dist_fit.P, 0.05, $"{dist},p");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_p25 = [
                0.000000000000000000e+00,
                2.500000000000000000e-01,
                1.875000000000000000e-01,
                1.406250000000000000e-01,
                1.054687500000000000e-01,
                7.910156250000000000e-02,
                5.932617187500000000e-02,
                4.449462890625000000e-02,
                3.337097167968750000e-02,
                2.502822875976562500e-02,
                1.877117156982421875e-02,
                1.407837867736816406e-02,
                1.055878400802612305e-02,
                7.919088006019592285e-03,
                5.939316004514694214e-03,
                4.454487003386020660e-03,
            ];
            ddouble[] expected_dist_p50 = [
                0.000000000000000000e+00,
                5.000000000000000000e-01,
                2.500000000000000000e-01,
                1.250000000000000000e-01,
                6.250000000000000000e-02,
                3.125000000000000000e-02,
                1.562500000000000000e-02,
                7.812500000000000000e-03,
                3.906250000000000000e-03,
                1.953125000000000000e-03,
                9.765625000000000000e-04,
                4.882812500000000000e-04,
                2.441406250000000000e-04,
                1.220703125000000000e-04,
                6.103515625000000000e-05,
                3.051757812500000000e-05,
            ];
            ddouble[] expected_dist_p875 = [
                0.000000000000000000e+00,
                8.750000000000000000e-01,
                1.093750000000000000e-01,
                1.367187500000000000e-02,
                1.708984375000000000e-03,
                2.136230468750000000e-04,
                2.670288085937500000e-05,
                3.337860107421875000e-06,
                4.172325134277343750e-07,
                5.215406417846679688e-08,
                6.519258022308349609e-09,
                8.149072527885437012e-10,
                1.018634065985679626e-10,
                1.273292582482099533e-11,
                1.591615728102624416e-12,
                1.989519660128280520e-13,
            ];

            foreach ((GeometricDistribution dist, ddouble[] expecteds) in new[]{
                (dist_p25, expected_dist_p25),
                (dist_p50, expected_dist_p50),
                (dist_p875, expected_dist_p875),
            }) {
                for ((int x, int i) = (-1, 0); i < expecteds.Length; x++, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PMF(x);

                    Console.WriteLine($"{dist} pmf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pmf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}