using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class BernoulliDistributionTests {
        readonly BernoulliDistribution dist_p25 = new(p: 0.25);
        readonly BernoulliDistribution dist_p50 = new(p: 0.50);
        readonly BernoulliDistribution dist_p875 = new(p: 0.875);

        BernoulliDistribution[] Dists => [
            dist_p25,
            dist_p50,
            dist_p875,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BernoulliDistribution dist in Dists) {
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
            Assert.AreEqual(0.25, (double)dist_p25.Mean, 1e-10);
            Assert.AreEqual(0.5, (double)dist_p50.Mean, 1e-10);
            Assert.AreEqual(0.875, (double)dist_p875.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(0.1875, (double)dist_p25.Variance, 1e-10);
            Assert.AreEqual(0.25, (double)dist_p50.Variance, 1e-10);
            Assert.AreEqual(0.109375, (double)dist_p875.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(1.1547005383792515, (double)dist_p25.Skewness, 1e-10);
            Assert.AreEqual(0, (double)dist_p50.Skewness, 1e-10);
            Assert.AreEqual(-2.2677868380553634, (double)dist_p875.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(-0.6666666666666665, (double)dist_p25.Kurtosis, 1e-10);
            Assert.AreEqual(-2, (double)dist_p50.Kurtosis, 1e-10);
            Assert.AreEqual(3.1428571428571432, (double)dist_p875.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(0.5623351446188083, (double)dist_p25.Entropy, 1e-10);
            Assert.AreEqual(0.6931471805599453, (double)dist_p50.Entropy, 1e-10);
            Assert.AreEqual(0.37677016125643675, (double)dist_p875.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (BernoulliDistribution dist in Dists) {
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

            foreach (BernoulliDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 100000).ToArray();

                for (int i = -1; i < 5; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, 0.01, $"{dist},{i}");
                }
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_p25 = [
                0.000000000000000000e+00,
                7.500000000000000000e-01,
                2.500000000000000000e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_p50 = [
                0.000000000000000000e+00,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_p875 = [
                0.000000000000000000e+00,
                1.2500000000000000000e-01,
                8.750000000000000000e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];

            foreach ((BernoulliDistribution dist, ddouble[] expecteds) in new[]{
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

                    Assert.AreEqual(expected, actual, $"{dist} pmf({x})\n{expected}\n{actual}");
                }
            }
        }
    }
}