using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class BinaryDistributionTests {
        readonly BinaryDistribution dist = new();

        BinaryDistribution[] Dists => [
            dist,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BinaryDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
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
            Assert.AreEqual(0.5, (double)dist.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(0.25, (double)dist.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0, (double)dist.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(-2, (double)dist.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(0.6931471805599453, (double)dist.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (BinaryDistribution dist in Dists) {
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

            foreach (BinaryDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 100000).ToArray();

                for (int i = -1; i < 5; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, 0.01, $"{dist},{i}");
                }

                Assert.AreEqual((double)dist.Mean, samples.Average(), 0.05, $"{dist},mean");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist = [
                0.000000000000000000e+00,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];

            foreach ((BinaryDistribution dist, ddouble[] expecteds) in new[]{
                (dist, expected_dist),
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