using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class CategoricalDistributionTests {
        static readonly BinomialDistribution dist_n8p25 = new(n: 8, p: 0.25);
        static readonly CategoricalDistribution dist;

        static CategoricalDistributionTests() {
            dist = new CategoricalDistribution(new ddouble[9].Select((_, k) => dist_n8p25.PMF(k)));
        }

        CategoricalDistribution[] Dists => [
            dist,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (CategoricalDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"P={dist.Probs}");
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
            Assert.AreEqual(2.0, (double)dist.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(1.5, (double)dist.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0.4082482904638631, (double)dist.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(-0.08333333333333333, (double)dist.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.5958372941143697, (double)dist.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (CategoricalDistribution dist in Dists) {
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

            foreach (CategoricalDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 1000000).ToArray();

                for (int i = -1; i <= dist.N + 1; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, (double)dist.PMF(i) * 0.1 + 1e-5, $"{dist},{i}");
                }

                Assert.AreEqual((double)dist.Mean, samples.Average(), 0.05, $"{dist},mean");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist = [
                0.000000000000000000e+00,
                1.001129150390625139e-01,
                2.669677734375001665e-01,
                3.114624023437500555e-01,
                2.076416015624999722e-01,
                8.651733398437498612e-02,
                2.307128906249999653e-02,
                3.845214843750001735e-03,
                3.662109375000002168e-04,
                1.525878906250000000e-05,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];

            foreach ((CategoricalDistribution dist, ddouble[] expecteds) in new[]{
                (dist, expected_dist),
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