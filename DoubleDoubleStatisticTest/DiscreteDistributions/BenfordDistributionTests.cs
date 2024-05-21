using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class BenfordDistributionTests {
        readonly BenfordDistribution dist_n5 = new(n: 5);
        readonly BenfordDistribution dist_n10 = new(n: 10);

        BenfordDistribution[] Dists => [
            dist_n5,
            dist_n10,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BenfordDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"N={dist.N}");
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
            Assert.AreEqual(2.0253641312938355, (double)dist_n5.Mean, 1e-10);
            Assert.AreEqual(3.440236967123206, (double)dist_n10.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(1.1633676759907585, (double)dist_n5.Variance, 1e-10);
            Assert.AreEqual(6.056512631375666, (double)dist_n10.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0.6126107364114097, (double)dist_n5.Skewness, 1e-10);
            Assert.AreEqual(0.7956043243647771, (double)dist_n10.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(-0.9715952643744011, (double)dist_n5.Kurtosis, 1e-10);
            Assert.AreEqual(-0.548225322238441, (double)dist_n10.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.2918182285868438, (double)dist_n5.Entropy, 1e-10);
            Assert.AreEqual(1.9934331507912042, (double)dist_n10.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (BenfordDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int x = -1; x <= 12; x++) {
                    ddouble pdf = dist.PMF(x);

                    Console.WriteLine($"pmf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (BenfordDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 100000).ToArray();

                for (int i = -1; i < 5; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, 0.01, $"{dist},{i}");
                }

                Assert.AreEqual((double)dist.Mean, samples.Average(), 0.05, $"{dist},mean");
                Assert.AreEqual((double)dist.Variance, samples.Select(n => (double)n).Variance(), 0.05, $"{dist},var");
                Assert.AreEqual((double)dist.Skewness, samples.Select(n => (double)n).Skewness(), 0.05, $"{dist},skew");
                Assert.AreEqual((double)dist.Kurtosis, samples.Select(n => (double)n).Kurtosis(), 0.05, $"{dist},kurt");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_n5 = [
                0, 
                0,
                4.306765580733931e-1,
                2.519296364125923e-1,
                1.787469216608008e-1,
                1.386468838532139e-1,
                0
            ];
            ddouble[] expected_dist_n10 = [
                0,
                0,
                3.010299956639812e-1,
                1.760912590556813e-1,
                1.249387366083000e-1,
                9.691001300805641e-2,
                7.918124604762483e-2,
                6.694678963061320e-2,
                5.799194697768675e-2,
                5.115252244738129e-2,
                4.575749056067512e-2,
                0
            ];

            foreach ((BenfordDistribution dist, ddouble[] expecteds) in new[]{
                (dist_n5, expected_dist_n5),
                (dist_n10, expected_dist_n10),
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