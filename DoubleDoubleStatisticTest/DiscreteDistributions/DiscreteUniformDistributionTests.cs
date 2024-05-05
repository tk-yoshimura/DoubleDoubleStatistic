using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class DiscreteUniformDistributionTests {
        readonly DiscreteUniformDistribution dist_a0b1 = new(a: 0, b: 1);
        readonly DiscreteUniformDistribution dist_a1b3 = new(a: 1, b: 3);
        readonly DiscreteUniformDistribution dist_a2b6 = new(a: 2, b: 6);

        DiscreteUniformDistribution[] Dists => [
            dist_a0b1,
            dist_a1b3,
            dist_a2b6,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (DiscreteUniformDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Min={dist.A}");
                Console.WriteLine($"Max={dist.B}");
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
            Assert.AreEqual(0, (double)dist_a0b1.Mean, 1e-10);
            Assert.AreEqual(1.5, (double)dist_a1b3.Mean, 1e-10);
            Assert.AreEqual(3.5, (double)dist_a2b6.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(0, (double)dist_a0b1.Variance, 1e-10);
            Assert.AreEqual(0.25, (double)dist_a1b3.Variance, 1e-10);
            Assert.AreEqual(1.25, (double)dist_a2b6.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0, (double)dist_a0b1.Skewness, 1e-10);
            Assert.AreEqual(0, (double)dist_a1b3.Skewness, 1e-10);
            Assert.AreEqual(0, (double)dist_a2b6.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.IsTrue(ddouble.IsNegativeInfinity(dist_a0b1.Kurtosis));
            Assert.AreEqual(-2, (double)dist_a1b3.Kurtosis, 1e-10);
            Assert.AreEqual(-1.36, (double)dist_a2b6.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(0, (double)dist_a0b1.Entropy, 1e-10);
            Assert.AreEqual(0.6931471805599453, (double)dist_a1b3.Entropy, 1e-10);
            Assert.AreEqual(1.3862943611198906, (double)dist_a2b6.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (DiscreteUniformDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int x = -1; x <= 5; x++) {
                    ddouble pdf = dist.PMF(x);

                    Console.WriteLine($"pmf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist_a0b1 = [
                0.000000000000000000e+00,
                1.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a1b3 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_a2b6 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                2.500000000000000000e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];

            foreach ((DiscreteUniformDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a0b1, expected_dist_a0b1),
                (dist_a1b3, expected_dist_a1b3),
                (dist_a2b6, expected_dist_a2b6),
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