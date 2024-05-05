using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class BinomialDistributionTests {
        readonly BinomialDistribution dist_n8p25 = new(n: 8, p: 0.25);
        readonly BinomialDistribution dist_n6p50 = new(n: 6, p: 0.50);
        readonly BinomialDistribution dist_n10p875 = new(n: 10, p: 0.875);

        BinomialDistribution[] Dists => [
            dist_n8p25,
            dist_n6p50,
            dist_n10p875,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BinomialDistribution dist in Dists) {
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
            Assert.AreEqual(2.0, (double)dist_n8p25.Mean, 1e-10);
            Assert.AreEqual(3.0, (double)dist_n6p50.Mean, 1e-10);
            Assert.AreEqual(8.75, (double)dist_n10p875.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(1.5, (double)dist_n8p25.Variance, 1e-10);
            Assert.AreEqual(1.5, (double)dist_n6p50.Variance, 1e-10);
            Assert.AreEqual(1.09375, (double)dist_n10p875.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0.4082482904638631, (double)dist_n8p25.Skewness, 1e-10);
            Assert.AreEqual(0, (double)dist_n6p50.Skewness, 1e-10);
            Assert.AreEqual(-0.7171371656006361, (double)dist_n10p875.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(-0.08333333333333333, (double)dist_n8p25.Kurtosis, 1e-10);
            Assert.AreEqual(-0.3333333333333333, (double)dist_n6p50.Kurtosis, 1e-10);
            Assert.AreEqual(0.3142857142857143, (double)dist_n10p875.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.5958372941143697, (double)dist_n8p25.Entropy, 1e-10);
            Assert.AreEqual(1.6173633156271288, (double)dist_n6p50.Entropy, 1e-10);
            Assert.AreEqual(1.3944699335258353, (double)dist_n10p875.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (BinomialDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int x = -1; x <= 5; x++) {
                    ddouble pdf = dist.PMF(x);

                    Console.WriteLine($"pmf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist_n8p25 = [
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
            ddouble[] expected_dist_n6p50 = [
                0.000000000000000000e+00,
                1.562500000000012143e-02,
                9.375000000000002776e-02,
                2.343750000000000555e-01,
                3.124999999999999445e-01,
                2.343750000000000278e-01,
                9.375000000000002776e-02,
                1.562500000000000000e-02,
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
            ddouble[] expected_dist_n10p875 = [
                0.000000000000000000e+00,
                9.313225746154772749e-10,
                6.519258022308344315e-08,
                2.053566277027131398e-06,
                3.833323717117307537e-05,
                4.695821553468702597e-04,
                3.944490104913712415e-03,
                2.300952561199662988e-02,
                9.203810244798654727e-02,
                2.416000189259648601e-01,
                3.758222516626119614e-01,
                2.630755761638283730e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];

            foreach ((BinomialDistribution dist, ddouble[] expecteds) in new[]{
                (dist_n8p25, expected_dist_n8p25),
                (dist_n6p50, expected_dist_n6p50),
                (dist_n10p875, expected_dist_n10p875),
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