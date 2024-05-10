using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class ZipfDistributionTests {
        readonly ZipfDistribution dist_rho1 = new(s: 1);
        readonly ZipfDistribution dist_rho2 = new(s: 2);
        readonly ZipfDistribution dist_rho4 = new(s: 4);
        readonly ZipfDistribution dist_rho8 = new(s: 8);

        ZipfDistribution[] Dists => [
            dist_rho1,
            dist_rho2,
            dist_rho4,
            dist_rho8,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (ZipfDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"S={dist.S}");
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
            Assert.AreEqual(1.3684327776202063, (double)dist_rho2.Mean, 1e-10);
            Assert.AreEqual(1.0437788248434832, (double)dist_rho4.Mean, 1e-10);
            Assert.AreEqual(1.0020648164093982, (double)dist_rho8.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(0.06977422469107308, (double)dist_rho4.Variance, 1e-10);
            Assert.AreEqual(0.002194278808777872, (double)dist_rho8.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(12.516969344281573, (double)dist_rho4.Skewness, 1e-10);
            Assert.AreEqual(24.324655746524662, (double)dist_rho8.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(687.096606711965, (double)dist_rho8.Kurtosis, 1e-5);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.6279134826784538, (double)dist_rho1.Entropy, 1e-2);
            Assert.AreEqual(0.1740432179156, (double)dist_rho4.Entropy, 1e-2);
            Assert.AreEqual(0.014724675410817834, (double)dist_rho8.Entropy, 1e-2);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (ZipfDistribution dist in Dists) {
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

            foreach (ZipfDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 1000000).ToArray();

                for (int i = -1; i <= 50; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, (double)dist.PMF(i) * 0.2 + 1e-5, $"{dist},{i}");
                }
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_rho1 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                6.079271018540265414e-01,
                1.519817754635066354e-01,
                6.754745576155850306e-02,
                3.799544386587665884e-02,
                2.431708407416106194e-02,
                1.688686394038962577e-02,
                1.240667554804135778e-02,
                9.498860966469164710e-03,
                7.505272862395389807e-03,
                6.079271018540265484e-03,
                5.024190924413442921e-03,
                4.221715985097406441e-03,
                3.597201786118500114e-03,
                3.101668887010339444e-03,
                2.701898230462340070e-03,
            ];
            ddouble[] expected_dist_rho2 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                8.319073725807076825e-01,
                1.039884215725884603e-01,
                3.081138416965583868e-02,
                1.299855269657355754e-02,
                6.655258980645661730e-03,
                3.851423021206979835e-03,
                2.425385925891275930e-03,
                1.624819087071694692e-03,
                1.141162376653920112e-03,
                8.319073725807077163e-04,
                6.250243219990290686e-04,
                4.814278776508724793e-04,
                3.786560639875774449e-04,
                3.031732407364094912e-04,
                2.464910733572467046e-04,
            ];

            foreach ((ZipfDistribution dist, ddouble[] expecteds) in new[]{
                (dist_rho1, expected_dist_rho1),
                (dist_rho2, expected_dist_rho2),
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