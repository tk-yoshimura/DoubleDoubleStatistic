using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class PoissonDistributionTests {
        readonly PoissonDistribution dist_lambda1 = new(lambda: 1);
        readonly PoissonDistribution dist_lambda2 = new(lambda: 2);
        readonly PoissonDistribution dist_lambda4 = new(lambda: 4);

        PoissonDistribution[] Dists => [
            dist_lambda1,
            dist_lambda2,
            dist_lambda4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (PoissonDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Lambda={dist.Lambda}");
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
            Assert.AreEqual(1, (double)dist_lambda1.Mean, 1e-10);
            Assert.AreEqual(2.0, (double)dist_lambda2.Mean, 1e-10);
            Assert.AreEqual(4.0, (double)dist_lambda4.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(1, (double)dist_lambda1.Variance, 1e-10);
            Assert.AreEqual(2.0, (double)dist_lambda2.Variance, 1e-10);
            Assert.AreEqual(4.0, (double)dist_lambda4.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(1, (double)dist_lambda1.Skewness, 1e-10);
            Assert.AreEqual(0.7071067811865476, (double)dist_lambda2.Skewness, 1e-10);
            Assert.AreEqual(0.5, (double)dist_lambda4.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(1, (double)dist_lambda1.Kurtosis, 1e-10);
            Assert.AreEqual(0.5, (double)dist_lambda2.Kurtosis, 1e-10);
            Assert.AreEqual(0.25, (double)dist_lambda4.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.3048422422562516, (double)dist_lambda1.Entropy, 1e-10);
            Assert.AreEqual(1.7048826439329838, (double)dist_lambda2.Entropy, 1e-10);
            Assert.AreEqual(2.086672699880963, (double)dist_lambda4.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (PoissonDistribution dist in Dists) {
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

            foreach (PoissonDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 1000000).ToArray();

                for (int i = -1; i <= 50; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, (double)dist.PMF(i) * 0.2 + 1e-5, $"{dist},{i}");
                }

                Assert.AreEqual((double)dist.Mean, samples.Average(), 0.05, $"{dist},mean");
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (PoissonDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 100000).ToArray();

                PoissonDistribution? dist_fit = PoissonDistribution.Fit(samples);

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);

                Assert.AreEqual((double)dist.Lambda, (double)dist_fit.Lambda, 0.05, $"{dist},p");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_lambda1 = [
                0.000000000000000000e+00,
                3.678794411714423340e-01,
                3.678794411714423340e-01,
                1.839397205857211393e-01,
                6.131324019524039132e-02,
                1.532831004881010130e-02,
                3.065662009762020000e-03,
                5.109436682936697831e-04,
                7.299195261338138597e-05,
                9.123994076672671552e-06,
                1.013777119630298016e-06,
                1.013777119630298677e-07,
                9.216155633002697802e-09,
                7.680129694168931034e-10,
                5.907792072437641444e-11,
                4.219851480312585337e-12,
                2.813234320208388979e-13,
            ];
            ddouble[] expected_dist_lambda2 = [
                0.000000000000000000e+00,
                1.353352832366127023e-01,
                2.706705664732254046e-01,
                2.706705664732254046e-01,
                1.804470443154835568e-01,
                9.022352215774177842e-02,
                3.608940886309672247e-02,
                1.202980295436556490e-02,
                3.437086558390161587e-03,
                8.592716395975401799e-04,
                1.909492532438982274e-04,
                3.818985064877954112e-05,
                6.943609208869006937e-06,
                1.157268201478170505e-06,
                1.780412617658724668e-07,
                2.543446596655316133e-08,
                3.391262128873747593e-09,
            ];
            ddouble[] expected_dist_lambda4 = [
                0.000000000000000000e+00,
                1.831563888873417867e-02,
                7.326255555493672855e-02,
                1.465251111098734293e-01,
                1.953668148131645355e-01,
                1.953668148131645355e-01,
                1.562934518505317005e-01,
                1.041956345670210227e-01,
                5.954036260972629668e-02,
                2.977018130486314487e-02,
                1.323119169105030531e-02,
                5.292476676420116918e-03,
                1.924536973243678919e-03,
                6.415123244145600010e-04,
                1.973884075121724456e-04,
                5.639668786062049352e-05,
                1.503911676283212788e-05,
            ];

            foreach ((PoissonDistribution dist, ddouble[] expecteds) in new[]{
                (dist_lambda1, expected_dist_lambda1),
                (dist_lambda2, expected_dist_lambda2),
                (dist_lambda4, expected_dist_lambda4),
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