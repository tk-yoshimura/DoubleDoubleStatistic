using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class YuleSimonDistributionTests {
        readonly YuleSimonDistribution dist_rho5 = new(rho: 5);
        readonly YuleSimonDistribution dist_rho6 = new(rho: 6);
        readonly YuleSimonDistribution dist_rho8 = new(rho: 8);

        YuleSimonDistribution[] Dists => [
            dist_rho5,
            dist_rho6,
            dist_rho8,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (YuleSimonDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Rho={dist.Rho}");
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
            Assert.AreEqual(1.25, (double)dist_rho5.Mean, 1e-10);
            Assert.AreEqual(1.2, (double)dist_rho6.Mean, 1e-10);
            Assert.AreEqual(1.1428571428571428, (double)dist_rho8.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(0.5208333333333334, (double)dist_rho5.Variance, 1e-10);
            Assert.AreEqual(0.36, (double)dist_rho6.Variance, 1e-10);
            Assert.AreEqual(0.21768707482993196, (double)dist_rho8.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(6.235382907247958, (double)dist_rho5.Skewness, 1e-10);
            Assert.AreEqual(5.444444444444445, (double)dist_rho6.Skewness, 1e-10);
            Assert.AreEqual(4.960216729135935, (double)dist_rho8.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(118.8, (double)dist_rho5.Kurtosis, 1e-10);
            Assert.AreEqual(66.222222222222222, (double)dist_rho6.Kurtosis, 1e-10);
            Assert.AreEqual(43.6125, (double)dist_rho8.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(0.6061305592752751, (double)dist_rho5.Entropy, 1e-5);
            Assert.AreEqual(0.527723375738105, (double)dist_rho6.Entropy, 1e-5);
            Assert.AreEqual(0.42364668385976184, (double)dist_rho8.Entropy, 1e-5);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (YuleSimonDistribution dist in Dists) {
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

            foreach (YuleSimonDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 1000000).ToArray();

                for (int i = -1; i <= 50; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, (double)dist.PMF(i) * 0.2 + 1e-5, $"{dist},{i}");
                }

                Assert.AreEqual((double)dist.Mean, samples.Average(), 0.05, $"{dist},mean");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_rho5 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                8.333333333333332593e-01,
                1.190476190476190410e-01,
                2.976190476190476025e-02,
                9.920634920634920084e-03,
                3.968253968253968034e-03,
                1.803751803751803750e-03,
                9.018759018759018751e-04,
                4.856254856254856251e-04,
                2.775002775002775000e-04,
                1.665001665001665000e-04,
                1.040626040626040625e-04,
                6.733462615815557465e-05,
                4.488975077210371643e-05,
                3.071404000196569965e-05,
                2.149982800137598840e-05,
                1.535702000098284982e-05,
                1.116874181889661928e-05,
                8.255156996575761786e-06,
                6.191367747431820492e-06,
                4.705439488048184286e-06,
                3.619568836960140716e-06,
                2.815220206524554360e-06,
                2.211958733697864382e-06,
                1.754312099139685508e-06,
                1.403449679311748025e-06,
                1.131814257509474282e-06,
                9.195990842264480264e-07,
                7.523992507307292024e-07,
                6.196229123664868711e-07,
                5.134018416750859333e-07,
                4.278348680625702346e-07,
            ];
            ddouble[] expected_dist_rho6 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                8.571428571428570953e-01,
                1.071428571428571369e-01,
                2.380952380952380820e-02,
                7.142857142857143501e-03,
                2.597402597402597400e-03,
                1.082251082251082250e-03,
                4.995004995004995001e-04,
                2.497502497502497500e-04,
                1.332001332001332000e-04,
                7.492507492507492501e-05,
                4.407357348533818680e-05,
                2.693385046326222715e-05,
                1.701085292416561946e-05,
                1.105705440070765214e-05,
                7.371369600471767813e-06,
                5.025933818503478170e-06,
                3.496301786785027813e-06,
                2.476547098972728282e-06,
                1.783113911260364464e-06,
                1.303044781305650547e-06,
                9.652183565227042191e-07,
                7.239137673920282702e-07,
                5.491759614698146703e-07,
                4.210349037935244500e-07,
                3.259625061627285552e-07,
                2.546582079396317400e-07,
                2.006398001948611383e-07,
                1.593316060370966686e-07,
                1.274652848296765090e-07,
                1.026803683350168320e-07,
                8.325435270406869349e-08,
            ];
            ddouble[] expected_dist_rho8 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                8.888888888888888395e-01,
                8.888888888888889228e-02,
                1.616161616161616160e-02,
                4.040404040404040401e-03,
                1.243201243201243200e-03,
                4.440004440004440001e-04,
                1.776001776001776000e-04,
                7.770007770007770001e-05,
                3.656474244709539143e-05,
                1.828237122354769571e-05,
                9.622300643972470983e-06,
                5.292265354184858617e-06,
                3.024151630962776474e-06,
                1.786998691023458806e-06,
                1.087738333666453159e-06,
                6.798364585415333036e-07,
                4.350953334665813058e-07,
                2.844854103435338534e-07,
                1.896569402290226307e-07,
                1.286957808696939365e-07,
                8.875571094461650060e-08,
                6.212899766123153454e-08,
                4.409154672732560857e-08,
                3.169079921026528282e-08,
                2.304785397110199409e-08,
                1.694695144933981199e-08,
                1.258916393379521375e-08,
                9.441872950346376398e-09,
                7.145201151613554804e-09,
                5.452916668336629618e-09,
                4.194551283335861492e-09,
            ];

            foreach ((YuleSimonDistribution dist, ddouble[] expecteds) in new[]{
                (dist_rho5, expected_dist_rho5),
                (dist_rho6, expected_dist_rho6),
                (dist_rho8, expected_dist_rho8),
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