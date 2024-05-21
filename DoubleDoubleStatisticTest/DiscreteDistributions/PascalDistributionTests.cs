using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class PascalDistributionTests {
        readonly PascalDistribution dist_n8p25 = new(n: 8, p: 0.25);
        readonly PascalDistribution dist_n6p50 = new(n: 6, p: 0.50);
        readonly PascalDistribution dist_n10p875 = new(n: 10, p: 0.875);

        PascalDistribution[] Dists => [
            dist_n8p25,
            dist_n6p50,
            dist_n10p875,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (PascalDistribution dist in Dists) {
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
            Assert.AreEqual(32.0, (double)dist_n8p25.Mean, 1e-10);
            Assert.AreEqual(12.0, (double)dist_n6p50.Mean, 1e-10);
            Assert.AreEqual(11.4285714285714286, (double)dist_n10p875.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(96.0, (double)dist_n8p25.Variance, 1e-10);
            Assert.AreEqual(12.0, (double)dist_n6p50.Variance, 1e-10);
            Assert.AreEqual(1.6326530612244898, (double)dist_n10p875.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0.7144345083117604, (double)dist_n8p25.Skewness, 1e-10);
            Assert.AreEqual(0.8660254037844387, (double)dist_n6p50.Skewness, 1e-10);
            Assert.AreEqual(1.0062305898749053, (double)dist_n10p875.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(0.7604166666666666, (double)dist_n8p25.Kurtosis, 1e-10);
            Assert.AreEqual(1.0833333333333333, (double)dist_n6p50.Kurtosis, 1e-10);
            Assert.AreEqual(1.2125, (double)dist_n10p875.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(3.6569514456676746, (double)dist_n8p25.Entropy, 1e-10);
            Assert.AreEqual(2.5912916959514245, (double)dist_n6p50.Entropy, 1e-10);
            Assert.AreEqual(1.5476613049845014, (double)dist_n10p875.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (PascalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int x = -1; x <= 25; x++) {
                    ddouble pdf = dist.PMF(x);

                    Console.WriteLine($"pmf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (PascalDistribution dist in Dists) {
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

            foreach (PascalDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 100000).ToArray();

                PascalDistribution? dist_fit = PascalDistribution.Fit(samples);

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);

                Assert.AreEqual(dist.N, dist_fit.N, $"{dist},N");
                Assert.AreEqual((double)dist.P, (double)dist_fit.P, 0.05, $"{dist},p");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_n8p25 = [
                0, 0, 0, 0, 0, 0, 0, 0,
                0.000000000000000000e+00,
                1.525878906249999831e-05,
                9.155273437500004066e-05,
                3.089904785156251084e-04,
                7.724761962890621747e-04,
                1.593232154846191840e-03,
                2.867817878723146700e-03,
                4.660204052925114200e-03,
                6.990306079387660458e-03,
                9.830117924138900148e-03,
                1.310682389885186629e-02,
                1.671120047103613615e-02,
                2.050920057808979305e-02,
                2.435467568648164854e-02,
                2.810154886901724061e-02,
                3.161424247764444556e-02,
                3.477566672540885473e-02,
            ];
            ddouble[] expected_dist_n6p50 = [
                0, 0, 0, 0, 0, 0,
                0.000000000000000000e+00,
                1.562499999999998265e-02,
                4.687500000000001388e-02,
                8.203125000000002776e-02,
                1.093749999999999306e-01,
                1.230468749999999584e-01,
                1.230468750000000000e-01,
                1.127929687499999722e-01,
                9.667968750000004163e-02,
                7.855224609375005551e-02,
                6.109619140624993061e-02,
                4.582214355468752082e-02,
                3.332519531250002082e-02,
                2.360534667968750347e-02,
                1.634216308593752082e-02,
                1.108932495117187673e-02,
                7.392883300781274286e-03,
            ];
            ddouble[] expected_dist_n10p875 = [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0.000000000000000000e+00,
                2.630755761638282619e-01,
                3.288444702047854662e-01,
                2.260805732657900913e-01,
                1.130402866328948930e-01,
                4.592261644461353987e-02,
                1.607291575561475491e-02,
                5.022786173629616331e-03,
                1.435081763894173988e-03,
                3.811935935343899249e-04,
                9.529839838359758966e-05,
                2.263336961610441789e-05,
                5.143947640023734726e-06,
                1.125238546255191442e-06,
                2.380312309385978272e-07,
                4.888141349631930821e-08,
                9.776282699263842452e-09,
            ];

            foreach ((PascalDistribution dist, ddouble[] expecteds) in new[]{
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