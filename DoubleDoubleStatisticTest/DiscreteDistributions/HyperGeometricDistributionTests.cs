using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class HyperGeometricDistributionTests {
        readonly HyperGeometricDistribution dist_n50m5r10 = new(50, 5, 10);
        readonly HyperGeometricDistribution dist_n40m10r20 = new(40, 10, 20);
        readonly HyperGeometricDistribution dist_n30m15r5 = new(30, 15, 5);
        readonly HyperGeometricDistribution dist_n36m12r10 = new(36, 12, 10);

        HyperGeometricDistribution[] Dists => [
            dist_n50m5r10,
            dist_n40m10r20,
            dist_n30m15r5,
            dist_n36m12r10,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (HyperGeometricDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"N={dist.N}");
                Console.WriteLine($"M={dist.M}");
                Console.WriteLine($"R={dist.R}");
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
            Assert.AreEqual(1.0, (double)dist_n50m5r10.Mean, 1e-10);
            Assert.AreEqual(5.0, (double)dist_n40m10r20.Mean, 1e-10);
            Assert.AreEqual(2.5, (double)dist_n30m15r5.Mean, 1e-10);
            Assert.AreEqual(3.333333333333333, (double)dist_n36m12r10.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(0.7346938775510204, (double)dist_n50m5r10.Variance, 1e-10);
            Assert.AreEqual(1.9230769230769231, (double)dist_n40m10r20.Variance, 1e-10);
            Assert.AreEqual(1.07758620689655187, (double)dist_n30m15r5.Variance, 1e-10);
            Assert.AreEqual(1.6507936507936507, (double)dist_n36m12r10.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(0.5833333333333334, (double)dist_n50m5r10.Skewness, 1e-10);
            Assert.AreEqual(0, (double)dist_n40m10r20.Skewness, 1e-10);
            Assert.AreEqual(0, (double)dist_n30m15r5.Skewness, 1e-10);
            Assert.AreEqual(0.12208812274418136, (double)dist_n36m12r10.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(-0.07505910165484633, (double)dist_n50m5r10.Kurtosis, 1e-10);
            Assert.AreEqual(-0.11891891891891893, (double)dist_n40m10r20.Kurtosis, 1e-10);
            Assert.AreEqual(-0.29333333333333333, (double)dist_n30m15r5.Kurtosis, 1e-10);
            Assert.AreEqual(-0.1255656108597285, (double)dist_n36m12r10.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(1.2143173878789477, (double)dist_n50m5r10.Entropy, 1e-10);
            Assert.AreEqual(1.745551479912629, (double)dist_n40m10r20.Entropy, 1e-10);
            Assert.AreEqual(1.4529702623936718, (double)dist_n30m15r5.Entropy, 1e-10);
            Assert.AreEqual(1.6673295785189002, (double)dist_n36m12r10.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (HyperGeometricDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int x = -1; x <= 5; x++) {
                    ddouble pdf = dist.PMF(x);

                    Console.WriteLine($"pmf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_n50m5r10 = [
                0.000000000000000000e+00,
                3.105627820045686049e-01,
                4.313371972285675593e-01,
                2.098397175706545048e-01,
                4.417678264645358288e-02,
                3.964583058015064797e-03,
                1.189374917404519499e-04,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_n40m10r20 = [
                0.000000000000000000e+00,
                2.179598953792502537e-04,
                3.962907188713640089e-03,
                2.823571371958468954e-02,
                1.042549429646203868e-01,
                2.215417537998183306e-01,
                2.835734448637674543e-01,
                2.215417537998183584e-01,
                1.042549429646204007e-01,
                2.823571371958469300e-02,
                3.962907188713640089e-03,
                2.179598953792502537e-04,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_n30m15r5 = [
                0.000000000000000000e+00,
                2.107279693486589431e-02,
                1.436781609195402210e-01,
                3.352490421455938674e-01,
                3.352490421455938119e-01,
                1.436781609195402210e-01,
                2.107279693486589431e-02,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_dist_n36m12r10 = [
                0.000000000000000000e+00,
                7.715804156293587165e-03,
                6.172643325034869732e-02,
                1.909661528682662923e-01,
                2.995547495972804564e-01,
                2.621104058976204132e-01,
                1.324347314009029497e-01,
                3.862679665859669598e-02,
                6.306415780995378723e-03,
                5.374786176984697184e-04,
                2.077212049076211295e-05,
                2.596515061345264437e-07,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];

            foreach ((HyperGeometricDistribution dist, ddouble[] expecteds) in new[]{
                (dist_n50m5r10, expected_dist_n50m5r10),
                (dist_n40m10r20, expected_dist_n40m10r20),
                (dist_n30m15r5, expected_dist_n30m15r5),
                (dist_n36m12r10, expected_dist_n36m12r10),
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