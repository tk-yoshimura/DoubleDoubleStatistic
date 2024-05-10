using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class LogarithmicDistributionTests {
        readonly LogarithmicDistribution dist_p25 = new(p: 0.25);
        readonly LogarithmicDistribution dist_p50 = new(p: 0.50);
        readonly LogarithmicDistribution dist_p875 = new(p: 0.875);

        LogarithmicDistribution[] Dists => [
            dist_p25,
            dist_p50,
            dist_p875,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (LogarithmicDistribution dist in Dists) {
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
            Assert.AreEqual(1.1586864989274024, (double)dist_p25.Mean, 1e-10);
            Assert.AreEqual(1.4426950408889634, (double)dist_p50.Mean, 1e-10);
            Assert.AreEqual(3.3662884287409147, (double)dist_p875.Mean, 1e-10);
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.AreEqual(0.202360929106562, (double)dist_p25.Variance, 1e-10);
            Assert.AreEqual(0.804021100772319, (double)dist_p50.Variance, 1e-10);
            Assert.AreEqual(15.59840964445234, (double)dist_p875.Variance, 1e-10);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.AreEqual(3.4695883757908987, (double)dist_p25.Skewness, 1e-10);
            Assert.AreEqual(3.01482443189054, (double)dist_p50.Skewness, 1e-10);
            Assert.AreEqual(3.3808906654430424, (double)dist_p875.Skewness, 1e-10);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.AreEqual(15.762245862037844, (double)dist_p25.Kurtosis, 1e-10);
            Assert.AreEqual(13.388420819538247, (double)dist_p50.Kurtosis, 1e-10);
            Assert.AreEqual(17.886884561151653, (double)dist_p875.Kurtosis, 1e-10);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(0.46169444526612724, (double)dist_p25.Entropy, 1e-10);
            Assert.AreEqual(0.8829244358028678, (double)dist_p50.Entropy, 1e-10);
            Assert.AreEqual(1.9797737654713417, (double)dist_p875.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (LogarithmicDistribution dist in Dists) {
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

            foreach (LogarithmicDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 1000000).ToArray();

                for (int i = -1; i <= 50; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, (double)dist.PMF(i) * 0.2 + 1e-5, $"{dist},{i}");
                }

                Assert.AreEqual((double)dist.Mean, samples.Average(), 0.05, $"{dist},mean");
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist_p25 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                8.690148741955517897e-01,
                1.086268592744439737e-01,
                1.810447654574066229e-02,
                3.394589352326374179e-03,
                6.789178704652748574e-04,
                1.414412230135989241e-04,
                3.030883350291405275e-05,
                6.630057328762449567e-06,
                1.473346073058322126e-06,
                3.315028664381224890e-07,
                7.534156055411874268e-08,
                1.726577429365221242e-08,
                3.984409452381280234e-09,
                9.249521943027970198e-10,
                2.158221786706526449e-10,
                5.058332312593421606e-11,
                1.190195838257275710e-11,
                2.810184618107456448e-12,
                6.655700411307133905e-13,
                1.580728847685444302e-13,
                3.763640113536771938e-14,
                8.981413907303660235e-15,
                2.147729412616092716e-15,
                5.145601717726055507e-16,
                1.234944412254253312e-16,
                2.968616375611186202e-17,
                7.146669052397298802e-18,
                1.722857717988634519e-18,
                4.158622077903601061e-19,
                1.005000335493370168e-19,
                2.431452424580734472e-20,
            ];
            ddouble[] expected_dist_p50 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                7.213475204444816935e-01,
                1.803368801111204234e-01,
                6.011229337037347215e-02,
                2.254211001389005292e-02,
                9.016844005556022210e-03,
                3.757018335648342009e-03,
                1.610150715277860892e-03,
                7.044409379340641538e-04,
                3.130848613040285188e-04,
                1.408881875868128470e-04,
                6.404008526673310366e-05,
                2.935170574725267195e-05,
                1.354694111411662004e-05,
                6.289651231554144109e-06,
                2.935170574725267364e-06,
                1.375861206902469050e-06,
                6.474640973658677697e-07,
                3.057469348672153504e-07,
                1.448274954634177864e-07,
                6.879306034512346046e-08,
                3.275860016434450278e-08,
                1.563478644207351164e-08,
                7.477506559252549403e-09,
                3.582971892975179681e-09,
                1.719826508628086429e-09,
                8.268396676096569847e-10,
                3.981079881083533036e-10,
                1.919449228379560580e-10,
                9.266306619763396529e-11,
                4.478714866218974860e-11,
                2.167120096557568626e-11,
            ];
            ddouble[] expected_dist_p875 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                4.207860535926143397e-01,
                1.840938984467687944e-01,
                1.073881074272817898e-01,
                7.047344549915367240e-02,
                4.933141184940757068e-02,
                3.597082114019302390e-02,
                2.697811585514476793e-02,
                2.065511995159521191e-02,
                1.606509329568516289e-02,
                1.265126097035206704e-02,
                1.006350304459823451e-02,
                8.071768067021501195e-03,
                6.519504977209674042e-03,
                5.297097793982860105e-03,
                4.325963198419335159e-03,
                3.548641686203361442e-03,
                2.922410800402768221e-03,
                2.415047814221732165e-03,
                2.001947530210120181e-03,
                1.664118884487162135e-03,
                1.386765737072635221e-03,
                1.158264564486803433e-03,
                9.694170811465636223e-04,
                8.128966149197746228e-04,
                6.828331565326107525e-04,
                5.744990499673408284e-04,
                4.840686439539631115e-04,
                4.084329183361563279e-04,
                3.450553965253734447e-04,
                2.918593562277117153e-04,
                2.471389709992720266e-04,
            ];

            foreach ((LogarithmicDistribution dist, ddouble[] expecteds) in new[]{
                (dist_p25, expected_dist_p25),
                (dist_p50, expected_dist_p50),
                (dist_p875, expected_dist_p875),
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