using DoubleDouble;
using DoubleDoubleStatistic.DiscreteDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.DiscreteDistributions {
    [TestClass()]
    public class GaussKuzminDistributionTests {
        readonly GaussKuzminDistribution dist = new();

        GaussKuzminDistribution[] Dists => [
            dist,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (GaussKuzminDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(2.379246769061239570532, (double)dist.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PMFTest() {
            foreach (GaussKuzminDistribution dist in Dists) {
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

            foreach (GaussKuzminDistribution dist in Dists) {
                int[] samples = dist.Sample(random, count: 1000000).ToArray();

                for (int i = -1; i <= 50; i++) {
                    Assert.AreEqual((double)dist.PMF(i), samples.Count(c => c == i) / (double)samples.Length, (double)dist.PMF(i) * 0.2 + 1e-5, $"{dist},{i}");
                }
            }
        }

        [TestMethod()]
        public void PMFExpectedTest() {
            ddouble[] expected_dist = [
                0,
                4.150374992788438133e-01,
                1.699250014423124289e-01,
                9.310940439148146508e-02,
                5.889368905356856532e-02,
                4.064198449734592739e-02,
                2.974734339405206776e-02,
                2.272007650008352891e-02,
                1.792190799726245745e-02,
                1.449956969511508910e-02,
                1.197264166607594355e-02,
                1.005366466392291006e-02,
                8.562013503424075259e-03,
                7.379530365597452950e-03,
                6.426269159432992928e-03,
                5.646563141142062724e-03,
                5.000681058366530235e-03,
                4.459648190699844021e-03,
                4.001930557496238307e-03,
                3.611253552378835847e-03,
                3.275132032861030780e-03,
                2.983858438821421297e-03,
                2.729792755572332981e-03,
                2.506855610574866382e-03,
                2.310160687201115795e-03,
                2.135744343991011207e-03,
                1.980364128240758373e-03,
                1.841346824167611976e-03,
                1.716472589021663395e-03,
                1.603885702589659184e-03,
                1.502025165231954135e-03,
                1.409570254671353631e-03,
                1.325397466567533817e-03,
                1.248546197258862574e-03,
                1.178191197281153527e-03,
                1.113620310708413486e-03,
                1.054216372001670317e-03,
                9.994423959729344716e-04,
                9.488293935489002743e-04,
                9.019662943926414801e-04,
                8.584915700446877403e-04,
                8.180862373389199884e-04,
                7.804679881383638513e-04,
                7.453862428218402262e-04,
                7.126179650393154111e-04,
                6.819641067135880173e-04,
                6.532465771057086121e-04,
                6.263056494667936739e-04,
                6.009977345355729155e-04,
                5.771934627455232569e-04,
                5.547760271743495649e-04,
                5.336397474895746094e-04,
                5.136888218375962323e-04,
                4.948362390782897347e-04,
                4.770028282465739976e-04,
                4.601164258069911721e-04,
                4.441111443071180847e-04,
                4.289267285613220502e-04,
                4.145079875918718823e-04,
                4.008042923095860347e-04,
                3.877691303788703995e-04,
                3.753597109477273069e-04,
                3.635366129577204120e-04,
                3.522634716290213848e-04,
                3.415066984555677893e-04,
            ];

            foreach ((GaussKuzminDistribution dist, ddouble[] expecteds) in new[]{
                (dist, expected_dist),
            }) {
                for ((int x, int i) = (0, 0); i < expecteds.Length; x++, i++) {
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