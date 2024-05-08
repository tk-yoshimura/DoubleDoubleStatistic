using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class BenktanderDistributionTests {
        readonly BenktanderDistribution dist_a1b1 = new(alpha: 1, beta: 1);
        readonly BenktanderDistribution dist_a2b1 = new(alpha: 2, beta: 1);
        readonly BenktanderDistribution dist_a2b2 = new(alpha: 2, beta: 2);
        readonly BenktanderDistribution dist_a3b4 = new(alpha: 3, beta: 4);

        BenktanderDistribution[] Dists => [
            dist_a1b1,
            dist_a2b1,
            dist_a2b2,
            dist_a3b4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Beta={dist.Beta}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode - 1e-4), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode + 1e-4), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (BenktanderDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 5) * 0.1, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (BenktanderDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.NegativeInfinity)) && dist.PDF(ddouble.NegativeInfinity) >= 0d, "pdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MinValue)) && dist.PDF(ddouble.MinValue) >= 0d, "pdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MinValue / 2)) && dist.PDF(ddouble.MinValue / 2) >= 0d, "pdf(-lval / 2)");

                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.PositiveInfinity)) && dist.PDF(ddouble.PositiveInfinity) >= 0d, "pdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MaxValue)) && dist.PDF(ddouble.MaxValue) >= 0d, "pdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.PDF(ddouble.MaxValue / 2)) && dist.PDF(ddouble.MaxValue / 2) >= 0d, "pdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.PDF(ddouble.NaN)), "pdf(NaN)");

                Assert.IsTrue(dist.CDF(ddouble.NegativeInfinity, Interval.Lower) == 0d, "cdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue, Interval.Lower)) && dist.CDF(ddouble.MinValue, Interval.Lower) >= 0d, "cdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue / 2, Interval.Lower)) && dist.CDF(ddouble.MinValue / 2, Interval.Lower) >= 0d, "cdf(-lval / 2)");

                Assert.IsTrue(dist.CDF(ddouble.PositiveInfinity, Interval.Lower) == 1d, "cdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue, Interval.Lower)) && dist.CDF(ddouble.MaxValue, Interval.Lower) <= 1d, "cdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue / 2, Interval.Lower)) && dist.CDF(ddouble.MaxValue / 2, Interval.Lower) <= 1d, "cdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.CDF(ddouble.NaN, Interval.Lower)), "cdf(NaN)");

                Assert.IsTrue(dist.CDF(ddouble.NegativeInfinity, Interval.Upper) == 1d, "ccdf(-inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue, Interval.Upper)) && dist.CDF(ddouble.MinValue, Interval.Upper) <= 1d, "ccdf(-lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MinValue / 2, Interval.Upper)) && dist.CDF(ddouble.MinValue / 2, Interval.Upper) <= 1d, "ccdf(-lval / 2)");

                Assert.IsTrue(dist.CDF(ddouble.PositiveInfinity, Interval.Upper) == 0d, "cdf(+inf)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue, Interval.Upper)) && dist.CDF(ddouble.MaxValue, Interval.Upper) >= 0d, "ccdf(+lval)");
                Assert.IsTrue(ddouble.IsFinite(dist.CDF(ddouble.MaxValue / 2, Interval.Upper)) && dist.CDF(ddouble.MaxValue / 2, Interval.Upper) >= 0d, "ccdf(+lval / 2)");

                Assert.IsTrue(ddouble.IsNaN(dist.CDF(ddouble.NaN, Interval.Upper)), "ccdf(NaN)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(0d, Interval.Lower)), "quantile(0)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(1d, Interval.Lower)), "quantile(1)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(0d, Interval.Upper)), "cquantile(0)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(1d, Interval.Upper)), "cquantile(1)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.Epsilon, Interval.Lower)), "quantile(0+eps)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)), "quantile(1-eps)");

                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.Epsilon, Interval.Upper)), "cquantile(0+eps)");
                Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)), "cquantile(1-eps)");
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist_a1b1 = [
                0.0,
                0.527935109400371,
                0.7492226660119179,
                0.8051141944648212,
                0.7768583547908648,
                0.7098993281332316,
                0.6289470533161791,
                0.5469169653977595,
                0.4701161279382973,
                0.4012167745178745,
                0.3409473393163515,
                0.2890457455580052,
                0.2447900186317338,
                0.2072875769631228,
                0.1756274518742331,
                0.1489553469161867,
                0.1265059330076899,
                0.1076120692056761,
                0.09170214834402787,
                0.07829186231420654,
                0.06697385298818576,
                0.05740709064012241,
                0.04930689767854474,
                0.04243601813882314,
                0.03659685151139837,
                0.03162482301186688,
                0.02738279281109533,
                0.02375637988848176,
                0.02065007223055612,
                0.01798400268241064,
                0.01569128243671761,
                0.01371579837876193,
                0.01201039444321453,
            ];
            ddouble[] expected_dist_a2b1 = [
                2.0,
                1.611062726660252,
                1.253017462289316,
                0.959331619053212,
                0.7300235650649433,
                0.5550355283186131,
                0.4228565601719388,
                0.3233577365207783,
                0.2484311374766015,
                0.1918610296975989,
                0.1489821081793818,
                0.1163278095256574,
                0.0913323116651013,
                0.0720970058205608,
                0.05721427931755407,
                0.04563709610061779,
                0.03658340198362015,
                0.02946634289634588,
                0.02384336303898903,
                0.01937902305458611,
                0.01581776525591041,
                0.0129638922997455,
                0.01066678697462909,
                0.008809951513662722,
                0.007302840799469102,
                0.006074747790967127,
                0.005070203133217545,
                0.004245497159992087,
                0.003566037807205997,
                0.003004334036105571,
                0.002538449539195614,
                0.002150811679368316,
                0.001827289996254191,
            ];
            ddouble[] expected_dist_a2b2 = [
                1.0,
                1.389794033203485,
                1.345831059131862,
                1.141068372356057,
                0.9056650030245392,
                0.6928379096247764,
                0.5187351876645658,
                0.3835094850635507,
                0.2815323369200114,
                0.2059523910782962,
                0.150501483826436,
                0.1100443696686589,
                0.0806018609448777,
                0.05918593981540956,
                0.04359418735973387,
                0.03222108897509605,
                0.02390364326501894,
                0.01780209992321863,
                0.01331077775891655,
                0.009992655696103102,
                0.007531958820859629,
                0.005700025965958773,
                0.004330830690009089,
                0.003303453546756722,
                0.002529532787392908,
                0.001944270543995673,
                0.001499975727940043,
                0.001161417320656689,
                9.024712801212231e-4,
                7.03693559044272e-4,
                5.505577254751596e-4,
                4.321708367252325e-4,
                3.403345091746621e-4,

            ];
            ddouble[] expected_dist_a3b4 = [
                1.333333333333333,
                2.009550049277416,
                1.761635459122323,
                1.280398460487352,
                0.8466307523781363,
                0.5302353696921698,
                0.3213811119059761,
                0.1909553176231714,
                0.1121421872866319,
                0.06545048073805115,
                0.0381066423220927,
                0.02219125179268764,
                0.01295010553085571,
                0.007583333852250254,
                0.004460275040305798,
                0.002636796283775047,
                0.001567527890240154,
                9.373938658473222e-4,
                5.640167239485485e-4,
                3.414935222072014e-4,
                2.080764013969504e-4,
                1.275912221072229e-4,
                7.873498820643093e-5,
                4.889255468892458e-5,
                3.055028959464566e-5,
                1.920649976643374e-5,
                1.214783586646449e-5,
                7.728962286290448e-6,
                4.94613808302185e-6,
                3.183343768651895e-6,
                2.060257756469235e-6,
                1.340689676847106e-6,
                8.771049913097381e-7,
            ];

            foreach ((BenktanderDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a2b2, expected_dist_a2b2), (dist_a3b4, expected_dist_a3b4),
            }) {
                for ((ddouble x, int i) = (1, 0); i < expecteds.Length; x += 1d / 8, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist_a1b1 = [
                0.0,
                0.03720011488456565,
                0.1193372179364306,
                0.2176956917758091,
                0.3171578270416288,
                0.4103254543557374,
                0.4940653899579578,
                0.5675261365434028,
                0.6310173624822851,
                0.6853869803298225,
                0.731682445082772,
                0.7709730759196431,
                0.8042626746133638,
                0.8324519808565118,
                0.8563279132832206,
                0.8765665508205519,
                0.8937425572548463,
                0.9083410484531992,
                0.9207697839529065,
                0.9313706300927096,
                0.9404298357237713,
                0.9481869838766239,
                0.954842649774187,
                0.9605648748013422,
                0.9654945966502264,
                0.9697501804816961,
                0.973431187519317,
                0.9766215033071645,
                0.979391931906996,
                0.9818023466702338,
                0.9839034738954638,
                0.985738373065672,
                0.9873436665168877,
            ];
            ddouble[] expected_dist_a2b1 = [
                0.0,
                0.2257609936037011,
                0.4041697525721573,
                0.5417380605171219,
                0.6466969148990087,
                0.7265084652744347,
                0.7872372731785982,
                0.8335813633612104,
                0.8690978945079728,
                0.8964518495757697,
                0.9176319631272833,
                0.9341221625777831,
                0.9470321150765568,
                0.9571945741012946,
                0.9652372289069214,
                0.9716355110844086,
                0.9767513802764922,
                0.9808618565258787,
                0.9841800721476498,
                0.9868708641739375,
                0.9890623740492032,
                0.9908547178667418,
                0.9923264988422067,
                0.9935397234769843,
                0.9945435312789114,
                0.995377038434097,
                0.9960715165446865,
                0.9966520699346377,
                0.9971389329914255,
                0.9975484782121086,
                0.9978940029535898,
                0.9981863461262303,
                0.9984343736201682,
            ];
            ddouble[] expected_dist_a2b2 = [
                0.0,
                0.1559685917746482,
                0.3296914092104671,
                0.4859208222950503,
                0.6137840412111077,
                0.7133261632503842,
                0.7886275351087502,
                0.8446369886954438,
                0.8858915404504959,
                0.916118956850219,
                0.9382166967703878,
                0.9543674369017697,
                0.966185349795011,
                0.974851255880775,
                0.9812236854513928,
                0.9859249713050054,
                0.9894058325910783,
                0.991992890203031,
                0.9939232188172118,
                0.9953693017507278,
                0.9964569893719568,
                0.9972784040950852,
                0.9979012147566693,
                0.9983753095485418,
                0.9987376065105515,
                0.9990155300515965,
                0.9992295306938842,
                0.9993949171691953,
                0.9995231930362468,
                0.9996230352502677,
                0.9997010131809446,
                0.9997621188525554,
                0.9998101594075233,
            ];
            ddouble[] expected_dist_a3b4 = [
                0.0,
                0.2239051643069938,
                0.4646537963216442,
                0.6551713135506335,
                0.7870087597185158,
                0.8718182303967116,
                0.9240692072516227,
                0.955425735143329,
                0.9739476156271725,
                0.9847902032490535,
                0.9911097082739417,
                0.9947887391737589,
                0.9969329718338453,
                0.9981861232846546,
                0.9989213711704249,
                0.9993548064043095,
                0.9996116870686121,
                0.9997648078842029,
                0.9998566312989309,
                0.999912038486925,
                0.9999456831562782,
                0.9999662431431469,
                0.9999788871777097,
                0.9999867122826552,
                0.9999915854435651,
                0.9999946390539072,
                0.9999965641785498,
                0.9999977851403442,
                0.9999985640654783,
                0.9999990638583381,
                0.9999993863614786,
                0.9999995956163334,
                0.9999997321257723,
            ];

            foreach ((BenktanderDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a2b2, expected_dist_a2b2), (dist_a3b4, expected_dist_a3b4),
            }) {
                for ((ddouble x, int i) = (1, 0); i < expecteds.Length; x += 1d / 8, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}