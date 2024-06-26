﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class SnedecorFDistributionTests {
        readonly SnedecorFDistribution dist_n_1_m_1 = new(n: 1, m: 1);
        readonly SnedecorFDistribution dist_n_2_m_1 = new(n: 2, m: 1);
        readonly SnedecorFDistribution dist_n_1_m_2 = new(n: 1, m: 2);
        readonly SnedecorFDistribution dist_n_2_m_2 = new(n: 2, m: 2);
        readonly SnedecorFDistribution dist_n_3_m_4 = new(n: 3, m: 4);
        readonly SnedecorFDistribution dist_n_4_m_2 = new(n: 4, m: 2);
        readonly SnedecorFDistribution dist_n_2_m_4 = new(n: 2, m: 4);
        readonly SnedecorFDistribution dist_n_4_m_4 = new(n: 4, m: 4);
        readonly SnedecorFDistribution dist_n_6_m_8 = new(n: 6, m: 8);
        readonly SnedecorFDistribution dist_n_8_m_6 = new(n: 8, m: 6);
        readonly SnedecorFDistribution dist_n_8_m_8 = new(n: 8, m: 8);
        readonly SnedecorFDistribution dist_n_10_m_10 = new(n: 10, m: 10);

        SnedecorFDistribution[] Dists => [
            dist_n_1_m_1,
            dist_n_2_m_1,
            dist_n_1_m_2,
            dist_n_2_m_2,
            dist_n_3_m_4,
            dist_n_4_m_2,
            dist_n_2_m_4,
            dist_n_4_m_4,
            dist_n_6_m_8,
            dist_n_8_m_6,
            dist_n_8_m_8,
            dist_n_10_m_10
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"D1={dist.N}");
                Console.WriteLine($"D2={dist.M}");
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
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (SnedecorFDistribution dist in Dists) {
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
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
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
            foreach (SnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 10; i++) {
                    ddouble p = (ddouble)i / 10;
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

            foreach (SnedecorFDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.04, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (SnedecorFDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 20000).ToArray();

                (SnedecorFDistribution? dist_fit, ddouble error) = SnedecorFDistribution.Fit(xs, (0.05, 0.80));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (SnedecorFDistribution dist in Dists) {
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
            ddouble[] expected_dist_n_1_m_1 = [
                "+inf",
                1.198343100927211,
                0.800281169917427,
                0.6190359527542925,
                0.5092958178940648,
                0.4338362169399933,
                0.378034812719358,
                0.334775174636754,
                0.3001054387190351,
                0.2716244362101679,
                0.2477745826682437,
                0.2274940249497353,
                0.2100300553987778,
                0.1948320621534664,
                0.1814867103882289,
                0.1696769222362837,
                0.1591549430918952,
                0.1497240963894878,
                0.1412260888089577,
                0.1335319700792855,
                0.1265355632741647,
                0.1201486042020076,
                0.1142970875915064,
                0.1089184810768201,
                0.1039595734978234,
                0.09937479373542729,
                0.0951248832832613,
                0.0911758380025702,
                0.08749805700733343,
                0.08406565256533029,
                0.08085588634887607,
                0.07784870569688242,
                0.07502635967975878,
                0.07237307932003163,
                0.06987480975026314,
                0.06751898469154644,
                0.0652943356274442,
                0.0631907295853871,
                0.06119903063313276,
                0.05931098113456614,
                0.05751909954798514,
                0.05581659213659687,
                0.05419727642947428,
                0.05265551464760641,
                0.05118615561369044,
                0.04978448391114679,
                0.04844617525923173,
                0.04716725723621016,
                0.04594407461848265,
                0.04477325871596174,
                0.04365170017731421,
                0.04257652481646656,
                0.0415450720768421,
                0.04055487580442758,
                0.03960364704678987,
                0.03868925863406004,
                0.03780973133088102,
                0.03696322137636136,
                0.03614800925299703,
                0.03536248954598202,
                0.03460516177187364,
                0.03387462207066411,
                0.03316955566831391,
                0.03248873002803941,
                0.03183098861837905,
            ];
            ddouble[] expected_dist_n_2_m_1 = [
                1,
                0.8380524814062784,
                0.7155417527999323,
                0.6202202657617755,
                0.5443310539518172,
                0.4827474095709144,
                0.4319593977248309,
                0.3894915964481179,
                0.3535533905932735,
                0.3228208661506047,
                0.2962962962962962,
                0.2732150250623369,
                0.2529822128134702,
                0.2351289142279495,
                0.2192809778747371,
                0.2051366346502916,
                0.1924500897298752,
                0.1810193359837561,
                0.1706769834539166,
                0.1612832752449828,
                0.1527207096642424,
                0.1448898569486593,
                0.1377060745318192,
                0.1310969052905192,
                0.1249999999999999,
                0.1193614457981408,
                0.1141344117818037,
                0.1092780442820017,
                0.1047565601757848,
                0.1005384983730329,
                0.09659609847181541,
                0.09290478228878728,
                0.08944271909999153,
                0.08619045937548993,
                0.08313062485180658,
                0.0802476451746518,
                0.07752753322022193,
                0.07495769268547267,
                0.07252675271550425,
                0.0702244252771944,
                0.06804138174397716,
                0.06596914576667497,
                0.06399999999999997,
                0.06212690465735994,
                0.0603434261963643,
                0.05864367470836308,
                0.05702224880885189,
                0.05547418701065262,
                0.05399492471560386,
                0.05258025608879683,
                0.05122630018677291,
                0.04992947080126498,
                0.04868644955601474,
                0.04749416185837035,
                0.04634975536174835,
                0.04525058064125701,
                0.04419417382415919,
                0.04317824095049989,
                0.04220064386804794,
                0.041259387490453,
                0.04035260826882558,
                0.03947856374533325,
                0.0386356230733036,
                0.03782225840210428,
                0.03703703703703702,
            ];
            ddouble[] expected_dist_n_1_m_2 = [
                "inf",
                1.350420603777534,
                0.9130752942544295,
                0.713801196213784,
                0.5925925925925923,
                0.5086890358270668,
                0.4461582676429621,
                0.3972775599056752,
                0.3577708763999661,
                0.3250445775892938,
                0.2974171650890767,
                0.2737416869933858,
                0.2532038632082877,
                0.23520660781361,
                0.2193002869534888,
                0.2051386763314952,
                0.1924500897298752,
                0.1810178425927149,
                0.1706666666666666,
                0.161253036223557,
                0.1526581348790258,
                0.1447826473854021,
                0.1375428392736377,
                0.1308675623832513,
                0.1246959372777099,
                0.1189755380381232,
                0.1136609550619817,
                0.1087126458580362,
                0.1040960077811768,
                0.09978062360455846,
                0.0957396429956476,
                0.09194927181235203,
                0.08838834764831839,
                0.0850379849033906,
                0.08188127629978158,
                0.0789030405315248,
                0.0760896078545358,
                0.07342863706227747,
                0.07090895856726684,
                0.06852043930909624,
                0.06625386599999376,
                0.06410084384746183,
                0.06205370839646919,
                0.06010544853849588,
                0.05824963906240625,
                0.05648038138875847,
                0.05479225134721455,
                0.0531802530358954,
                0.0516397779494322,
                0.05016656868508718,
                0.04875668663839406,
                0.04740648318508286,
                0.04611257391762302,
                0.04487181556497755,
                0.04368128527506969,
                0.04253826198261752,
                0.0414402096216865,
                0.04038476197360965,
                0.03936970896769892,
                0.03839298427514017,
                0.0374526540562202,
                0.03654690673807482,
                0.03567404371487631,
                0.03483247087514862,
                0.03402069087198858,
            ];
            ddouble[] expected_dist_n_2_m_2 = [
                1,
                0.8858131487889269,
                0.7901234567901231,
                0.7091412742382268,
                0.6399999999999997,
                0.5804988662131517,
                0.5289256198347105,
                0.4839319470699431,
                0.4444444444444442,
                0.4095999999999999,
                0.3786982248520708,
                0.3511659807956102,
                0.3265306122448978,
                0.304399524375743,
                0.2844444444444443,
                0.2663891779396461,
                0.2499999999999999,
                0.2350780532598713,
                0.2214532871972317,
                0.2089795918367346,
                0.1975308641975308,
                0.1869978086194302,
                0.1772853185595567,
                0.1683103221564759,
                0.1599999999999999,
                0.1522903033908387,
                0.1451247165532879,
                0.1384532179556517,
                0.1322314049586776,
                0.1264197530864197,
                0.1209829867674858,
                0.1158895427795382,
                0.111111111111111,
                0.1066222407330279,
                0.1024,
                0.09842368319876967,
                0.09467455621301771,
                0.09113563545745813,
                0.08779149519890256,
                0.08462809917355368,
                0.08163265306122444,
                0.07879347491535853,
                0.07609988109393576,
                0.073542085607584,
                0.07111111111111108,
                0.06879871002418701,
                0.06659729448491153,
                0.06449987402368353,
                0.06249999999999997,
                0.06059171597633133,
                0.05876951331496783,
                0.05702829137892624,
                0.05536332179930793,
                0.05377021634110479,
                0.05224489795918365,
                0.05078357468756197,
                0.04938271604938269,
                0.04803903171326701,
                0.04674945215485754,
                0.04551111111111109,
                0.04432132963988918,
                0.04317760161916004,
                0.04207758053911898,
                0.04101906745713826,
                0.03999999999999998,
            ];
            ddouble[] expected_dist_n_3_m_4 = [
                0,
                0.5187179430251337,
                0.6293104423762348,
                0.6654626687632927,
                0.6673835582470904,
                0.6516032819743296,
                0.6264987694166048,
                0.5967272359245283,
                0.5650003559193203,
                0.5329291330975444,
                0.5014707343549458,
                0.4711816330883287,
                0.4423680000000004,
                0.4151781315771823,
                0.3896606672245072,
                0.3658019086462543,
                0.3435500312753235,
                0.3228309022472129,
                0.3035584365814747,
                0.2856413551477751,
                0.2689875505600926,
                0.2535068532499524,
                0.2391127244147413,
                0.2257232294021309,
                0.2132615306649621,
                0.2016560629167905,
                0.1908405014679057,
                0.1807535995667191,
                0.1713389464855033,
                0.1625446815025117,
                0.1543231874680197,
                0.1466307797012701,
                0.1394274004634671,
                0.1326763254459722,
                0.1263438860905969,
                0.1203992097650164,
                0.1148139786029416,
                0.1095622070120911,
                0.104620037330952,
                0.09996555279162257,
                0.0955786067597303,
                0.09144066713041234,
                0.08753467473112153,
                0.08384491459605942,
                0.08035689901816681,
                0.07705726134221298,
                0.07393365952924455,
                0.07097468859342142,
                0.06816980108366502,
                0.06550923485232246,
                0.06298394741979299,
                0.06058555630691751,
                0.05830628476546532,
                0.05613891239109615,
                0.05407673015274769,
                0.05211349941764841,
                0.05024341459228914,
                0.04846106903697258,
                0.04676142394527352,
                0.04513977991017121,
                0.04359175092604076,
                0.04211324060038237,
                0.04070042037138346,
                0.03934970954738051,
                0.03805775700224587,
            ];
            ddouble[] expected_dist_n_4_m_2 = [
                0,
                0.3511659807956106,
                0.5120000000000003,
                0.577009767092412,
                0.592592592592593,
                0.5826126536185711,
                0.5597667638483969,
                0.5309629629629633,
                0.5000000000000003,
                0.4689599023000207,
                0.4389574759945133,
                0.4105554745589739,
                0.3840000000000003,
                0.3593564409890943,
                0.3365890308039071,
                0.3156077915673545,
                0.2962962962962965,
                0.2785280000000002,
                0.262175694128357,
                0.2471168013006149,
                0.2332361516034987,
                0.2204272417893314,
                0.2085925925925927,
                0.1976435836326409,
                0.1875000000000001,
                0.1780894342877815,
                0.1693466313861186,
                0.1612128279883383,
                0.1536351165980797,
                0.146565849999013,
                0.1399620935996502,
                0.1337851278679682,
                0.1280000000000001,
                0.1225751222414069,
                0.1174819134002808,
                0.1126944797313445,
                0.1081893313298273,
                0.1039451303155008,
                0.09994246732966228,
                0.09616366315748925,
                0.09259259259259264,
                0.08921452796028871,
                0.08601600000000005,
                0.08298467406955098,
                0.08010923987255354,
                0.077379313124257,
                0.07478534776202819,
                0.07231855747558232,
                0.0699708454810496,
                0.06773474159390479,
                0.06560334577063434,
                0.06357027738960658,
                0.06162962962962967,
                0.05977592838167076,
                0.0580040951965359,
                0.05630941383019996,
                0.05468750000000003,
                0.05313427401001369,
                0.05164593594345664,
                0.05021894315457689,
                0.04884998982291882,
                0.04753598835952747,
                0.04627405247813414,
                0.04506148176501982,
                0.04389574759945133,
            ];
            ddouble[] expected_dist_n_2_m_4 = [
                1,
                0.9118179035534413,
                0.8337064929778145,
                0.7642682215743446,
                0.7023319615912212,
                0.646911337926678,
                0.5971715993585075,
                0.5524031086161267,
                0.5120000000000003,
                0.4754428983909116,
                0.4422848504481161,
                0.412139811588917,
                0.3846731780616081,
                0.3595939643347053,
                0.3366483110051782,
                0.3156140739527852,
                0.2962962962962965,
                0.2785234043638282,
                0.2621440000000002,
                0.2470241460675006,
                0.2330450614474285,
                0.2201011573312199,
                0.2080983589899915,
                0.1969526671675433,
                0.1865889212827989,
                0.1769397331432615,
                0.1679445651728239,
                0.15954893148764,
                0.1517037037037038,
                0.1443645062802614,
                0.137491188614011,
                0.1310473630957381,
                0.1250000000000001,
                0.1193190714610834,
                0.1139772379441802,
                0.1089495715895905,
                0.1042133116222268,
                0.09974764770523797,
                0.09553352769679307,
                0.09155348676067519,
                0.08779149519890266,
                0.08423282273011211,
                0.08086391724083475,
                0.07767229629629635,
                0.07464644991981344,
                0.07177575334094144,
                0.06905038857701584,
                0.06646127385460383,
                0.06400000000000004,
                0.06165877303407156,
                0.05943036229886395,
                0.05730805352342747,
                0.05528560630601451,
                0.05335721555058013,
                0.05151747644861462,
                0.04976135264379968,
                0.04808414725770101,
                0.04648147649045564,
                0.04494924554183816,
                0.04348362662575925,
                0.04208103887564727,
                0.04073812995970697,
                0.03945175924409815,
                0.03821898235894448,
                0.03703703703703706,
            ];
            ddouble[] expected_dist_n_4_m_4 = [
                0,
                0.2942493504627582,
                0.4682213077274811,
                0.5657415151817442,
                0.6144000000000007,
                0.631835500640166,
                0.6294652004644499,
                0.61474908966163,
                0.5925925925925932,
                0.5662310400000007,
                0.5377962956479121,
                0.5086848775310906,
                0.4798000832986262,
                0.4517129683958715,
                0.4247703703703709,
                0.3991679669439032,
                0.3750000000000004,
                0.3522932809183753,
                0.331030519270603,
                0.311166347355269,
                0.2926383173296757,
                0.2753744208741942,
                0.2592981944582994,
                0.2443321441955946,
                0.2304000000000003,
                0.2174281547519414,
                0.205346537708054,
                0.1940890973180366,
                0.1835940168021312,
                0.173803749428441,
                0.1646649347307937,
                0.1561282387106864,
                0.1481481481481483,
                0.1406827399592806,
                0.1336934400000002,
                0.1271447810641548,
                0.1210041665207802,
                0.1152416436970067,
                0.10982968946694,
                0.1047430093572845,
                0.09995835068721379,
                0.0954543297220227,
                0.09121127246455098,
                0.08721106848265074,
                0.08343703703703713,
                0.07987380470424306,
                0.07650719366424812,
                0.0733241198273773,
                0.07031250000000008,
                0.0674611673260741,
                0.06475979428646604,
                0.06219882258659463,
                0.05976939931274776,
                0.05746331878671321,
                0.05527296959600173,
                0.0531912863222233,
                0.05121170553269325,
                0.04932812563989444,
                0.04753487026995018,
                0.04582665481481487,
                0.04419855587357377,
                0.04264598331620874,
                0.04116465472860562,
                0.03975057202063333,
                0.03840000000000005,
            ];
            ddouble[] expected_dist_n_6_m_8 = [
                0,
                0.07175154784432666,
                0.2112168708713445,
                0.3542724508140492,
                0.4750942014064494,
                0.5661157827729545,
                0.6280011329288305,
                0.6646890023090353,
                0.6810141913364425,
                0.6816323061698375,
                0.6705943931053011,
                0.651235585322031,
                0.6262062317567987,
                0.5975593205917207,
                0.5668522928388978,
                0.5352441454769685,
                0.5035802623542416,
                0.4724630545540623,
                0.4423090772149864,
                0.4133943177713291,
                0.3858896125670295,
                0.3598880580102785,
                0.3354260505319665,
                0.3124993212559247,
                0.2910750749539997,
                0.2711011180740026,
                0.2525126725369495,
                0.2352374190455452,
                0.2191991915512896,
                0.2043206482820888,
                0.1905251695223156,
                0.1777381739279657,
                0.1658879999999997,
                0.1549064645378399,
                0.1447291831424281,
                0.1352957173123768,
                0.1265495969533121,
                0.1184382550928878,
                0.1109129024061149,
                0.1039283621473202,
                0.09744288074806537,
                0.09141792528411963,
                0.08581797593917202,
                0.08061031926818135,
                0.07576484631250577,
                0.07125385830674509,
                0.06705188173967108,
                0.06313549380902851,
                0.05948315878126637,
                0.05607507538557314,
                0.05289303510105644,
                0.04992029100905371,
                0.04714143675836775,
                0.04454229511363286,
                0.04210981551378833,
                0.0398319800494474,
                0.03769771726773016,
                0.03569682322553443,
                0.03381988923322566,
                0.03205823575732872,
                0.03040385198076745,
                0.02884934055088186,
                0.02738786707766766,
                0.02601311397658181,
                0.02471923828124995,
            ];
            ddouble[] expected_dist_n_8_m_6 = [
                0,
                0.02643695945834061,
                0.1258950655885604,
                0.2621439999999995,
                0.3955078124999992,
                0.5053386717951384,
                0.5852766346593495,
                0.6365519964156724,
                0.6635519999999987,
                0.6714403498056554,
                0.6650529212269947,
                0.648482371378375,
                0.6249999999999987,
                0.5971246055423988,
                0.5667423183881768,
                0.5352324048096475,
                0.5035802623542416,
                0.4724716452458119,
                0.4423679999999992,
                0.4135651332183211,
                0.3862380981445305,
                0.3604751043648145,
                0.3363028860796646,
                0.3137055258656797,
                0.2926383173296748,
                0.2730378960450335,
                0.2548295801977464,
                0.2379326351250655,
                0.2222639999999995,
                0.2077408804118008,
                0.1942825086243216,
                0.1818112964330691,
                0.1702535478341106,
                0.1595398555098305,
                0.1496052728485141,
                0.140389329122997,
                0.1318359374999997,
                0.1238932321883831,
                0.1165133611007998,
                0.1096522530305105,
                0.103269372884143,
                0.09732747447671616,
                0.09179235742485455,
                0.08663263249693846,
                0.08181949819256537,
                0.07732653017683096,
                0.07312948437443763,
                0.06920611394991129,
                0.06553599999999987,
                0.06210039551536997,
                0.05888208199570152,
                0.05586523799886849,
                0.05303531885147084,
                0.050378946729976,
                0.04788381032834362,
                0.04553857335128727,
                0.04333279110643311,
                0.04125683450931977,
                0.03930182085938431,
                0.03745955079061558,
                0.03572245084590756,
                0.0340835211682738,
                0.03253628784431037,
                0.03107475947519993,
                0.02969338758790309,
            ];
            ddouble[] expected_dist_n_8_m_8 = [
                0,
                0.02104440144910709,
                0.1065707188243211,
                0.2333796077084316,
                0.3670015999999999,
                0.4851584547021628,
                0.577830222952134,
                0.6429828374566072,
                0.6828227404359088,
                0.7013509795676157,
                0.7029770796139966,
                0.6918244255701088,
                0.6714403498056565,
                0.6447242753245772,
                0.6139627437280899,
                0.580909094186927,
                0.5468749999999998,
                0.5128179214530457,
                0.4794177705124708,
                0.4471409004032574,
                0.4162918704075043,
                0.3870544283248003,
                0.3595234465971325,
                0.3337295023139544,
                0.3096575999999999,
                0.2872612997845819,
                0.2664732812451629,
                0.2472131661413542,
                0.229393248232155,
                0.2129226351744096,
                0.1977101928221619,
                0.1836665917454694,
                0.1707056851089772,
                0.1587453922738219,
                0.1477082202767359,
                0.1375215229584616,
                0.1281175727596509,
                0.119433501335932,
                0.1114111508026244,
                0.1039968665457327,
                0.09714125431216095,
                0.09079891808990735,
                0.08492819062615085,
                0.07949086493859578,
                0.07445193257125435,
                0.06977933241487623,
                0.06544371249048464,
                0.0614182060570555,
                0.05767822265624998,
                0.05420125417603536,
                0.05096669564614187,
                0.04795568022998761,
                0.04515092771844557,
                0.04253660573680029,
                0.04009820282948186,
                0.03782241257416324,
                0.03569702788744349,
                0.0337108447110742,
                0.03185357430492643,
                0.03011576341653424,
                0.02848872164409566,
                0.02696445535811496,
                0.0255356075949161,
                0.02419540338200315,
                0.02293759999999999,
            ];
            ddouble[] expected_dist_n_10_m_10 = [
                0,
                0.005242895862753338,
                0.04736476392192063,
                0.1396398760803917,
                0.2642411520000007,
                0.3960477181242157,
                0.5157492899076914,
                0.6125961249114008,
                0.682822740435911,
                0.7271606956157063,
                0.748732984204259,
                0.7516117216070343,
                0.7399546712143994,
                0.7175528200996504,
                0.687638272975463,
                0.6528426864951957,
                0.6152343750000017,
                0.5763903910546648,
                0.5374787461800727,
                0.4993377565319657,
                0.4625465226750062,
                0.4274845110497443,
                0.3943802904500414,
                0.3633504640577977,
                0.3344302080000009,
                0.3075968706795057,
                0.2827879719336432,
                0.2599147620220628,
                0.2388723080764597,
                0.2195468949353919,
                0.2018213688543814,
                0.1855789193190986,
                0.1707056851089777,
                0.1570924831497718,
                0.1446358892949803,
                0.1332388457729043,
                0.1228109277341036,
                0.1132683686575022,
                0.1045339192715985,
                0.09653659545798951,
                0.08921135600096443,
                0.08249873998196593,
                0.07634448527154605,
                0.07069914330649045,
                0.06551770066270403,
                0.06075921446498246,
                0.05638646612083398,
                0.05236563600329466,
                0.04866600036621107,
                0.04525965082439132,
                0.04212123607119176,
                0.03922772506226667,
                0.03655819060939895,
                0.03409361215955269,
                0.03181669644918898,
                0.02971171469821204,
                0.02776435502356725,
                0.02596158879606848,
                0.02429154972559773,
                0.02274342453216673,
                0.02130735413824336,
                0.01997434439751076,
                0.01873618545425803,
                0.01758537890414463,
                0.01651507200000004,
            ];

            foreach ((SnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_n_1_m_1,   expected_dist_n_1_m_1),
                 (dist_n_2_m_1,   expected_dist_n_2_m_1),
                 (dist_n_1_m_2,   expected_dist_n_1_m_2),
                 (dist_n_2_m_2,   expected_dist_n_2_m_2),
                 (dist_n_3_m_4,   expected_dist_n_3_m_4),
                 (dist_n_4_m_2,   expected_dist_n_4_m_2),
                 (dist_n_2_m_4,   expected_dist_n_2_m_4),
                 (dist_n_4_m_4,   expected_dist_n_4_m_4),
                 (dist_n_6_m_8,   expected_dist_n_6_m_8),
                 (dist_n_8_m_6,   expected_dist_n_8_m_6),
                 (dist_n_8_m_8,   expected_dist_n_8_m_8),
                 (dist_n_10_m_10, expected_dist_n_10_m_10),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (ddouble.IsInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else if (expected > 0) {
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
            ddouble[] expected_dist_n_1_m_1 = [
                0.0,
                0.1559582607547385,
                0.2163468959387853,
                0.2601469382930059,
                0.2951672353008664,
                0.3245103583044378,
                0.3498017122810425,
                0.3720239771892146,
                0.3918265520306068,
                0.4096655293982663,
                0.4258757566828432,
                0.4407115039679876,
                0.4543710516570099,
                0.4670123752688586,
                0.4787635903929203,
                0.4897301575265077,
                0.4999999999999994,
                0.5096472309774334,
                0.5187349260190778,
                0.5273172240862714,
                0.5354409456024599,
                0.543146856491735,
                0.5504706682632844,
                0.557443838299598,
                0.5640942168489756,
                0.5704465749545546,
                0.5765230388779268,
                0.5823434503447235,
                0.5879256673992912,
                0.5932858173032707,
                0.5984385104078173,
                0.6033970220363094,
                0.608173447969393,
                0.6127788380105946,
                0.6172233112449614,
                0.6215161559244841,
                0.6256659163780024,
                0.6296804689167759,
                0.6335670883653126,
                0.6373325065717252,
                0.6409829640286244,
                0.6445242555535369,
                0.6479617708286441,
                0.6513005304766762,
                0.6545452182480777,
                0.6577002098099427,
                0.6607695985565596,
                0.6637572188021679,
                0.6666666666666669,
                0.6695013189228645,
                0.6722643500381511,
                0.6749587476130693,
                0.6775873263933165,
                0.6801527410094983,
                0.6826574975798751,
                0.6851039642949148,
                0.6874943810882701,
                0.689830868486515,
                0.6921154357193064,
                0.6943499881623495,
                0.6965363341774535,
                0.6986761914068774,
                0.7007711925729545,
                0.702822890828557,
                0.7048327646991338,
            ];
            ddouble[] expected_dist_n_2_m_1 = [
                0.0,
                0.05719095841793662,
                0.1055728090000841,
                0.1471971345775582,
                0.1835034190722737,
                0.2155354594472638,
                0.2440710539815453,
                0.2697032566597785,
                0.2928932188134525,
                0.3140056594299645,
                0.333333333333333,
                0.3511143154769501,
                0.3675444679663244,
                0.3827866001516325,
                0.3969773108444731,
                0.4102321753804118,
                0.4226497308103745,
                0.4343145750507623,
                0.4452998037747711,
                0.4556689460481829,
                0.4654775161751513,
                0.4747742685611102,
                0.4836022205056779,
                0.4919994919992382,
                0.5000000000000002,
                0.5076340360826692,
                0.5149287499273343,
                0.5219085562662427,
                0.5285954792089685,
                0.535009445024723,
                0.5411685322588767,
                0.5470891863421619,
                0.5527864045000423,
                0.5582738957006141,
                0.5635642195280155,
                0.5686689071862465,
                0.5735985672887793,
                0.5783629786442164,
                0.5829711718858508,
                0.5874315014964829,
                0.5917517095361371,
                0.5959389821791159,
                0.6000000000000002,
                0.6039409828093305,
                0.6077677297236324,
                0.6114856550570945,
                0.6150998205402498,
                0.6186149643017632,
                0.6220355269907729,
                0.6253656753673226,
                0.6286093236458963,
                0.6317701528406708,
                0.6348516283298894,
                0.637857015829926,
                0.6407893959464502,
                0.643651677450101,
                0.6464466094067264,
                0.6491767922771884,
                0.6518446880886044,
                0.6544526297674561,
                0.6570028297149825,
                0.6594973876965008,
                0.6619382981085935,
                0.6643274566813245,
                0.6666666666666666,
            ];
            ddouble[] expected_dist_n_1_m_2 = [
                0.0,
                0.1740776559556977,
                0.2425356250363329,
                0.2927700218845598,
                0.3333333333333331,
                0.3676073110469036,
                0.397359707119513,
                0.4236592728681615,
                0.4472135954999578,
                0.468521285665818,
                0.4879500364742662,
                0.5057805388588729,
                0.5222329678670932,
                0.5374838498865698,
                0.5516772843673703,
                0.5649326828660317,
                0.5773502691896255,
                0.5890150893739512,
                0.5999999999999999,
                0.6103679378930735,
                0.6201736729460421,
                0.6294651817966893,
                0.6382847385042251,
                0.6466697906828628,
                0.6546536707079766,
                0.662266178532522,
                0.6695340634119866,
                0.6764814252025458,
                0.683130051063973,
                0.6894996998299376,
                0.6956083436402523,
                0.7014723744122018,
                0.7071067811865475,
                0.7125253031944253,
                0.7177405625652735,
                0.7227641798688512,
                0.7276068751089989,
                0.7322785563281032,
                0.7367883976130075,
                0.7411449079956548,
                0.7453559924999299,
                0.7494290063884895,
                0.7533708035008843,
                0.7571877794400367,
                0.7608859102526825,
                0.7644707871564386,
                0.7679476477883046,
                0.7713214043839051,
                0.7745966692414835,
                0.7777777777777778,
                0.7808688094430304,
                0.7838736067283434,
                0.7867957924694432,
                0.7896387856258748,
                0.7924058156930617,
                0.7950999358860352,
                0.7977240352174659,
                0.800280849578601,
                0.8027729719194866,
                0.8052028616141702,
                0.8075728530872482,
                0.8098851637699161,
                0.8121419014464816,
                0.8143450710459555,
                0.8164965809277264,
            ];
            ddouble[] expected_dist_n_2_m_2 = [
                0.0,
                0.05882352941176468,
                0.111111111111111,
                0.1578947368421051,
                0.1999999999999999,
                0.2380952380952379,
                0.2727272727272726,
                0.3043478260869564,
                0.3333333333333331,
                0.3599999999999998,
                0.3846153846153845,
                0.4074074074074072,
                0.4285714285714284,
                0.4482758620689653,
                0.4666666666666665,
                0.4838709677419352,
                0.4999999999999998,
                0.5151515151515154,
                0.5294117647058827,
                0.542857142857143,
                0.5555555555555558,
                0.5675675675675678,
                0.5789473684210529,
                0.5897435897435901,
                0.6000000000000003,
                0.6097560975609758,
                0.6190476190476193,
                0.6279069767441862,
                0.6363636363636366,
                0.6444444444444446,
                0.6521739130434784,
                0.6595744680851064,
                0.6666666666666669,
                0.6734693877551021,
                0.6800000000000002,
                0.6862745098039217,
                0.6923076923076924,
                0.6981132075471699,
                0.7037037037037038,
                0.7090909090909093,
                0.7142857142857143,
                0.7192982456140352,
                0.7241379310344829,
                0.7288135593220341,
                0.7333333333333334,
                0.737704918032787,
                0.7419354838709679,
                0.7460317460317462,
                0.7500000000000001,
                0.7538461538461539,
                0.7575757575757577,
                0.7611940298507464,
                0.7647058823529412,
                0.7681159420289856,
                0.7714285714285715,
                0.7746478873239437,
                0.7777777777777778,
                0.7808219178082193,
                0.7837837837837839,
                0.7866666666666668,
                0.7894736842105264,
                0.7922077922077924,
                0.7948717948717949,
                0.7974683544303799,
                0.8000000000000002,
            ];
            ddouble[] expected_dist_n_3_m_4 = [
                0.0,
                0.02305061301992792,
                0.05950998812184121,
                0.1002174431705656,
                0.1419928039291337,
                0.1832766799579897,
                0.2232513681251164,
                0.261493081351286,
                0.2978022709324751,
                0.3321113912812833,
                0.3644314500154354,
                0.3948197627096052,
                0.4233600000000004,
                0.4501496586840313,
                0.4752921504417158,
                0.4988918242615086,
                0.5210508807675741,
                0.5418675184386659,
                0.5614348857828648,
                0.5798405606290171,
                0.5971663720116628,
                0.6134884416045314,
                0.6288773622750321,
                0.6433984584623859,
                0.6571120913614144,
                0.6700739842912747,
                0.682335552074555,
                0.6939442240201532,
                0.7049437540428283,
                0.7153745141356442,
                0.7252737692283632,
                0.7346759326820409,
                0.7436128024718239,
                0.7521137786240062,
                0.7602060627888022,
                0.7679148410091915,
                0.7752634508322054,
                0.7822735339326629,
                0.7889651754020064,
                0.7953570308114802,
                0.8014664420998208,
                0.8073095432680274,
                0.8129013567925771,
                0.8182558815968799,
                0.8233861733509911,
                0.8283044178029679,
                0.833021997782513,
                0.8375495544591554,
                0.8418970433832618,
                0.846073785788637,
                0.8500885155902226,
                0.8539494224691936,
                0.85766419140034,
                0.8612400389427338,
                0.864683746584027,
                0.8680016914010269,
                0.8711998742742127,
                0.8742839458713043,
                0.8772592305946668,
                0.8801307486690185,
                0.8829032365293789,
                0.8855811656543152,
                0.8881687599761116,
                0.8906700119873793,
                0.8930886976527022,
            ];
            ddouble[] expected_dist_n_4_m_2 = [
                0.0,
                0.01234567901234569,
                0.04000000000000003,
                0.07438016528925623,
                0.1111111111111112,
                0.1479289940828403,
                0.1836734693877552,
                0.2177777777777778,
                0.2500000000000002,
                0.2802768166089967,
                0.3086419753086422,
                0.3351800554016622,
                0.3599999999999999,
                0.3832199546485259,
                0.4049586776859503,
                0.425330812854442,
                0.4444444444444441,
                0.4623999999999999,
                0.4792899408284021,
                0.4951989026063097,
                0.5102040816326527,
                0.5243757431629011,
                0.5377777777777776,
                0.5504682622268467,
                0.5624999999999998,
                0.5739210284664829,
                0.58477508650519,
                0.595102040816326,
                0.6049382716049378,
                0.6143170197224249,
                0.6232686980609417,
                0.6318211702827085,
                0.64,
                0.6478286734086852,
                0.6553287981859409,
                0.6625202812330987,
                0.6694214876033057,
                0.6760493827160493,
                0.6824196597353496,
                0.6885468537799907,
                0.6944444444444441,
                0.7001249479383589,
                0.7055999999999999,
                0.710880430603614,
                0.7159763313609467,
                0.7208971164115342,
                0.7256515775034293,
                0.7302479338842974,
                0.7346938775510202,
                0.7389966143428746,
                0.7431629013079666,
                0.7471990807239299,
                0.7511111111111111,
                0.7549045955388337,
                0.7585848074921955,
                0.7621567145376668,
                0.7656249999999998,
                0.7689940828402366,
                0.7722681359044995,
                0.7754511026954778,
                0.7785467128027681,
                0.7815584961142618,
                0.7844897959183672,
                0.7873437809958341,
                0.7901234567901232,
            ];
            ddouble[] expected_dist_n_2_m_4 = [
                0.0,
                0.05968778696051429,
                0.1141868512110727,
                0.1640816326530613,
                0.2098765432098767,
                0.2520087655222792,
                0.2908587257617729,
                0.326758711374096,
                0.3600000000000003,
                0.3908387864366452,
                0.4195011337868481,
                0.4461871281773934,
                0.4710743801652895,
                0.4943209876543212,
                0.5160680529300571,
                0.5364418288818473,
                0.555555555555556,
                0.5735110370678889,
                0.5904000000000005,
                0.6063052672049214,
                0.6213017751479294,
                0.635457458170168,
                0.6488340192043894,
                0.6614876033057849,
                0.6734693877551018,
                0.6848261003385656,
                0.6956004756242566,
                0.7058316575696637,
                0.7155555555555555,
                0.7248051599032518,
                0.7336108220603537,
                0.7420005039052657,
                0.7499999999999999,
                0.7576331360946744,
                0.7649219467401285,
                0.7718868344842948,
                0.7785467128027681,
                0.7849191346355806,
                0.7910204081632652,
                0.796865701249752,
                0.802469135802469,
                0.8078438731469317,
                0.8130021913805696,
                0.8179555555555555,
                0.8227146814404431,
                0.8272895935233597,
                0.8316896778435239,
                0.8359237301714467,
                0.84,
                0.8439262307575064,
                0.8477096966091611,
                0.8513572361736099,
                0.8548752834467119,
                0.8582698961937716,
                0.8615467820443482,
                0.8647113224996696,
                0.8677685950413223,
                0.8707233935109202,
                0.8735802469135802,
                0.8763434367829972,
                0.8790170132325141,
                0.8816048098046017,
                0.8841104572204616,
                0.8865373961218835,
                0.8888888888888887,
            ];
            ddouble[] expected_dist_n_4_m_4 = [
                0.0,
                0.00997353958884593,
                0.03429355281207137,
                0.06691937600233276,
                0.1040000000000002,
                0.1430731022567758,
                0.1825694966190836,
                0.2215007808005263,
                0.2592592592592596,
                0.2954880000000003,
                0.3299954483386441,
                0.3626987755931517,
                0.3935860058309043,
                0.4226905572184187,
                0.4500740740740747,
                0.4758148434090839,
                0.5000000000000004,
                0.5227203161087451,
                0.5440667616527574,
                0.5641282798833812,
                0.5829903978052121,
                0.6007344086233781,
                0.6174369441609561,
                0.6331698106846031,
                0.6479999999999999,
                0.6619898144252111,
                0.6751970629521646,
                0.6876752990302738,
                0.6994740796393686,
                0.7106392318244167,
                0.7212131174488366,
                0.7312348901495813,
                0.7407407407407405,
                0.7497641288918729,
                0.7583359999999999,
                0.7664849869205657,
                0.7742375967228037,
                0.7816183829604301,
                0.7886501041507897,
                0.7953538692712244,
                0.801749271137026,
                0.8078545085397393,
                0.8136864980113983,
                0.8192609760491576,
                0.8245925925925923,
                0.8296949964975039,
                0.8345809136987679,
                0.8392622187028836,
                0.8437499999999999,
                0.8480546199362766,
                0.8521857695411413,
                0.8561525187606186,
                0.8599633625076325,
                0.8636262629029948,
                0.867148688046647,
                0.8705376476274826,
                0.8737997256515773,
                0.8769411105427268,
                0.8799676228456359,
                0.8828847407407407,
                0.8856976235602856,
                0.888411133477755,
                0.8910298555268967,
                0.893558116092236,
                0.8959999999999998,
            ];
            ddouble[] expected_dist_n_6_m_8 = [
                0.0,
                0.001620949186622189,
                0.01032835774209723,
                0.02807278777034187,
                0.05413730056410238,
                0.08683340150054124,
                0.1242926931800777,
                0.1648080738536674,
                0.206954770397406,
                0.2496080609181597,
                0.2919164826900167,
                0.3332594877599288,
                0.3732032102399989,
                0.411460244551413,
                0.4478554818353038,
                0.48229822441511,
                0.5147600064598922,
                0.5452572743648657,
                0.5738380523823595,
                0.6005717984191429,
                0.6255417716966075,
                0.6488393545929863,
                0.6705598807659474,
                0.6907996154354564,
                0.7096536109332309,
                0.7272142227620843,
                0.7435701206387811,
                0.7588056675845101,
                0.7730005701476963,
                0.7862297260675103,
                0.7985632135758549,
                0.8100663802697065,
                0.8208000000000004,
                0.8308204742522153,
                0.8401800606099317,
                0.8489271155396491,
                0.8571063422559013,
                0.8647590370842954,
                0.8719233297396002,
                0.8786344144323892,
                0.8849247698306458,
                0.8908243667238681,
                0.8963608628368551,
                0.9015597846722743,
                0.9064446965664839,
                0.9110373573533446,
                0.9153578651699359,
                0.9194247910245136,
                0.9232553017946227,
                0.9268652733426388,
                0.9302693944351376,
                0.9334812621374579,
                0.9365134693301559,
                0.9393776849631611,
                0.942084727628823,
                0.9446446329985452,
                0.947066715630674,
                0.9493596256206677,
                0.9515314005290139,
                0.9535895129882881,
                0.9555409143584583,
                0.9573920747691856,
                0.9591490198595174,
                0.9608173644990168,
                0.9624023437500001,
            ];
            ddouble[] expected_dist_n_8_m_6 = [
                0.0,
                4.626244792366958e-4,
                0.004904419077085219,
                0.01695999999999997,
                0.03759765624999992,
                0.06589830980907797,
                0.1001371742112481,
                0.1384587313818182,
                0.1791999999999996,
                0.221004853419918,
                0.2628331736812896,
                0.3039249016162554,
                0.3437499999999992,
                0.3819588198399991,
                0.4183386995424923,
                0.4527783893226143,
                0.4852399935401074,
                0.5157373999463624,
                0.5443200000000008,
                0.5710605824518273,
                0.5960464477539071,
                0.6193729710690182,
                0.6411390061691808,
                0.6614436610595932,
                0.6803840877914956,
                0.698054016483408,
                0.7145428310716514,
                0.7299350357555069,
                0.7443100000000009,
                0.7577418992056184,
                0.7702997900534644,
                0.7820477758799267,
                0.7930452296025939,
                0.8033470507544586,
                0.8130039398756883,
                0.8220626784549137,
                0.8305664062500003,
                0.8385548904845677,
                0.8460647833600003,
                0.8531298657292291,
                0.8597812757869644,
                0.8660477223448874,
                0.8719556827568818,
                0.8775295858962803,
                0.8827919808073169,
                0.8877636917884483,
                0.8924639607396969,
                0.8969105776373268,
                0.9011200000000003,
                0.9051074621907096,
                0.9088870753652684,
                0.9124719188361103,
                0.9158741235733034,
                0.9191049485156758,
                0.9221748503156258,
                0.9250935470928411,
                0.9278700767256223,
                0.9305128501643275,
                0.9330297002099467,
                0.9354279261621218,
                0.9377143347050755,
                0.939895277366868,
                0.9419766848570488,
                0.9439640985600002,
                0.9458626994358975,
            ];
            ddouble[] expected_dist_n_8_m_8 = [
                0.0,
                3.627467011865095e-4,
                0.004039541130205942,
                0.01454611711356498,
                0.03334400000000003,
                0.06009039674413202,
                0.0934524051746659,
                0.1317430305847922,
                0.1732967535436672,
                0.2166516523007998,
                0.2606137289268526,
                0.3042579803220433,
                0.3468999190084789,
                0.3880565363536832,
                0.4274065441243711,
                0.4647544331270105,
                0.4999999999999998,
                0.533113520952927,
                0.5641160734562304,
                0.5930642554134999,
                0.6200385158256304,
                0.6451343750761104,
                0.668455914792044,
                0.6901110254515995,
                0.7102080000000003,
                0.7288531481046395,
                0.7461491772380334,
                0.7621941445409169,
                0.7770808292286244,
                0.7908964111822594,
                0.8037223692072578,
                0.8156345338944682,
                0.8267032464563329,
                0.8369935874466867,
                0.8465656487936002,
                0.8554748297773941,
                0.863772143013356,
                0.871504520570951,
                0.8787151134020843,
                0.8854435795085907,
                0.8917263579436653,
                0.8975969269585451,
                0.9030860454849886,
                0.9082219777701771,
                0.9130307014174669,
                0.9175360993824988,
                0.9217601366662078,
                0.9257230225620321,
                0.9294433593750001,
                0.932938278551125,
                0.9362235651483829,
                0.9393137715543226,
                0.942222321316519,
                0.9449616039054411,
                0.9475430611783476,
                0.949977266260055,
                0.9522739955036297,
                0.9544422941424641,
                0.9564905361956212,
                0.9584264791413174,
                0.9602573138292272,
                0.9619897100611166,
                0.9636298582311298,
                0.965183507381862,
                0.9666560000000001,
            ];
            ddouble[] expected_dist_n_10_m_10 = [
                0.0,
                7.261249003154098e-5,
                0.001449280603225923,
                0.007061080353179447,
                0.01958144000000006,
                0.04023644212462403,
                0.06882895817386492,
                0.1042240640360383,
                0.1448458060255045,
                0.1890359574803256,
                0.2352659976907717,
                0.2822381056655883,
                0.3289149096386856,
                0.3745089896200322,
                0.4184529207783381,
                0.4603624792899131,
                0.5000000000000016,
                0.5372413167221974,
                0.5720476174536777,
                0.6044423943969748,
                0.6344930946592232,
                0.6622968391699028,
                0.6879695229132687,
                0.7116376479966895,
                0.7334323199999994,
                0.753484928764521,
                0.7719241225965681,
                0.7888737631252122,
                0.8044516145290513,
                0.8187685755783538,
                0.8319283070419405,
                0.8440271419934184,
                0.8551541939744953,
                0.8653915992660786,
                0.8748148459215255,
                0.8834931547673989,
                0.8914898871200109,
                0.898862961177678,
                0.9056652644642929,
                0.9119450537391365,
                0.9177463367772798,
                0.9231092326170303,
                0.9280703084601383,
                0.9326628925460534,
                0.9369173631174108,
                0.940861414136959,
                0.9445202987722735,
                0.9479170518837408,
                0.9510726928710936,
                0.954006410282531,
                0.9567357295893016,
                0.9592766654933448,
                0.9616438600777619,
                0.9638507080381098,
                0.9659094701529228,
                0.9678313760689414,
                0.969626717393359,
                0.9713049320041464,
                0.9728746804115841,
                0.9743439149304312,
                0.9757199423531606,
                0.9770094807506509,
                0.9782187109676483,
                0.9793533233261086,
                0.9804185599999999,
            ];

            foreach ((SnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_n_1_m_1,   expected_dist_n_1_m_1),
                 (dist_n_2_m_1,   expected_dist_n_2_m_1),
                 (dist_n_1_m_2,   expected_dist_n_1_m_2),
                 (dist_n_2_m_2,   expected_dist_n_2_m_2),
                 (dist_n_3_m_4,   expected_dist_n_3_m_4),
                 (dist_n_4_m_2,   expected_dist_n_4_m_2),
                 (dist_n_2_m_4,   expected_dist_n_2_m_4),
                 (dist_n_4_m_4,   expected_dist_n_4_m_4),
                 (dist_n_6_m_8,   expected_dist_n_6_m_8),
                 (dist_n_8_m_6,   expected_dist_n_8_m_6),
                 (dist_n_8_m_8,   expected_dist_n_8_m_8),
                 (dist_n_10_m_10, expected_dist_n_10_m_10),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
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