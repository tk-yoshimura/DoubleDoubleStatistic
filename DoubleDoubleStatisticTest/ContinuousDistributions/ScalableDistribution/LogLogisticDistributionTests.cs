﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class LogLogisticDistributionTests {
        readonly LogLogisticDistribution dist_gamma1sigma1 = new(gamma: 1, sigma: 1);
        readonly LogLogisticDistribution dist_gamma1sigma2 = new(gamma: 1, sigma: 2);
        readonly LogLogisticDistribution dist_gamma2sigma1 = new(gamma: 2, sigma: 1);
        readonly LogLogisticDistribution dist_gamma2sigma2 = new(gamma: 2, sigma: 2);
        readonly LogLogisticDistribution dist_gamma4sigma3 = new(gamma: 4, sigma: 3);
        readonly LogLogisticDistribution dist_gamma5sigma4 = new(gamma: 5, sigma: 4);
        readonly LogLogisticDistribution dist_gamma6sigma5 = new(gamma: 6, sigma: 5);

        LogLogisticDistribution[] Dists => [
            dist_gamma1sigma1,
            dist_gamma1sigma2,
            dist_gamma2sigma1,
            dist_gamma2sigma2,
            dist_gamma4sigma3,
            dist_gamma5sigma4,
            dist_gamma6sigma5,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Sigma={dist.Sigma}");
                Console.WriteLine($"Gamma={dist.Gamma}");
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
            foreach (LogLogisticDistribution dist in Dists) {
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
            foreach (LogLogisticDistribution dist in Dists) {
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
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (LogLogisticDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (LogLogisticDistribution dist in Dists) {
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
            foreach (LogLogisticDistribution dist in Dists) {
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
            foreach (LogLogisticDistribution dist in Dists) {
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

            foreach (LogLogisticDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.02, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (LogLogisticDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                (LogLogisticDistribution? dist_fit, ddouble error) = LogLogisticDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (LogLogisticDistribution dist in Dists) {
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
            ddouble[] expected_dist_gamma1sigma1 = [
                1.0,
                0.8858131487889274,
                0.7901234567901234,
                0.7091412742382271,
                0.64,
                0.5804988662131519,
                0.5289256198347108,
                0.4839319470699433,
                0.4444444444444444,
                0.4096,
                0.378698224852071,
                0.3511659807956104,
                0.3265306122448979,
                0.3043995243757431,
                0.2844444444444444,
                0.2663891779396462,
                0.25,
                0.2350780532598714,
                0.2214532871972318,
                0.2089795918367347,
                0.1975308641975309,
                0.1869978086194302,
                0.1772853185595568,
                0.168310322156476,
                0.16,
                0.1522903033908388,
                0.145124716553288,
                0.1384532179556517,
                0.1322314049586777,
                0.1264197530864198,
                0.1209829867674858,
                0.1158895427795383,
                0.1111111111111111,
                0.1066222407330279,
                0.1024,
                0.09842368319876971,
                0.09467455621301775,
                0.09113563545745818,
                0.0877914951989026,
                0.08462809917355373,
                0.08163265306122448,
                0.07879347491535857,
                0.07609988109393578,
                0.07354208560758403,
                0.07111111111111111,
                0.06879871002418704,
                0.06659729448491156,
                0.06449987402368355,
                0.0625,
                0.06059171597633136,
                0.05876951331496786,
                0.05702829137892627,
                0.05536332179930796,
                0.05377021634110481,
                0.05224489795918368,
                0.05078357468756199,
                0.04938271604938271,
                0.04803903171326703,
                0.04674945215485756,
                0.04551111111111111,
                0.0443213296398892,
                0.04317760161916006,
                0.042077580539119,
                0.04101906745713828,
                0.04,
            ];
            ddouble[] expected_dist_gamma1sigma2 = [
                0.5,
                0.4701561065197429,
                0.4429065743944637,
                0.4179591836734694,
                0.3950617283950617,
                0.3739956172388605,
                0.3545706371191136,
                0.336620644312952,
                0.32,
                0.3045806067816776,
                0.290249433106576,
                0.2769064359113034,
                0.2644628099173554,
                0.2528395061728395,
                0.2419659735349717,
                0.2317790855590765,
                0.2222222222222222,
                0.2132444814660558,
                0.2048,
                0.1968473663975394,
                0.1893491124260355,
                0.1822712709149164,
                0.1755829903978052,
                0.1692561983471075,
                0.163265306122449,
                0.1575869498307171,
                0.1521997621878716,
                0.1470841712151681,
                0.1422222222222222,
                0.1375974200483741,
                0.1331945889698231,
                0.1289997480473671,
                0.125,
                0.1211834319526627,
                0.1175390266299357,
                0.1140565827578525,
                0.1107266435986159,
                0.1075404326822096,
                0.1044897959183674,
                0.101567149375124,
                0.09876543209876543,
                0.09607806342653406,
                0.09349890430971512,
                0.09102222222222223,
                0.0886426592797784,
                0.08635520323832012,
                0.084155161078238,
                0.08203813491427656,
                0.08,
                0.07803688462124676,
                0.0761451516954194,
                0.07432138191319494,
                0.07256235827664399,
                0.07086505190311419,
                0.06922660897782586,
                0.06764433875016515,
                0.06611570247933884,
                0.06463830324453983,
                0.06320987654320988,
                0.06182828160850139,
                0.06049149338374291,
                0.05919759509769915,
                0.05794477138976913,
                0.05673130193905817,
                0.05555555555555555,
            ];
            ddouble[] expected_dist_gamma2sigma1 = [
                0,
                0.1240291298884162,
                0.2423668639053254,
                0.3499608401566394,
                0.4429065743944637,
                0.5187370980610682,
                0.5764683805592044,
                0.616436441816716,
                0.64,
                0.6491912405674083,
                0.6463830324453983,
                0.6340155773979976,
                0.6144000000000001,
                0.58959723183391,
                0.5613595426423369,
                0.5311180363155416,
                0.5,
                0.4688628903290969,
                0.4383353151010702,
                0.4088586746661974,
                0.380725758477097,
                0.3541144770887324,
                0.3291161431701972,
                0.3057584486186052,
                0.2840236686390533,
                0.2638627810467158,
                0.2452062112030061,
                0.2279718621969131,
                0.2120710059171598,
                0.1974125172738446,
                0.1839058440392237,
                0.1714630248418562,
                0.16,
                0.1494374041265322,
                0.1397009846800793,
                0.1307217553334813,
                0.1224359655648847,
                0.1147849467455621,
                0.1077148788927336,
                0.1011765100805041,
                0.09512485136741973,
                0.08951886329551231,
                0.0843211449857857,
                0.07949763316614102,
                0.07501731578666951,
                0.07085196294955892,
                0.06697587651322767,
                0.06336565877662528,
                0.06,
                0.05685948409049799,
                0.05392641151328886,
                0.05118463833587199,
                0.04861943024105186,
                0.04621733032794084,
                0.04396603954208403,
                0.04185430862220495,
                0.03987184051263795,
                0.03800920226047566,
                0.03625774549006441,
                0.03460953462145714,
                0.03305728207158967,
                0.03159428974587282,
                0.03021439619274242,
                0.02891192885403172,
                0.02768166089965398,
            ];
            ddouble[] expected_dist_gamma2sigma2 = [
                0,
                0.03118905413444378,
                0.06201456494420809,
                0.09212352484188291,
                0.1211834319526627,
                0.1488911769436778,
                0.1749804200783197,
                0.1992271540107128,
                0.2214532871972318,
                0.241528224237833,
                0.2593685490305341,
                0.274936023340516,
                0.2882341902796022,
                0.2993039166020844,
                0.308218220908358,
                0.3150767211046659,
                0.32,
                0.3231241397032081,
                0.3245956202837041,
                0.3245667218392003,
                0.3231915162226991,
                0.3206224883225198,
                0.3170077886989988,
                0.3124890901393933,
                0.3072,
                0.3012649671723422,
                0.294798615916955,
                0.2879054373091651,
                0.2806797713211684,
                0.2732060174370548,
                0.2655590181577708,
                0.2578045669980775,
                0.25,
                0.2421948369557305,
                0.2344314451645484,
                0.2267457044789038,
                0.2191676575505351,
                0.2117221335595585,
                0.2044293373330987,
                0.1973053986977046,
                0.1903628792385485,
                0.1836112354406333,
                0.1770572385443662,
                0.1707053524363382,
                0.1645580715850986,
                0.1586162214818401,
                0.1528792243093026,
                0.1473453326826306,
                0.1420118343195266,
                0.136875230433161,
                0.1319313905233579,
                0.1271756860879905,
                0.1226031056015031,
                0.118208352921806,
                0.1139859310984565,
                0.1099302133697962,
                0.1060355029585799,
                0.1022960831072427,
                0.09870625863692228,
                0.09526039016951199,
                0.09195292201961183,
                0.08877840464312164,
                0.08573151242092811,
                0.08280705745905273,
                0.08,
            ];
            ddouble[] expected_dist_gamma4sigma3 = [
                0,
                1.205632261815059e-5,
                9.645003586648087e-5,
                3.255108994532398e-4,
                7.715305218821825e-4,
                0.001506686087058777,
                0.002602895566421209,
                0.004131581936447937,
                0.006163324533396266,
                0.008767376819269222,
                0.01201103174677983,
                0.01595881887055681,
                0.02067152164806936,
                0.02620500929136825,
                0.03260888531528643,
                0.03992496474309808,
                0.04818560380725757,
                0.05741191975951572,
                0.06761195366924072,
                0.07877884512610826,
                0.09088910352273932,
                0.1039010746881262,
                0.1177537123998442,
                0.13236576986814,
                0.1476355247981546,
                0.1634411414684366,
                0.1796417532966499,
                0.1960793192820127,
                0.2125812682601149,
                0.2289638980361508,
                0.2450364453575083,
                0.2606056915670161,
                0.2754809225209905,
                0.2894790249362655,
                0.3024294791586047,
                0.3141790036084425,
                0.3245956202837041,
                0.3335719430389362,
                0.3410275381977395,
                0.3469102659173124,
                0.3511965749662421,
                0.3538907871427561,
                0.3550234647711611,
                0.3546490009568229,
                0.3528426044905796,
                0.3496968681694126,
                0.3453181112436396,
                0.3398226755368858,
                0.3333333333333333,
                0.3259759367173535,
                0.3178764060655656,
                0.3091581228925036,
                0.2999397617279073,
                0.2903335689467773,
                0.280444074559488,
                0.2703672063303911,
                0.2601897641384953,
                0.2499892057533781,
                0.2398336924684581,
                0.2297823434973679,
                0.2198856508722632,
                0.2101860110201296,
                0.2007183345724436,
                0.1915107017463013,
                0.1825850364095836,
            ];
            ddouble[] expected_dist_gamma5sigma4 = [
                0,
                7.45058058304604e-8,
                1.192092824453542e-6,
                6.034967551944254e-6,
                1.907344994838897e-5,
                4.656585768141409e-5,
                9.655812599005835e-5,
                1.788828400756735e-4,
                3.051571556511243e-4,
                4.887788320055939e-4,
                7.449193011989432e-4,
                0.001090512347825994,
                0.001544236578046494,
                0.002126489417355512,
                0.00285934989873469,
                0.003766526993666184,
                0.004873289708506841,
                0.006206374632517117,
                0.007793866101123726,
                0.009665043654609913,
                0.01185019106136511,
                0.01438036087711176,
                0.01728708837682898,
                0.0206020487828626,
                0.02435665208763037,
                0.02858157050576695,
                0.03330619476572648,
                0.03855801714140762,
                0.0443619414000013,
                0.05073952275738945,
                0.05770814451564217,
                0.06528014229896427,
                0.07346189164370982,
                0.08225288000847644,
                0.09164478984817204,
                0.1016206249519847,
                0.1121539174022507,
                0.1232080568113771,
                0.1347357864172471,
                0.1466789116143347,
                0.1589682650356235,
                0.1715239679235985,
                0.1842560199264155,
                0.1970652385326057,
                0.2098445552973305,
                0.2224806593204605,
                0.2348559599528351,
                0.2468508215833914,
                0.258346004987264,
                0.2692252336093213,
                0.2793777908104795,
                0.2887010468151301,
                0.297102812810899,
                0.3045034248194429,
                0.3108374714570588,
                0.3160550968239891,
                0.3201228312568074,
                0.323023926899319,
                0.3247582001170471,
                0.3253414067983099,
                0.3248041978274817,
                0.3231907191048215,
                0.3205569325049039,
                0.3169687407078312,
                0.3125,
            ];

            foreach ((LogLogisticDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1sigma1, expected_dist_gamma1sigma1),
                (dist_gamma1sigma2, expected_dist_gamma1sigma2),
                (dist_gamma2sigma1, expected_dist_gamma2sigma1),
                (dist_gamma2sigma2, expected_dist_gamma2sigma2),
                (dist_gamma4sigma3, expected_dist_gamma4sigma3),
                (dist_gamma5sigma4, expected_dist_gamma5sigma4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
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
            ddouble[] expected_dist_gamma1sigma1 = [
                0,
                0.05882352941176471,
                0.1111111111111111,
                0.1578947368421053,
                0.2,
                0.2380952380952381,
                0.2727272727272728,
                0.3043478260869565,
                0.3333333333333333,
                0.36,
                0.3846153846153846,
                0.4074074074074074,
                0.4285714285714286,
                0.4482758620689655,
                0.4666666666666667,
                0.4838709677419356,
                0.5,
                0.5151515151515151,
                0.5294117647058824,
                0.5428571428571429,
                0.5555555555555556,
                0.5675675675675675,
                0.5789473684210527,
                0.5897435897435898,
                0.6000000000000001,
                0.6097560975609756,
                0.6190476190476191,
                0.627906976744186,
                0.6363636363636364,
                0.6444444444444445,
                0.6521739130434783,
                0.6595744680851063,
                0.6666666666666666,
                0.673469387755102,
                0.6799999999999999,
                0.6862745098039216,
                0.6923076923076923,
                0.6981132075471698,
                0.7037037037037037,
                0.7090909090909091,
                0.7142857142857143,
                0.7192982456140351,
                0.7241379310344828,
                0.728813559322034,
                0.7333333333333333,
                0.7377049180327868,
                0.7419354838709677,
                0.746031746031746,
                0.75,
                0.7538461538461539,
                0.7575757575757576,
                0.7611940298507462,
                0.7647058823529411,
                0.7681159420289856,
                0.7714285714285715,
                0.7746478873239436,
                0.7777777777777779,
                0.7808219178082192,
                0.7837837837837838,
                0.7866666666666667,
                0.7894736842105263,
                0.7922077922077921,
                0.7948717948717948,
                0.7974683544303798,
                0.8,
            ];
            ddouble[] expected_dist_gamma1sigma2 = [
                0,
                0.0303030303030303,
                0.05882352941176471,
                0.08571428571428572,
                0.1111111111111111,
                0.1351351351351351,
                0.1578947368421053,
                0.1794871794871795,
                0.2,
                0.2195121951219512,
                0.2380952380952381,
                0.2558139534883721,
                0.2727272727272728,
                0.2888888888888889,
                0.3043478260869565,
                0.3191489361702128,
                0.3333333333333333,
                0.3469387755102041,
                0.36,
                0.3725490196078431,
                0.3846153846153846,
                0.3962264150943396,
                0.4074074074074074,
                0.4181818181818182,
                0.4285714285714286,
                0.4385964912280702,
                0.4482758620689655,
                0.4576271186440678,
                0.4666666666666667,
                0.4754098360655738,
                0.4838709677419356,
                0.4920634920634921,
                0.5,
                0.5076923076923077,
                0.5151515151515151,
                0.5223880597014926,
                0.5294117647058824,
                0.536231884057971,
                0.5428571428571429,
                0.5492957746478874,
                0.5555555555555556,
                0.5616438356164384,
                0.5675675675675675,
                0.5733333333333334,
                0.5789473684210527,
                0.5844155844155844,
                0.5897435897435898,
                0.5949367088607594,
                0.6000000000000001,
                0.6049382716049383,
                0.6097560975609756,
                0.6144578313253012,
                0.6190476190476191,
                0.6235294117647059,
                0.627906976744186,
                0.632183908045977,
                0.6363636363636364,
                0.6404494382022472,
                0.6444444444444445,
                0.6483516483516484,
                0.6521739130434783,
                0.6559139784946236,
                0.6595744680851063,
                0.6631578947368422,
                0.6666666666666666,
            ];
            ddouble[] expected_dist_gamma2sigma1 = [
                0,
                0.003891050583657588,
                0.01538461538461539,
                0.0339622641509434,
                0.05882352941176471,
                0.08896797153024912,
                0.1232876712328767,
                0.160655737704918,
                0.2,
                0.2403560830860534,
                0.2808988764044944,
                0.3209549071618037,
                0.36,
                0.3976470588235295,
                0.4336283185840709,
                0.4677754677754677,
                0.5,
                0.5302752293577981,
                0.5586206896551724,
                0.5850891410048623,
                0.6097560975609756,
                0.6327116212338595,
                0.654054054054054,
                0.6738853503184714,
                0.6923076923076923,
                0.7094211123723042,
                0.7253218884120172,
                0.7401015228426395,
                0.7538461538461539,
                0.7666362807657248,
                0.7785467128027681,
                0.7896466721446178,
                0.8,
                0.8096654275092936,
                0.8186968838526911,
                0.8271438217420661,
                0.8350515463917526,
                0.8424615384615385,
                0.8494117647058823,
                0.8559369724254361,
                0.8620689655172414,
                0.8678368611254518,
                0.8732673267326732,
                0.8783847980997624,
                0.8832116788321168,
                0.8877685225778168,
                0.8920741989881956,
                0.8961460446247463,
                0.8999999999999999,
                0.9036507339104253,
                0.9071117561683599,
                0.9103955197759888,
                0.9135135135135135,
                0.9164763458401305,
                0.9192938209331653,
                0.9219750076196281,
                0.9245283018867924,
                0.9269614835948645,
                0.9292817679558012,
                0.9314958522879314,
                0.9336099585062241,
                0.9356298717626351,
                0.937560975609756,
                0.9394082840236687,
                0.9411764705882353,
            ];
            ddouble[] expected_dist_gamma2sigma2 = [
                0,
                9.75609756097561e-4,
                0.003891050583657588,
                0.008712487899322363,
                0.01538461538461539,
                0.02383222116301239,
                0.0339622641509434,
                0.0456663560111836,
                0.05882352941176471,
                0.07330316742081448,
                0.08896797153024912,
                0.1056768558951965,
                0.1232876712328767,
                0.1416596814752724,
                0.160655737704918,
                0.1801441152922338,
                0.2,
                0.2201066260472201,
                0.2403560830860534,
                0.2606498194945848,
                0.2808988764044944,
                0.3010238907849829,
                0.3209549071618037,
                0.3406310367031552,
                0.36,
                0.3790175864160097,
                0.3976470588235295,
                0.4158585282373075,
                0.4336283185840709,
                0.4509383378016086,
                0.4677754677754677,
                0.4841309823677581,
                0.5,
                0.5153809749171794,
                0.5302752293577981,
                0.5446865273454868,
                0.5586206896551724,
                0.572085248641872,
                0.5850891410048623,
                0.5976424361493123,
                0.6097560975609756,
                0.621441774491682,
                0.6327116212338595,
                0.6435781413156979,
                0.654054054054054,
                0.6641521810429649,
                0.6738853503184714,
                0.6832663161150635,
                0.6923076923076923,
                0.701021897810219,
                0.7094211123723042,
                0.7175172413793104,
                0.7253218884120172,
                0.7328463344638664,
                0.7401015228426395,
                0.7470980489009632,
                0.7538461538461539,
                0.760355721975193,
                0.7666362807657248,
                0.7726970033296338,
                0.7785467128027681,
                0.7841938883034774,
                0.7896466721446178,
                0.7949128780292409,
                0.8,
            ];
            ddouble[] expected_dist_gamma4sigma3 = [
                0,
                1.883800763956562e-7,
                3.014072705461801e-6,
                1.525855623540901e-5,
                4.822298307373293e-5,
                1.177237094232273e-4,
                2.44081034903588e-4,
                4.520961652416189e-4,
                7.710100231303005e-4,
                0.00123443619793651,
                0.001880259084659793,
                0.002750487173066153,
                0.003891050583657588,
                0.005351531400641223,
                0.007184815232646169,
                0.009446652861957951,
                0.01219512195121951,
                0.01548998068783073,
                0.01939190806799138,
                0.02396162932680878,
                0.0292589298253827,
                0.03534156645127103,
                0.0422640921201904,
                0.05007661601819671,
                0.05882352941176471,
                0.06854223368457953,
                0.07926191315352016,
                0.09100239954505734,
                0.1037731771621212,
                0.1175725772092577,
                0.1323872061004024,
                0.1481916457114376,
                0.1649484536082474,
                0.182608478740786,
                0.2011114937020975,
                0.220387129406329,
                0.2403560830860534,
                0.2609315570163744,
                0.2820208744051573,
                0.303527211281975,
                0.3253513794898489,
                0.3473935961806585,
                0.3695551793135293,
                0.391740116007199,
                0.4138564604121322,
                0.4358175291190675,
                0.4575428740535334,
                0.4789590244380281,
                0.5,
                0.5206076066241635,
                0.5407315327636588,
                0.5603292700108012,
                0.5793658843337323,
                0.597813665793437,
                0.6156516843389322,
                0.6328652778562561,
                0.6494454963483904,
                0.6653885232593167,
                0.6806950917995837,
                0.6953699109084202,
                0.7094211123723042,
                0.7228597277357195,
                0.7356992010655644,
                0.7479549414076772,
                0.7596439169139466,
            ];
            ddouble[] expected_dist_gamma5sigma4 = [
                0,
                9.313225737481168e-10,
                2.980232149951692e-8,
                2.263113344147296e-7,
                9.536734069124155e-7,
                2.910374575368549e-6,
                7.241911894542264e-6,
                1.565249350717443e-5,
                3.051664683084623e-5,
                5.499064257139973e-5,
                9.312358465188632e-5,
                1.499679382085647e-4,
                2.316891665768831e-4,
                3.456740208489605e-4,
                5.006368695537272e-4,
                7.067232690914846e-4,
                9.75609756097561e-4,
                0.001320598590045919,
                0.001756705888330935,
                0.002300741278399988,
                0.002971376845700441,
                0.003789202788010504,
                0.004776766804734787,
                0.005958593865932118,
                0.007361182636090999,
                0.00901297449225012,
                0.01094429079997553,
                0.01318723392378572,
                0.01577554738530651,
                0.01874443068111271,
                0.02213030457389936,
                0.02597052321701413,
                0.0303030303030303,
                0.03516595757571341,
                0.04059716553494483,
                0.04663372800404889,
                0.05331136440582327,
                0.06066382606853545,
                0.06872224558666648,
                0.07751446109007257,
                0.08706433009221853,
                0.09739105022471771,
                0.1085085064287708,
                0.1204246658583774,
                0.1331410426519911,
                0.1466522546665382,
                0.1609456931004713,
                0.1760013235833246,
                0.191791633780584,
                0.2082817379485976,
                0.225429643363581,
                0.2431866774184902,
                0.2614980677794924,
                0.2803036617120825,
                0.2995387649282653,
                0.3191350754471518,
                0.3390216843166919,
                0.3591261128376333,
                0.3793753552745507,
                0.3996968969261347,
                0.4200196797368955,
                0.440274991146308,
                0.4603972562983384,
                0.4803247187419613,
                0.5,

            ];

            foreach ((LogLogisticDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1sigma1, expected_dist_gamma1sigma1),
                (dist_gamma1sigma2, expected_dist_gamma1sigma2),
                (dist_gamma2sigma1, expected_dist_gamma2sigma1),
                (dist_gamma2sigma2, expected_dist_gamma2sigma2),
                (dist_gamma4sigma3, expected_dist_gamma4sigma3),
                (dist_gamma5sigma4, expected_dist_gamma5sigma4),
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