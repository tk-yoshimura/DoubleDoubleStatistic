﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class InverseGammaDistributionTests {
        readonly InverseGammaDistribution dist_alpha1beta1 = new(alpha: 1, beta: 1);
        readonly InverseGammaDistribution dist_alpha2beta1 = new(alpha: 2, beta: 1);
        readonly InverseGammaDistribution dist_alpha1beta2 = new(alpha: 1, beta: 2);
        readonly InverseGammaDistribution dist_alpha2beta2 = new(alpha: 2, beta: 2);
        readonly InverseGammaDistribution dist_alpha3beta4 = new(alpha: 3, beta: 4);
        readonly InverseGammaDistribution dist_alpha5beta6 = new(alpha: 5, beta: 6);
        readonly InverseGammaDistribution dist_alpha7beta8 = new(alpha: 7, beta: 8);

        InverseGammaDistribution[] Dists => [
            dist_alpha1beta1,
            dist_alpha2beta1,
            dist_alpha1beta2,
            dist_alpha2beta2,
            dist_alpha3beta4,
            dist_alpha5beta6,
            dist_alpha7beta8,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (InverseGammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (InverseGammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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
            foreach (InverseGammaDistribution dist in Dists) {
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

            foreach (InverseGammaDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
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

            foreach (InverseGammaDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 20000).ToArray();

                (InverseGammaDistribution? dist_fit, ddouble error) = InverseGammaDistribution.Fit(xs, (0.05, 0.90));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (InverseGammaDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha1beta1 = [
                0,
                2.880900472813034e-5,
                0.02146960818576076,
                0.1373283553800943,
                0.2930502222197469,
                0.4174049687384701,
                0.4941045420288109,
                0.5313378863241238,
                0.5413411329464508,
                0.5341655400488013,
                0.5168550860663178,
                0.4940302367543291,
                0.4686171344279587,
                0.4424222654733848,
                0.4165310136476348,
                0.3915705308335359,
                0.3678794411714423,
                0.3456164260087822,
                0.3248294641044444,
                0.3054999154862698,
                0.2875705370350218,
                0.2709632183737651,
                0.2555901255880069,
                0.2413606089133221,
                0.2281853862367075,
                0.2159789768880327,
                0.204661016420936,
                0.1941568619031186,
                0.1843977541249826,
                0.1753207100486841,
                0.1668682579939646,
                0.1589880881683073,
                0.1516326649281584,
                0.144758829907436,
                0.1383274138132095,
                0.1323028672656575,
                0.1266529162330775,
                0.1213482445111999,
                0.1163622037399827,
                0.1116705502539737,
                0.1072512073657023,
                0.1030840513122371,
                0.09915071893695394,
                0.0954344351540591,
                0.0919198583006782,
                0.08859294158601898,
                0.08544080897667222,
                0.08245164399692298,
                0.07961459006375435,
                0.07691966111237715,
                0.07435766139634596,
                0.07192011346493676,
                0.06959919342879854,
                0.06738767272293118,
                0.06527886566416735,
                0.06326658217914784,
                0.06134508514900293,
                0.05950905187934784,
                0.05775353925953209,
                0.05607395222407881,
                0.05446601517259723,
                0.05292574604277619,
                0.05144943276494369,
                0.0500336118566202,
                0.0486750489419628,
            ];
            ddouble[] expected_dist_alpha2beta1 = [
                0,
                4.609440756500853e-4,
                0.1717568654860861,
                0.7324178953605031,
                1.172200888878987,
                1.335695899963104,
                1.317612112076829,
                1.214486597312283,
                1.082682265892902,
                0.9496276267534247,
                0.8269681377061087,
                0.7185894352790243,
                0.6248228459039449,
                0.5445197113518583,
                0.4760354441687255,
                0.4176752328891049,
                0.3678794411714423,
                0.3252860480082656,
                0.2887373014261728,
                0.2572630867252798,
                0.2300564296280175,
                0.2064481663800115,
                0.1858837277003686,
                0.1679030322875284,
                0.1521235908244717,
                0.1382265452083409,
                0.1259452408744222,
                0.115055918164811,
                0.1053701452142758,
                0.09672866761306709,
                0.08899640426344778,
                0.0820583680868683,
                0.07581633246407918,
                0.07018609934905987,
                0.06509525355915745,
                0.06048131075001484,
                0.05629018499247886,
                0.05247491654538374,
                0.04899461210104536,
                0.04581355907855329,
                0.04290048294628093,
                0.04022792246331203,
                0.03777170245217292,
                0.03551048749918478,
                0.03342540301842843,
                0.03149971256391786,
                0.02971854225275555,
                0.02806864476490995,
                0.02653819668791812,
                0.02511662403669458,
                0.02379445164683071,
                0.02256317285174487,
                0.02141513643963032,
                0.02034344836918677,
                0.01934188612271625,
                0.01840482390666119,
                0.01752716718542941,
                0.01670429526437834,
                0.01593201083021575,
                0.01520649551839426,
                0.0145242707126926,
                0.01388216289646589,
                0.01327727297159837,
                0.01270694904295116,
                0.0121687622354907,
            ];
            ddouble[] expected_dist_alpha1beta2 = [
                0,
                6.484052761136218e-12,
                1.440450236406517e-5,
                0.001326028865020417,
                0.01073480409288038,
                0.03402869295460217,
                0.06866417769004715,
                0.1080756056462968,
                0.1465251111098734,
                0.1805621777986394,
                0.208702484369235,
                0.2307185222979075,
                0.2470522710144055,
                0.2584346164588731,
                0.2656689431620619,
                0.2695209620225222,
                0.2706705664732254,
                0.2696973150384945,
                0.2670827700244007,
                0.263220325068161,
                0.2584275430331589,
                0.2529585154590446,
                0.2470151183771646,
                0.240756758828282,
                0.2343085672139793,
                0.2277681565312566,
                0.2212111327366924,
                0.2146955518792547,
                0.2082655068238174,
                0.201954003936178,
                0.1957852654167679,
                0.1897765695657527,
                0.1839397205857212,
                0.1782822219733523,
                0.1728082130043911,
                0.1675192159470693,
                0.1624147320522222,
                0.1574927166758241,
                0.1527499577431349,
                0.1481823768649407,
                0.1437852685175109,
                0.139553489596412,
                0.1354816091868826,
                0.1315640264293695,
                0.1277950627940034,
                0.1241690338296807,
                0.1206803044566611,
                0.1173233310744522,
                0.1140926931183538,
                0.1109831161860056,
                0.1079894884440163,
                0.1051068716939468,
                0.102330508210468,
                0.09965582424961866,
                0.09707843095155928,
                0.09459412322194045,
                0.0921988770624913,
                0.08988884572949914,
                0.08766035502434205,
                0.08550989795981578,
                0.0834341289969823,
                0.08142985800750233,
                0.07949404408415367,
                0.07762378929607594,
                0.07581633246407918,
            ];
            ddouble[] expected_dist_alpha2beta2 = [
                0,
                2.07489688356359e-10,
                2.304720378250427e-4,
                0.01414430789355111,
                0.08587843274304303,
                0.2177836349094539,
                0.3662089476802515,
                0.4940599115259283,
                0.5861004444394937,
                0.6419988543951622,
                0.6678479499815522,
                0.6711811557757309,
                0.6588060560384146,
                0.6361467482064569,
                0.6072432986561415,
                0.5749780523147141,
                0.5413411329464508,
                0.5076655341901073,
                0.4748138133767124,
                0.4433184422200605,
                0.4134840688530543,
                0.3854605949852108,
                0.3592947176395122,
                0.3349659253263054,
                0.3124114229519724,
                0.2915432403600085,
                0.2722598556759291,
                0.25445398741245,
                0.2380177220843628,
                0.2228457974468171,
                0.2088376164445525,
                0.1958983943904544,
                0.1839397205857212,
                0.1728797303984022,
                0.1626430240041328,
                0.1531604260087491,
                0.1443686507130864,
                0.136209917125037,
                0.1286315433626399,
                0.1215855399917462,
                0.1150282148140087,
                0.1089197967581752,
                0.1032240831900058,
                0.09790811269162385,
                0.09294186385018431,
                0.08829797961221739,
                0.08395151614376421,
                0.07987971477409508,
                0.07606179541223584,
                0.0724787697541261,
                0.06911327260417047,
                0.06594940969031957,
                0.06297262043721108,
                0.0601695542639207,
                0.0575279590824055,
                0.05503658078367445,
                0.05268507260713788,
                0.05046391339199952,
                0.04836433380653354,
                0.04637824974091703,
                0.04449820213172389,
                0.0427173025613127,
                0.04102918404343415,
                0.03942795646784809,
                0.03790816623203959,
            ];
            ddouble[] expected_dist_alpha3beta4 = [
                0,
                3.363435216735857e-22,
                1.659917506850872e-9,
                1.406682040693868e-5,
                9.218881513001707e-4,
                0.009263615553588731,
                0.03771815438280296,
                0.09344269228651089,
                0.1717568654860861,
                0.2608216004110998,
                0.3484538158551262,
                0.4258482922506402,
                0.4882785969070019,
                0.5343076078739588,
                0.564639898886775,
                0.5811323917563668,
                0.5861004444394937,
                0.5818931339117835,
                0.5706656483512553,
                0.5542795162319066,
                0.5342783599852418,
                0.5119040843459497,
                0.4881317496550771,
                0.4637105353719962,
                0.439204037358943,
                0.4150266650371762,
                0.3914749219732043,
                0.3687534399739019,
                0.3469961706606522,
                0.3262833576468303,
                0.3066549612345141,
                0.2881211708491598,
                0.2706705664732254,
                0.2542764053740453,
                0.2389014278541681,
                0.2245015016921668,
                0.211028361500761,
                0.1984316464474512,
                0.186660396724236,
                0.1756641345067466,
                0.1653936275412217,
                0.155801411668287,
                0.1468421314229374,
                0.1384727444024628,
                0.1306526245961862,
                0.1233435916975711,
                0.1165098870700193,
                0.1101181121152441,
                0.1041371409839908,
                0.09853801662685135,
                0.09329383691520274,
                0.08837963581830237,
                0.08377226328490128,
                0.0794502664549671,
                0.07539377404813331,
                0.07158438518502122,
                0.06800506345267505,
                0.06464003669265356,
                0.06147470274394953,
                0.05849554119278487,
                0.05569003105188065,
                0.05304657420097593,
                0.05055442435882695,
                0.04820362131745275,
                0.04598493014643029,
            ];

            foreach ((InverseGammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1beta1, expected_dist_alpha1beta1), (dist_alpha2beta1, expected_dist_alpha2beta1),
                (dist_alpha1beta2, expected_dist_alpha1beta2), (dist_alpha2beta2, expected_dist_alpha2beta2),
                (dist_alpha3beta4, expected_dist_alpha3beta4),
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
            ddouble[] expected_dist_alpha1beta1 = [
                0,
                1.125351747192591e-7,
                3.354626279025119e-4,
                0.004827949993831442,
                0.01831563888873418,
                0.04076220397836622,
                0.06948345122280156,
                0.1017013923042268,
                0.1353352832366128,
                0.1690133154060661,
                0.2018965179946556,
                0.2335064790909132,
                0.2635971381157268,
                0.292067823691414,
                0.3189065573239707,
                0.3441537868654124,
                0.3678794411714422,
                0.3901685434239768,
                0.4111122905071876,
                0.4308026151974352,
                0.4493289641172217,
                0.4667764816516813,
                0.4832250811898253,
                0.4987490707622945,
                0.5134171190325922,
                0.5272924240430487,
                0.540432996486534,
                0.5528920012788026,
                0.5647181220077593,
                0.5759559263708725,
                0.5866462195100317,
                0.59682637785056,
                0.6065306597126334,
                0.6157904912859286,
                0.6246347280002745,
                0.6330898921891812,
                0.6411803884299545,
                0.6489286981868463,
                0.6563555554708403,
                0.6634801052198983,
                0.6703200460356393,
                0.6768917588119943,
                0.683210422674948,
                0.6892901195306848,
                0.6951439283988786,
                0.7007840105925329,
                0.7062216866978064,
                0.7114675062078237,
                0.7165313105737892,
                0.7214222903547562,
                0.726149037073691,
                0.7307195903214863,
                0.7351414805916845,
                0.7394217682762254,
                0.7435670792059064,
                0.7475836370778213,
                0.751477293075286,
                0.7552535529531295,
                0.7589176018322887,
                0.7624743269219467,
                0.7659283383646487,
                0.769283988379571,
                0.7725453888611077,
                0.7757164275739283,
                0.7788007830714049,
            ];
            ddouble[] expected_dist_alpha2beta1 = [
                0,
                1.913097970227405e-6,
                0.003019163651122607,
                0.03057701662759913,
                0.0915781944436709,
                0.1712012567091381,
                0.2547726544836055,
                0.3341617175710309,
                0.4060058497098382,
                0.4694814316835169,
                0.5249309467861041,
                0.5731522668595147,
                0.6150599889366959,
                0.6515359143885393,
                0.6833711942656508,
                0.7112511595218522,
                0.7357588823428847,
                0.757385996058308,
                0.7765454376246874,
                0.7935837648373807,
                0.8087921354109989,
                0.8224157057672481,
                0.8346615038733349,
                0.8457049460751951,
                0.8556951983876534,
                0.8647595754305997,
                0.873007148170555,
                0.8805317057403151,
                0.8874141917264787,
                0.8937247133341124,
                0.8995242032487154,
                0.9048657986766555,
                0.9097959895689501,
                0.9143555779700152,
                0.9185804823533448,
                0.9225024143328071,
                0.9261494499543789,
                0.929546513618996,
                0.9327157893532992,
                0.9356770714639591,
                0.938448064449895,
                0.9410446402996017,
                0.9434810598844523,
                0.9457701640072186,
                0.9479235387257438,
                0.9499516588032112,
                0.951864012505739,
                0.953669210448785,
                0.9553750807650523,
                0.9569887525114112,
                0.958516728937272,
                0.9599649519909722,
                0.9613388592352797,
                0.962643434170935,
                0.9638832508224712,
                0.965062513318642,
                0.9661850910967962,
                0.9672545502733062,
                0.9682741816480926,
                0.9692470257482374,
                0.9701758952618883,
                0.9710633951676552,
                0.9719119408252644,
                0.9727237742593704,
                0.9735009788392561,
            ];
            ddouble[] expected_dist_alpha1beta2 = [
                0,
                1.266416554909418e-14,
                1.125351747192591e-7,
                2.330910114293702e-5,
                3.354626279025119e-4,
                0.001661557273173934,
                0.004827949993831442,
                0.01034317319661825,
                0.01831563888873418,
                0.02856550078455037,
                0.04076220397836622,
                0.05452527577743517,
                0.06948345122280156,
                0.08530361363583897,
                0.1017013923042268,
                0.1184418290138037,
                0.1353352832366128,
                0.1522314922775879,
                0.1690133154060661,
                0.1855908932609494,
                0.2018965179946556,
                0.2178802838231223,
                0.2335064790909132,
                0.2487506355862523,
                0.2635971381157268,
                0.2780373004531941,
                0.292067823691414,
                0.3056895650780793,
                0.3189065573239707,
                0.3317252291217299,
                0.3441537868654124,
                0.3562017252982195,
                0.3678794411714422,
                0.3791979291581653,
                0.3901685434239768,
                0.4008028115921093,
                0.4111122905071876,
                0.4211084553304749,
                0.4308026151974352,
                0.4402058500226073,
                0.4493289641172217,
                0.4581824531475951,
                0.4667764816516813,
                0.4751208688826257,
                0.4832250811898253,
                0.4910982295021552,
                0.4987490707622945,
                0.5061860123895796,
                0.5134171190325922,
                0.5204501210207022,
                0.5272924240430487,
                0.5339511196796007,
                0.540432996486534,
                0.5467445514007401,
                0.5528920012788026,
                0.5588812944265036,
                0.5647181220077593,
                0.5704079292483255,
                0.5759559263708725,
                0.5813670992150757,
                0.5866462195100317,
                0.5917978547771798,
                0.59682637785056,
                0.6017359760080574,
                0.6065306597126334,
            ];
            ddouble[] expected_dist_alpha2beta2 = [
                0,
                4.179174631201078e-13,
                1.913097970227405e-6,
                2.719395133342652e-4,
                0.003019163651122607,
                0.01229552382148711,
                0.03057701662759913,
                0.0576262506668731,
                0.0915781944436709,
                0.130131725796285,
                0.1712012567091381,
                0.2131442598572464,
                0.2547726544836055,
                0.2952817395086733,
                0.3341617175710309,
                0.3711177309099181,
                0.4060058497098382,
                0.4387848895059879,
                0.4694814316835169,
                0.4981650292793908,
                0.5249309467861041,
                0.5498883353631184,
                0.5731522668595147,
                0.5948384764019078,
                0.6150599889366959,
                0.6339250450332825,
                0.6515359143885393,
                0.6679883088743217,
                0.6833711942656508,
                0.6977668612560525,
                0.7112511595218522,
                0.7238938288318655,
                0.7357588823428847,
                0.7469050119782044,
                0.757385996058308,
                0.7672510964763234,
                0.7765454376246874,
                0.7853103626433182,
                0.7935837648373807,
                0.8014003936309002,
                0.8087921354109989,
                0.8157882702384007,
                0.8224157057672481,
                0.8286991899115563,
                0.8346615038733349,
                0.8403236371481321,
                0.8457049460751951,
                0.8508232974207829,
                0.8556951983876534,
                0.8603359143403443,
                0.8647595754305997,
                0.8689792732040562,
                0.873007148170555,
                0.8768544692276019,
                0.8805317057403151,
                0.884048593001924,
                0.8874141917264787,
                0.8906369421596663,
                0.8937247133341124,
                0.8966848479418964,
                0.8995242032487154,
                0.9022491884307825,
                0.9048657986766555,
                0.9073796463613566,
                0.9097959895689501,

            ];
            ddouble[] expected_dist_alpha3beta4 = [
                0,
                3.388852411729272e-25,
                6.901970224256324e-12,
                1.357681807789067e-7,
                1.631760033429257e-5,
                2.642611505954007e-4,
                0.001597968378354682,
                0.005556474180979418,
                0.01375396774400298,
                0.02724996412185383,
                0.04632421677608928,
                0.07059045762060492,
                0.09924119431764629,
                0.1312821190323359,
                0.1657018563131699,
                0.2015734815947414,
                0.2381033055535439,
                0.2746450287219193,
                0.3106939035949244,
                0.3458709869046112,
                0.3799037410783731,
                0.4126065276558063,
                0.4438627821551542,
                0.4736096065685852,
                0.5018249254980112,
                0.5285170494627516,
                0.5537163559675466,
                0.5774687605491114,
                0.599830660733093,
                0.620865069122682,
                0.6406386929324406,
                0.6592197584218555,
                0.6766764161830634,
                0.6930755956741848,
                0.7084822045444825,
                0.7229585905387028,
                0.7365642017079178,
                0.7493553949652927,
                0.7613853543475516,
                0.772704089244165,
                0.7833584898192629,
                0.7933924222757448,
                0.8028468508221631,
                0.8117599764617073,
                0.8201673852366793,
                0.8281022004953644,
                0.8355952352301896,
                0.8426751416676529,
                0.8493685561506752,
                0.8557002380043588,
                0.8616932015645393,
                0.8673688409122559,
                0.8727470471252315,
                0.8778463180520769,
                0.8826838607535763,
                0.8872756868515179,
                0.8916367010894684,
                0.8957807834496924,
                0.8997208651922305,
                0.9034689991907323,
                0.9070364249386202,
                0.9104336285913062,
                0.9136703983976181,
                0.9167558758578966,
                0.9196986029286058,
            ];

            foreach ((InverseGammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1beta1, expected_dist_alpha1beta1), (dist_alpha2beta1, expected_dist_alpha2beta1),
                (dist_alpha1beta2, expected_dist_alpha1beta2), (dist_alpha2beta2, expected_dist_alpha2beta2),
                (dist_alpha3beta4, expected_dist_alpha3beta4),
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