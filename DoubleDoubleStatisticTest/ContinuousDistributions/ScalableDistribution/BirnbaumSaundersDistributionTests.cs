﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class BirnbaumSaundersDistributionTests {
        readonly BirnbaumSaundersDistribution dist_alpha1theta1 = new(alpha: 1, theta: 1);
        readonly BirnbaumSaundersDistribution dist_alpha2theta1 = new(alpha: 2, theta: 1);
        readonly BirnbaumSaundersDistribution dist_alpha1theta2 = new(alpha: 1, theta: 2);
        readonly BirnbaumSaundersDistribution dist_alpha2theta2 = new(alpha: 2, theta: 2);
        readonly BirnbaumSaundersDistribution dist_alpha3theta4 = new(alpha: 3, theta: 4);
        readonly BirnbaumSaundersDistribution dist_alpha4theta5 = new(alpha: 4, theta: 5);
        readonly BirnbaumSaundersDistribution dist_alpha5theta7 = new(alpha: 5, theta: 7);
        readonly BirnbaumSaundersDistribution dist_alpha6theta3 = new(alpha: 6, theta: 3);

        BirnbaumSaundersDistribution[] Dists => [
            dist_alpha1theta1,
            dist_alpha2theta1,
            dist_alpha1theta2,
            dist_alpha2theta2,
            dist_alpha3theta4,
            dist_alpha4theta5,
            dist_alpha5theta7,
            dist_alpha6theta3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BirnbaumSaundersDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Theta={dist.Theta}");
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
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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
            foreach (BirnbaumSaundersDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode - 1e-15), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) > dist.PDF(dist.Mode + 1e-15), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (BirnbaumSaundersDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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
            foreach (BirnbaumSaundersDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (BirnbaumSaundersDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (BirnbaumSaundersDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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

            foreach (BirnbaumSaundersDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
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

            foreach (BirnbaumSaundersDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                (BirnbaumSaundersDistribution? dist_fit, ddouble error) = BirnbaumSaundersDistribution.Fit(xs, (0.15, 0.85));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-1);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (BirnbaumSaundersDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha1theta1 = [
                0,
                3.542921657540213e-12,
                1.121413859877021e-5,
                0.001273515423683824,
                0.01198825320049314,
                0.0427769002388211,
                0.09517895074891927,
                0.1628467686328756,
                0.2374874816856767,
                0.312269196933268,
                0.382670214921155,
                0.4461824290554608,
                0.501732849554182,
                0.549160929723652,
                0.588834202050792,
                0.6213924360040674,
                0.6475879783294587,
                0.6681914967990588,
                0.6839398256383017,
                0.6955098819185134,
                0.7035081739121458,
                0.7084692686729321,
                0.7108591213205049,
                0.7110807912743663,
                0.709481089182417,
                0.7063573276913794,
                0.7019637323879333,
                0.6965172984906686,
                0.6902030128417632,
                0.6831784363173025,
                0.675577681972508,
                0.6675148433489467,
                0.6590869342015836,
                0.6503764008400458,
                0.6414532645783437,
                0.6323769463900142,
                0.6231978198980616,
                0.6139585329211287,
                0.604695132278312,
                0.5954380215753158,
                0.5862127773009127,
                0.5770408447438719,
                0.5679401319557152,
                0.5589255171796713,
                0.5500092827830521,
                0.5412014867126865,
                0.532510280788929,
                0.523942183716592,
                0.5155023154801236,
                0.507194598770214,
                0.4990219322296308,
                0.4909863395818636,
                0.4830890980954587,
                0.4753308493215654,
                0.4677116946068868,
                0.460231277516125,
                0.4528888549864022,
                0.4456833587720241,
                0.4386134485138333,
                0.431677557576954,
                0.4248739326387049,
                0.4182006678704205,
                0.4116557344391664,
                0.4052370059547541,
                0.3989422804014327,
            ];
            ddouble[] expected_dist_alpha2theta1 = [
                0,
                0.02229579445531405,
                0.4361523668375034,
                0.9125218871251315,
                1.169366689053762,
                1.264140770135086,
                1.271265170995535,
                1.235522902070302,
                1.180676351902834,
                1.118920306852991,
                1.056397833988747,
                0.9961043572022626,
                0.9393922515621947,
                0.8867561877280106,
                0.8382483669474211,
                0.793700565686822,
                0.7528435803870112,
                0.7153712971585767,
                0.6809744876052147,
                0.6493579196827782,
                0.6202482585127667,
                0.5933969217060171,
                0.56858022995123,
                0.5455981736684339,
                0.524272539133992,
                0.5044448076789149,
                0.4859740522014361,
                0.4687349464605572,
                0.4526159404951146,
                0.437517620556978,
                0.4233512528465828,
                0.4100375004718443,
                0.3975052985160428,
                0.3856908705900682,
                0.3745368703959067,
                0.3639916328389826,
                0.3540085206295525,
                0.3445453538410398,
                0.3355639113989233,
                0.3270294948816232,
                0.3189105462904595,
                0.3111783125787912,
                0.3038065507236617,
                0.2967712679862759,
                0.2900504927530162,
                0.2836240719899601,
                0.2774734918941429,
                0.2715817187963037,
                0.2659330577735701,
                0.2605130267761513,
                0.2553082443680902,
                0.2503063294357505,
                0.2454958114352487,
                0.2408660499368198,
                0.2364071623846813,
                0.2321099591291961,
                0.2279658849073081,
                0.2239669660501311,
                0.2201057627855544,
                0.2163753260808093,
                0.2127691585368117,
                0.2092811789042025,
                0.2059056898415889,
                0.2026373485805795,
                0.1994711402007164,
            ];
            ddouble[] expected_dist_alpha1theta2 = [
                0,
                6.321147826859193e-26,
                1.771460828770107e-12,
                4.152413736644897e-8,
                5.607069299385105e-6,
                9.878579362259047e-5,
                6.367577118419119e-4,
                0.002327386999002704,
                0.005994126600246572,
                0.01226088624394835,
                0.02138845011941055,
                0.03328204200415862,
                0.04758947537445964,
                0.06381671826062005,
                0.0814233843164378,
                0.09988768737481453,
                0.1187437408428383,
                0.1375983742795632,
                0.156134598466634,
                0.1741073097838625,
                0.1913351074605775,
                0.2076906822628985,
                0.2230912145277304,
                0.2374895408280922,
                0.250866424777091,
                0.2632240184265117,
                0.274580464861826,
                0.2849655267440575,
                0.294417101025396,
                0.3029784785993434,
                0.3106962180020337,
                0.3176185177964216,
                0.3237939891647293,
                0.3292707464935488,
                0.3340957483995294,
                0.3383143343561201,
                0.3419699128191508,
                0.3451037656544055,
                0.3477549409592567,
                0.3499602122773884,
                0.3517540869560729,
                0.353168850192606,
                0.3542346343364661,
                0.3549795054053766,
                0.3554295606602524,
                0.3556090325672615,
                0.3555403956371831,
                0.355244473539175,
                0.3547405445912085,
                0.3540464442752581,
                0.3531786638456897,
                0.3521524444210556,
                0.3509818661939667,
                0.3496799325777321,
                0.3482586492453343,
                0.3467290981165481,
                0.3451015064208816,
                0.3433853110140154,
                0.3415892181586512,
                0.3397212590011596,
                0.337788840986254,
                0.3357987954555806,
                0.3337574216744734,
                0.3316705275256876,
                0.3295434671007918,
            ];
            ddouble[] expected_dist_alpha2theta2 = [
                0,
                1.050636702689661e-5,
                0.01114789722765703,
                0.08791814449679355,
                0.2180761834187517,
                0.3495685861270489,
                0.4562609435625657,
                0.533388295172528,
                0.5846833445268812,
                0.6158279087580286,
                0.6320703850675428,
                0.6376144038458009,
                0.6356325854977677,
                0.6284532532447158,
                0.6177614510351509,
                0.604767786568179,
                0.5903381759514168,
                0.5750897003880914,
                0.5594601534264957,
                0.5437581133682673,
                0.5281989169943736,
                0.5129305269816801,
                0.4980521786011313,
                0.4836278652665228,
                0.4696961257810974,
                0.4562771717516527,
                0.4433780938640053,
                0.4309966743193601,
                0.4191241834737106,
                0.4077474330154584,
                0.396850282843411,
                0.3864147450941715,
                0.3764217901935056,
                0.366851931959602,
                0.3576856485792884,
                0.3489036815405221,
                0.3404872438026074,
                0.3324181605311908,
                0.3246789598413891,
                0.3172529266221656,
                0.3101241292563833,
                0.3032774266124698,
                0.2966984608530085,
                0.2903736402274,
                0.284290114975615,
                0.2784357486834876,
                0.2727990868342169,
                0.267369323848998,
                0.262136269566996,
                0.2570903158549213,
                0.2522224038394574,
                0.2475239921067487,
                0.2429870261007181,
                0.2386039088675126,
                0.2343674732302786,
                0.230270955431711,
                0.2263079702475573,
                0.2224724875495157,
                0.218758810278489,
                0.2151615537771895,
                0.2116756264232914,
                0.2082962114996679,
                0.2050187502359221,
                0.2018389259548369,
                0.1987526492580214,
            ];
            ddouble[] expected_dist_alpha3theta4 = [
                0,
                5.084871625975336e-5,
                0.02211279552420058,
                0.1292800000384727,
                0.2756906012646811,
                0.4031451914296517,
                0.4944745061820324,
                0.5525182300465604,
                0.5850724294058501,
                0.5995357886490398,
                0.6016579274382067,
                0.5955897268618552,
                0.5842333078294426,
                0.5695871599478478,
                0.5530163145790561,
                0.5354475351571393,
                0.5175057419355621,
                0.4996080967475628,
                0.4820286304373316,
                0.4649426674189299,
                0.4484574545179039,
                0.4326333616455811,
                0.4174986164584204,
                0.4030595828740504,
                0.3893079518784005,
                0.3762257810544309,
                0.3637890274054536,
                0.3519700198970125,
                0.3407391828165107,
                0.3300662280392729,
                0.3199209699437968,
                0.3102738719147563,
                0.3010964019789289,
                0.2923612529907436,
                0.2840424670985903,
                0.276115493046598,
                0.2685571968652959,
                0.2613458407521088,
                0.2544610407916796,
                0.2478837111614217,
                0.2415960002875581,
                0.2355812228322557,
                0.2298237902396099,
                0.2243091417295445,
                0.2190236770192627,
                0.21395469161014,
                0.209090315159099,
                0.2044194532248091,
                0.1999317325165346,
                0.195617449659491,
                0.1914675234121806,
                0.1874734502187877,
                0.1836272629462514,
                0.1799214926358571,
                0.1763491330892362,
                0.1729036081056738,
                0.1695787411894865,
                0.1663687275513919,
                0.1632681082351129,
                0.1602717462091111,
                0.1573748042727226,
                0.1545727246356512,
                0.1518612100394493,
                0.1492362062990831,
                0.1466938861517915,
            ];

            foreach ((BirnbaumSaundersDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1theta1, expected_dist_alpha1theta1), (dist_alpha2theta1, expected_dist_alpha2theta1),
                (dist_alpha1theta2, expected_dist_alpha1theta2), (dist_alpha2theta2, expected_dist_alpha2theta2),
                (dist_alpha3theta4, expected_dist_alpha3theta4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 64, i++) {
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
            ddouble[] expected_dist_alpha1theta1 = [
                0,
                1.70371429163288e-15,
                2.125697204124573e-8,
                5.355569218449527e-6,
                8.841728520080404e-5,
                4.865335295668477e-4,
                0.001539193558335859,
                0.003540646687975497,
                0.006664164390408786,
                0.0109624426491709,
                0.01639940693277282,
                0.02288476040686283,
                0.03030098485600309,
                0.0385215191244514,
                0.04742182899680732,
                0.05688551222772088,
                0.06680720126885814,
                0.07709351782149178,
                0.08766290881043286,
                0.09844488793755406,
                0.109379002266696,
                0.1204137115947362,
                0.1315052859634286,
                0.1426167762841981,
                0.1537170829636977,
                0.1647801299803106,
                0.1757841421865406,
                0.1867110186267617,
                0.1975457924520751,
                0.2082761674058426,
                0.2188921211263572,
                0.2293855662349518,
                0.2397500610934768,
                0.2499805630782042,
                0.2600732181476975,
                0.2700251813396976,
                0.2798344635997059,
                0.2894998010182215,
                0.2990205431382895,
                0.3083965574975897,
                0.3176281479986242,
                0.3267159850657604,
                0.3356610458576623,
                0.34446456306591,
                0.3531279810524264,
                0.3616529182658222,
                0.3700411350352769,
                0.3782945059745655,
                0.3864149963422239,
                0.3944046417998098,
                0.4022655310915618,
                0.4099997912377387,
                0.4176095748925146,
                0.4250970495671153,
                0.4324643884612825,
                0.4397137626822947,
                0.4468473346616164,
                0.4538672526056051,
                0.4607756458392525,
                0.4675746209212606,
                0.4742662584253222,
                0.4808526102967101,
                0.4873356977055178,
                0.4937175093284321,
                0.5,
            ];
            ddouble[] expected_dist_alpha2theta1 = [
                0,
                4.11674659715994e-5,
                0.003071596802393406,
                0.01386277756691388,
                0.03039636176526126,
                0.04956356116266086,
                0.06944996078719838,
                0.08907091236582491,
                0.1079624694700702,
                0.1259320600066064,
                0.1429252565771281,
                0.1589563730553997,
                0.1740723124765312,
                0.1883336833773987,
                0.2018049808831361,
                0.2145495863977466,
                0.2266273523768683,
                0.2380935771554052,
                0.2489987208946829,
                0.2593885064177917,
                0.2693042083406119,
                0.2787830222520503,
                0.2878584551319938,
                0.2965607060429858,
                0.3049170218367299,
                0.3129520214178876,
                0.3206879869549681,
                0.328145123003716,
                0.3353417857817041,
                0.3422946853666827,
                0.349019063713827,
                0.3555288512898281,
                0.3618368049158816,
                0.3679546291607841,
                0.3738930833652396,
                0.3796620761282188,
                0.3852707488552772,
                0.3907275497609426,
                0.3960402995333206,
                0.4012162497079351,
                0.4062621346576843,
                0.4111842179844662,
                0.4159883339933134,
                0.4206799248396442,
                0.4252640738625542,
                0.4297455355502034,
                0.4341287625257796,
                0.4384179298929083,
                0.4426169572366009,
                0.4467295285388821,
                0.4507591102363087,
                0.4547089676189277,
                0.4585821797462482,
                0.4623816530349687,
                0.4661101336550808,
                0.4697702188551832,
                0.473364367324052,
                0.4768949086834643,
                0.4803640521967158,
                0.4837738947680115,
                0.487126428299772,
                0.4904235464677346,
                0.4936670509674097,
                0.4968586572798813,
                0.5,
            ];
            ddouble[] expected_dist_alpha1theta2 = [
                0,
                1.531375814024154e-29,
                1.70371429163288e-15,
                8.919413949475407e-11,
                2.125697204124573e-8,
                5.810225289367713e-7,
                5.355569218449527e-6,
                2.646175940940346e-5,
                8.841728520080404e-5,
                2.273868403253242e-4,
                4.865335295668477e-4,
                9.102272215671691e-4,
                0.001539193558335859,
                0.002407404712469713,
                0.003540646687975497,
                0.004956342205234481,
                0.006664164390408786,
                0.00866706332993168,
                0.0109624426491709,
                0.0135433234832981,
                0.01639940693277282,
                0.01951799505727069,
                0.02288476040686283,
                0.02648437086510061,
                0.03030098485600309,
                0.03431863508079232,
                0.0385215191244514,
                0.04289421391633379,
                0.04742182899680732,
                0.05209011132514627,
                0.05688551222772088,
                0.06179522515797345,
                0.06680720126885814,
                0.07191014838841218,
                0.07709351782149178,
                0.08234748244712133,
                0.08766290881043286,
                0.09303132529090609,
                0.09844488793755406,
                0.1038963451734504,
                0.109379002266696,
                0.1148866862261746,
                0.1204137115947362,
                0.1259548474688501,
                0.1315052859634286,
                0.1370606122564169,
                0.1426167762841981,
                0.1481700661113899,
                0.1537170829636977,
                0.159254717887269,
                0.1647801299803106,
                0.1702907261307921,
                0.1757841421865406,
                0.1812582254798158,
                0.1867110186267617,
                0.1921407445222452,
                0.1975457924520751,
                0.2029247052469854,
                0.2082761674058426,
                0.2135989941189889,
                0.2188921211263572,
                0.2241545953488037,
                0.2293855662349518,
                0.2345842777696064,
                0.2397500610934768,
            ];
            ddouble[] expected_dist_alpha2theta2 = [
                0,
                9.961987590979348e-9,
                4.11674659715994e-5,
                7.127442037770395e-4,
                0.003071596802393406,
                0.00752858751409391,
                0.01386277756691388,
                0.02163143989123384,
                0.03039636176526126,
                0.0397979515476222,
                0.04956356116266086,
                0.05949461170716999,
                0.06944996078719838,
                0.07933117053047795,
                0.08907091236582491,
                0.09862428614685768,
                0.1079624694700702,
                0.1170681199525334,
                0.1259320600066064,
                0.134550888161284,
                0.1429252565771281,
                0.15105862728863,
                0.1589563730553997,
                0.1666251270009396,
                0.1740723124765312,
                0.181305803938775,
                0.1883336833773987,
                0.1951640666246246,
                0.2018049808831361,
                0.2082642798436197,
                0.2145495863977466,
                0.220668255590194,
                0.2266273523768683,
                0.232433640165911,
                0.2380935771554052,
                0.2436133182483978,
                0.2489987208946829,
                0.2542553536322829,
                0.2593885064177917,
                0.2644032020714771,
                0.2693042083406119,
                0.2740960502178898,
                0.2787830222520503,
                0.2833692006632084,
                0.2878584551319938,
                0.2922544601740296,
                0.2965607060429858,
                0.300780509129009,
                0.3049170218367299,
                0.3089732419397754,
                0.3129520214178876,
                0.3168560747892352,
                0.3206879869549681,
                0.3244502205759908,
                0.328145123003716,
                0.3317749327874925,
                0.3353417857817041,
                0.3388477208753807,
                0.3422946853666827,
                0.3456845400039106,
                0.349019063713827,
                0.35229995803712,
                0.3555288512898281,
                0.358707302468506,
                0.3618368049158816,
            ];
            ddouble[] expected_dist_alpha3theta4 = [
                0,
                5.406571334852242e-8,
                9.136312729735371e-5,
                0.001170706226257928,
                0.004332448363012565,
                0.009679252881425104,
                0.016739427803989,
                0.0249575570834591,
                0.03387300692034461,
                0.04314715540929641,
                0.05254450502164954,
                0.06190652372111661,
                0.07112929634065022,
                0.08014678897727456,
                0.0889189455012411,
                0.09742340277568703,
                0.1056497736668554,
                0.113595711589412,
                0.1212642007296284,
                0.1286616882487395,
                0.1357967960009852,
                0.142679432526038,
                0.1493201826382396,
                0.1557298902921164,
                0.1619193764668066,
                0.1678992515976764,
                0.1736797942872152,
                0.1792708764518079,
                0.1846819209126646,
                0.189921881528574,
                0.1949992388425585,
                0.1999220062458291,
                0.2046977431049594,
                0.2093335723269906,
                0.213836200573566,
                0.2182119398637879,
                0.2224667296856972,
                0.2266061590100274,
                0.2306354877969715,
                0.2345596677283751,
                0.2383833619992901,
                0.2421109640751807,
                0.2457466153722282,
                0.2492942218538995,
                0.2527574695614466,
                0.2561398391123633,
                0.2594446192112878,
                0.2626749192240664,
                0.2658336808688819,
                0.2689236890794098,
                0.2719475820945233,
                0.2749078608276295,
                0.2778068975666015,
                0.2806469440527589,
                0.2834301389845986,
                0.2861585149891342,
                0.2888340051008468,
                0.2914584487854626,
                0.2940335975430683,
                0.2965611201225106,
                0.2990426073765928,
                0.3014795767853069,
                0.3038734766722049,
                0.3062256901370424,
                0.308537538725987,
            ];

            foreach ((BirnbaumSaundersDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1theta1, expected_dist_alpha1theta1), (dist_alpha2theta1, expected_dist_alpha2theta1),
                (dist_alpha1theta2, expected_dist_alpha1theta2), (dist_alpha2theta2, expected_dist_alpha2theta2),
                (dist_alpha3theta4, expected_dist_alpha3theta4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 64, i++) {
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