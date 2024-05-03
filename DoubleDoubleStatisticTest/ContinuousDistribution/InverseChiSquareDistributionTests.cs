using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.InternalUtils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class InverseChiSquareDistributionTests {
        readonly InverseChiSquareDistribution dist_nu1 = new(nu: 1);
        readonly InverseChiSquareDistribution dist_nu2 = new(nu: 2);
        readonly InverseChiSquareDistribution dist_nu3 = new(nu: 3);
        readonly InverseChiSquareDistribution dist_nu4 = new(nu: 4);
        readonly InverseChiSquareDistribution dist_nu5 = new(nu: 5);
        readonly InverseChiSquareDistribution dist_nu8 = new(nu: 8);
        readonly InverseChiSquareDistribution dist_nu16 = new(nu: 16);
        readonly InverseChiSquareDistribution dist_nu32 = new(nu: 32);
        readonly InverseChiSquareDistribution dist_nu64 = new(nu: 64);
        readonly InverseChiSquareDistribution dist_nu128 = new(nu: 128);


        InverseChiSquareDistribution[] Dists => [
            dist_nu1,
            dist_nu2,
            dist_nu3,
            dist_nu4,
            dist_nu5,
            dist_nu8,
            dist_nu16,
            dist_nu32,
            dist_nu64,
            dist_nu128,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (InverseChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.Nu}");
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
            foreach (InverseChiSquareDistribution dist in Dists) {
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
            foreach (InverseChiSquareDistribution dist in Dists) {
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
            foreach (InverseChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (InverseChiSquareDistribution dist in Dists) {
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
            foreach (InverseChiSquareDistribution dist in Dists) {
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
            foreach (InverseChiSquareDistribution dist in Dists) {
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
            foreach (InverseChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (InverseChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (InverseChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);

                Console.WriteLine($"cdf(median)={cdf_median}");

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (InverseChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);
                ddouble ccdf_median = dist.CDF(dist.Median, Interval.Upper);

                Console.WriteLine($"ccdf(median)={ccdf_median}");

                Assert.IsTrue(ddouble.Abs(cdf_median + ccdf_median - 1) < 1e-28);

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (InverseChiSquareDistribution dist in Dists) {
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
            foreach (InverseChiSquareDistribution dist in Dists) {
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
        public void IrregularValueTest() {
            foreach (InverseChiSquareDistribution dist in Dists) {
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
            ddouble[] expected_nu1 = [
                0,
                2.586762794770889e-12,
                8.126870189994599e-6,
                9.162704297893385e-4,
                0.008565134448952662,
                0.0303557058707318,
                0.06709897939421094,
                0.1140740245479289,
                0.1653358828327365,
                0.2161016608957144,
                0.2632920128692105,
                0.3052733340606918,
                0.3414204569963355,
                0.3717341220586473,
                0.3965654475692404,
                0.4164359781350884,
                0.4319277321055044,
                0.4436200325765856,
                0.4520559503785433,
                0.4577267433189473,
                0.4610668118035667,
                0.4624545077172937,
                0.46221596392128,
                0.460630267808355,
                0.4579350180044978,
                0.45433173884705,
                0.44999088770528,
                0.4450563425508125,
                0.4396493435467702,
                0.4338719089958032,
                0.4278097689282868,
                0.4215348688204982,
                0.4151074974205948,
                0.4085780900403695,
                0.4019887540808652,
                0.3953745582364216,
                0.388764621454207,
                0.3821830326762757,
                0.375649627826573,
                0.3691806464824544,
                0.3627892871839332,
                0.3564861773470029,
                0.350279771208785,
                0.3441766870866545,
                0.338181993427016,
                0.3322994516023839,
                0.3265317221435988,
                0.3208805400287553,
                0.3153468637585577,
                0.3099310022011155,
                0.3046327225638241,
                0.2994513423259197,
                0.294385807525805,
                0.2894347594283207,
                0.2845965912871562,
                0.2798694966568057,
                0.2752515104888594,
                0.2707405440622088,
                0.2663344146403772,
                0.2620308706169872,
                0.2578276127984735,
                0.2537223123783056,
                0.2497126260764819,
                0.2457962088496677,
                0.2419707245191433,
            ];
            ddouble[] expected_nu2 = [
                0,
                2.593621104454487e-11,
                5.761800945626067e-5,
                0.005304115460081667,
                0.04293921637152152,
                0.1361147718184087,
                0.2746567107601887,
                0.4323024225851873,
                0.5861004444394937,
                0.7222487111945576,
                0.8348099374769402,
                0.9228740891916299,
                0.9882090840576219,
                1.033738465835492,
                1.062675772648248,
                1.078083848090089,
                1.082682265892902,
                1.078789260153978,
                1.068331080097603,
                1.052881300272644,
                1.033710172132636,
                1.011834061836179,
                0.9880604735086583,
                0.9630270353131279,
                0.9372342688559173,
                0.9110726261250266,
                0.8848445309467696,
                0.8587822075170187,
                0.8330620272952697,
                0.8078160157447118,
                0.7831410616670718,
                0.7591062782630108,
                0.7357588823428847,
                0.7131288878934092,
                0.6912328520175645,
                0.6700768637882772,
                0.6496589282088888,
                0.6299708667032964,
                0.6109998309725396,
                0.5927295074597629,
                0.5751410740700437,
                0.5582139583856481,
                0.5419264367475303,
                0.5262561057174782,
                0.5111802511760137,
                0.4966761353187228,
                0.4827212178266442,
                0.4692933242978086,
                0.4563707724734151,
                0.4439324647440224,
                0.4319579537760654,
                0.4204274867757872,
                0.409322032841872,
                0.3986232969984747,
                0.3883137238062371,
                0.3783764928877618,
                0.3687955082499652,
                0.3595553829179965,
                0.3506414200973682,
                0.3420395918392631,
                0.3337365159879292,
                0.3257194320300093,
                0.3179761763366147,
                0.3104951571843038,
                0.3032653298563167,
            ];
            ddouble[] expected_nu4 = [
                0,
                8.299587534254359e-10,
                9.218881513001707e-4,
                0.05657723157420444,
                0.3435137309721721,
                0.8711345396378156,
                1.464835790721006,
                1.976239646103713,
                2.344401777757975,
                2.567995417580649,
                2.671391799926209,
                2.684724623102924,
                2.635224224153658,
                2.544586992825828,
                2.428973194624566,
                2.299912209258856,
                2.165364531785803,
                2.030662136760429,
                1.899255253506849,
                1.773273768880242,
                1.653936275412217,
                1.541842379940843,
                1.437178870558048,
                1.339863701305221,
                1.24964569180789,
                1.166172961440034,
                1.089039422703717,
                1.0178159496498,
                0.9520708883374509,
                0.8913831897872683,
                0.8353504657782098,
                0.7835935775618177,
                0.7357588823428847,
                0.691518921593609,
                0.6505720960165313,
                0.6126417040349963,
                0.5774746028523454,
                0.5448396685001482,
                0.5145261734505596,
                0.486342159966985,
                0.460112859256035,
                0.435679187032701,
                0.412896332760023,
                0.3916324507664954,
                0.3717674554007372,
                0.3531919184488696,
                0.3358060645750569,
                0.3195188590963803,
                0.3042471816489433,
                0.2899150790165044,
                0.2764530904166819,
                0.2637976387612783,
                0.2518904817488443,
                0.2406782170556828,
                0.230111836329622,
                0.2201463231346978,
                0.2107402904285515,
                0.2018556535679981,
                0.1934573352261342,
                0.1855129989636682,
                0.1779928085268956,
                0.1708692102452508,
                0.1641167361737366,
                0.1577118258713924,
                0.1516326649281584,
            ];

            foreach ((InverseChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
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
            ddouble[] expected_nu1 = [
                0,
                1.244192114854357e-15,
                1.541725790028002e-8,
                3.859616436928451e-6,
                6.334248366623986e-5,
                3.46619351134667e-4,
                0.001090835176125297,
                0.00249690891514155,
                0.004677734981047264,
                0.007660761135179466,
                0.01141203638600162,
                0.01586133273977301,
                0.02092133533779397,
                0.02650028060249181,
                0.03250944464571953,
                0.03886710381241729,
                0.04550026389635835,
                0.0523450632731632,
                0.05934643879191983,
                0.06645742001693121,
                0.07363827012030258,
                0.0808555983700521,
                0.08808151166219025,
                0.09529283802345645,
                0.1024704348597497,
                0.1095985833991157,
                0.1166644647810232,
                0.1236577104028335,
                0.1305700181157362,
                0.1373948258558012,
                0.144127034816016,
                0.1507627750266455,
                0.157299207050285,
                0.1637343543245922,
                0.1700669614539048,
                0.176296374440511,
                0.182422439451736,
                0.1884454172427062,
                0.1943659108031092,
                0.2001848041775792,
                0.2059032107320683,
                0.2115224294107291,
                0.217043907756942,
                0.2224692106646891,
                0.2277999939882293,
                0.2330379822739055,
                0.2381849499920227,
                0.2432427057426569,
                0.2482130789899236,
                0.2530979089471155,
                0.2578990352923395,
                0.2626182904425206,
                0.2672574931543877,
                0.2718184432554954,
                0.2763029173374835,
                0.2807126652684959,
                0.2850494074026128,
                0.2893148323819879,
                0.2935105954425006,
                0.2976383171466433,
                0.3016995824783477,
                0.3056959402438416,
                0.309628902730633,
                0.3134999455835882,
                0.317310507862914,
            ];
            ddouble[] expected_nu2 = [
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
            ddouble[] expected_nu4 = [
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

            foreach ((InverseChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
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