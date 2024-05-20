using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class BeniniDistributionTests {
        readonly BeniniDistribution dist_alpha0beta1sigma1 = new(alpha: 0, beta: 1, sigma: 1);
        readonly BeniniDistribution dist_alpha0beta2sigma2 = new(alpha: 0, beta: 2, sigma: 2);
        readonly BeniniDistribution dist_alpha1beta1sigma1 = new(alpha: 1, beta: 1, sigma: 1);
        readonly BeniniDistribution dist_alpha2beta2sigma1 = new(alpha: 2, beta: 2, sigma: 1);
        readonly BeniniDistribution dist_alpha1beta3sigma2 = new(alpha: 1, beta: 3, sigma: 2);
        readonly BeniniDistribution dist_alpha2beta1sigma2 = new(alpha: 2, beta: 1, sigma: 2);
        readonly BeniniDistribution dist_alpha3beta2sigma4 = new(alpha: 3, beta: 2, sigma: 4);
        readonly BeniniDistribution dist_alpha4beta3sigma5 = new(alpha: 4, beta: 3, sigma: 5);
        readonly BeniniDistribution dist_alpha5beta1sigma7 = new(alpha: 5, beta: 1, sigma: 7);
        readonly BeniniDistribution dist_alpha6beta2sigma3 = new(alpha: 6, beta: 2, sigma: 3);

        BeniniDistribution[] Dists => [
            dist_alpha0beta1sigma1,
            dist_alpha0beta2sigma2,
            dist_alpha1beta1sigma1,
            dist_alpha2beta2sigma1,
            dist_alpha1beta3sigma2,
            dist_alpha2beta1sigma2,
            dist_alpha3beta2sigma4,
            dist_alpha4beta3sigma5,
            dist_alpha5beta1sigma7,
            dist_alpha6beta2sigma3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BeniniDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Beta={dist.Beta}");
                Console.WriteLine($"Sigma={dist.Sigma}");
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
            foreach (BeniniDistribution dist in Dists) {
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
            foreach (BeniniDistribution dist in Dists) {
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
            foreach (BeniniDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (BeniniDistribution dist in Dists) {
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
            foreach (BeniniDistribution dist in Dists) {
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
            foreach (BeniniDistribution dist in Dists) {
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
            foreach (BeniniDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (BeniniDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (BeniniDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (BeniniDistribution dist in Dists) {
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
            foreach (BeniniDistribution dist in Dists) {
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
            foreach (BeniniDistribution dist in Dists) {
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

            foreach (BeniniDistribution dist in Dists) {

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

            foreach (BeniniDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                (BeniniDistribution? dist_fit, ddouble error) = BeniniDistribution.Fit(xs, (0.15, 0.85));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-1);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (BeniniDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha0beta1sigma1 = [
                0, 0, 0, 0,
                0,
                0.339687431946915,
                0.4586626359206531,
                0.467599433366062,
                0.4287137061346417,
                0.37345262379127,
                0.3165880653901416,
                0.2644095641594739,
                0.2190695014477827,
                0.1808004284252094,
                0.14901931413908,
                0.1228629066306601,
                0.101436227723353,
                0.08391844126260738,
                0.06959939898976976,
                0.05788388127383388,
                0.04828201583953105,
                0.04039490396115618,
                0.03389968297745798,
                0.02853584254346752,
                0.02409342726071066,
                0.02040320065280995,
                0.01732860476932321,
                0.01475926562936285,
                0.01260578471376748,
                0.01079557880294206,
                0.009269563789175855,
                0.00797951276445083,
                0.006885950391840917,
            ];
            ddouble[] expected_dist_alpha0beta2sigma2 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                0,
                0.2036621939029456,
                0.3231875589253868,
                0.3781704746006501,
                0.389129809290942,
                0.37293310948674,
                0.3418739793413527,
                0.3042195752631103,
                0.2651607724628187,
                0.2277209989682512,
                0.1934817779089898,
                0.163111685589876,
                0.1367306244362925,
                0.1141502930972116,
                0.09502717239215198,
                0.07895620634322802,
                0.06552554567193393,
                0.05434643952836809,
                0.04506771863507551,
                0.03738104650749978,
                0.03102089058636384,
                0.02576168540997147,
                0.02141369462277574,
                0.01781845976540501,
                0.01484433405100474,
            ];
            ddouble[] expected_dist_alpha1beta1sigma1 = [
                0, 0, 0, 0,
                1.0,
                0.8806627820635692,
                0.6828421729583712,
                0.5059346100420423,
                0.3689826375177149,
                0.268317554917228,
                0.1957373253866363,
                0.1436720867167794,
                0.1062574427451536,
                0.0792302160470934,
                0.05957016427622874,
                0.04515735022581306,
                0.03450540334977357,
                0.02656881248068304,
                0.02060806809300393,
                0.01609652610453613,
                0.01265633348311233,
                0.01001429430098892,
                0.007971346906413813,
                0.006381333900025148,
                0.005136137692138618,
                0.004155198867637859,
                0.003378070858599063,
                0.002759092356482137,
                0.002263547235443572,
                0.001864876813017722,
                0.001542641919393464,
                0.001281022996002698,
                0.001067708934024205,
            ];
            ddouble[] expected_dist_alpha2beta2sigma1 = [
                0, 0, 0, 0,
                2.0,
                1.340617181579065,
                0.7724319175777845,
                0.4227449297824998,
                0.2282169190990082,
                0.1235666064592245,
                0.067629300409978,
                0.03755262909721428,
                0.02118833481784338,
                0.01215356236557637,
                0.007086021256086286,
                0.004197570486661366,
                0.002524786978896961,
                0.00154093861223158,
                9.536139275063332e-4,
                5.979736381108507e-4,
                3.796811849534689e-4,
                2.439510355871991e-4,
                1.585138204647863e-4,
                1.041028781746988e-4,
                6.906460493821202e-5,
                4.62622578727489e-5,
                3.127320585144776e-5,
                2.132561955604595e-5,
                1.46634669402889e-5,
                1.016275096667693e-5,
                7.096926020700475e-6,
                4.991930006142649e-6,
                3.535656379489612e-6,
            ];
            ddouble[] expected_dist_alpha1beta3sigma2 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                0.5,
                0.6467657230198793,
                0.6445845834812304,
                0.5678535737694549,
                0.4658415482550439,
                0.3653132763356537,
                0.278052575848181,
                0.2073974097454759,
                0.152577825176144,
                0.111208080665884,
                0.08056189809483673,
                0.0581417576324068,
                0.04187611508408082,
                0.03013915996207544,
                0.02169760176085379,
                0.0156363418864218,
                0.01128629010188817,
                0.008163017309951125,
                0.00591805538622622,
                0.004301759130380682,
                0.003135697353024648,
                0.002292469436600155,
                0.001681114934612257,
                0.001236644899483714,
                9.125669078436865e-4,
            ];
            ddouble[] expected_dist_alpha2beta1sigma2 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                1.0,
                0.7742390063962991,
                0.5958302474278426,
                0.4582619394828781,
                0.3533030851009912,
                0.2734915347255653,
                0.2127627268214018,
                0.1664186366387898,
                0.1309021054920272,
                0.1035481504242303,
                0.08236803687271665,
                0.06587783742221688,
                0.05296788492344319,
                0.04280542589870534,
                0.03476277109307855,
                0.02836448891559142,
                0.02324861972350773,
                0.0191381434741214,
                0.01581992785235019,
                0.01312913582606243,
                0.01093762595079678,
                0.009145282133258105,
                0.007673501157793345,
                0.006460276523015739,
                0.00545646872108861,
            ];
            ddouble[] expected_dist_alpha3beta2sigma4 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0.75,
                0.6314103143864337,
                0.5269279884471111,
                0.4369928135111838,
                0.3608171469312701,
                0.29702643983212,
                0.2440388146242752,
                0.2002780531719002,
                0.1642835490765985,
                0.1347586218260343,
                0.1105834467227162,
                0.09080884501544986,
                0.07464072498411689,
                0.06142088878801354,
                0.05060739544517711,
                0.04175613930040694,
                0.03450439819160096,
            ];

            foreach ((BeniniDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha0beta1sigma1, expected_dist_alpha0beta1sigma1),
                (dist_alpha0beta2sigma2, expected_dist_alpha0beta2sigma2),
                (dist_alpha1beta1sigma1, expected_dist_alpha1beta1sigma1),
                (dist_alpha2beta2sigma1, expected_dist_alpha2beta2sigma1),
                (dist_alpha1beta3sigma2, expected_dist_alpha1beta3sigma2),
                (dist_alpha2beta1sigma2, expected_dist_alpha2beta1sigma2),
                (dist_alpha3beta2sigma4, expected_dist_alpha3beta2sigma4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 4, i++) {
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
            ddouble[] expected_dist_alpha0beta1sigma1 = [
                0, 0, 0, 0,
                0,
                0.04857369295931657,
                0.1515990647246444,
                0.2688742651368541,
                0.381496862198424,
                0.4819107817618907,
                0.5681118798088772,
                0.6406061456429085,
                0.7008915196369652,
                0.7507322353845062,
                0.791833087102978,
                0.8257106623144792,
                0.853658457297035,
                0.876754199933744,
                0.8958839165706342,
                0.9117705658171145,
                0.9250017421198471,
                0.9360542591250628,
                0.9453150124570011,
                0.9530982425553567,
                0.9596596066472736,
                0.9652075483129584,
                0.9699124372247905,
                0.9739138975059819,
                0.9773266784596373,
                0.980245358837086,
                0.9827481204529366,
                0.9848997802270819,
                0.9867542313571782,
            ];
            ddouble[] expected_dist_alpha0beta2sigma2 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                0.0,
                0.02736431072620149,
                0.0947879822709271,
                0.1835793527218156,
                0.2802158530239018,
                0.3758943901446522,
                0.4654551598208249,
                0.5462902912103687,
                0.6174538685296047,
                0.6790095349269214,
                0.7315835619454247,
                0.7760743347393126,
                0.8134726516377783,
                0.8447567863955276,
                0.8708360574503538,
                0.8925247296977203,
                0.910534116974916,
                0.925475024741697,
                0.9378655815235948,
                0.9481414350088556,
                0.9566665363749236,
                0.9637435265545696,
                0.9696232267691425,
                0.9745130291702582,
                0.9785841528793163,
            ];
            ddouble[] expected_dist_alpha1beta1sigma1 = [
                0, 0, 0, 0,
                0.0,
                0.2388589543674533,
                0.4343993764830962,
                0.582213865792488,
                0.690748431099212,
                0.769738125227507,
                0.8272447519235508,
                0.8693113256883304,
                0.9002971732123217,
                0.9233022262721557,
                0.9405237391722794,
                0.9535228432838612,
                0.9634146143242588,
                0.9710009882197045,
                0.976863092571252,
                0.9814253822772873,
                0.9850003484239694,
                0.9878198588809644,
                0.990057274992182,
                0.9918431726183229,
                0.993276601107879,
                0.9944332077300734,
                0.9953711441884293,
                0.9961353922231084,
                0.9967609540656625,
                0.9972752219085635,
                0.9976997493937249,
                0.99805158454543,
                0.9983442789196473,
            ];
            ddouble[] expected_dist_alpha2beta2sigma1 = [
                0, 0, 0, 0,
                0.0,
                0.4206643086533934,
                0.6800959346772897,
                0.8254547460639428,
                0.9043634671324011,
                0.9469794690262567,
                0.9701556242620445,
                0.9829204704066583,
                0.9900593463305463,
                0.9941174515051924,
                0.9964625743979529,
                0.9978398739035834,
                0.9986615095549572,
                0.9991590573157663,
                0.9994646835146336,
                0.9996549835764551,
                0.9997750104525976,
                0.9998516441623204,
                0.999901142219419,
                0.9999334661670656,
                0.9999547959073374,
                0.9999690108238235,
                0.9999785736938757,
                0.9999850648067308,
                0.9999895085814352,
                0.9999925755843524,
                0.9999947088471484,
                0.9999962036772164,
                0.999997258587704,
            ];
            ddouble[] expected_dist_alpha1beta3sigma2 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                0.0,
                0.1473461067046404,
                0.3110059783065461,
                0.4635017976996315,
                0.5928896376727362,
                0.6965872937270612,
                0.7766745776609701,
                0.8370081695752125,
                0.8816970086658531,
                0.9144186242524552,
                0.9381939277537987,
                0.9553838393003249,
                0.9677764216606421,
                0.9766981915086592,
                0.9831197355793181,
                0.9877446711279396,
                0.9910799985613433,
                0.9934896739550316,
                0.9952344284309855,
                0.9965008950768649,
                0.9974226876148664,
                0.9980955466977492,
                0.9985881739500189,
                0.998949960525418,
                0.9992164929735171,
            ];
            ddouble[] expected_dist_alpha2beta1sigma2 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                0.0,
                0.2207621207237613,
                0.3910871634939627,
                0.5220839746462762,
                0.6229329176553975,
                0.7008270577861754,
                0.7612650661671361,
                0.8084039760611358,
                0.845374215549606,
                0.8745333803671989,
                0.8976613889900031,
                0.9161071963248266,
                0.9308979007694204,
                0.942819533175285,
                0.9524768457048474,
                0.9603376379147974,
                0.9667657244041072,
                0.972045554833542,
                0.9764006850068172,
                0.980007703081369,
                0.9830067826206512,
                0.9855097206572548,
                0.987606091542363,
                0.9893679821452858,
                0.9908536535810647,
            ];
            ddouble[] expected_dist_alpha3beta2sigma4 = [
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0.0,
                0.1723993560270805,
                0.316886868438704,
                0.4370788941638101,
                0.5365314469227147,
                0.6185200699683757,
                0.6859448749763859,
                0.7413080418174625,
                0.7867306231181932,
                0.8239884504810371,
                0.8545552698015758,
                0.8796464564161901,
                0.9002598548936815,
                0.9172121867379575,
                0.931170556770284,
                0.9426791789050019,
                0.9521817335662006,
            ];

            foreach ((BeniniDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha0beta1sigma1, expected_dist_alpha0beta1sigma1),
                (dist_alpha0beta2sigma2, expected_dist_alpha0beta2sigma2),
                (dist_alpha1beta1sigma1, expected_dist_alpha1beta1sigma1),
                (dist_alpha2beta2sigma1, expected_dist_alpha2beta2sigma1),
                (dist_alpha1beta3sigma2, expected_dist_alpha1beta3sigma2),
                (dist_alpha2beta1sigma2, expected_dist_alpha2beta1sigma2),
                (dist_alpha3beta2sigma4, expected_dist_alpha3beta2sigma4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 4, i++) {
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