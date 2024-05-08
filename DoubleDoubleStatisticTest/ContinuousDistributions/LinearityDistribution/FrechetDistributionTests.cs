
using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.LinearityDistribution {
    [TestClass()]
    public class FrechetDistributionTests {
        readonly FrechetDistribution dist_alpha1mu0theta1 = new(alpha: 1, mu: 0, theta: 1);
        readonly FrechetDistribution dist_alpha2mu1theta1 = new(alpha: 2, mu: 1, theta: 1);
        readonly FrechetDistribution dist_alpha1mu0theta2 = new(alpha: 1, mu: 0, theta: 2);
        readonly FrechetDistribution dist_alpha3mu0theta4 = new(alpha: 3, mu: 0, theta: 4);
        readonly FrechetDistribution dist_alpha4mu0theta4 = new(alpha: 4, mu: 0, theta: 4);
        readonly FrechetDistribution dist_alpha5mu0theta6 = new(alpha: 5, mu: 0, theta: 6);

        FrechetDistribution[] Dists => [
            dist_alpha1mu0theta1,
            dist_alpha2mu1theta1,
            dist_alpha1mu0theta2,
            dist_alpha3mu0theta4,
            dist_alpha4mu0theta4,
            dist_alpha5mu0theta6,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (FrechetDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Mu={dist.Mu}");
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
            foreach (FrechetDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (FrechetDistribution dist in Dists) {
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
            foreach (FrechetDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (FrechetDistribution dist in Dists) {
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
            foreach (FrechetDistribution dist in Dists) {
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
            foreach (FrechetDistribution dist in Dists) {
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
            foreach (FrechetDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (FrechetDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (FrechetDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (FrechetDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (FrechetDistribution dist in Dists) {
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
            foreach (FrechetDistribution dist in Dists) {
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

            foreach (FrechetDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
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
            foreach (FrechetDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha1mu0theta1 = [
                0,
                6.569209407687221e-25,
                1.296810552227244e-11,
                2.47268327465719e-7,
                2.880900472813034e-5,
                4.523249782025746e-4,
                0.002652057730040834,
                0.00894275766023249,
                0.02146960818576076,
                0.04126279225253728,
                0.06805738590920433,
                0.1006399284420458,
                0.1373283553800943,
                0.1763632533802715,
                0.2161512112925936,
                0.2553804455960597,
                0.2930502222197469,
                0.3284513978525497,
                0.3611243555972788,
                0.3908103620307,
                0.41740496873847,
                0.4409173851495388,
                0.461437044595815,
                0.4791071742417697,
                0.4941045420288109,
                0.5066243469692092,
                0.5168692329177461,
                0.5250415190253408,
                0.5313378863241238,
                0.535945905822235,
                0.5390419240450445,
                0.5407899320039895,
                0.5413411329464508,
                0.5408339950240924,
                0.539394630076989,
                0.5371373819783287,
                0.5341655400488015,
                0.5305721171612513,
                0.526440650136322,
                0.5218459933296126,
                0.5168550860663178,
                0.5115276816687313,
                0.5059170309180893,
                0.5000705164065503,
                0.4940302367543292,
                0.4878335413819951,
                0.481513517656564,
                0.4750994329347152,
                0.4686171344279587,
                0.4620894100020902,
                0.4555363130625133,
                0.4489754546160243,
                0.4424222654733848,
                0.4358902313906301,
                0.4293911037585094,
                0.4229350882513462,
                0.4165310136476348,
                0.4101864828406864,
                0.4039080078723559,
                0.3977011306486018,
                0.3915705308335359,
                0.3855201222692018,
                0.3795531391315055,
                0.3736722129081446,
                0.3678794411714423,
            ];
            ddouble[] expected_dist_alpha2mu1theta1 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.0,
                4.32475933767926e-194,
                5.420041057656411e-108,
                2.936648652185994e-68,
                9.376809962577158e-47,
                7.599803643498192e-34,
                1.642302351921805e-25,
                7.861074673028946e-20,
                8.528407410524779e-16,
                7.834030681931388e-13,
                1.345119648071087e-10,
                7.110181739988021e-9,
                1.60441773539919e-7,
                1.928441070156118e-6,
                1.440450236406517e-5,
                7.463850865341053e-5,
                2.906780816707425e-4,
                9.030064202402137e-4,
                0.002340477314114205,
                0.005238328909347808,
                0.0103980099012845,
                0.01869316424647374,
                0.03094709418940295,
                0.04781427499693126,
                0.069692345977235,
                0.09667795650059999,
                0.1285667703725734,
                0.1648888861284311,
                0.2049668899230828,
                0.2479837572675829,
                0.2930502222197469,
                0.3392645454346194,
                0.3857608201976285,
                0.4317445335348731,
                0.4765158886799499,
                0.5194824532686129,
                0.5601631796355283,
                0.5981859248382492,
                0.6332804337115114,
                0.6652684578813839,
                0.6940523499859765,
                0.7196031469426074,
                0.7419488670991105,
                0.7611635055255935,
                0.7773570214287091,
                0.7906664678684259,
                0.8012483100732022,
                0.8092719073045349,
                0.8149140873066081,
                0.8183547156048132,
                0.8197731490370573,
                0.8193454597100533,
                0.8172423188189989,
                0.8136274370174086,
                0.8086564675036577,
                0.802476288455511,
                0.7952245920488015,
                0.7870297174962653,
                0.7780106750101671,
                0.7682773161548931,
                0.757930613646292,
                0.7470630202764424,
                0.7357588823428847,
            ];
            ddouble[] expected_dist_alpha1mu0theta2 = [
                0,
                2.107153918068667e-52,
                3.28460470384361e-25,
                2.686887851119367e-16,
                6.484052761136218e-12,
                2.497532786937984e-9,
                1.236341637328595e-7,
                1.913414459918309e-6,
                1.440450236406517e-5,
                6.733987303835941e-5,
                2.261624891012873e-4,
                5.984061615308756e-4,
                0.001326028865020417,
                0.002566687264231289,
                0.004471378830116245,
                0.007165192235541721,
                0.01073480409288038,
                0.01522334604351754,
                0.02063139612626864,
                0.02692212832249736,
                0.03402869295460217,
                0.04186229979123735,
                0.05031996422102292,
                0.05929131301408926,
                0.06866417769004717,
                0.0783289272894219,
                0.08818162669013575,
                0.0981261752903412,
                0.1080756056462968,
                0.1179527196029126,
                0.1276902227980298,
                0.1372304952564264,
                0.1465251111098734,
                0.1555341971089751,
                0.1642256989262748,
                0.1725746067969741,
                0.1805621777986394,
                0.1881751807596599,
                0.19540518101535,
                0.2022478755796351,
                0.208702484369235,
                0.2147711995505877,
                0.2204586925747694,
                0.225771676772452,
                0.2307185222979075,
                0.235308919578628,
                0.2395535871208848,
                0.2434640194459607,
                0.2470522710144055,
                0.2503307721823377,
                0.2533121734846046,
                0.2560092148255754,
                0.2584346164588731,
                0.2606009889381697,
                0.2625207595126704,
                0.2642061127169806,
                0.2656689431620619,
                0.2669208187701348,
                0.2679729529111175,
                0.2688361840918225,
                0.2695209620225222,
                0.2700373390398765,
                0.2703949660019948,
                0.2706030918920845,
                0.2706705664732254,
            ];
            ddouble[] expected_dist_alpha3mu0theta4 = [
                0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                8.069511631208691e-296,
                5.47442713777133e-267,
                9.195946867460515e-242,
                1.344765246582697e-219,
                4.822292507963397e-200,
                1.000348187880722e-182,
                2.45472943022143e-167,
                1.297587595842781e-153,
                2.447945615691748e-141,
                2.526629710002853e-130,
                2.051560388589519e-120,
                1.786851640792088e-111,
                2.177636340036034e-103,
                4.667738673544714e-96,
                2.144093785817623e-89,
                2.504689364740962e-83,
                8.635704203923243e-78,
                1.000653794060498e-72,
                4.365871206068317e-68,
                7.9245529632541e-64,
                6.533018308942095e-60,
                2.64315397455981e-56,
                5.620228473311283e-53,
                6.674111052483511e-50,
                4.671904799007846e-47,
                2.02278607883108e-44,
                5.654791735479906e-42,
                1.060681761389934e-39,
                1.381718073414179e-37,
                1.289322065464924e-35,
                8.861518781583158e-34,
                4.600091230847113e-32,
                1.84496041707192e-30,
                5.835554503063314e-29,
                1.482954535722172e-27,
                3.079316909853385e-26,
            ];

            foreach ((FrechetDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1mu0theta1, expected_dist_alpha1mu0theta1),
                (dist_alpha2mu1theta1, expected_dist_alpha2mu1theta1),
                (dist_alpha1mu0theta2, expected_dist_alpha1mu0theta2),
                (dist_alpha3mu0theta4, expected_dist_alpha3mu0theta4),
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
            ddouble[] expected_dist_alpha1mu0theta1 = [
                0,
                1.603810890548638e-28,
                1.266416554909418e-14,
                5.433141960916678e-10,
                1.125351747192591e-7,
                2.760772572037198e-6,
                2.330910114293702e-5,
                1.069812317752422e-4,
                3.354626279025119e-4,
                8.159878350721482e-4,
                0.001661557273173934,
                0.002973005698605358,
                0.004827949993831441,
                0.007276706499332492,
                0.01034317319661825,
                0.01402846686013511,
                0.01831563888873418,
                0.02317442724106124,
                0.02856550078455038,
                0.03444397966139714,
                0.04076220397836621,
                0.04747181807884439,
                0.05452527577743516,
                0.06187687870456449,
                0.06948345122280154,
                0.07730474044329974,
                0.08530361363583897,
                0.09344611019762536,
                0.1017013923042268,
                0.1100416276358642,
                0.1184418290138037,
                0.1268796691054282,
                0.1353352832366127,
                0.143791069477841,
                0.1522314922775877,
                0.1606428937801398,
                0.1690133154060661,
                0.1773323311508186,
                0.1855908932609495,
                0.1937811903941261,
                0.2018965179946554,
                0.209931160372348,
                0.2178802838231224,
                0.2257398400477811,
                0.2335064790909134,
                0.2411774710201514,
                0.2487506355862523,
                0.2562242791388638,
                0.2635971381157268,
                0.2708683284704635,
                0.2780373004531941,
                0.2851037982070994,
                0.2920678236914142,
                0.2989296044863964,
                0.3056895650780794,
                0.3123483012598443,
                0.3189065573239704,
                0.3253652057493628,
                0.3317252291217299,
                0.3379877040497516,
                0.3441537868654124,
                0.3502247009188721,
                0.3562017252982195,
                0.3620861848223695,
                0.3678794411714423,
            ];
            ddouble[] expected_dist_alpha2mu1theta1 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.0,
                2.22718242869072e-198,
                6.616261056709486e-112,
                7.001515989747034e-72,
                3.863126663049062e-50,
                4.971947955550917e-37,
                1.603810890548638e-28,
                1.093048751189824e-22,
                1.626664621453243e-18,
                1.988810508279929e-15,
                4.433377746328046e-13,
                2.979482513952957e-11,
                8.39714482485843e-10,
                1.241395685534839e-8,
                1.125351747192591e-7,
                6.994228229793662e-7,
                3.233403343780081e-6,
                1.181358535085225e-5,
                3.571284964163521e-5,
                9.252960973638543e-5,
                2.111778439118907e-4,
                4.338068568932457e-4,
                8.159878350721482e-4,
                0.001424976438192465,
                0.002336335511962666,
                0.00362951701698553,
                0.005383105741918052,
                0.007670354926655402,
                0.0105554695661988,
                0.01409088919212067,
                0.01831563888873418,
                0.02325468057495864,
                0.0289191117802574,
                0.03530701994954621,
                0.04240479526949261,
                0.05018872204859743,
                0.05862669752685682,
                0.06767996001335164,
                0.07730474044329974,
                0.08745377995613644,
                0.09807767964508253,
                0.1091260669783895,
                0.1205485769175923,
                0.1322956551380533,
                0.1443191967731186,
                0.1565730375166008,
                0.1690133154060661,
                0.1815987217378068,
                0.1942906587854882,
                0.2070533206552393,
                0.2198537119670917,
                0.2326616172890732,
                0.2454495324907586,
                0.2581925675082614,
                0.2708683284704635,
                0.2834567857512311,
                0.2959401332928195,
                0.3083026434891996,
                0.3205305210155412,
                0.3326117582285953,
                0.3445359941274518,
                0.3562943783398888,
                0.3678794411714423,
            ];
            ddouble[] expected_dist_alpha1mu0theta2 = [
                0,
                2.572209372642415e-56,
                1.603810890548638e-28,
                2.951903156747352e-19,
                1.266416554909418e-14,
                7.62186519451289e-12,
                5.433141960916678e-10,
                1.144498395214809e-8,
                1.125351747192591e-7,
                6.658361469857315e-7,
                2.760772572037198e-6,
                8.838762883939932e-6,
                2.330910114293702e-5,
                5.295045747742773e-5,
                1.069812317752422e-4,
                1.967978824459091e-4,
                3.354626279025119e-4,
                5.370540779512414e-4,
                8.159878350721482e-4,
                0.00118638773491474,
                0.001661557273173934,
                0.002253573511710897,
                0.002973005698605358,
                0.003828748118219387,
                0.004827949993831441,
                0.005976022895005943,
                0.007276706499332492,
                0.008732175511066741,
                0.01034317319661825,
                0.01210915981275018,
                0.01402846686013511,
                0.01609845043230295,
                0.01831563888873418,
                0.02067587166158129,
                0.02317442724106124,
                0.02580613932205729,
                0.02856550078455038,
                0.03144675567138359,
                0.03444397966139714,
                0.03755114975056457,
                0.04076220397836621,
                0.0440710920952805,
                0.04747181807884439,
                0.05095847538479781,
                0.05452527577743516,
                0.05816657252767598,
                0.06187687870456449,
                0.06565088122023037,
                0.06948345122280154,
                0.0733696513683829,
                0.07730474044329974,
                0.08128417575211445,
                0.08530361363583897,
                0.0893589084383934,
                0.09344611019762536,
                0.09756146129991042,
                0.1017013923042268,
                0.1058625171123252,
                0.1100416276358642,
                0.1142356880888225,
                0.1184418290138037,
                0.1226573411337134,
                0.1268796691054282,
                0.1311064052392192,
                0.1353352832366127,
            ];
            ddouble[] expected_dist_alpha3mu0theta4 = [
                0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.771813958893504e-299,
                1.376583545653391e-270,
                2.636465569021802e-245,
                4.377491037053051e-223,
                1.775367170987173e-203,
                4.14997741579939e-186,
                1.143548746663466e-170,
                6.765899830122378e-157,
                1.424254291706552e-144,
                1.635516238397766e-133,
                1.473403199556881e-123,
                1.420062097543151e-114,
                1.9102915630541e-106,
                4.50902424737637e-99,
                2.275602885536562e-92,
                2.914361714456274e-86,
                1.099326478665469e-80,
                1.390894047338186e-75,
                6.613650288649728e-71,
                1.305926086473759e-66,
                1.16916840524747e-62,
                5.128393676442037e-59,
                1.180357463902296e-55,
                1.514905975421123e-52,
                1.14439601855912e-49,
                5.339540325204869e-47,
                1.606372452790408e-44,
                3.238292578502949e-42,
                4.527909483971119e-40,
                4.529520868182664e-38,
                3.33345874165257e-36,
                1.850760925305963e-34,
                7.930220597135365e-33,
                2.67687297375179e-31,
                7.252185594912769e-30,
                1.603810890548638e-28,
            ];

            foreach ((FrechetDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1mu0theta1, expected_dist_alpha1mu0theta1),
                (dist_alpha2mu1theta1, expected_dist_alpha2mu1theta1),
                (dist_alpha1mu0theta2, expected_dist_alpha1mu0theta2),
                (dist_alpha3mu0theta4, expected_dist_alpha3mu0theta4),
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