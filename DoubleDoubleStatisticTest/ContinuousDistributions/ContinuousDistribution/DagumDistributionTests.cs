using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class DagumDistributionTests {
        readonly DagumDistribution dist_a1b1 = new(a: 1, p: 1);
        readonly DagumDistribution dist_a2b1 = new(a: 2, p: 1);
        readonly DagumDistribution dist_a1b2 = new(a: 1, p: 2);
        readonly DagumDistribution dist_a2b2 = new(a: 2, p: 2);
        readonly DagumDistribution dist_a3b4 = new(a: 3, p: 4);

        DagumDistribution[] Dists => [
            dist_a1b1,
            dist_a2b1,
            dist_a1b2,
            dist_a2b2,
            dist_a3b4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"A={dist.A}");
                Console.WriteLine($"P={dist.P}");
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
            foreach (DagumDistribution dist in Dists) {
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
            foreach (DagumDistribution dist in Dists) {
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
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (DagumDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (DagumDistribution dist in Dists) {
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
            foreach (DagumDistribution dist in Dists) {
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

            foreach (DagumDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.08, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (DagumDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 20000).ToArray();

                (DagumDistribution? dist_fit, ddouble error) = DagumDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (DagumDistribution dist in Dists) {
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
                1.0,
                0.7901234567901234,
                0.64,
                0.5289256198347108,
                0.4444444444444444,
                0.378698224852071,
                0.3265306122448979,
                0.2844444444444444,
                0.25,
                0.2214532871972318,
                0.1975308641975309,
                0.1772853185595568,
                0.16,
                0.145124716553288,
                0.1322314049586777,
                0.1209829867674858,
                0.1111111111111111,
                0.1024,
                0.09467455621301775,
                0.0877914951989026,
                0.08163265306122448,
                0.07609988109393578,
                0.07111111111111111,
                0.06659729448491156,
                0.0625,
                0.05876951331496786,
                0.05536332179930796,
                0.05224489795918368,
                0.04938271604938271,
                0.04674945215485756,
                0.0443213296398892,
                0.042077580539119,
                0.04,
            ];
            ddouble[] expected_dist_a2b1 = [
                0.0,
                0.2423668639053254,
                0.4429065743944637,
                0.5764683805592044,
                0.64,
                0.6463830324453983,
                0.6144000000000001,
                0.5613595426423369,
                0.5,
                0.4383353151010702,
                0.380725758477097,
                0.3291161431701972,
                0.2840236686390533,
                0.2452062112030061,
                0.2120710059171598,
                0.1839058440392237,
                0.16,
                0.1397009846800793,
                0.1224359655648847,
                0.1077148788927336,
                0.09512485136741973,
                0.0843211449857857,
                0.07501731578666951,
                0.06697587651322767,
                0.06,
                0.05392641151328886,
                0.04861943024105186,
                0.04396603954208403,
                0.03987184051263795,
                0.03625774549006441,
                0.03305728207158967,
                0.03021439619274242,
                0.02768166089965398,
            ];
            ddouble[] expected_dist_a1b2 = [
                0.0,
                0.1755829903978052,
                0.256,
                0.2885048835462058,
                0.2962962962962963,
                0.2913063268092854,
                0.2798833819241983,
                0.2654814814814815,
                0.25,
                0.2344799511500102,
                0.2194787379972565,
                0.2052777372794868,
                0.192,
                0.179678220494547,
                0.1682945154019534,
                0.1578038957836772,
                0.1481481481481481,
                0.139264,
                0.1310878470641784,
                0.1235584006503074,
                0.1166180758017493,
                0.1102136208946656,
                0.1042962962962963,
                0.09882179181632036,
                0.09375,
                0.0890447171438907,
                0.08467331569305923,
                0.0806064139941691,
                0.07681755829903977,
                0.07328292499950644,
                0.06998104679982504,
                0.06689256393398406,
                0.064,
            ];
            ddouble[] expected_dist_a2b2 = [
                0.0,
                0.007457441966317706,
                0.05210665581111337,
                0.1421428883570641,
                0.256,
                0.3631365350816844,
                0.442368,
                0.4868427891942391,
                0.5,
                0.4897263520439542,
                0.4642997054598743,
                0.4305194953902039,
                0.3932635411925353,
                0.3557068643202407,
                0.3197378243058717,
                0.2863585806839123,
                0.256,
                0.2287455216574669,
                0.2044806847578486,
                0.1829885707307144,
                0.1640083644265858,
                0.1472698017375505,
                0.1325123388348469,
                0.1194949027841398,
                0.108,
                0.09783456370335425,
                0.08882901308905691,
                0.08083541696388212,
                0.07372529000450036,
                0.06738732366219706,
                0.06172521548637077,
                0.05665567754385456,
                0.05210665581111337,
            ];
            ddouble[] expected_dist_a3b4 = [
                0.0,
                1.383421014580496e-9,
                2.647612532420487e-6,
                1.913952564503686e-4,
                0.003251536859218615,
                0.02288305898344285,
                0.08720633837104555,
                0.2127105646901902,
                0.375,
                0.5240224464629127,
                0.6219883326255994,
                0.6595292400852505,
                0.6475820598560125,
                0.603772105004783,
                0.543910703631041,
                0.4791509355714292,
                0.4161967179799827,
                0.3584929476662869,
                0.3074134597020288,
                0.2631471970922004,
                0.2252750379712653,
                0.1931140812337905,
                0.1659089099913359,
                0.1429296667089772,
                0.1235164756500268,
                0.1070945312793911,
                0.093174114708599,
                0.0813435856332996,
                0.07125971589488087,
                0.06263761735836203,
                0.05524133708224029,
                0.04887555444281098,
                0.04337848373117727,
            ];

            foreach ((DagumDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a1b2, expected_dist_a1b2), (dist_a2b2, expected_dist_a2b2),
                (dist_a3b4, expected_dist_a3b4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 8, i++) {
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
                0,
                0.1111111111111111,
                0.2,
                0.2727272727272728,
                0.3333333333333333,
                0.3846153846153846,
                0.4285714285714286,
                0.4666666666666667,
                0.5,
                0.5294117647058824,
                0.5555555555555556,
                0.5789473684210527,
                0.6000000000000001,
                0.6190476190476191,
                0.6363636363636364,
                0.6521739130434783,
                0.6666666666666666,
                0.6799999999999999,
                0.6923076923076923,
                0.7037037037037037,
                0.7142857142857143,
                0.7241379310344828,
                0.7333333333333333,
                0.7419354838709677,
                0.75,
                0.7575757575757576,
                0.7647058823529411,
                0.7714285714285715,
                0.7777777777777779,
                0.7837837837837838,
                0.7894736842105263,
                0.7948717948717948,
                0.8,
            ];
            ddouble[] expected_dist_a2b1 = [
                0,
                0.01538461538461539,
                0.05882352941176471,
                0.1232876712328767,
                0.2,
                0.2808988764044944,
                0.36,
                0.4336283185840709,
                0.5,
                0.5586206896551724,
                0.6097560975609756,
                0.654054054054054,
                0.6923076923076923,
                0.7253218884120172,
                0.7538461538461539,
                0.7785467128027681,
                0.8,
                0.8186968838526911,
                0.8350515463917526,
                0.8494117647058823,
                0.8620689655172414,
                0.8732673267326732,
                0.8832116788321168,
                0.8920741989881956,
                0.8999999999999999,
                0.9071117561683599,
                0.9135135135135135,
                0.9192938209331653,
                0.9245283018867924,
                0.9292817679558012,
                0.9336099585062241,
                0.937560975609756,
                0.9411764705882353,
            ];
            ddouble[] expected_dist_a1b2 = [
                0,
                0.01234567901234568,
                0.04,
                0.0743801652892562,
                0.1111111111111111,
                0.1479289940828402,
                0.1836734693877551,
                0.2177777777777778,
                0.25,
                0.2802768166089966,
                0.308641975308642,
                0.3351800554016621,
                0.36,
                0.383219954648526,
                0.4049586776859505,
                0.4253308128544424,
                0.4444444444444444,
                0.4623999999999999,
                0.4792899408284024,
                0.4951989026063101,
                0.5102040816326531,
                0.5243757431629014,
                0.5377777777777777,
                0.5504682622268471,
                0.5625000000000001,
                0.573921028466483,
                0.5847750865051903,
                0.5951020408163266,
                0.6049382716049384,
                0.614317019722425,
                0.6232686980609419,
                0.6318211702827088,
                0.64,
            ];
            ddouble[] expected_dist_a2b2 = [
                0,
                2.366863905325444e-4,
                0.003460207612456748,
                0.0151998498780259,
                0.04,
                0.0789041787653074,
                0.1296,
                0.1880335186780485,
                0.25,
                0.3120570749108205,
                0.37180249851279,
                0.4277867056245435,
                0.4792899408284024,
                0.5260918418095747,
                0.5682840236686392,
                0.606134984015996,
                0.64,
                0.6702645876301069,
                0.6973110851312573,
                0.7215003460207612,
                0.7431629013079668,
                0.7625958239388294,
                0.7800628696254462,
                0.7957963765004308,
                0.8099999999999999,
                0.822851738178846,
                0.8345069393718042,
                0.8451011292058985,
                0.8547525809896759,
                0.8635646042550594,
                0.8716275546219934,
                0.8790205829863176,
                0.8858131487889274,
            ];
            ddouble[] expected_dist_a3b4 = [
                0,
                1.443878134114165e-11,
                5.602044746332411e-8,
                6.296511427413863e-6,
                1.524157902758726e-4,
                0.001482799130498565,
                0.007749782023208152,
                0.02590073099102353,
                0.0625,
                0.1190756572376357,
                0.1913343015401014,
                0.2720262915903377,
                0.3541464389837568,
                0.4325977014557587,
                0.504427602423449,
                0.5683775898944536,
                0.6242950769699741,
                0.6726483581242534,
                0.7141959185850848,
                0.7497864136874506,
                0.7802494804733925,
                0.8063432257864281,
                0.8287344478571125,
                0.8479964318170533,
                0.8646153295501876,
                0.8790000667600177,
                0.8914930838638833,
                0.9023805790271222,
                0.9119016768423038,
                0.9202563473587361,
                0.9276121007120325,
                0.9341095726459869,
                0.9398671475088407,
            ];

            foreach ((DagumDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a1b2, expected_dist_a1b2), (dist_a2b2, expected_dist_a2b2),
                (dist_a3b4, expected_dist_a3b4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 8, i++) {
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