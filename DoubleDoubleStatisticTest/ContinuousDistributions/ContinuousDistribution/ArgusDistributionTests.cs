using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class ArgusDistributionTests {
        readonly ArgusDistribution dist_a1 = new(alpha: 1);
        readonly ArgusDistribution dist_a2 = new(alpha: 2);
        readonly ArgusDistribution dist_a3 = new(alpha: 3);

        ArgusDistribution[] Dists => [
            dist_a1,
            dist_a2,
            dist_a3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
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
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (ArgusDistribution dist in Dists) {
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
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (ArgusDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (ArgusDistribution dist in Dists) {
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
            foreach (ArgusDistribution dist in Dists) {
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
            foreach (ArgusDistribution dist in Dists) {
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

            foreach (ArgusDistribution dist in Dists) {

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
            foreach (ArgusDistribution dist in Dists) {
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
            ddouble[] expected_dist1_a1 = [
                0.0,
                0.07609215385801377,
                0.1521837619708911,
                0.2282720816212929,
                0.3043499156261317,
                0.3804032313652878,
                0.4564085883385007,
                0.5323302982622434,
                0.6081172301875762,
                0.6836991571769357,
                0.758982519442446,
                0.8338454496655971,
                0.9081318667927001,
                0.981644390978912,
                1.054135758655976,
                1.125298314028752,
                1.194751007931254,
                1.26202312521267,
                1.326533652561861,
                1.387564731391194,
                1.444226914361304,
                1.495412779534782,
                1.539733518370517,
                1.575429749302947,
                1.600241667217864,
                1.611211737436153,
                1.604368331527387,
                1.574182025877242,
                1.512540496548527,
                1.40654684005021,
                1.232747604323243,
                0.9354262812899866,
                0.0,
            ];
            ddouble[] expected_dist1_a2 = [
                0.0,
                0.03660627878095898,
                0.07353473670527115,
                0.111111214096793,
                0.1496688750852653,
                0.1895518678171996,
                0.2311189671396426,
                0.2747471660730665,
                0.3208351542369129,
                0.3698065802096646,
                0.4221129355105346,
                0.4782358130662516,
                0.5386881718660567,
                0.6040140661796135,
                0.6747860488206565,
                0.7515990984988878,
                0.8350593982287975,
                0.9257655230871523,
                1.024278452115589,
                1.131075092131814,
                1.246477343679383,
                1.370544554673818,
                1.502910427619266,
                1.642534067388801,
                1.787314898047045,
                1.933484177028422,
                2.074612248519293,
                2.199910459486697,
                2.291113541344451,
                2.316091395389702,
                2.213145163229877,
                1.836334899702237,
                0.0,
            ];
            ddouble[] expected_dist1_a3 = [
                0.0,
                0.007734565666430715,
                0.01565142043050036,
                0.02393980206554285,
                0.03280319771457822,
                0.04246736928514048,
                0.05318951867780745,
                0.06526907445562752,
                0.07906067537632926,
                0.09499005162722683,
                0.1135736664862875,
                0.1354431841549912,
                0.1613760765133034,
                0.1923339703117352,
                0.2295106533636894,
                0.274391967300232,
                0.3288300346647916,
                0.3951342345259016,
                0.4761807299406002,
                0.5755405418625463,
                0.697621975893433,
                0.8478143870598708,
                1.03260247018106,
                1.259584955096584,
                1.537261495940612,
                1.874310206901225,
                2.277783999807646,
                2.749007570095817,
                3.274412977025323,
                3.804337355552066,
                4.198462668538827,
                4.043057596142619,
                0.0,
            ];

            foreach ((ArgusDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1, expected_dist1_a1),
                (dist_a2, expected_dist1_a2),
                (dist_a3, expected_dist1_a3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 32, i++) {
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
            ddouble[] expected_dist1_a1 = [
                0.0,
                0.001188940093146074,
                0.004755754691631364,
                0.01070039250312749,
                0.01902265293303862,
                0.02972201008680075,
                0.0427973595610196,
                0.0582466808777462,
                0.0760666058715388,
                0.09625188036781818,
                0.1187947029367229,
                0.1436839201640275,
                0.1709040524745671,
                0.2004341176958946,
                0.2322462107315592,
                0.2663037861675353,
                0.3025595752580474,
                0.3409530478813269,
                0.381407301206006,
                0.4238252160208369,
                0.4680846625879214,
                0.5140324498657082,
                0.5614765766352401,
                0.6101761272500241,
                0.6598277953628332,
                0.7100473890211734,
                0.760343495068251,
                0.8100781100819731,
                0.858403759910964,
                0.9041531397072684,
                0.9456153776771437,
                0.9799518108582926,
                1.0,
            ];
            ddouble[] expected_dist1_a2 = [
                0.0,
                5.715544757124436e-4,
                0.002291247185402523,
                0.005174223179557291,
                0.009245914914416287,
                0.014542328330723,
                0.02111044309980381,
                0.02900872602246041,
                0.03830775510243423,
                0.04909094927998403,
                0.06145539472082218,
                0.07551275223662446,
                0.0913902209095262,
                0.109231518937409,
                0.1291978221376269,
                0.1514685705732163,
                0.1762420101721118,
                0.2037352726964612,
                0.2341837044480199,
                0.2678390169589802,
                0.3049656284292662,
                0.3458342552959465,
                0.390711336316009,
                0.4398421178912141,
                0.4934240004975684,
                0.5515646588871075,
                0.6142157094930538,
                0.6810655093391718,
                0.7513594647420306,
                0.8235794828603886,
                0.8948060568501798,
                0.9591434972096481,
                1.0,
            ];
            ddouble[] expected_dist1_a3 = [
                0.0,
                1.206168745463421e-4,
                4.85306800605545e-4,
                0.001102695805016984,
                0.001987526471419643,
                0.003161226227797131,
                0.004652743082039001,
                0.00649969174355769,
                0.008749869508721142,
                0.01146322114098353,
                0.01471435630585261,
                0.01859575314089079,
                0.02322181858978445,
                0.02873402156627303,
                0.03530736997382433,
                0.04315856748790181,
                0.05255625945245823,
                0.06383385419218901,
                0.07740547406372089,
                0.09378562296281534,
                0.113613099839909,
                0.13767943526544,
                0.1669614750321499,
                0.202656282476378,
                0.2462134947685822,
                0.2993540636932566,
                0.3640516118849164,
                0.4424260294665008,
                0.5364403948576169,
                0.6471517964648982,
                0.7728736595239356,
                0.9041344744030125,
                1.0,
            ];

            foreach ((ArgusDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1, expected_dist1_a1),
                (dist_a2, expected_dist1_a2),
                (dist_a3, expected_dist1_a3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 32, i++) {
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