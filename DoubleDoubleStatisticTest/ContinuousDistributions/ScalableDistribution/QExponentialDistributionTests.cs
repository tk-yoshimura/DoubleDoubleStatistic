using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class QExponentialDistributionTests {
        readonly QExponentialDistribution dist_q025theta1 = new(q: 0.25, theta: 1);
        readonly QExponentialDistribution dist_q050theta2 = new(q: 0.50, theta: 2);
        readonly QExponentialDistribution dist_q075theta2 = new(q: 0.75, theta: 2);
        readonly QExponentialDistribution dist_q100theta3 = new(q: 1.00, theta: 3);
        readonly QExponentialDistribution dist_q125theta6 = new(q: 1.25, theta: 6);
        readonly QExponentialDistribution dist_q150theta1 = new(q: 1.50, theta: 1);
        readonly QExponentialDistribution dist_q175theta4 = new(q: 1.75, theta: 4);

        QExponentialDistribution[] Dists => [
            dist_q025theta1,
            dist_q050theta2,
            dist_q075theta2,
            dist_q100theta3,
            dist_q125theta6,
            dist_q150theta1,
            dist_q175theta4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (QExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Q={dist.Q}");
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
            foreach (QExponentialDistribution dist in Dists) {
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
            foreach (QExponentialDistribution dist in Dists) {
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
            foreach (QExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (QExponentialDistribution dist in Dists) {
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
            foreach (QExponentialDistribution dist in Dists) {
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
            foreach (QExponentialDistribution dist in Dists) {
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
            foreach (QExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (QExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (QExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (QExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (QExponentialDistribution dist in Dists) {
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
            foreach (QExponentialDistribution dist in Dists) {
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

            foreach (QExponentialDistribution dist in Dists) {
                // ignore
                if (dist.Q > 1.5) {
                    continue;
                }

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

            foreach (QExponentialDistribution dist in Dists) {
                // ignore
                if (dist.Q > 1.5) {
                    continue;
                }

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                (QExponentialDistribution? dist_fit, ddouble error) = QExponentialDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 5e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (QExponentialDistribution dist in Dists) {
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
            ddouble[] expected_dist_q025theta1 = [
                1.75,
                1.534742105561412,
                1.326791075274881,
                1.126700600065561,
                0.9351430958388186,
                0.7529566491250401,
                0.581222107898076,
                0.4214004462118454,
                0.2756077296645035,
                0.1472758088888983,
                0.04340549751475546,
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

            ];
            ddouble[] expected_dist_q050theta2 = [
                0.75,
                0.703857421875,
                0.6591796875,
                0.615966796875,
                0.57421875,
                0.533935546875,
                0.4951171875,
                0.457763671875,
                0.421875,
                0.387451171875,
                0.3544921875,
                0.322998046875,
                0.29296875,
                0.264404296875,
                0.2373046875,
                0.211669921875,
                0.1875,
                0.164794921875,
                0.1435546875,
                0.123779296875,
                0.10546875,
                0.088623046875,
                0.0732421875,
                0.059326171875,
                0.046875,
                0.035888671875,
                0.0263671875,
                0.018310546875,
                0.01171875,
                0.006591796875,
                0.0029296875,
                7.32421875e-4,
                0.0
            ];
            ddouble[] expected_dist_q075theta2 = [
                0.625,
                0.5868435278534889,
                0.5504614114761353,
                0.5157977715134621,
                0.4827976226806641,
                0.4514068737626076,
                0.4215723276138306,
                0.3932416811585426,
                0.366363525390625,
                0.3408873453736305,
                0.3167635202407837,
                0.2939433231949806,
                0.2723789215087891,
                0.2520233765244484,
                0.2328306436538696,
                0.2147555723786354,
                0.19775390625,
                0.1817822828888893,
                0.1667982339859009,
                0.1527601853013039,
                0.1396274566650391,
                0.1273602619767189,
                0.1159197092056274,
                0.1052678003907204,
                0.095367431640625,
                0.08618239313364029,
                0.07767736911773682,
                0.0698179379105568,
                0.06257057189941406,
                0.0559026375412941,
                0.049782395362854,
                0.04417899996042252,
                0.0390625,
            ];
            ddouble[] expected_dist_q125theta6 = [
                0.125,
                0.122429391733529,
                0.119924525310445,
                0.1174833934386051,
                0.115104059966684,
                0.1127846570172055,
                0.1105233822489019,
                0.1083184962419518,
                0.10616832,
                0.1040712325631909,
                0.1020256687267619,
                0.1000301168600346,
                0.09808311682091927,
                0.09618325796130674,
                0.09432917721896816,
                0.09251955729181485,
                0.09075312489058508,
                0.08902864906623432,
                0.08734493960849646,
                0.085700845512266,
                0.08409525350862637,
                0.08252708665750999,
                0.08099530299912952,
                0.07949889426146767,
                0.07803688462124676,
                0.07660832951593102,
                0.0752123145044373,
                0.07384795417434481,
                0.07251439109350456,
                0.07121079480405385,
                0.06993636085693812,
                0.06869030988513558,
                0.0674718867138692,
            ];
            ddouble[] expected_dist_q150theta1 = [
                0.5,
                0.4429065743944637,
                0.3950617283950617,
                0.3545706371191136,
                0.32,
                0.290249433106576,
                0.2644628099173554,
                0.2419659735349717,
                0.2222222222222222,
                0.2048,
                0.1893491124260355,
                0.1755829903978052,
                0.163265306122449,
                0.1521997621878716,
                0.1422222222222222,
                0.1331945889698231,
                0.125,
                0.1175390266299357,
                0.1107266435986159,
                0.1044897959183674,
                0.09876543209876543,
                0.09349890430971512,
                0.0886426592797784,
                0.084155161078238,
                0.08,
                0.0761451516954194,
                0.07256235827664399,
                0.06922660897782586,
                0.06611570247933884,
                0.06320987654320988,
                0.06049149338374291,
                0.05794477138976913,
                0.05555555555555555,
            ];
            ddouble[] expected_dist_q175theta4 = [
                0.0625,
                0.06059892443782244,
                0.05879678493110772,
                0.05708638655836507,
                0.05546120007320198,
                0.05391528768875981,
                0.05244323848649584,
                0.05104011203316278,
                0.04970138902196668,
                0.04842292794422125,
                0.04720092695453278,
                0.04603189022209832,
                0.04491259816820047,
                0.04384008107951341,
                0.042811595661665,
                0.04182460416025485,
                0.04087675572932811,
                0.03996586977187009,
                0.03908992101461425,
                0.03824702611149019,
                0.03743543159731375,
                0.03665350303660951,
                0.0358997152323914,
                0.03517264337683682,
                0.03447095504051014,
                0.03379340290948635,
                0.03313881819069948,
                0.03250610461534786,
                0.03189423297844265,
                0.03130223615976842,
                0.03072920457778631,
                0.03017428203348052,
                0.02963666190593727,
            ];

            foreach ((QExponentialDistribution dist, ddouble[] expecteds) in new[]{
                (dist_q025theta1, expected_dist_q025theta1),
                (dist_q050theta2, expected_dist_q050theta2),
                (dist_q075theta2, expected_dist_q075theta2),
                (dist_q125theta6, expected_dist_q125theta6),
                (dist_q150theta1, expected_dist_q150theta1),
                (dist_q175theta4, expected_dist_q175theta4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
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
        public void CDFIntegrateTest() {
            foreach (QExponentialDistribution dist in new[]{
                dist_q025theta1,
                dist_q050theta2,
                dist_q075theta2,
                dist_q100theta3,
                dist_q125theta6,
                dist_q150theta1,
                dist_q175theta4,
            }) {
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble expected = GaussKronrodIntegral.AdaptiveIntegrate(dist.PDF, 0, x, 1e-28, 65536).value;
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-5, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}