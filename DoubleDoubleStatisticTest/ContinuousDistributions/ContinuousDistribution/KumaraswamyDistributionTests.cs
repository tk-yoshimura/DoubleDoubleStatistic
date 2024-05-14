using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class KumaraswamyDistributionTests {
        readonly KumaraswamyDistribution dist_a1b1 = new(alpha: 1, beta: 1);
        readonly KumaraswamyDistribution dist_a2b1 = new(alpha: 2, beta: 1);
        readonly KumaraswamyDistribution dist_a1b2 = new(alpha: 1, beta: 2);
        readonly KumaraswamyDistribution dist_a2b2 = new(alpha: 2, beta: 2);
        readonly KumaraswamyDistribution dist_a3b4 = new(alpha: 3, beta: 4);

        KumaraswamyDistribution[] Dists => [
            dist_a1b1,
            dist_a2b1,
            dist_a1b2,
            dist_a2b2,
            dist_a3b4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (KumaraswamyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (KumaraswamyDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
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
            foreach (KumaraswamyDistribution dist in Dists) {
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

            foreach (KumaraswamyDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.01, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (KumaraswamyDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                (KumaraswamyDistribution? dist_fit, ddouble error) = KumaraswamyDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (KumaraswamyDistribution dist in Dists) {
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
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ];
            ddouble[] expected_dist_a2b1 = [
                0.0,
                0.03125,
                0.0625,
                0.09375,
                0.125,
                0.15625,
                0.1875,
                0.21875,
                0.25,
                0.28125,
                0.3125,
                0.34375,
                0.375,
                0.40625,
                0.4375,
                0.46875,
                0.5,
                0.53125,
                0.5625,
                0.59375,
                0.625,
                0.65625,
                0.6875,
                0.71875,
                0.75,
                0.78125,
                0.8125,
                0.84375,
                0.875,
                0.90625,
                0.9375,
                0.96875,
                1.0,
                1.03125,
                1.0625,
                1.09375,
                1.125,
                1.15625,
                1.1875,
                1.21875,
                1.25,
                1.28125,
                1.3125,
                1.34375,
                1.375,
                1.40625,
                1.4375,
                1.46875,
                1.5,
                1.53125,
                1.5625,
                1.59375,
                1.625,
                1.65625,
                1.6875,
                1.71875,
                1.75,
                1.78125,
                1.8125,
                1.84375,
                1.875,
                1.90625,
                1.9375,
                1.96875,
                2.0,
            ];
            ddouble[] expected_dist_a1b2 = [
                2.0,
                1.96875,
                1.9375,
                1.90625,
                1.875,
                1.84375,
                1.8125,
                1.78125,
                1.75,
                1.71875,
                1.6875,
                1.65625,
                1.625,
                1.59375,
                1.5625,
                1.53125,
                1.5,
                1.46875,
                1.4375,
                1.40625,
                1.375,
                1.34375,
                1.3125,
                1.28125,
                1.25,
                1.21875,
                1.1875,
                1.15625,
                1.125,
                1.09375,
                1.0625,
                1.03125,
                1.0,
                0.96875,
                0.9375,
                0.90625,
                0.875,
                0.84375,
                0.8125,
                0.78125,
                0.75,
                0.71875,
                0.6875,
                0.65625,
                0.625,
                0.59375,
                0.5625,
                0.53125,
                0.5,
                0.46875,
                0.4375,
                0.40625,
                0.375,
                0.34375,
                0.3125,
                0.28125,
                0.25,
                0.21875,
                0.1875,
                0.15625,
                0.125,
                0.09375,
                0.0625,
                0.03125,
                0.0,
            ];
            ddouble[] expected_dist_a2b2 = [
                0.0,
                0.0624847412109375,
                0.1248779296875,
                0.1870880126953125,
                0.2490234375,
                0.3105926513671875,
                0.3717041015625,
                0.4322662353515625,
                0.4921875,
                0.5513763427734375,
                0.6097412109375,
                0.6671905517578125,
                0.7236328125,
                0.7789764404296875,
                0.8331298828125,
                0.8860015869140625,
                0.9375,
                0.9875335693359375,
                1.0360107421875,
                1.082839965820313,
                1.1279296875,
                1.171188354492188,
                1.2125244140625,
                1.251846313476563,
                1.2890625,
                1.324081420898438,
                1.3568115234375,
                1.387161254882813,
                1.4150390625,
                1.440353393554688,
                1.4630126953125,
                1.482925415039063,
                1.5,
                1.514144897460938,
                1.5252685546875,
                1.533279418945313,
                1.5380859375,
                1.539596557617188,
                1.5377197265625,
                1.532363891601563,
                1.5234375,
                1.510848999023438,
                1.4945068359375,
                1.474319458007813,
                1.4501953125,
                1.422042846679688,
                1.3897705078125,
                1.353286743164063,
                1.3125,
                1.267318725585938,
                1.2176513671875,
                1.163406372070313,
                1.1044921875,
                1.040817260742188,
                0.9722900390625,
                0.8988189697265625,
                0.8203125,
                0.7366790771484375,
                0.6478271484375,
                0.5536651611328125,
                0.4541015625,
                0.3490447998046875,
                0.2384033203125,
                0.1220855712890625,
                0.0,
            ];
            ddouble[] expected_dist_a3b4 = [
                0.0,
                0.002929653972515212,
                0.01171767714913552,
                0.02635904112922521,
                0.04684067610583043,
                0.07313746366245094,
                0.1052082540437612,
                0.1429919258967557,
                0.1864035115577281,
                0.235330415971266,
                0.2896287622375954,
                0.3491199015551849,
                0.4135871299110931,
                0.4827726582201555,
                0.5563748866673629,
                0.6340460377007406,
                0.7153902053833008,
                0.7999618815643089,
                0.8872650214877839,
                0.9767527129289524,
                1.067827513639941,
                1.159842521690428,
                1.252103242096951,
                1.34387031082921,
                1.43436313373968,
                1.522764493054351,
                1.608226167651036,
                1.689875605294744,
                1.766823675147634,
                1.838173517068435,
                1.903030490300436,
                1.960513207949643,
                2.009765625,
                2.049970126319103,
                2.08036153698711,
                2.100241950139106,
                2.108996237145448,
                2.106108071158178,
                2.091176257609973,
                2.063931123944725,
                2.024250675458461,
                1.972176174402235,
                1.907926745204601,
                1.831912549563741,
                1.744746010984954,
                1.647250498838497,
                1.540465806919376,
                1.425649680532176,
                1.304274559020996,
                1.178018607131651,
                1.048750008338164,
                0.9185033859897396,
                0.789447103532666,
                0.6638400728224534,
                0.5439765693467241,
                0.4321174147045366,
                0.3304057396017015,
                0.2407653845867826,
                0.1647798304247435,
                0.1035493750340926,
                0.05752408894181826,
                0.02630988587488794,
                0.008444839037562324,
                0.001142656443741005,
                0.0,
            ];

            foreach ((KumaraswamyDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a1b2, expected_dist_a1b2), (dist_a2b2, expected_dist_a2b2),
                (dist_a3b4, expected_dist_a3b4),
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
            ddouble[] expected_dist_a1b1 = [
                0.0,
                0.015625,
                0.03125,
                0.046875,
                0.0625,
                0.078125,
                0.09375,
                0.109375,
                0.125,
                0.140625,
                0.15625,
                0.171875,
                0.1875,
                0.203125,
                0.21875,
                0.234375,
                0.25,
                0.265625,
                0.28125,
                0.296875,
                0.3125,
                0.328125,
                0.34375,
                0.359375,
                0.375,
                0.390625,
                0.40625,
                0.421875,
                0.4375,
                0.453125,
                0.46875,
                0.484375,
                0.5,
                0.515625,
                0.53125,
                0.546875,
                0.5625,
                0.578125,
                0.59375,
                0.609375,
                0.625,
                0.640625,
                0.65625,
                0.671875,
                0.6875,
                0.703125,
                0.71875,
                0.734375,
                0.75,
                0.765625,
                0.78125,
                0.796875,
                0.8125,
                0.828125,
                0.84375,
                0.859375,
                0.875,
                0.890625,
                0.90625,
                0.921875,
                0.9375,
                0.953125,
                0.96875,
                0.984375,
                1.0,
            ];
            ddouble[] expected_dist_a2b1 = [
                0.0,
                2.44140625e-4,
                9.765625e-4,
                0.002197265625,
                0.00390625,
                0.006103515625,
                0.0087890625,
                0.011962890625,
                0.015625,
                0.019775390625,
                0.0244140625,
                0.029541015625,
                0.03515625,
                0.041259765625,
                0.0478515625,
                0.054931640625,
                0.0625,
                0.070556640625,
                0.0791015625,
                0.088134765625,
                0.09765625,
                0.107666015625,
                0.1181640625,
                0.129150390625,
                0.140625,
                0.152587890625,
                0.1650390625,
                0.177978515625,
                0.19140625,
                0.205322265625,
                0.2197265625,
                0.234619140625,
                0.25,
                0.265869140625,
                0.2822265625,
                0.299072265625,
                0.31640625,
                0.334228515625,
                0.3525390625,
                0.371337890625,
                0.390625,
                0.410400390625,
                0.4306640625,
                0.451416015625,
                0.47265625,
                0.494384765625,
                0.5166015625,
                0.539306640625,
                0.5625,
                0.586181640625,
                0.6103515625,
                0.635009765625,
                0.66015625,
                0.685791015625,
                0.7119140625,
                0.738525390625,
                0.765625,
                0.793212890625,
                0.8212890625,
                0.849853515625,
                0.87890625,
                0.908447265625,
                0.9384765625,
                0.968994140625,
                1.0,
            ];
            ddouble[] expected_dist_a1b2 = [
                0.0,
                0.031005859375,
                0.0615234375,
                0.091552734375,
                0.12109375,
                0.150146484375,
                0.1787109375,
                0.206787109375,
                0.234375,
                0.261474609375,
                0.2880859375,
                0.314208984375,
                0.33984375,
                0.364990234375,
                0.3896484375,
                0.413818359375,
                0.4375,
                0.460693359375,
                0.4833984375,
                0.505615234375,
                0.52734375,
                0.548583984375,
                0.5693359375,
                0.589599609375,
                0.609375,
                0.628662109375,
                0.6474609375,
                0.665771484375,
                0.68359375,
                0.700927734375,
                0.7177734375,
                0.734130859375,
                0.75,
                0.765380859375,
                0.7802734375,
                0.794677734375,
                0.80859375,
                0.822021484375,
                0.8349609375,
                0.847412109375,
                0.859375,
                0.870849609375,
                0.8818359375,
                0.892333984375,
                0.90234375,
                0.911865234375,
                0.9208984375,
                0.929443359375,
                0.9375,
                0.945068359375,
                0.9521484375,
                0.958740234375,
                0.96484375,
                0.970458984375,
                0.9755859375,
                0.980224609375,
                0.984375,
                0.988037109375,
                0.9912109375,
                0.993896484375,
                0.99609375,
                0.997802734375,
                0.9990234375,
                0.999755859375,
                1.0,
            ];
            ddouble[] expected_dist_a2b2 = [
                0.0,
                4.882216453552246e-4,
                0.001952171325683594,
                0.004389703273773193,
                0.0077972412109375,
                0.01216977834701538,
                0.01750087738037109,
                0.02378267049789429,
                0.031005859375,
                0.03915971517562866,
                0.04823207855224609,
                0.05820935964584351,
                0.0690765380859375,
                0.08081716299057007,
                0.0934133529663086,
                0.1068457961082458,
                0.12109375,
                0.1361350417137146,
                0.1519460678100586,
                0.1685017943382263,
                0.1857757568359375,
                0.2037400603294373,
                0.2223653793334961,
                0.2416209578514099,
                0.261474609375,
                0.281892716884613,
                0.3028402328491211,
                0.3242806792259216,
                0.3461761474609375,
                0.3684872984886169,
                0.3911733627319336,
                0.4141921401023865,
                0.4375,
                0.461051881313324,
                0.4848012924194336,
                0.5087003111839294,
                0.5326995849609375,
                0.5567483305931091,
                0.5807943344116211,
                0.6047839522361755,
                0.628662109375,
                0.6523723006248474,
                0.6758565902709961,
                0.6990556120872498,
                0.7219085693359375,
                0.7443532347679138,
                0.7663259506225586,
                0.7877616286277771,
                0.80859375,
                0.8287543654441833,
                0.8481740951538086,
                0.8667821288108826,
                0.8845062255859375,
                0.901272714138031,
                0.9170064926147461,
                0.9316310286521912,
                0.945068359375,
                0.9572390913963318,
                0.9680624008178711,
                0.9774560332298279,
                0.9853363037109375,
                0.9916180968284607,
                0.9962148666381836,
                0.9990386366844177,
                1.0,
            ];
            ddouble[] expected_dist_a3b4 = [
                0.0,
                1.525870175123067e-5,
                1.220647246782391e-4,
                4.119236589806841e-4,
                9.762049303354559e-4,
                0.001905984824388973,
                0.003291827069795694,
                0.005223501496159022,
                0.007789641604176722,
                0.01107734228509982,
                0.01517169940397645,
                0.02015529351419965,
                0.02610762059708449,
                0.03310447341627565,
                0.04121727783858198,
                0.05051238929523749,
                0.06105035543441772,
                0.07288515193974932,
                0.08606339945193198,
                0.1006235705214552,
                0.1165951965282908,
                0.1339980855163442,
                0.1528415628916612,
                0.173123747907428,
                0.1948308797873324,
                0.2179367082015302,
                0.2424019635838907,
                0.2681739234407705,
                0.2951860913234405,
                0.323358005489178,
                0.3525951944282323,
                0.3827892963510189,
                0.413818359375,
                0.4455473384839574,
                0.4778288043110586,
                0.5105038773755828,
                0.5434033995326182,
                0.5763493520234944,
                0.6091565265869057,
                0.6416344525479324,
                0.6735895785823232,
                0.7048277028906691,
                0.7351566397420366,
                0.764389103685922,
                0.7923457851078943,
                0.8188585821367922,
                0.8437739441145534,
                0.8669562708241765,
                0.8882912993431091,
                0.907689396650328,
                0.9250886608627452,
                0.9404577171030126,
                0.9537980753942144,
                0.9651458975204418,
                0.9745729973640336,
                0.9821868747034957,
                0.9881295556988334,
                0.9925749841660665,
                0.9957246761075328,
                0.9978013156713041,
                0.9990399336081701,
                0.9996762692191073,
                0.9999318735751203,
                0.999995465274169,
                1.0,
            ];

            foreach ((KumaraswamyDistribution dist, ddouble[] expecteds) in new[]{
                (dist_a1b1, expected_dist_a1b1), (dist_a2b1, expected_dist_a2b1),
                (dist_a1b2, expected_dist_a1b2), (dist_a2b2, expected_dist_a2b2),
                (dist_a3b4, expected_dist_a3b4),
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