using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class TukeyLambdaDistributionTests {
        readonly TukeyLambdaDistribution dist_lambda_m20 = new(lambda: -2.0);
        readonly TukeyLambdaDistribution dist_lambda_m15 = new(lambda: -1.5);
        readonly TukeyLambdaDistribution dist_lambda_m10 = new(lambda: -1.0);
        readonly TukeyLambdaDistribution dist_lambda_m05 = new(lambda: -0.5);
        readonly TukeyLambdaDistribution dist_lambda_z00 = new(lambda: 0.0);
        readonly TukeyLambdaDistribution dist_lambda_p05 = new(lambda: 0.5);
        readonly TukeyLambdaDistribution dist_lambda_p10 = new(lambda: 1.0);
        readonly TukeyLambdaDistribution dist_lambda_p15 = new(lambda: 1.5);
        readonly TukeyLambdaDistribution dist_lambda_p20 = new(lambda: 2.0);
        readonly TukeyLambdaDistribution dist_lambda_p25 = new(lambda: 2.5);


        TukeyLambdaDistribution[] Dists => [
            dist_lambda_m20,
            dist_lambda_m15,
            dist_lambda_m10,
            dist_lambda_m05,
            dist_lambda_z00,
            dist_lambda_p05,
            dist_lambda_p10,
            dist_lambda_p15,
            dist_lambda_p20,
            dist_lambda_p25,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Lambda={dist.Lambda}");
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
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m20.Mean));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m15.Mean));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m10.Mean));
            Assert.IsTrue(dist_lambda_m05.Mean == 0d);
            Assert.IsTrue(dist_lambda_z00.Mean == 0d);
            Assert.IsTrue(dist_lambda_p05.Mean == 0d);
            Assert.IsTrue(dist_lambda_p10.Mean == 0d);
            Assert.IsTrue(dist_lambda_p15.Mean == 0d);
            Assert.IsTrue(dist_lambda_p20.Mean == 0d);
            Assert.IsTrue(dist_lambda_p25.Mean == 0d);
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Mode)) {
                    continue;
                }

                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode - 1e-4), $"{dist}\n{dist.Mode}");
                Assert.IsTrue(dist.PDF(dist.Mode) >= dist.PDF(dist.Mode + 1e-4), $"{dist}\n{dist.Mode}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m20.Variance));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m15.Variance));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m10.Variance));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m05.Variance));
            Assert.IsTrue(ddouble.Abs(dist_lambda_z00.Variance - ddouble.PI * ddouble.PI / 3) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p05.Variance - 2 * (2 - ddouble.PI / 2)) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p10.Variance - ddouble.One / 3) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p15.Variance - "0.1567723752724348630875838184038992454681145264019") < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p20.Variance - ddouble.One / 12) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p25.Variance - "0.0484245948120992813982354530469591100767752561468") < 1e-30);
        }

        [TestMethod()]
        public void SkewnessTest() {
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m20.Skewness));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m15.Skewness));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m10.Skewness));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m05.Skewness));
            Assert.IsTrue(dist_lambda_z00.Skewness == 0d);
            Assert.IsTrue(dist_lambda_p05.Skewness == 0d);
            Assert.IsTrue(dist_lambda_p10.Skewness == 0d);
            Assert.IsTrue(dist_lambda_p15.Skewness == 0d);
            Assert.IsTrue(dist_lambda_p20.Skewness == 0d);
            Assert.IsTrue(dist_lambda_p25.Skewness == 0d);
        }

        [TestMethod()]
        public void KurtosisTest() {
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m20.Kurtosis));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m15.Kurtosis));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m10.Kurtosis));
            Assert.IsTrue(ddouble.IsNaN(dist_lambda_m05.Kurtosis));
            Assert.IsTrue(ddouble.Abs(dist_lambda_z00.Kurtosis - (ddouble)6 / 5) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p05.Kurtosis - (ddouble)("-0.91830356643745397905044652979341519486")) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p10.Kurtosis - (ddouble)(-6) / 5) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p15.Kurtosis - (ddouble)("-1.24692314260983033024627372517968272710")) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p20.Kurtosis - (ddouble)(-6) / 5) < 1e-30);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p25.Kurtosis - (ddouble)("-1.09348935542751987712395481688704467492")) < 1e-30);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.IsTrue(ddouble.Abs(dist_lambda_m20.Entropy - 5.209199576156147) < 1e-13);
            Assert.IsTrue(ddouble.Abs(dist_lambda_m15.Entropy - 4.385460155603392) < 1e-13);
            Assert.IsTrue(ddouble.Abs(dist_lambda_m10.Entropy - 3.570796326794898) < 1e-13);
            Assert.IsTrue(ddouble.Abs(dist_lambda_m05.Entropy - 2.771404238276070) < 1e-13);
            Assert.AreEqual(dist_lambda_z00.Entropy, 2.0);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p05.Entropy - 1.285398163397444) < 1e-13);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p10.Entropy - 0.693147180559945) < 1e-13);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p15.Entropy - 0.2853981633974442) < 1e-13);
            Assert.AreEqual(dist_lambda_p20.Entropy, 0.0);
            Assert.IsTrue(ddouble.Abs(dist_lambda_p25.Entropy - -0.2285957617423647) < 1e-10);

            ddouble actual = dist_lambda_p25.Entropy;
            ddouble expected = IntegrationStatistics.Entropy(dist_lambda_p25, eps: 1e-28, discontinue_eval_points: 131072);

            Console.WriteLine(actual);
            Console.WriteLine(expected);

            Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist_lambda_p25}\n{expected}\n{actual}");
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);

                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);

                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-30);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 100; i++) {
                    ddouble p = (ddouble)i / 100;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-30);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 100; i++) {
                    ddouble p = (ddouble)i / 100;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-30);
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateTest() {
            Random random = new(1234);

            foreach (TukeyLambdaDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 25; i <= 75; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.1, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (TukeyLambdaDistribution dist in Dists) {
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
            ddouble[] expected_lambda_m20 = [
                2.830805273851736900e-02,
                3.004460489242231913e-02,
                3.192501892670308178e-02,
                3.395747379621803091e-02,
                3.614764274709745401e-02,
                3.849697739436401533e-02,
                4.100034266066432398e-02,
                4.364290482036356500e-02,
                4.639628594618375179e-02,
                4.921423900955655684e-02,
                5.202851925763357643e-02,
                5.474622784848030138e-02,
                5.725053669109826660e-02,
                5.940697580176809017e-02,
                6.107675866993042546e-02,
                6.213645025641784114e-02,
                6.250000000000000000e-02,
                6.213645025641784114e-02,
                6.107675866993042546e-02,
                5.940697580176809017e-02,
                5.725053669109826660e-02,
                5.474622784848030138e-02,
                5.202851925763357643e-02,
                4.921423900955655684e-02,
                4.639628594618375179e-02,
                4.364290482036356500e-02,
                4.100034266066432398e-02,
                3.849697739436401533e-02,
                3.614764274709745401e-02,
                3.395747379621803091e-02,
                3.192501892670308178e-02,
                3.004460489242231913e-02,
                2.830805273851736900e-02,
            ];
            ddouble[] expected_lambda_m15 = [
                3.195169948908938712e-02,
                3.430047267889300977e-02,
                3.689343404654390340e-02,
                3.975589546410227687e-02,
                4.291268213134602494e-02,
                4.638578194769487489e-02,
                5.019054753407113656e-02,
                5.432986735207233353e-02,
                5.878569760065602973e-02,
                6.350758661779863645e-02,
                6.839863081709525450e-02,
                7.330106538386459147e-02,
                7.798660723699346597e-02,
                8.216001399819397588e-02,
                8.548552768769503152e-02,
                8.764070664891701612e-02,
                8.838834764831843271e-02,
                8.764070664891701612e-02,
                8.548552768769503152e-02,
                8.216001399819397588e-02,
                7.798660723699346597e-02,
                7.330106538386459147e-02,
                6.839863081709525450e-02,
                6.350758661779863645e-02,
                5.878569760065602973e-02,
                5.432986735207233353e-02,
                5.019054753407113656e-02,
                4.638578194769487489e-02,
                4.291268213134602494e-02,
                3.975589546410227687e-02,
                3.689343404654390340e-02,
                3.430047267889300977e-02,
                3.195169948908938712e-02,
            ];
            ddouble[] expected_lambda_m10 = [
                3.454915028125172399e-02,
                3.764705882352792354e-02,
                4.113151523617866734e-02,
                4.505586502586052949e-02,
                4.947775597497328787e-02,
                5.445663501817691460e-02,
                6.004879239129268242e-02,
                6.629850097186698599e-02,
                7.322330470336367580e-02,
                8.079128336101698560e-02,
                8.888888888888976658e-02,
                9.728108543674525432e-02,
                1.055728090000846581e-01,
                1.131925732105883986e-01,
                1.194299994186717034e-01,
                1.235539725813169648e-01,
                1.250000000000000000e-01,
                1.235539725813169648e-01,
                1.194299994186717034e-01,
                1.131925732105883986e-01,
                1.055728090000846581e-01,
                9.728108543674525432e-02,
                8.888888888888976658e-02,
                8.079128336101698560e-02,
                7.322330470336367580e-02,
                6.629850097186698599e-02,
                6.004879239129268242e-02,
                5.445663501817691460e-02,
                4.947775597497328787e-02,
                4.505586502586052949e-02,
                4.113151523617866734e-02,
                3.764705882352792354e-02,
                3.454915028125172399e-02,
            ];
            ddouble[] expected_lambda_m05 = [
                3.357267187562316263e-02,
                3.762868992327970952e-02,
                4.230488183820344666e-02,
                4.770422319074930917e-02,
                5.394232598274689011e-02,
                6.114440274875233144e-02,
                6.943790176747545451e-02,
                7.893789707208610784e-02,
                8.972114427017276184e-02,
                1.017839259315299616e-01,
                1.149798518690761662e-01,
                1.289395213487233016e-01,
                1.429883459748274577e-01,
                1.561033095710896157e-01,
                1.669736777654491078e-01,
                1.742233443506109591e-01,
                1.767766952966368654e-01,
                1.742233443506109591e-01,
                1.669736777654491078e-01,
                1.561033095710896157e-01,
                1.429883459748274577e-01,
                1.289395213487233016e-01,
                1.149798518690761662e-01,
                1.017839259315299616e-01,
                8.972114427017276184e-02,
                7.893789707208610784e-02,
                6.943790176747545451e-02,
                6.114440274875233144e-02,
                5.394232598274689011e-02,
                4.770422319074930917e-02,
                4.230488183820344666e-02,
                3.762868992327970952e-02,
                3.357267187562316263e-02,
            ];
            ddouble[] expected_lambda_z00 = [
                1.766270621329111071e-02,
                2.244941038204345887e-02,
                2.845302387973555613e-02,
                3.593359082532813359e-02,
                4.517665973091214426e-02,
                5.647624464487404489e-02,
                7.010371654510816342e-02,
                8.625794444256298932e-02,
                1.049935854035065202e-01,
                1.261292251866552028e-01,
                1.491464520703328356e-01,
                1.731047869924970117e-01,
                1.966119332414818510e-01,
                2.178949937618140098e-01,
                2.350037122015944946e-01,
                2.461340827375984031e-01,
                2.500000000000000000e-01,
                2.461340827375983475e-01,
                2.350037122015944946e-01,
                2.178949937618140098e-01,
                1.966119332414818510e-01,
                1.731047869924970672e-01,
                1.491464520703328633e-01,
                1.261292251866551750e-01,
                1.049935854035066174e-01,
                8.625794444256293381e-02,
                7.010371654510806627e-02,
                5.647624464487403795e-02,
                4.517665973091200549e-02,
                3.593359082532804338e-02,
                2.845302387973560124e-02,
                2.244941038204348316e-02,
                1.766270621329110377e-02,
            ];
            ddouble[] expected_lambda_p05 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.054769907118275341e-01,
                1.824501122999327984e-01,
                2.401740115525483055e-01,
                2.834733547569234791e-01,
                3.151151083534468045e-01,
                3.367599413002040909e-01,
                3.493966255825947664e-01,
                3.535533905932737309e-01,
                3.493966255825947664e-01,
                3.367599413002040909e-01,
                3.151151083534468045e-01,
                2.834733547569234791e-01,
                2.401740115525483055e-01,
                1.824501122999327984e-01,
                1.054769907118275341e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_lambda_p10 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                5.000000000000000000e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_lambda_p15 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                7.695095231009843051e-01,
                7.189432841846109534e-01,
                7.071067811865474617e-01,
                7.189432841846109534e-01,
                7.695095231009843051e-01,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_lambda_p20 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];
            ddouble[] expected_lambda_p25 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.205333313346505619e+00,
                1.414213562373094923e+00,
                1.205333313346505619e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
            ];

            foreach ((TukeyLambdaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_lambda_m20, expected_lambda_m20),
                (dist_lambda_m15, expected_lambda_m15),
                (dist_lambda_m10, expected_lambda_m10),
                (dist_lambda_m05, expected_lambda_m05),
                (dist_lambda_z00, expected_lambda_z00),
                (dist_lambda_p05, expected_lambda_p05),
                (dist_lambda_p10, expected_lambda_p10),
                (dist_lambda_p15, expected_lambda_p15),
                (dist_lambda_p20, expected_lambda_p20),
                (dist_lambda_p25, expected_lambda_p25),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.25, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else if (ddouble.IsPositiveInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-12, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_lambda_m20 = [
                3.142465129999649776e-01,
                3.215376937078744390e-01,
                3.292808074393107631e-01,
                3.375128800304665333e-01,
                3.462726986882316282e-01,
                3.555999832562903862e-01,
                3.655340499567287793e-01,
                3.761117917707892389e-01,
                3.873647869505347785e-01,
                3.993153772948190294e-01,
                4.119716699693114492e-01,
                4.253216576252967229e-01,
                4.393270520257956946e-01,
                4.539179510965638542e-01,
                4.689899435404285555e-01,
                4.844053842195776838e-01,
                5.000000000000000000e-01,
                5.155946157804223162e-01,
                5.310100564595714445e-01,
                5.460820489034361458e-01,
                5.606729479742043054e-01,
                5.746783423747032771e-01,
                5.880283300306885508e-01,
                6.006846227051809706e-01,
                6.126352130494652215e-01,
                6.238882082292107611e-01,
                6.344659500432712207e-01,
                6.444000167437096138e-01,
                6.537273013117683718e-01,
                6.624871199695334667e-01,
                6.707191925606892369e-01,
                6.784623062921255610e-01,
                6.857534870000350224e-01,
            ];
            ddouble[] expected_lambda_m15 = [
                2.594023118178725440e-01,
                2.676790057111659848e-01,
                2.765728945102026159e-01,
                2.861481821064373321e-01,
                2.964753767556871367e-01,
                3.076309003636410466e-01,
                3.196959402841130782e-01,
                3.327541107634814921e-01,
                3.468873389152449249e-01,
                3.621692598880557057e-01,
                3.786553989281813415e-01,
                3.963697246444795041e-01,
                4.152880510217258347e-01,
                4.353204726857811124e-01,
                4.562973702470571880e-01,
                4.779654696738688813e-01,
                5.000000000000000000e-01,
                5.220345303261311187e-01,
                5.437026297529428120e-01,
                5.646795273142188876e-01,
                5.847119489782741653e-01,
                6.036302753555204959e-01,
                6.213446010718186585e-01,
                6.378307401119442943e-01,
                6.531126610847550751e-01,
                6.672458892365185079e-01,
                6.803040597158869218e-01,
                6.923690996363589534e-01,
                7.035246232443128633e-01,
                7.138518178935626679e-01,
                7.234271054897973841e-01,
                7.323209942888340152e-01,
                7.405976881821274560e-01,
            ];
            ddouble[] expected_lambda_m10 = [
                1.909830056250498842e-01,
                1.999999999999957367e-01,
                2.098387322643944231e-01,
                2.206024029817754695e-01,
                2.324081207559984819e-01,
                2.453877041483920607e-01,
                2.596875762567165680e-01,
                2.754669678448706804e-01,
                2.928932188134538706e-01,
                3.121324419475186573e-01,
                3.333333333333357018e-01,
                3.566018867943441251e-01,
                3.819660112501068738e-01,
                4.093327091137481943e-01,
                4.384471871911657104e-01,
                4.688711258507183288e-01,
                5.000000000000000000e-01,
                5.311288741492816712e-01,
                5.615528128088342896e-01,
                5.906672908862518057e-01,
                6.180339887498931262e-01,
                6.433981132056558749e-01,
                6.666666666666642982e-01,
                6.878675580524813427e-01,
                7.071067811865461294e-01,
                7.245330321551293196e-01,
                7.403124237432834320e-01,
                7.546122958516079393e-01,
                7.675918792440015181e-01,
                7.793975970182245305e-01,
                7.901612677356055769e-01,
                8.000000000000042633e-01,
                8.090169943749501158e-01,
            ];
            ddouble[] expected_lambda_m05 = [
                1.069243111212827557e-01,
                1.158125457540322145e-01,
                1.257903000529623228e-01,
                1.370252132593421379e-01,
                1.497122803873836006e-01,
                1.640766988723569852e-01,
                1.803754642924744189e-01,
                1.988962921649388704e-01,
                2.199515671420542162e-01,
                2.438638817857921026e-01,
                2.709385763545739678e-01,
                3.014183121193170223e-01,
                3.354167687866222991e-01,
                3.728355198891009081e-01,
                4.132816105975720689e-01,
                4.560198031513706951e-01,
                5.000000000000000000e-01,
                5.439801968486293049e-01,
                5.867183894024279311e-01,
                6.271644801108990919e-01,
                6.645832312133777009e-01,
                6.985816878806829777e-01,
                7.290614236454260322e-01,
                7.561361182142078974e-01,
                7.800484328579457838e-01,
                8.011037078350611296e-01,
                8.196245357075255811e-01,
                8.359233011276430148e-01,
                8.502877196126163994e-01,
                8.629747867406578621e-01,
                8.742096999470376772e-01,
                8.841874542459677855e-01,
                8.930756888787172443e-01,
            ];
            ddouble[] expected_lambda_z00 = [
                1.798620996209155526e-02,
                2.297736991002561138e-02,
                2.931223075135631559e-02,
                3.732688734412946407e-02,
                4.742587317756678800e-02,
                6.008665017400761921e-02,
                7.585818002124355974e-02,
                9.534946489910950396e-02,
                1.192029220221175467e-01,
                1.480471980316895031e-01,
                1.824255238063563211e-01,
                2.227001388253088410e-01,
                2.689414213699951040e-01,
                3.208213008246069697e-01,
                3.775406687981454623e-01,
                4.378234991142018750e-01,
                5.000000000000000000e-01,
                5.621765008857980694e-01,
                6.224593312018545932e-01,
                6.791786991753929748e-01,
                7.310585786300048960e-01,
                7.772998611746910758e-01,
                8.175744761936436511e-01,
                8.519528019683105802e-01,
                8.807970779778823145e-01,
                9.046505351008905516e-01,
                9.241418199787565513e-01,
                9.399133498259923947e-01,
                9.525741268224333647e-01,
                9.626731126558706331e-01,
                9.706877692486436393e-01,
                9.770226300899743643e-01,
                9.820137900379084517e-01,
            ];
            ddouble[] expected_lambda_p05 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.392686780305751881e-02,
                5.039079468942730955e-02,
                1.035589232385873970e-01,
                1.692810861169320447e-01,
                2.443270598132230020e-01,
                3.260073636615672399e-01,
                4.119575963313693023e-01,
                5.000000000000000000e-01,
                5.880424036686306977e-01,
                6.739926363384327601e-01,
                7.556729401867769980e-01,
                8.307189138830679553e-01,
                8.964410767614126030e-01,
                9.496092053105726905e-01,
                9.860731321969424812e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,

            ];
            ddouble[] expected_lambda_p10 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.250000000000000000e-01,
                2.500000000000000000e-01,
                3.750000000000000000e-01,
                5.000000000000000000e-01,
                6.250000000000000000e-01,
                7.500000000000000000e-01,
                8.750000000000000000e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_lambda_p15 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.375128042058619826e-01,
                3.222642331250895609e-01,
                5.000000000000000000e-01,
                6.777357668749104391e-01,
                8.624871957941380174e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_lambda_p20 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                2.500000000000000000e-01,
                5.000000000000000000e-01,
                7.500000000000000000e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_lambda_p25 = [
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                0.000000000000000000e+00,
                1.655079272592203665e-01,
                5.000000000000000000e-01,
                8.344920727407796335e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];

            foreach ((TukeyLambdaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_lambda_m20, expected_lambda_m20),
                (dist_lambda_m15, expected_lambda_m15),
                (dist_lambda_m10, expected_lambda_m10),
                (dist_lambda_m05, expected_lambda_m05),
                (dist_lambda_z00, expected_lambda_z00),
                (dist_lambda_p05, expected_lambda_p05),
                (dist_lambda_p10, expected_lambda_p10),
                (dist_lambda_p15, expected_lambda_p15),
                (dist_lambda_p20, expected_lambda_p20),
                (dist_lambda_p25, expected_lambda_p25),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.25, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-12, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }
    }
}