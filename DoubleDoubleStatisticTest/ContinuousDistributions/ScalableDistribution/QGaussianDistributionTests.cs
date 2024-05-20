using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class QGaussianDistributionTests {
        readonly QGaussianDistribution dist_q025sigma1 = new(q: 0.25, sigma: 1);
        readonly QGaussianDistribution dist_q050sigma2 = new(q: 0.50, sigma: 2);
        readonly QGaussianDistribution dist_q075sigma2 = new(q: 0.75, sigma: 2);
        readonly QGaussianDistribution dist_q100sigma3 = new(q: 1.00, sigma: 3);
        readonly QGaussianDistribution dist_q125sigma6 = new(q: 1.25, sigma: 6);
        readonly QGaussianDistribution dist_q150sigma1 = new(q: 1.50, sigma: 1);
        readonly QGaussianDistribution dist_q175sigma4 = new(q: 1.75, sigma: 4);
        readonly QGaussianDistribution dist_q200sigma4 = new(q: 2.00, sigma: 4);
        readonly QGaussianDistribution dist_q225sigma5 = new(q: 2.25, sigma: 5);
        readonly QGaussianDistribution dist_q250sigma2 = new(q: 2.50, sigma: 2);

        QGaussianDistribution[] Dists => [
            dist_q025sigma1,
            dist_q050sigma2,
            dist_q075sigma2,
            dist_q100sigma3,
            dist_q125sigma6,
            dist_q150sigma1,
            dist_q175sigma4,
            dist_q200sigma4,
            dist_q225sigma5,
            dist_q250sigma2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (QGaussianDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Q={dist.Q}");
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
            foreach (QGaussianDistribution dist in Dists) {
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
            foreach (QGaussianDistribution dist in Dists) {
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
            foreach (QGaussianDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (QGaussianDistribution dist in Dists) {
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
            foreach (QGaussianDistribution dist in Dists) {
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
            foreach (QGaussianDistribution dist in Dists) {
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
            foreach (QGaussianDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (QGaussianDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (QGaussianDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (QGaussianDistribution dist in Dists) {
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
            foreach (QGaussianDistribution dist in Dists) {
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
            foreach (QGaussianDistribution dist in Dists) {
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

            foreach (QGaussianDistribution dist in Dists) {
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

            foreach (QGaussianDistribution dist in Dists) {
                // ignore
                if (dist.Q > 1.5) {
                    continue;
                }

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                (QGaussianDistribution? dist_fit, ddouble error) = QGaussianDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 5e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (QGaussianDistribution dist in Dists) {
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
            ddouble[] expected_dist_q025sigma1 = [
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
                0.01330818551100658,
                0.111455372677331,
                0.2266886139112359,
                0.3408343286901829,
                0.4457579966184176,
                0.5365524179591407,
                0.6099698273970459,
                0.6638156364533219,
                0.6966613835852355,
                0.7076975033420114,
                0.6966613835852355,
                0.6638156364533219,
                0.6099698273970459,
                0.5365524179591407,
                0.4457579966184176,
                0.3408343286901829,
                0.2266886139112359,
                0.111455372677331,
                0.01330818551100658,
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
            ];
            ddouble[] expected_dist_q050sigma2 = [
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
                9.912938769762276e-4,
                0.006373867900302504,
                0.01586070203161964,
                0.02882970878690684,
                0.04468914641307116,
                0.06287761888794748,
                0.08286407592029854,
                0.1041478129498149,
                0.126258471147115,
                0.1487560374137452,
                0.1712308443821794,
                0.1933035704158197,
                0.2146252396089959,
                0.2348772217869656,
                0.2537712325059143,
                0.2710493330529553,
                0.2864839304461298,
                0.2998777774344068,
                0.3110639724976832,
                0.3199059598467837,
                0.3262975294234607,
                0.3301628169003948,
                0.3314563036811942,
                0.3301628169003948,
                0.3262975294234607,
                0.3199059598467837,
                0.3110639724976832,
                0.2998777774344068,
                0.2864839304461298,
                0.2710493330529553,
                0.2537712325059143,
                0.2348772217869656,
                0.2146252396089959,
                0.1933035704158197,
                0.1712308443821794,
                0.1487560374137452,
                0.126258471147115,
                0.1041478129498149,
                0.08286407592029854,
                0.06287761888794748,
                0.04468914641307116,
                0.02882970878690684,
                0.01586070203161964,
                0.006373867900302504,
                9.912938769762276e-4,
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
            ddouble[] expected_dist_q075sigma2 = [
                0.0,
                4.407303602071977e-6,
                6.614507469748788e-5,
                3.137726817170617e-4,
                9.282302926294508e-4,
                0.002118843275833449,
                0.004103247810007816,
                0.007090916206489163,
                0.01126997172832489,
                0.01679699397072465,
                0.02378952714821026,
                0.03232101491434048,
                0.04241789662046359,
                0.05405861120052703,
                0.06717426614955004,
                0.0816507403439415,
                0.09733200073242189,
                0.1140244242068844,
                0.1315019272431073,
                0.149511717181805,
                0.1677804903010838,
                0.1860209131119417,
                0.2039382345890317,
                0.2212368883294813,
                0.2376269549131394,
                0.2528303660181975,
                0.2665867431267089,
                0.2786587749351055,
                0.2888370488653891,
                0.2969442633532494,
                0.3028387588699389,
                0.3064173169153107,
                0.3076171875000001,
                0.3064173169153107,
                0.3028387588699389,
                0.2969442633532494,
                0.2888370488653891,
                0.2786587749351055,
                0.2665867431267089,
                0.2528303660181975,
                0.2376269549131394,
                0.2212368883294813,
                0.2039382345890317,
                0.1860209131119417,
                0.1677804903010838,
                0.149511717181805,
                0.1315019272431073,
                0.1140244242068844,
                0.09733200073242189,
                0.0816507403439415,
                0.06717426614955004,
                0.05405861120052703,
                0.04241789662046359,
                0.03232101491434048,
                0.02378952714821026,
                0.01679699397072465,
                0.01126997172832489,
                0.007090916206489163,
                0.004103247810007816,
                0.002118843275833449,
                9.282302926294508e-4,
                3.137726817170617e-4,
                6.614507469748788e-5,
                4.407303602071977e-6,
                0.0,
            ];
            ddouble[] expected_dist_q125sigma6 = [
                0.05569149768671602,
                0.05708337281510048,
                0.05847273947407457,
                0.05985699682075867,
                0.06123344927480449,
                0.06259931401009491,
                0.06395172937866175,
                0.0652877642539597,
                0.06660442827054563,
                0.06789868292672101,
                0.06916745350592474,
                0.07040764176176334,
                0.07161613930068153,
                0.07278984158558219,
                0.07392566247337204,
                0.07502054918961479,
                0.07607149763440445,
                0.07707556790540022,
                0.0780298999168686,
                0.07893172898771883,
                0.07977840126704616,
                0.08056738886274209,
                0.08129630453740296,
                0.08196291583615428,
                0.08256515851315949,
                0.0831011491275313,
                0.08356919668510018,
                0.08396781320998067,
                0.0842957231390396,
                0.08455187144310652,
                0.0847354303909363,
                0.08484580488537247,
                0.08488263631567752,
                0.08484580488537247,
                0.0847354303909363,
                0.08455187144310652,
                0.0842957231390396,
                0.08396781320998067,
                0.08356919668510018,
                0.0831011491275313,
                0.08256515851315949,
                0.08196291583615428,
                0.08129630453740296,
                0.08056738886274209,
                0.07977840126704616,
                0.07893172898771883,
                0.0780298999168686,
                0.07707556790540022,
                0.07607149763440445,
                0.07502054918961479,
                0.07392566247337204,
                0.07278984158558219,
                0.07161613930068153,
                0.07040764176176334,
                0.06916745350592474,
                0.06789868292672101,
                0.06660442827054563,
                0.0652877642539597,
                0.06395172937866175,
                0.06259931401009491,
                0.06123344927480449,
                0.05985699682075867,
                0.05847273947407457,
                0.05708337281510048,
                0.05569149768671602,
            ];
            ddouble[] expected_dist_q150sigma1 = [
                0.00555750812442658,
                0.006219125272222191,
                0.006979090582331879,
                0.007854843025385519,
                0.008867381384126622,
                0.01004207407452255,
                0.01140966693577977,
                0.01300753826122515,
                0.0148812614240844,
                0.01708654777830833,
                0.01969165508447342,
                0.02278035730665217,
                0.02645557586503893,
                0.03084376220390101,
                0.03610008253367047,
                0.04241436354202401,
                0.05001757311983922,
                0.05918827100738319,
                0.07025788048657801,
                0.08361268421543168,
                0.09968900386514669,
                0.1189560049347432,
                0.1418781021460259,
                0.1688466670167582,
                0.2000702924793569,
                0.2354173852328199,
                0.274218889870576,
                0.3150664813515747,
                0.3556805199633011,
                0.3929560052191919,
                0.4232892138406228,
                0.4432060129775262,
                0.450158158078553,
                0.4432060129775262,
                0.4232892138406228,
                0.3929560052191919,
                0.3556805199633011,
                0.3150664813515747,
                0.274218889870576,
                0.2354173852328199,
                0.2000702924793569,
                0.1688466670167582,
                0.1418781021460259,
                0.1189560049347432,
                0.09968900386514669,
                0.08361268421543168,
                0.07025788048657801,
                0.05918827100738319,
                0.05001757311983922,
                0.04241436354202401,
                0.03610008253367047,
                0.03084376220390101,
                0.02645557586503893,
                0.02278035730665217,
                0.01969165508447342,
                0.01708654777830833,
                0.0148812614240844,
                0.01300753826122515,
                0.01140966693577977,
                0.01004207407452255,
                0.008867381384126622,
                0.007854843025385519,
                0.006979090582331879,
                0.006219125272222191,
                0.00555750812442658,
            ];
            ddouble[] expected_dist_q175sigma4 = [
                0.04582204383298002,
                0.04748402455169609,
                0.04919647186553652,
                0.05095834059157775,
                0.05276813522557096,
                0.05462385823907096,
                0.0565229572647805,
                0.05846227216144162,
                0.06043798319272065,
                0.06244556181320611,
                0.06447972582010235,
                0.06653440088931949,
                0.06860269075365882,
                0.07067685847923943,
                0.07274832143143277,
                0.07480766256802589,
                0.07684466062854069,
                0.07884834157870137,
                0.0808070532952982,
                0.08270856492275287,
                0.0845401915918099,
                0.08628894426923206,
                0.08794170342737033,
                0.08948541402408519,
                0.09090729802507806,
                0.09219507945729057,
                0.09333721584141948,
                0.09432312890796811,
                0.09514342684684568,
                0.09579011005594017,
                0.09625675249841094,
                0.09653865138049042,
                0.0966329389136607,
                0.09653865138049042,
                0.09625675249841094,
                0.09579011005594017,
                0.09514342684684568,
                0.09432312890796811,
                0.09333721584141948,
                0.09219507945729057,
                0.09090729802507806,
                0.08948541402408519,
                0.08794170342737033,
                0.08628894426923206,
                0.0845401915918099,
                0.08270856492275287,
                0.0808070532952982,
                0.07884834157870137,
                0.07684466062854069,
                0.07480766256802589,
                0.07274832143143277,
                0.07067685847923943,
                0.06860269075365882,
                0.06653440088931949,
                0.06447972582010235,
                0.06244556181320611,
                0.06043798319272065,
                0.05846227216144162,
                0.0565229572647805,
                0.05462385823907096,
                0.05276813522557096,
                0.05095834059157775,
                0.04919647186553652,
                0.04748402455169609,
                0.04582204383298002,
            ];
            ddouble[] expected_dist_q200sigma4 = [
                0.03978873577297384,
                0.0410515520720657,
                0.04235308256915302,
                0.04369293880056323,
                0.04507042636230665,
                0.04648450134800366,
                0.04793372403708848,
                0.0494162103475139,
                0.05092958178940651,
                0.05247091491503568,
                0.05403669155374696,
                0.05562275144235523,
                0.05722424920158035,
                0.05883561795166094,
                0.06045054218327182,
                0.06206194277460047,
                0.06366197723675814,
                0.06524205833711001,
                0.06679289415004133,
                0.06830455227414117,
                0.06976655039644727,
                0.07116797455288246,
                0.07249762532299858,
                0.07374419082628997,
                0.07489644380795075,
                0.0759434583998606,
                0.07687484043684002,
                0.07768096364447132,
                0.0783532027529331,
                0.07888415378804493,
                0.07926783157884282,
                0.07949983498834187,
                0.07957747154594767,
                0.07949983498834187,
                0.07926783157884282,
                0.07888415378804493,
                0.0783532027529331,
                0.07768096364447132,
                0.07687484043684002,
                0.0759434583998606,
                0.07489644380795075,
                0.07374419082628997,
                0.07249762532299858,
                0.07116797455288246,
                0.06976655039644727,
                0.06830455227414117,
                0.06679289415004133,
                0.06524205833711001,
                0.06366197723675814,
                0.06206194277460047,
                0.06045054218327182,
                0.05883561795166094,
                0.05722424920158035,
                0.05562275144235523,
                0.05403669155374696,
                0.05247091491503568,
                0.05092958178940651,
                0.0494162103475139,
                0.04793372403708848,
                0.04648450134800366,
                0.04507042636230665,
                0.04369293880056323,
                0.04235308256915302,
                0.0410515520720657,
                0.03978873577297384,
            ];
            ddouble[] expected_dist_q225sigma5 = [
                0.0306783302693022,
                0.03136636653393082,
                0.03206656965754934,
                0.03277820473529708,
                0.03350039624975531,
                0.03423211630090735,
                0.03497217299287687,
                0.03571919925002438,
                0.036471642386684,
                0.03722775480802028,
                0.03798558627202481,
                0.0387429781917318,
                0.03949756049877878,
                0.04024675162026111,
                0.0409877621356239,
                0.0417176026738935,
                0.04243309657855939,
                0.04313089780286366,
                0.04380751439796643,
                0.04445933781771517,
                0.04508267808597841,
                0.04567380465793242,
                0.0462289925608997,
                0.0467445731326243,
                0.04721698839826913,
                0.04764284785832155,
                0.04801898621687054,
                0.04834252038338201,
                0.04861090395053418,
                0.04882197730261885,
                0.04897401155551338,
                0.04906574467581735,
                0.04909640837119979,
                0.04906574467581735,
                0.04897401155551338,
                0.04882197730261885,
                0.04861090395053418,
                0.04834252038338201,
                0.04801898621687054,
                0.04764284785832155,
                0.04721698839826913,
                0.0467445731326243,
                0.0462289925608997,
                0.04567380465793242,
                0.04508267808597841,
                0.04445933781771517,
                0.04380751439796643,
                0.04313089780286366,
                0.04243309657855939,
                0.0417176026738935,
                0.0409877621356239,
                0.04024675162026111,
                0.03949756049877878,
                0.0387429781917318,
                0.03798558627202481,
                0.03722775480802028,
                0.036471642386684,
                0.03571919925002438,
                0.03497217299287687,
                0.03423211630090735,
                0.03350039624975531,
                0.03277820473529708,
                0.03206656965754934,
                0.03136636653393082,
                0.0306783302693022,
            ];
            ddouble[] expected_dist_q250sigma2 = [
                0.02296839445808278,
                0.02381311008640864,
                0.0247092149659419,
                0.0256607990168721,
                0.02667230997271809,
                0.02774857694552555,
                0.02889483115359334,
                0.03011672149063021,
                0.03142032164568311,
                0.0328121241704548,
                0.03429901513076639,
                0.03588822064764702,
                0.03758721359178691,
                0.03940356480523261,
                0.04134471838589381,
                0.04341766479053071,
                0.04562847902985673,
                0.0479816847271291,
                0.05047939974033257,
                0.05312021809331236,
                0.05589779061791582,
                0.0587990898126919,
                0.06180239222598165,
                0.06487509463075287,
                0.06797160620784935,
                0.07103172512436691,
                0.0739800883728714,
                0.07672741579544806,
                0.07917424927172735,
                0.08121759566230725,
                0.0827602474987702,
                0.08372165780861819,
                0.08404837699052227,
                0.08372165780861819,
                0.0827602474987702,
                0.08121759566230725,
                0.07917424927172735,
                0.07672741579544806,
                0.0739800883728714,
                0.07103172512436691,
                0.06797160620784935,
                0.06487509463075287,
                0.06180239222598165,
                0.0587990898126919,
                0.05589779061791582,
                0.05312021809331236,
                0.05047939974033257,
                0.0479816847271291,
                0.04562847902985673,
                0.04341766479053071,
                0.04134471838589381,
                0.03940356480523261,
                0.03758721359178691,
                0.03588822064764702,
                0.03429901513076639,
                0.0328121241704548,
                0.03142032164568311,
                0.03011672149063021,
                0.02889483115359334,
                0.02774857694552555,
                0.02667230997271809,
                0.0256607990168721,
                0.0247092149659419,
                0.02381311008640864,
                0.02296839445808278,
            ];

            foreach ((QGaussianDistribution dist, ddouble[] expecteds) in new[]{
                (dist_q025sigma1, expected_dist_q025sigma1),
                (dist_q050sigma2, expected_dist_q050sigma2),
                (dist_q075sigma2, expected_dist_q075sigma2),
                (dist_q125sigma6, expected_dist_q125sigma6),
                (dist_q150sigma1, expected_dist_q150sigma1),
                (dist_q175sigma4, expected_dist_q175sigma4),
                (dist_q200sigma4, expected_dist_q200sigma4),
                (dist_q225sigma5, expected_dist_q225sigma5),
                (dist_q250sigma2, expected_dist_q250sigma2),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.125, i++) {
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
            foreach (QGaussianDistribution dist in new[]{
                dist_q025sigma1,
                dist_q050sigma2,
                dist_q075sigma2,
                dist_q100sigma3,
                dist_q125sigma6,
                dist_q150sigma1,
                dist_q175sigma4,
                dist_q200sigma4,
                dist_q225sigma5,
            }) {
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble expected = GaussKronrodIntegral.AdaptiveIntegrate(dist.PDF, ddouble.NegativeInfinity, x, 1e-28, 65536).value;
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