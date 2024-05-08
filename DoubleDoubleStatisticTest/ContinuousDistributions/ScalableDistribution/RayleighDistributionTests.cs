using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class RayleighDistributionTests {
        readonly RayleighDistribution dist1 = new(sigma: 1);
        readonly RayleighDistribution dist2 = new(sigma: 2);

        RayleighDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (RayleighDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
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
            foreach (RayleighDistribution dist in Dists) {
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
            foreach (RayleighDistribution dist in Dists) {
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
            foreach (RayleighDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (RayleighDistribution dist in Dists) {
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
            foreach (RayleighDistribution dist in Dists) {
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
            foreach (RayleighDistribution dist in Dists) {
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
            foreach (RayleighDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (RayleighDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (RayleighDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (RayleighDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (RayleighDistribution dist in Dists) {
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
            foreach (RayleighDistribution dist in Dists) {
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

            foreach (RayleighDistribution dist in Dists) {

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
            foreach (RayleighDistribution dist in Dists) {
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
            ddouble[] expected_dist1 = [
                0.0,
                0.06237804881921722,
                0.1240272422825304,
                0.1842329004296063,
                0.242308308619086,
                0.2976077499672426,
                0.3495384346348228,
                0.3975710184079682,
                0.4412484512922977,
                0.4801929532566456,
                0.5141109764991654,
                0.5427960791417917,
                0.5661297014917555,
                0.5840799003887053,
                0.5966981572915546,
                0.6041144295236206,
                0.6065306597126334,
                0.604212994576089,
                0.5974829899147633,
                0.5867080935688205,
                0.5722917022145179,
                0.5546630817304353,
                0.5342674253295007,
                0.511556299954324,
                0.4869787010375246,
                0.4609729002725669,
                0.433959232242808,
                0.4063339253344194,
                0.3784640419523028,
                0.3506835541627722,
                0.3232905448007866,
                0.2965454918641939,
                0.2706705664732254,
                0.2458498523304601,
                0.2222303777391281,
                0.1999238398413533,
                0.1790088946160123,
                0.1595338849336161,
                0.1415198820597169,
                0.1249639227784783,
                0.1098423340585186,
                0.09611404916225266,
                0.08372383257566239,
                0.07260534541571251,
                0.06268399742993394,
                0.0538795457919257,
                0.04610841416663283,
                0.03928571761875986,
                0.03332698961472692,
                0.02814961646638679,
                0.0236739920133123,
                0.01982441114569532,
                0.01652972500079128,
                0.01372378344108214,
                0.01134569189677137,
                0.009339910007049819,
                0.007656218913640098,
                0.00624958273785796,
                0.005079927893459547,
                0.004111861623228684,
                0.003314348651006438,
                0.002660362245168146,
                0.002126523404549189,
                0.001692739391146265,
                0.001341850511610047,
            ];
            ddouble[] expected_dist2 = [
                0.0,
                0.01561737246781077,
                0.03118902440960861,
                0.04666945830813365,
                0.06201362114126522,
                0.07717712283703128,
                0.09211645021480316,
                0.1067891749697384,
                0.121154154309543,
                0.1351717229165845,
                0.1488038749836213,
                0.1620144351571199,
                0.1747692173174114,
                0.1870361702287891,
                0.1987855092039841,
                0.2099898330451332,
                0.2206242256461489,
                0.2306663417680795,
                0.2400964766283228,
                0.2488976190751457,
                0.2570554882495827,
                0.2645585537661701,
                0.2713980395708959,
                0.2775679117580142,
                0.2830648507458777,
                0.2878882083246136,
                0.2920399501943526,
                0.2955245847109348,
                0.2983490786457773,
                0.3005227608472417,
                0.3020572147618103,
                0.3029661608342442,
                0.3032653298563167,
                0.3029723283735011,
                0.3021064972880445,
                0.3006887648152236,
                0.2987414949573817,
                0.2962883326578508,
                0.2933540467844102,
                0.2899643720699744,
                0.2861458511072589,
                0.281925677454867,
                0.2773315408652177,
                0.2723914755907616,
                0.2671337126647503,
                0.2615865369872712,
                0.255778149977162,
                0.2497365384766278,
                0.2434893505187623,
                0.2370637784895639,
                0.2304864501362835,
                0.2237833277938419,
                0.216979616121404,
                0.210099678562716,
                0.2031669626672097,
                0.1962039343347811,
                0.1892320209761514,
                0.1822715635133293,
                0.1753417770813861,
                0.1684607202339019,
                0.1616452724003933,
                0.1549111192950185,
                0.148272745932097,
                0.1417434368655761,
                0.1353352832366127,
            ];

            foreach ((RayleighDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
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
            ddouble[] expected_dist1 = [
                0.0,
                0.001951218892524476,
                0.007782061739756485,
                0.01742453104209962,
                0.03076676552365587,
                0.0476552001048236,
                0.06789750764047242,
                0.09126624363892988,
                0.1175030974154045,
                0.1463236386548523,
                0.1774224376013354,
                0.2104784303392121,
                0.2451603980109927,
                0.2811324302908242,
                0.3180592488096519,
                0.3556112751748047,
                0.3934693402873666,
                0.431328946281328,
                0.4689040089646548,
                0.5059300264683617,
                0.5421666382283857,
                0.5773995567768112,
                0.6114418724876358,
                0.6441347478578616,
                0.6753475326416503,
                0.7049773438255571,
                0.7329481647736567,
                0.7592095257277515,
                0.7837348331701127,
                0.8065194183929533,
                0.8275783761062472,
                0.8469442622636418,
                0.8646647167633873,
                0.8808000715973527,
                0.8954209987109986,
                0.9086062446439528,
                0.9204404912817723,
                0.9310123740827606,
                0.9404126812380139,
                0.9487327496293423,
                0.9560630663765926,
                0.9624920783757063,
                0.9681052066378429,
                0.9729840575197349,
                0.9772058191163877,
                0.9808428281628708,
                0.9839622907246495,
                0.9866261386829753,
                0.9888910034617577,
                0.9908082885007716,
                0.9924243225557401,
                0.9937805768954682,
                0.9949139307689873,
                0.9958569710366545,
                0.9966383135120678,
                0.9972829352706765,
                0.9978125088818172,
                0.998245731161303,
                0.9985986405811146,
                0.9988849188818363,
                0.9991161736930649,
                0.9993022000668411,
                0.9994512197665679,
                0.9995700979324073,
                0.9996645373720975,
            ];
            ddouble[] expected_dist2 = [
                0.0,
                4.881620601105974e-4,
                0.001951218892524476,
                0.004384889426482075,
                0.007782061739756485,
                0.01213282768599966,
                0.01742453104209962,
                0.02364182884810584,
                0.03076676552365587,
                0.03877885925984348,
                0.0476552001048236,
                0.05737055908584809,
                0.06789750764047242,
                0.07920654656596138,
                0.09126624363892988,
                0.1040433790074319,
                0.1175030974154045,
                0.1316090662848772,
                0.1463236386548523,
                0.1616080199574039,
                0.1774224376013354,
                0.193726312331672,
                0.2104784303392121,
                0.2276371151081346,
                0.2451603980109927,
                0.2630061866889891,
                0.2811324302908242,
                0.2994972806851915,
                0.3180592488096519,
                0.3367773553716046,
                0.3556112751748047,
                0.3745214744067218,
                0.3934693402873666,
                0.4124173025483615,
                0.431328946281328,
                0.4501691157664482,
                0.4689040089646548,
                0.4875012624296635,
                0.5059300264683617,
                0.5241610304492728,
                0.5421666382283857,
                0.5599208937289881,
                0.5773995567768112,
                0.5945801293532851,
                0.6114418724876358,
                0.6279658140625476,
                0.6441347478578616,
                0.6599332242020388,
                0.6753475326416503,
                0.6903656770748552,
                0.7049773438255571,
                0.719173863160669,
                0.7329481647736567,
                0.7462947277733241,
                0.7592095257277515,
                0.7716899673195274,
                0.7837348331701127,
                0.7953442093885424,
                0.8065194183929533,
                0.8172629475428861,
                0.8275783761062472,
                0.8374703010675216,
                0.8469442622636418,
                0.8560066673111608,
                0.8646647167633873,
            ];

            foreach ((RayleighDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
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