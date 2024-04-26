
using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.LinearityDistribution {
    [TestClass()]
    public class JohnsonSBDistributionTests {
        readonly JohnsonSBDistribution dist_gamma1delta1mu0sigma1 = new(gamma: 1, delta: 1, mu: 0, sigma: 1);
        readonly JohnsonSBDistribution dist_gamma2delta1mu1sigma1 = new(gamma: 2, delta: 1, mu: 1, sigma: 1);
        readonly JohnsonSBDistribution dist_gamma1delta1mu0sigma2 = new(gamma: 1, delta: 1, mu: 0, sigma: 2);
        readonly JohnsonSBDistribution dist_gamma2delta2mu2sigma2 = new(gamma: 2, delta: 2, mu: 2, sigma: 2);
        readonly JohnsonSBDistribution dist_gamma3delta2mu0sigma4 = new(gamma: 3, delta: 2, mu: 0, sigma: 4);

        JohnsonSBDistribution[] Dists => [
            dist_gamma1delta1mu0sigma1,
            dist_gamma2delta1mu1sigma1,
            dist_gamma1delta1mu0sigma2,
            dist_gamma2delta2mu2sigma2,
            dist_gamma3delta2mu0sigma4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Gamma={dist.Gamma}");
                Console.WriteLine($"Delta={dist.Delta}");
                Console.WriteLine($"Mu={dist.Mu}");
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            foreach (JohnsonSBDistribution dist in Dists) {
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
            ddouble[] expected_dist_gamma1delta1mu0sigma1 = [
                0,
                1.583236176729337,
                2.331848722396709,
                2.348888257446448,
                2.117372024464338,
                1.815808061979008,
                1.510209773266209,
                1.22487542167473,
                0.9678828980765735,
                0.7409740205192809,
                0.5436755183758354,
                0.3751669549543404,
                0.2352635582738152,
                0.1250887237693375,
                0.04758874943666752,
                0.007036605229908166,
                0,
            ];
            ddouble[] expected_dist_gamma2delta1mu1sigma1 = [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,
                5.29901612499607,
                3.642140452092024,
                2.271133856570766,
                1.417348676737599,
                0.8913553960700091,
                0.5616222476110218,
                0.3513942630002032,
                0.2159638660527523,
                0.1285928403388946,
                0.07278624329038846,
                0.03805048306425488,
                0.01749813181157528,
                0.006441015454018838,
                0.001516926469009589,
                1.046719234567122e-4,
                0,
            ];
            ddouble[] expected_dist_gamma1delta1mu0sigma2 = [
                0,
                0.3407016986888479,
                0.7916180883646685,
                1.04990429161722,
                1.165924361198354,
                1.195470165174227,
                1.174444128723224,
                1.124507304110124,
                1.058686012232169,
                0.9848774909746337,
                0.9079040309895038,
                0.8307155850465431,
                0.7551048866331046,
                0.6821410600695512,
                0.6124377108373649,
                0.546321736088979,
                0.4839414490382867,
                0.425336991764776,
                0.3704870102596404,
                0.3193402746585988,
                0.2718377591879177,
                0.2279287659651513,
                0.1875834774771702,
                0.1508035477673825,
                0.1176317791369076,
                0.0881613726422337,
                0.06254436188466873,
                0.04099691924465798,
                0.02379437471833376,
                0.01123559884013672,
                0.003518302614954083,
                3.545283024857936e-4,
                0,
            ];
            ddouble[] expected_dist_gamma2delta2mu2sigma2 = [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,
                9.420461790050085e-5,
                0.01990703148319264,
                0.1877787816624732,
                0.6092920198403676,
                1.179359444017865,
                1.695091533914113,
                2.011190385025742,
                2.086710981615408,
                1.958528065441002,
                1.697921998023665,
                1.377638554113488,
                1.054757843469345,
                0.7657333144067501,
                0.5283713195522136,
                0.3466605133390498,
                0.2159638660527523,
                0.1273629404704437,
                0.07075929224262907,
                0.03677863044090418,
                0.0177158815001661,
                0.007807687305278805,
                0.0030941138402285,
                0.001076582179386316,
                3.180477033402531e-4,
                7.598317464968739e-5,
                1.363378289881825e-5,
                1.631158300790792e-6,
                1.056917697315775e-7,
                2.462819431717822e-9,
                7.767413595027165e-12,
                1.104532862523103e-16,
                0,
            ];
            ddouble[] expected_dist_gamma3delta2mu0sigma4 = [
                0,
                1.108656395712615e-5,
                0.003715603639714572,
                0.04606619729261412,
                0.183832723339175,
                0.4249735662477178,
                0.7201679716610009,
                1.003570363207189,
                1.225336497948653,
                1.361773634359102,
                1.411456533386827,
                1.386816342966606,
                1.306384892925227,
                1.189415024955719,
                1.052860714204314,
                0.9101073227382754,
                0.7707945062446888,
                0.6412181735187481,
                0.5249697173545501,
                0.423611415504514,
                0.3372848981884491,
                0.2652109264193784,
                0.2060735757810673,
                0.1582994592407999,
                0.1202497168468439,
                0.09034378523670177,
                0.06713234828232523,
                0.04933407906796535,
                0.03584775593986156,
                0.02574854050104021,
                0.01827483649500986,
                0.01281024854062774,
                0.008863696823876015,
                0.006049654929822183,
                0.004069692090798817,
                0.002695951404624962,
                0.00175682434513487,
                0.001124841955227757,
                7.066573090347017e-4,
                4.349137141297019e-4,
                2.617574595412118e-4,
                1.537471885277313e-4,
                8.792310000778556e-5,
                4.882041675444218e-5,
                2.623763742947755e-5,
                1.359755114801474e-5,
                6.765662555565809e-6,
                3.215278407777746e-6,
                1.450385849501046e-6,
                6.163604407478016e-7,
                2.444934884331147e-7,
                8.949762999129371e-8,
                2.979929424376413e-8,
                8.859735113915236e-9,
                2.295844082232052e-9,
                5.019526608637457e-10,
                8.852764054055068e-11,
                1.180935123499138e-11,
                1.081715081433543e-12,
                5.83146709907311e-14,
                1.416862520575667e-15,
                9.223152717461595e-18,
                4.717251742989597e-21,
                2.836027079166137e-27,
                0,
            ];

            foreach ((JohnsonSBDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1delta1mu0sigma1, expected_dist_gamma1delta1mu0sigma1),
                (dist_gamma2delta1mu1sigma1, expected_dist_gamma2delta1mu1sigma1),
                (dist_gamma1delta1mu0sigma2, expected_dist_gamma1delta1mu0sigma2),
                (dist_gamma2delta2mu2sigma2, expected_dist_gamma2delta2mu2sigma2),
                (dist_gamma3delta2mu0sigma4, expected_dist_gamma3delta2mu0sigma4),
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
            ddouble[] expected_dist_gamma1delta1mu0sigma1 = [
                0,
                0.04381351399555408,
                0.1720972082596352,
                0.320487124626068,
                0.4607230563177007,
                0.5837680683314519,
                0.6876408756174844,
                0.7729766290101729,
                0.8413447460685429,
                0.8945901082661456,
                0.9345835577820661,
                0.9631488754519609,
                0.9820744538996565,
                0.9931748593051336,
                0.9983899703851478,
                0.9998955693717297,
                1,
            ];
            ddouble[] expected_dist_gamma2delta1mu1sigma1 = [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,
                0.2394570415614924,
                0.5215682109078519,
                0.7032126199969604,
                0.816308893562055,
                0.8871562484608512,
                0.9317792717940867,
                0.9598273072906667,
                0.9772498680518208,
                0.9878171851977938,
                0.9939775410701952,
                0.997352014073915,
                0.9990278536212877,
                0.9997361993016903,
                0.9999602512884768,
                0.9999987495124183,
                1,
            ];
            ddouble[] expected_dist_gamma1delta1mu0sigma2 = [
                0,
                0.007466761491408705,
                0.04381351399555408,
                0.1022769774639053,
                0.1720972082596352,
                0.246230782088612,
                0.320487124626068,
                0.3924398003322776,
                0.4607230563177007,
                0.5246112194186608,
                0.5837680683314519,
                0.6380955258460423,
                0.6876408756174844,
                0.7325390882933769,
                0.7729766290101729,
                0.8091686601235483,
                0.8413447460685429,
                0.8697400458730423,
                0.8945901082661456,
                0.9161280780478223,
                0.9345835577820661,
                0.950182648337472,
                0.9631488754519609,
                0.9737048327293886,
                0.9820744538996565,
                0.9884858757701355,
                0.9931748593051336,
                0.9963886604531905,
                0.9983899703851478,
                0.9994597546300091,
                0.9998955693717297,
                0.9999953746940579,
                1,
            ];
            ddouble[] expected_dist_gamma2delta2mu2sigma2 = [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,
                5.637396298579248e-7,
                3.176241029613558e-4,
                0.005584488593645441,
                0.02925746393726292,
                0.0849075907848647,
                0.175494120990015,
                0.2925565632465681,
                0.4218258991971268,
                0.5491288625282047,
                0.6638834788809165,
                0.760137728626455,
                0.8360490673695286,
                0.8927009594209057,
                0.9328516437488955,
                0.9599126840247024,
                0.9772498680518208,
                0.9877858797453531,
                0.9938362629202276,
                0.9971008891166004,
                0.998742999945778,
                0.9995048257668838,
                0.9998261632322979,
                0.9999470236435457,
                0.999986489702077,
                0.99999726535691,
                0.9999995944429446,
                0.9999999612379,
                0.9999999980901778,
                0.9999999999686944,
                0.9999999999999397,
                1.0,
            ];
            ddouble[] expected_dist_gamma3delta2mu0sigma4 = [
                0,
                6.241798648555458e-8,
                5.487158935837243e-5,
                0.001245127059401179,
                0.00784386745671148,
                0.02642166260597611,
                0.06210172126578334,
                0.1161837916916424,
                0.1862446300559798,
                0.2675573943292653,
                0.3546494303059947,
                0.4424370434526402,
                0.526838855816584,
                0.6049756790177654,
                0.6751102911008358,
                0.7364567675497111,
                0.7889477249531487,
                0.8330108041443387,
                0.8693790635585421,
                0.8989432075158943,
                0.9226442824031176,
                0.9414011117625819,
                0.9560653327173981,
                0.9673970826303606,
                0.9760553116028413,
                0.9825978658628528,
                0.987487634736846,
                0.9911020573736932,
                0.9937441035874011,
                0.995653477860396,
                0.9970172683321324,
                0.9979796020019763,
                0.9986501019683699,
                0.9991110977411649,
                0.9994236366564808,
                0.9996324000043708,
                0.9997696545539952,
                0.999858378355257,
                0.9999146960053054,
                0.9999497479456885,
                0.9999711042095204,
                0.9999838176349906,
                0.9999911963440573,
                0.9999953611120459,
                0.9999976405814502,
                0.9999988462882141,
                0.9999994601840055,
                0.9999997596602858,
                0.9999998988570361,
                0.9999999600911379,
                0.9999999853803232,
                0.9999999950883691,
                0.9999999985097999,
                0.9999999995996613,
                0.9999999999071948,
                0.9999999999820665,
                0.9999999999972459,
                0.999999999999686,
                0.999999999999976,
                0.999999999999999,
                1.0,
                1.0,
                1.0,
                1.0,
                1,
            ];

            foreach ((JohnsonSBDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1delta1mu0sigma1, expected_dist_gamma1delta1mu0sigma1),
                (dist_gamma2delta1mu1sigma1, expected_dist_gamma2delta1mu1sigma1),
                (dist_gamma1delta1mu0sigma2, expected_dist_gamma1delta1mu0sigma2),
                (dist_gamma2delta2mu2sigma2, expected_dist_gamma2delta2mu2sigma2),
                (dist_gamma3delta2mu0sigma4, expected_dist_gamma3delta2mu0sigma4),
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