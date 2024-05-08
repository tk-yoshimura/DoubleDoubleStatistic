using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.LinearityDistribution {
    [TestClass()]
    public class WeibullDistributionTests {
        readonly WeibullDistribution dist_alpha1mu0theta1 = new(alpha: 1, mu: 0, theta: 1);
        readonly WeibullDistribution dist_alpha2mu1theta1 = new(alpha: 2, mu: 1, theta: 1);
        readonly WeibullDistribution dist_alpha1mu0theta2 = new(alpha: 1, mu: 0, theta: 2);
        readonly WeibullDistribution dist_alpha3mu0theta4 = new(alpha: 3, mu: 0, theta: 4);

        WeibullDistribution[] Dists => [
            dist_alpha1mu0theta1,
            dist_alpha2mu1theta1,
            dist_alpha1mu0theta2,
            dist_alpha3mu0theta4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (WeibullDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (WeibullDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
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
            foreach (WeibullDistribution dist in Dists) {
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

            foreach (WeibullDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
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
        public void IrregularValueTest() {
            foreach (WeibullDistribution dist in Dists) {
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
                1.0,
                0.9394130628134758,
                0.8824969025845955,
                0.8290291181804004,
                0.7788007830714049,
                0.7316156289466418,
                0.6872892787909722,
                0.645648526427892,
                0.6065306597126334,
                0.569782824730923,
                0.5352614285189903,
                0.5028315779709409,
                0.4723665527410147,
                0.4437473100810799,
                0.4168620196785084,
                0.391605626676799,
                0.3678794411714423,
                0.3455907525769745,
                0.3246524673583497,
                0.3049827687110593,
                0.2865047968601901,
                0.2691463487291839,
                0.2528395958047465,
                0.2375208190954581,
                0.2231301601484298,
                0.2096113871510978,
                0.1969116752041941,
                0.1849813999073043,
                0.1737739434504451,
                0.1632455124539584,
                0.1533549668449285,
                0.1440636591014533,
                0.1353352832366127,
                0.1271357329320356,
                0.1194329682667196,
                0.1121968905203437,
                0.1053992245618643,
                0.0990134083638263,
                0.09301448921066349,
                0.08737902619542039,
                0.0820849986238988,
                0.07711171996831671,
                0.07243975703425146,
                0.0680508540250102,
                0.06392786120670757,
                0.06005466789530794,
                0.05641613950377735,
                0.0529980584033558,
                0.04978706836786394,
                0.04677062238395898,
                0.04393693362340742,
                0.04127492938579755,
                0.03877420783172201,
                0.03642499733736423,
                0.03421811831166603,
                0.03214494732687607,
                0.0301973834223185,
                0.0283678164497131,
                0.02664909733635549,
                0.02503451014996015,
                0.02351774585600911,
                0.02209287766506244,
                0.02075433787369974,
                0.01949689610859799,
                0.01831563888873418,
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
                0.0,
                0.1245126711837647,
                0.2461241092513521,
                0.3620454570741892,
                0.4697065314067379,
                0.5668503861796147,
                0.6516112921971324,
                0.7225724099563381,
                0.7788007830714049,
                0.8198587461594277,
                0.8457923077021612,
                0.8570984248194972,
                0.8546742370963845,
                0.839752197016812,
                0.8138255792345985,
                0.7785690537784524,
                0.7357588823428847,
                0.6871968805921976,
                0.6346416413110848,
                0.5797497045206904,
                0.5240284678777445,
                0.4688017283576435,
                0.4151879007537652,
                0.3640902233362898,
                0.316197673685593,
                0.2719948989256984,
                0.2317792187677136,
                0.1956826771883586,
                0.1636971783438564,
                0.1357009160388789,
                0.1114845614480953,
                0.09077597805933081,
                0.07326255555493671,
                0.05861056959118451,
                0.04648126192007111,
                0.03654358101660443,
                0.02848371942368586,
                0.02201172794984845,
                0.01686558064690206,
                0.01281311343276757,
                0.009652270681138546,
                0.007210076445942573,
                0.005340708678977179,
                0.00392300367102134,
                0.002857660751851611,
                0.002064359434484778,
                0.001478946683103824,
                0.001050803478346171,
                7.404588245200774e-4,
                5.174863067459099e-4,
                3.586930546216797e-4,
                2.46592801426553e-4,
                1.681426514472518e-4,
                1.137160645661708e-4,
                7.628131829123758e-5,
                5.07542801104229e-5,
                3.349582174490306e-5,
                2.192689650376544e-5,
                1.423760960151633e-5,
                9.170118513128587e-6,
                5.858617056228368e-6,
                3.712801193713363e-6,
                2.333988020694438e-6,
                1.455424328298783e-6,
                9.002813977540729e-7,
            ];
            ddouble[] expected_dist_alpha1mu0theta2 = [
                0.5,
                0.4846166172381721,
                0.4697065314067379,
                0.4552551806900171,
                0.4412484512922977,
                0.4276726636537113,
                0.4145145590902002,
                0.4017612868445304,
                0.3894003915357024,
                0.3774198009945037,
                0.3658078144733209,
                0.3545530912186992,
                0.3436446393954861,
                0.3330718053517439,
                0.322824263213946,
                0.3128920048022956,
                0.3032653298563167,
                0.2939348365611733,
                0.2848914123654615,
                0.2761262250815102,
                0.2676307142594951,
                0.2593965828269447,
                0.2514157889854705,
                0.2436805383568096,
                0.2361832763705073,
                0.2289166808858071,
                0.2218736550405399,
                0.2150473203200311,
                0.2084310098392542,
                0.202018261831671,
                0.1958028133383995,
                0.1897785940915448,
                0.1839397205857212,
                0.1782804903319735,
                0.1727953762884873,
                0.1674790214626475,
                0.1623262336791749,
                0.1573319805092295,
                0.1524913843555296,
                0.1477997176886854,
                0.143252398430095,
                0.138844985476895,
                0.1345731743645919,
                0.1304327930631425,
                0.1264197979023732,
                0.1225302696227629,
                0.1187604095477291,
                0.1151065358736807,
                0.1115650800742149,
                0.1081325834149437,
                0.1048056935755489,
                0.1015811613757659,
                0.09845583760209703,
                0.09542666993215816,
                0.09249069995365214,
                0.08964506027505932,
                0.08688697172522257,
                0.08421374063909212,
                0.0816227562269792,
                0.07911148802474921,
                0.07667748342246423,
                0.07431836526906126,
                0.07203182955072664,
                0.06981564314069949,
                0.06766764161830635,
            ];
            ddouble[] expected_dist_alpha3mu0theta4 = [
                0.0,
                1.831047702594013e-4,
                7.323995235992663e-4,
                0.001647779493951484,
                0.002928972331567082,
                0.00457545445180008,
                0.006586367638828394,
                0.008960436091775464,
                0.01169588415360323,
                0.01479035501534695,
                0.01824083073942164,
                0.02204355399545109,
                0.02619395194989625,
                0.03068656279617422,
                0.03551496545439451,
                0.04067171300867993,
                0.04614827048462852,
                0.05193495759912278,
                0.05802089713868872,
                0.06439396964022638,
                0.07104077505844805,
                0.07794660210705687,
                0.08509540595488846,
                0.09246979494327763,
                0.1000510269662106,
                0.1078190161198724,
                0.1157523501825748,
                0.1238283194294483,
                0.1320229572185155,
                0.1403110927058044,
                0.1486664159571202,
                0.1570615556232912,
                0.1654681692346117,
                0.1738570460495128,
                0.1821982222631128,
                0.190461108244305,
                0.198614627326791,
                0.2066273655314567,
                0.2144677314464662,
                0.2221041253393288,
                0.2295051164240678,
                0.2366396270587297,
                0.2434771225061823,
                0.2499878047569186,
                0.2561428087889086,
                0.2619143995289569,
                0.2672761676850005,
                0.2722032225417365,
                0.2766723797551647,
                0.2806623421471554,
                0.28415387149086,
                0.2871299492932247,
                0.2895759246232709,
                0.291479647105005,
                0.2928315832922252,
                0.2936249147690485,
                0.2938556164741504,
                0.2935225139274373,
                0.2926273182435873,
                0.2911746380455146,
                0.2891719676397348,
                0.2866296510817592,
                0.2835608220395059,
                0.2799813196523486,
                0.2759095808785818,
            ];

            foreach ((WeibullDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1mu0theta1, expected_dist_alpha1mu0theta1),
                (dist_alpha2mu1theta1, expected_dist_alpha2mu1theta1),
                (dist_alpha1mu0theta2, expected_dist_alpha1mu0theta2),
                (dist_alpha3mu0theta4, expected_dist_alpha3mu0theta4),
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
            ddouble[] expected_dist_alpha1mu0theta1 = [
                0.0,
                0.06058693718652419,
                0.1175030974154045,
                0.1709708818195996,
                0.2211992169285951,
                0.2683843710533582,
                0.3127107212090278,
                0.3543514735721079,
                0.3934693402873666,
                0.430217175269077,
                0.4647385714810097,
                0.4971684220290591,
                0.5276334472589853,
                0.5562526899189202,
                0.5831379803214916,
                0.608394373323201,
                0.6321205588285577,
                0.6544092474230254,
                0.6753475326416503,
                0.6950172312889407,
                0.7134952031398099,
                0.7308536512708161,
                0.7471604041952535,
                0.7624791809045419,
                0.7768698398515702,
                0.7903886128489022,
                0.8030883247958059,
                0.8150186000926958,
                0.8262260565495548,
                0.8367544875460415,
                0.8466450331550716,
                0.8559363408985468,
                0.8646647167633873,
                0.8728642670679644,
                0.8805670317332803,
                0.8878031094796562,
                0.8946007754381357,
                0.9009865916361737,
                0.9069855107893365,
                0.9126209738045796,
                0.9179150013761012,
                0.9228882800316833,
                0.9275602429657486,
                0.9319491459749898,
                0.9360721387932924,
                0.939945332104692,
                0.9435838604962227,
                0.9470019415966442,
                0.950212931632136,
                0.953229377616041,
                0.9560630663765926,
                0.9587250706142024,
                0.961225792168278,
                0.9635750026626357,
                0.965781881688334,
                0.9678550526731239,
                0.9698026165776815,
                0.9716321835502869,
                0.9733509026636445,
                0.9749654898500398,
                0.9764822541439909,
                0.9779071223349376,
                0.9792456621263003,
                0.980503103891402,
                0.9816843611112658,
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
                0.0,
                0.003898630529882485,
                0.01550356299459155,
                0.03454544780216218,
                0.06058693718652419,
                0.09303938211261642,
                0.1311849437371568,
                0.1742029600498993,
                0.2211992169285951,
                0.2712366700805088,
                0.3233661538382711,
                0.3766556910403657,
                0.430217175269077,
                0.4832294172204233,
                0.5349568118659437,
                0.5847631713181587,
                0.6321205588285577,
                0.6766132326624952,
                0.7179370483061845,
                0.7558948612544462,
                0.7903886128489022,
                0.8214088653875644,
                0.8490225815440854,
                0.8733599223178122,
                0.8946007754381357,
                0.9129616323437765,
                0.928683317302242,
                0.9420199474997456,
                0.953229377616041,
                0.9625652645409989,
                0.9702707836138412,
                0.9765739411459792,
                0.9816843611112658,
                0.9857913770688037,
                0.989063232489395,
                0.991647181481919,
                0.9936702845725143,
                0.995240707470303,
                0.9964493514427575,
                0.9973716690394323,
                0.9980695458637723,
                0.9985931558154258,
                0.9989827221563853,
                0.999270138851903,
                0.9994804253178452,
                0.9996330027672027,
                0.9997427918811993,
                0.999821139833473,
                0.9998765901959134,
                0.999915512439715,
                0.9999426091112605,
                0.9999613187762468,
                0.9999741318997774,
                0.9999828353110088,
                0.9999886990639568,
                0.9999926175592566,
                0.9999952148826079,
                0.9999969225408416,
                0.9999980361917791,
                0.9999987565940999,
                0.9999992188510591,
                0.9999995130752533,
                0.9999996988402554,
                0.9999998151842123,
                0.9999998874648253,

            ];
            ddouble[] expected_dist_alpha1mu0theta2 = [
                0.0,
                0.03076676552365587,
                0.06058693718652419,
                0.08948963861996584,
                0.1175030974154045,
                0.1446546726925775,
                0.1709708818195996,
                0.1964774263109392,
                0.2211992169285951,
                0.2451603980109927,
                0.2683843710533582,
                0.2908938175626016,
                0.3127107212090278,
                0.3338563892965122,
                0.3543514735721079,
                0.3742159903954089,
                0.3934693402873666,
                0.4121303268776535,
                0.430217175269077,
                0.4477475498369796,
                0.4647385714810097,
                0.4812068343461107,
                0.4971684220290591,
                0.5126389232863808,
                0.5276334472589853,
                0.5421666382283857,
                0.5562526899189202,
                0.5699053593599377,
                0.5831379803214916,
                0.5959634763366579,
                0.608394373323201,
                0.6204428118169104,
                0.6321205588285577,
                0.6434390193360531,
                0.6544092474230254,
                0.665041957074705,
                0.6753475326416503,
                0.685336038981541,
                0.6950172312889407,
                0.7044005646226292,
                0.7134952031398099,
                0.7223100290462101,
                0.7308536512708161,
                0.739134413873715,
                0.7471604041952535,
                0.7549394607544742,
                0.7624791809045419,
                0.7697869282526385,
                0.7768698398515702,
                0.7837348331701127,
                0.7903886128489022,
                0.7968376772484682,
                0.8030883247958059,
                0.8091466601356837,
                0.8150186000926958,
                0.8207098794498814,
                0.8262260565495548,
                0.8315725187218157,
                0.8367544875460415,
                0.8417770239505016,
                0.8466450331550716,
                0.8513632694618775,
                0.8559363408985468,
                0.860368713718601,
                0.8646647167633873,
            ];
            ddouble[] expected_dist_alpha3mu0theta4 = [
                0.0,
                3.814689989667386e-6,
                3.051711246848665e-5,
                1.029915221808508e-4,
                2.441108251027835e-4,
                4.767234894332839e-4,
                8.236352355147636e-4,
                0.001307585526195809,
                0.001951218892524476,
                0.002777051146319209,
                0.003807430551052815,
                0.005064494045535328,
                0.006570118640971789,
                0.008345868140989476,
                0.01041293536598686,
                0.0127920800974638,
                0.01550356299459155,
                0.01856707577390126,
                0.0220016679832552,
                0.02582567074287268,
                0.03005661786865599,
                0.03471116483596404,
                0.03980500608478477,
                0.0453527912094136,
                0.05136804061666977,
                0.05786306127573881,
                0.06484886321927685,
                0.07233507748876056,
                0.0803298762465312,
                0.08883989580186291,
                0.09787016331797915,
                0.1074240279805745,
                0.1175030974154045,
                0.1281071801422657,
                0.1392342348446252,
                0.1508803272177709,
                0.1630395951331937,
                0.1757042228226475,
                0.1888644247417122,
                0.2025084397195784,
                0.2166225359391818,
                0.2311910272198638,
                0.2461963009937093,
                0.2616188582770229,
                0.2774373658406545,
                0.2936287206777893,
                0.3101681267562938,
                0.3270291839257868,
                0.3441839887284984,
                0.3616032467389985,
                0.3792563959324999,
                0.39711174045621,
                0.4151365940547744,
                0.4332974322809539,
                0.4515600525080227,
                0.4698897406527724,
                0.4882514434191666,
                0.5066099447843303,
                0.5249300453722816,
                0.543176743298122,
                0.5613154150176468,
                0.579311994685699,
                0.5971331505120305,
                0.6147464566066952,
                0.6321205588285577,
            ];

            foreach ((WeibullDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha1mu0theta1, expected_dist_alpha1mu0theta1),
                (dist_alpha2mu1theta1, expected_dist_alpha2mu1theta1),
                (dist_alpha1mu0theta2, expected_dist_alpha1mu0theta2),
                (dist_alpha3mu0theta4, expected_dist_alpha3mu0theta4),
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