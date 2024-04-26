using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class NakagamiDistributionTests {
        readonly NakagamiDistribution dist_m1omega1 = new(m: 1, omega: 1);
        readonly NakagamiDistribution dist_m2omega1 = new(m: 2, omega: 1);
        readonly NakagamiDistribution dist_m1omega2 = new(m: 1, omega: 2);
        readonly NakagamiDistribution dist_m2omega2 = new(m: 2, omega: 2);
        readonly NakagamiDistribution dist_m3omega4 = new(m: 3, omega: 4);

        NakagamiDistribution[] Dists => [
            dist_m1omega1,
            dist_m2omega1,
            dist_m1omega2,
            dist_m2omega2,
            dist_m3omega4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NakagamiDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"M={dist.M}");
                Console.WriteLine($"Omega={dist.Omega}");
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
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-15, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NakagamiDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (NakagamiDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
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
            foreach (NakagamiDistribution dist in Dists) {
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
            ddouble[] expected_dist_m1omega1 = [
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
            ddouble[] expected_dist_m2omega1 = [
                0.0,
                0.001937925660664538,
                0.01514426928869288,
                0.04915384237052196,
                0.1103121128230744,
                0.2008246001949865,
                0.3184479570891125,
                0.4568470266763465,
                0.6065306597126334,
                0.7561894091108724,
                0.8942057847101841,
                1.010099351013587,
                1.09570207733443,
                1.145923597641165,
                1.159046128478927,
                1.136568321565265,
                1.082682265892902,
                1.00350904947825,
                0.9062325289935622,
                0.7982605847430909,
                0.6865145878657409,
                0.5769095338416736,
                0.4740477305638754,
                0.3811148608460745,
                0.2999429065325422,
                0.2311913282550029,
                0.1745952203208579,
                0.1292345217616613,
                0.0937886816920912,
                0.06675342747499186,
                0.04660802790477803,
                0.03193107799643392,
                0.02146960818576076,
                0.01417019532968772,
                0.009182157766149574,
                0.005842520746637452,
                0.003650950224932666,
                0.002240887273938895,
                0.001351127100146497,
                8.003573947247227e-4,
                4.658316465098339e-4,
                2.664241620762214e-4,
                1.497466382669363e-4,
                8.272102319030211e-5,
                4.491423734970214e-5,
                2.3971386795448e-5,
                1.257687892591686e-5,
                6.487104206863419e-6,
                3.289675624857928e-6,
                1.640226475725823e-6,
                8.041294214614455e-7,
                3.876510619356454e-7,
                1.83767683032128e-7,
                8.567014963026606e-8,
                3.927716676168141e-8,
                1.770997902799997e-8,
                7.853790520564247e-9,
                3.425620130793742e-9,
                1.469644071947627e-9,
                6.201716673930744e-10,
                2.574254535864746e-10,
                1.051098068683002e-10,
                4.221812562577485e-11,
                1.66812973063062e-11,
                6.484052761136218e-12,
            ];
            ddouble[] expected_dist_m1omega2 = [
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
            ddouble[] expected_dist_m2omega2 = [
                0.0,
                4.863776218115808e-4,
                0.003845689207052377,
                0.01272816060026446,
                0.02935665821292112,
                0.055356483025353,
                0.09163283796522174,
                0.1383048753432053,
                0.1947001957678512,
                0.2594084314020064,
                0.3303876201961567,
                0.4051129273560905,
                0.4807542583667163,
                0.5543676613118799,
                0.6230852091014896,
                0.684289207422468,
                0.7357588823428847,
                0.7757808534810356,
                0.8032183272843416,
                0.8175376692655048,
                0.8187944810589758,
                0.8075842273660968,
                0.7849646248625873,
                0.752358313066005,
                0.7114447657925842,
                0.6640500462053185,
                0.6120419995584938,
                0.5572369987121618,
                0.5013226086780603,
                0.4457987124558482,
                0.3919379113409601,
                0.3407645113867848,
                0.2930502222197469,
                0.2493238682999998,
                0.2098919483578211,
                0.1748667450989861,
                0.1441988295824097,
                0.1177111545443068,
                0.09513241583643194,
                0.07612791223140421,
                0.06032669175711591,
                0.0473442910376151,
                0.03680082074107712,
                0.02833450698327522,
                0.02161105943587781,
                0.01632940568293623,
                0.01222441867753005,
                0.009067284701823015,
                0.006664129420680696,
                0.004853455556628632,
                0.003502861861539841,
                0.002505421392619001,
                0.001776006755911597,
                0.001247767286587397,
                8.68891891161128e-4,
                5.997331927110517e-4,
                4.103238163750624e-4,
                2.782831513309919e-4,
                1.870910886699255e-4,
                1.246921193132836e-4,
                8.238680235321143e-5,
                5.396614547581025e-5,
                3.504628887323992e-5,
                2.256476233991355e-5,
                1.440450236406517e-5,
            ];
            ddouble[] expected_dist_m3omega4 = [
                0.0,
                4.011543720343199e-7,
                1.272460960010234e-5,
                9.522237876663804e-5,
                3.931210324274915e-4,
                0.001168490593731596,
                0.002815370999154394,
                0.005857726202300761,
                0.01092958310101114,
                0.01873856406126449,
                0.03001574764377109,
                0.04545621827150435,
                0.0656556526176807,
                0.09104874753045546,
                0.1218551870739589,
                0.1580381988119083,
                0.1992796394376156,
                0.2449740981720289,
                0.2942428626903922,
                0.3459669183426167,
                0.3988366035876467,
                0.451414258669926,
                0.5022052834208746,
                0.5497325254507687,
                0.5926088695272574,
                0.6296032677368038,
                0.6596961765920599,
                0.6821213625393486,
                0.6963921957743655,
                0.7023117634095141,
                0.699967292350842,
                0.6897103903510876,
                0.6721254229661633,
                0.6479889023874231,
                0.6182230554636722,
                0.5838467712327376,
                0.5459269320516508,
                0.5055327512374558,
                0.463695227467616,
                0.4213732386654074,
                0.3794271901048886,
                0.3386005508181308,
                0.2995090972377834,
                0.2626372602069591,
                0.2283406561978159,
                0.1968536800250182,
                0.168300939398794,
                0.1427113089943382,
                0.1200334561251702,
                0.100151821882196,
                0.08290221059340197,
                0.06808632780477761,
                0.05548479633616167,
                0.04486835826817125,
                0.03600712851823486,
                0.02867789706903242,
                0.0226695792392785,
                0.01778698663676174,
                0.01385313766465323,
                0.01071034913023647,
                0.008220353898923945,
                0.006263678195930019,
                0.004738490502753536,
                0.003559106018343376,
                0.002654299736637787,
            ];

            foreach ((NakagamiDistribution dist, ddouble[] expecteds) in new[]{
                (dist_m1omega1, expected_dist_m1omega1),
                (dist_m2omega1, expected_dist_m2omega1),
                (dist_m1omega2, expected_dist_m1omega2),
                (dist_m2omega2, expected_dist_m2omega2),
                (dist_m3omega4, expected_dist_m3omega4),
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
            ddouble[] expected_dist_m1omega1 = [
                0.0,
                0.003898630529882485,
                0.01550356299459155,
                0.03454544780216218,
                0.06058693718652419,
                0.09303938211261642,
                0.1311849437371568,
                0.1742029600498993,
                0.2211992169285951,
                0.2712366700805089,
                0.3233661538382711,
                0.3766556910403656,
                0.4302171752690769,
                0.4832294172204233,
                0.5349568118659437,
                0.5847631713181586,
                0.6321205588285578,
                0.6766132326624952,
                0.717937048306184,
                0.7558948612544462,
                0.7903886128489017,
                0.8214088653875641,
                0.8490225815440855,
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
            ddouble[] expected_dist_m2omega1 = [
                0.0,
                3.035909709836471e-5,
                4.782269462701993e-4,
                0.002359051146443192,
                0.007190984592330141,
                0.01676275744534628,
                0.03286175995158436,
                0.05700380499459667,
                0.09020401043104986,
                0.1328198271376004,
                0.1844843243443122,
                0.2441330175736043,
                0.3101135068635068,
                0.3803562885763752,
                0.4525787964618478,
                0.5244934903555099,
                0.5939941502901618,
                0.6593012223631752,
                0.7190554848387585,
                0.7723578212921002,
                0.8187601488034444,
                0.858217676382286,
                0.8910153226502284,
                0.9176814453601148,
                0.9389005190396673,
                0.9554337100349396,
                0.9680531276427015,
                0.9774924584362661,
                0.984414125782947,
                0.9893912712742192,
                0.9929017699724279,
                0.9953310806702539,
                0.9969808363488774,
                0.9980805156004976,
                0.9988001332615779,
                0.9992625137660455,
                0.9995542735665034,
                0.9997350910509477,
                0.9998451689895443,
                0.9999110039384412,
                0.9999496901821769,
                0.9999720281882621,
                0.9999847035611921,
                0.9999917723238031,
                0.9999956469296634,
                0.9999977345230937,
                0.9999988402023379,
                0.9999994159150658,
                0.9999997106303848,
                0.9999998589658133,
                0.9999999323759322,
                0.9999999680997581,
                0.9999999851948658,
                0.9999999932397018,
                0.9999999969628691,
                0.9999999986575011,
                0.9999999994161176,
                0.9999999997501349,
                0.9999999998947887,
                0.9999999999564084,
                0.9999999999822281,
                0.9999999999928705,
                0.9999999999971856,
                0.9999999999989068,
                0.9999999999995821,
            ];
            ddouble[] expected_dist_m1omega2 = [
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
                0.3180592488096518,
                0.3556112751748048,
                0.3934693402873666,
                0.4313289462813279,
                0.4689040089646548,
                0.5059300264683616,
                0.5421666382283856,
                0.5773995567768111,
                0.6114418724876358,
                0.6441347478578616,
                0.6753475326416503,
                0.7049773438255571,
                0.7329481647736567,
                0.7592095257277514,
                0.7837348331701127,
                0.8065194183929532,
                0.8275783761062468,
                0.8469442622636415,
                0.8646647167633872,
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
                0.9729840575197348,
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
            ddouble[] expected_dist_m2omega2 = [
                0.0,
                7.60955538980923e-6,
                1.208061663821036e-4,
                6.036862014568989e-4,
                0.001873620760681982,
                0.004469009272051627,
                0.00900782645019449,
                0.01614024537195036,
                0.02649902116074387,
                0.04065139772316972,
                0.0590560576813457,
                0.0820281075086634,
                0.1097143363579328,
                0.1420800871823433,
                0.1789081209508068,
                0.2198089273595092,
                0.2642411176571153,
                0.3115398898478902,
                0.3609511250686992,
                0.4116684741952863,
                0.4628708204253118,
                0.5137577311528606,
                0.563580899775872,
                0.6116700742948542,
                0.6574525201739408,
                0.7004656175580745,
                0.740362702053475,
                0.776912688622068,
                0.8099943465651667,
                0.839586309380765,
                0.865754007256252,
                0.8886347124010024,
                0.9084218055563291,
                0.9253492271778947,
                0.9396768916993194,
                0.9516776397450079,
                0.9616261002208677,
                0.9697896470282907,
                0.9764214744245613,
                0.9817556870432468,
                0.9860042075123491,
                0.9893552453690619,
                0.9919730420152277,
                0.9939986026689681,
                0.9955511417840492,
                0.9967299973124585,
                0.9976168060242376,
                0.9982777722246521,
                0.9987659019591332,
                0.9991231115325103,
                0.9993821512134141,
                0.9995683114989733,
                0.9997009000911755,
                0.9997944930790712,
                0.9998599743393404,
                0.9999053836403169,
                0.9999365971945543,
                0.9999578652564443,
                0.9999722305243763,
                0.9999818491880913,
                0.9999882339440788,
                0.9999924355479778,
                0.9999951767384653,
                0.9999969498175659,
                0.9999980869020297,
            ];
            ddouble[] expected_dist_m3omega4 = [
                0.0e0,
                4.1817537345054746328390475557141e-9,
                2.6587450570180753031918637844933e-7,
                2.9954184378914347064865548164375e-6,
                1.6573810366004086631349679100396e-5,
                6.1989769550581799842283852969679e-5,
                1.8069988877626324396953472280955e-4,
                4.4290570831317064619597282705846e-4,
                9.5514469275974289055414206117026e-4,
                1.866109942385609958412637954311e-3,
                3.3697440841221416936168940286297e-3,
                5.7048827314497612790092954581088e-3,
                9.1510288687286724029719925286063e-3,
                1.402019094001733908638973308594e-2,
                2.0645079803218997707236336136503e-2,
                2.9364298234851545352887240419881e-2,
                4.0505439744813876125842943396487e-2,
                5.4367215761335695879153622422931e-2,
                7.120183499029593686291958660035e-2,
                9.1198858975921616968508143247462e-2,
                1.1447165665546218044816343477913e-1,
                1.4104739028051949833257919055661e-1,
                1.708612052539284118970397158514e-1,
                2.0375499215729606210622322285363e-1,
                2.3948076795922748948524352721738e-1,
                2.777084122893516425530494522051e-1,
                3.1803721819046939009026680612403e-1,
                3.6001049461425858810266499021131e-1,
                4.0313230365720934454217107307577e-1,
                4.4688533588396533044743048062199e-1,
                4.9074892202486175960184208080769e-1,
                5.3421624273244621420512843556551e-1,
                5.7680991887315648467558946697448e-1,
                6.1809532850689066212732057937274e-1,
                6.5769118694954119607833902714853e-1,
                6.9527712658910223178719274858914e-1,
                7.3059820810403112399810153194007e-1,
                7.6346647140619223434729583583715e-1,
                7.9375978307962807874576589516017e-1,
                8.2141835086065143946548807909925e-1,
                8.4643935184946092548777123450699e-1,
                8.6887015991697008744569530602467e-1,
                8.888006621299879509780621551118e-1,
                9.0635512895721015178123196629466e-1,
                9.216840548431236089062166869819e-1,
                9.3495632135181220532057106641933e-1,
                9.4635196135469050932241191509446e-1,
                9.5605572597464175487317824792671e-1,
                9.6425158157771945657098859282225e-1,
                9.7111819623117328617108644751399e-1,
                9.7682541667614761882980334582734e-1,
                9.8153168972123166412735062376851e-1,
                9.8538234666978576993582667929948e-1,
                9.8850864564867630921940291212761e-1,
                9.9102745350377593608390339073183e-1,
                9.9304144478493804631268746015899e-1,
                9.9463969851103798754200922796921e-1,
                9.9589858205433922192321353511889e-1,
                9.9688282385089869469211042695743e-1,
                9.9764669114390939512028724756148e-1,
                9.9823520425839310352790669980018e-1,
                9.9868533392513734338471031409171e-1,
                9.9902714213800145221464078264993e-1,
                9.9928483942460310058022171919681e-1,
                9.994777419499671021705120038048e-1,
            ];

            foreach ((NakagamiDistribution dist, ddouble[] expecteds) in new[]{
                (dist_m1omega1, expected_dist_m1omega1),
                (dist_m2omega1, expected_dist_m2omega1),
                (dist_m1omega2, expected_dist_m1omega2),
                (dist_m2omega2, expected_dist_m2omega2),
                (dist_m3omega4, expected_dist_m3omega4),
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