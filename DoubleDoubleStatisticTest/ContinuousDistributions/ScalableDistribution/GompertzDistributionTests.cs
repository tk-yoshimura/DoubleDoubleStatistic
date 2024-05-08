using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class GompertzDistributionTests {
        readonly GompertzDistribution dist_eta1theta1 = new(eta: 1, theta: 1);
        readonly GompertzDistribution dist_eta2theta1 = new(eta: 2, theta: 1);
        readonly GompertzDistribution dist_eta1theta2 = new(eta: 1, theta: 2);
        readonly GompertzDistribution dist_eta2theta2 = new(eta: 2, theta: 2);
        readonly GompertzDistribution dist_eta3theta4 = new(eta: 3, theta: 4);

        GompertzDistribution[] Dists => [
            dist_eta1theta1,
            dist_eta2theta1,
            dist_eta1theta2,
            dist_eta2theta2,
            dist_eta3theta4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (GompertzDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Eta={dist.Eta}");
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
            foreach (GompertzDistribution dist in Dists) {
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
            foreach (GompertzDistribution dist in Dists) {
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
            foreach (GompertzDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (GompertzDistribution dist in Dists) {
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
            foreach (GompertzDistribution dist in Dists) {
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
            foreach (GompertzDistribution dist in Dists) {
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
            foreach (GompertzDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (GompertzDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (GompertzDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (GompertzDistribution dist in Dists) {
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
            foreach (GompertzDistribution dist in Dists) {
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
            foreach (GompertzDistribution dist in Dists) {
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

            foreach (GompertzDistribution dist in Dists) {

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
            foreach (GompertzDistribution dist in Dists) {
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
            ddouble[] expected_dist_eta1theta1 = [
                1.0,
                0.999877298941705,
                0.9995067142063597,
                0.9988846204820263,
                0.9980075286937053,
                0.9968720946717472,
                0.9954751279125088,
                0.9938136004164351,
                0.9918846555877275,
                0.9896856171787258,
                0.9872139982610871,
                0.9844675102048033,
                0.981444071645045,
                0.9781418174157783,
                0.9745591074280632,
                0.970694535469914,
                0.9665469379035974,
                0.9621154022352542,
                0.9573992755307759,
                0.952398172650939,
                0.9471119842779214,
                0.9415408847044824,
                0.9356853393563084,
                0.9295461120173021,
                0.9231242717269351,
                0.9164211993182065,
                0.9094385935642532,
                0.9021784769012459,
                0.8946432006949029,
                0.886835450017749,
                0.8787582479041574,
                0.8704149590502526,
                0.8618092929259169,
                0.8529453062664449,
                0.8438274049118419,
                0.8344603449623726,
                0.8248492332197181,
                0.8149995268840468,
                0.8049170324783979,
                0.7946079039730736,
                0.7840786400841948,
                0.7733360807222509,
                0.7623874025683124,
                0.7512401137576392,
                0.7399020476526581,
                0.7283813556897354,
                0.7166864992868156,
                0.70482624080185,
                0.6928096335349705,
                0.6806460107706079,
                0.6683449738591679,
                0.6559163793414775,
                0.6433703251229922,
                0.6307171357086732,
                0.6179673465135394,
                0.6051316872681026,
                0.5922210645422418,
                0.579246543415504,
                0.5662193283263486,
                0.5531507431374415,
                0.5400522104587291,
                0.5269352302746738,
                0.5138113579266691,
                0.5006921815062507,
                0.487589298719261,
            ];
            ddouble[] expected_dist_eta2theta1 = [
                2.0,
                1.968509708636254,
                1.936554504700264,
                1.90415849585641,
                1.871346570167866,
                1.838144365366024,
                1.804578235511186,
                1.770675215056875,
                1.736462980340607,
                1.701969808534734,
                1.667224534102237,
                1.632256502813934,
                1.597095523395436,
                1.561771816884368,
                1.526315963790742,
                1.490758849165871,
                1.455131605697816,
                1.419465554964044,
                1.383792146984462,
                1.348142898230533,
                1.312549328258373,
                1.277042895145682,
                1.241654929923983,
                1.206416570208691,
                1.171358693240119,
                1.136511848558421,
                1.101906190544597,
                1.06757141106801,
                1.033536672488209,
                0.9998305412651931,
                0.9664809224374705,
                0.9335149952312376,
                0.9009591500667369,
                0.8688389272291629,
                0.8371789574713709,
                0.8060029048140233,
                0.7753334118055879,
                0.7451920474998048,
                0.7155992584017373,
                0.6865743226253739,
                0.6581353074958641,
                0.6302990308179209,
                0.6030810260186305,
                0.5764955113579964,
                0.5505553633839685,
                0.525272094790578,
                0.5006558368181505,
                0.4767153263135284,
                0.4534578975458534,
                0.4308894788499086,
                0.4090145941444019,
                0.387836369347039,
                0.3673565436819706,
                0.3475754858483466,
                0.3284922149915016,
                0.310104426390885,
                0.2924085217514908,
                0.275399643958413,
                0.2590717161275028,
                0.2434174847591555,
                0.2284285667772147,
                0.2140955002111085,
                0.2004077982568267,
                0.1873540064314292,
                0.1749217625166659,
            ];
            ddouble[] expected_dist_eta1theta2 = [
                0.5,
                0.4999847016308209,
                0.4999386494708525,
                0.4998616108539113,
                0.4997533571031799,
                0.4996136636645775,
                0.4994423102410132,
                0.4992390809274676,
                0.4990037643468527,
                0.4987361537865931,
                0.4984360473358736,
                0.498103248023491,
                0.4977375639562544,
                0.4973388084578653,
                0.4969068002082175,
                0.4964413633830471,
                0.4959423277938638,
                0.4954095290280943,
                0.4948428085893629,
                0.4942420140378344,
                0.4936069991305436,
                0.4929376239616306,
                0.4922337551024016,
                0.4914952657411314,
                0.4907220358225225,
                0.4899139521867338,
                0.4890709087078892,
                0.4881928064319747,
                0.4872795537140316,
                0.4863310663545483,
                0.485347267734957,
                0.4843280889521333,
                0.4832734689517987,
                0.4821833546607238,
                0.4810577011176271,
                0.4798964716026622,
                0.4786996377653879,
                0.4774671797511073,
                0.4761990863254695,
                0.4748953549972159,
                0.4735559921389607,
                0.4721810131058869,
                0.4707704423522412,
                0.4693243135455099,
                0.4678426696781542,
                0.4663255631767844,
                0.464773056008651,
                0.4631852197853286,
                0.4615621358634676,
                0.4599038954424903,
                0.4582105996591033,
                0.4564823596785001,
                0.4547192967821266,
                0.4529215424518815,
                0.4510892384506229,
                0.449222536898853,
                0.4473216003474514,
                0.4453866018463296,
                0.4434177250088745,
                0.4414151640720582,
                0.4393791239520787,
                0.4373098202954091,
                0.4352074795251263,
                0.4330723388823914,
                0.4309046464629585,
            ];
            ddouble[] expected_dist_eta2theta2 = [
                1.0,
                0.992157221923816,
                0.984254854318127,
                0.9762943663218145,
                0.968277252350132,
                0.9602050316655495,
                0.9520792479282052,
                0.9439014687259424,
                0.9356732850839329,
                0.9273963109539043,
                0.9190721826830122,
                0.9107025584624081,
                0.9022891177555928,
                0.8938335607066427,
                0.8853376075284374,
                0.8768029978710271,
                0.8682314901703035,
                0.859624860977163,
                0.8509849042673668,
                0.8423134307323338,
                0.8336122670511186,
                0.8248832551438554,
                0.8161282514069672,
                0.8073491259304696,
                0.7985477616977179,
                0.7897260537679726,
                0.7808859084421839,
                0.7720292424124173,
                0.7631579818953712,
                0.7542740617504556,
                0.7453794245829353,
                0.736476019832657,
                0.7275658028489081,
                0.7186507339519782,
                0.7097327774820219,
                0.700813900835841,
                0.6918960734922308,
                0.6829812660265584,
                0.6740714491152667,
                0.665168592531013,
                0.6562746641291863,
                0.6473916288265567,
                0.6385214475728408,
                0.6296660763159857,
                0.6208274649619915,
                0.6120075563301197,
                0.6032082851043454,
                0.5944315767819367,
                0.5856793466200596,
                0.5769534985813254,
                0.5682559242792103,
                0.5595885019242969,
                0.5509530952722984,
                0.5423515525748436,
                0.5337857055340052,
                0.525257368261577,
                0.5167683362441045,
                0.5083203853146898,
                0.4999152706325966,
                0.4915547256716879,
                0.4832404612187353,
                0.4749741643826408,
                0.4667574976156188,
                0.45859209774738,
                0.4504795750333684,
            ];
            ddouble[] expected_dist_eta3theta4 = [
                0.75,
                0.7441463991413859,
                0.73830455447297,
                0.7324747767683183,
                0.7266573753155671,
                0.7208526578733684,
                0.7150609306269667,
                0.7092824981444135,
                0.7035176633329379,
                0.6977667273954792,
                0.6920299897873929,
                0.686307748173345,
                0.6806002983843973,
                0.6749079343753085,
                0.6692309481820448,
                0.6635696298795252,
                0.6579242675396066,
                0.65229514718932,
                0.6466825527693707,
                0.6410867660929096,
                0.6355080668045927,
                0.6299467323399354,
                0.624403037884975,
                0.6188772563362483,
                0.6133696582611042,
                0.6078805118583511,
                0.6024100829192595,
                0.5969586347889255,
                0.5915264283280092,
                0.586113721874854,
                0.5807207712080071,
                0.5753478295091393,
                0.569995147326387,
                0.5646629725381201,
                0.5593515503171477,
                0.554061123095372,
                0.5487919305289058,
                0.5435442094636538,
                0.5383181939013805,
                0.5331141149662644,
                0.527932200871956,
                0.5227726768891449,
                0.5176357653136504,
                0.5125216854350403,
                0.5074306535057942,
                0.5023628827110153,
                0.4973185831387035,
                0.4922979617505983,
                0.4873012223536018,
                0.4823285655717897,
                0.4773801888190204,
                0.4724562862721516,
                0.4675570488448721,
                0.4626826641621588,
                0.4578333165353663,
                0.4530091869379594,
                0.4482104529818933,
                0.4434372888946548,
                0.438689865496968,
                0.4339683501811733,
                0.4292729068902911,
                0.4246036960977699,
                0.4199608747879359,
                0.4153445964371412,
                0.4107550109956275,
            ];

            foreach ((GompertzDistribution dist, ddouble[] expecteds) in new[]{
                (dist_eta1theta1, expected_dist_eta1theta1), (dist_eta2theta1, expected_dist_eta2theta1),
                (dist_eta1theta2, expected_dist_eta1theta2), (dist_eta2theta2, expected_dist_eta2theta2),
                (dist_eta3theta4, expected_dist_eta3theta4),
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
            ddouble[] expected_dist_eta1theta1 = [
                0.0,
                0.01562436174929982,
                0.03124487450894708,
                0.04685763660194764,
                0.0624586907589384,
                0.07804402631448248,
                0.0936095815395459,
                0.1091512461114873,
                0.1246648637226423,
                0.1401462348283387,
                0.1555911195349012,
                0.1709952406279164,
                0.1863542867407267,
                0.2016639156628011,
                0.2169197577872889,
                0.232117419696716,
                0.2472524878854099,
                0.2623205326168557,
                0.2773171119137852,
                0.2922377756783886,
                0.3070780699396066,
                0.3218335412240213,
                0.3364997410464066,
                0.3510722305155339,
                0.3655465850503533,
                0.3799183992011775,
                0.3941832915700066,
                0.4083369098236265,
                0.4223749357926031,
                0.436293090648788,
                0.4500871401534311,
                0.4637528999674844,
                0.4772862410151655,
                0.4906830948913428,
                0.5039394593027975,
                0.5170514035329219,
                0.5300150739189331,
                0.5428266993302033,
                0.5554825966358574,
                0.5679791761493481,
                0.5803129470373067,
                0.5924805226795714,
                0.6044786259669391,
                0.6163040945228451,
                0.6279538858348837,
                0.6394250822818114,
                0.6507148960414632,
                0.6618206738648165,
                0.6727399017013203,
                0.6834702091605093,
                0.6940093737948956,
                0.7043553251891527,
                0.7145061488406824,
                0.7244600898167994,
                0.7342155561739665,
                0.7437711221247834,
                0.7531255309387648,
                0.7622776975633466,
                0.7712267109520297,
                0.7799718360871175,
                0.7885125156851188,
                0.796848371573573,
                0.8049792057288185,
                0.8129050009650588,
                0.8206259212659828,
            ];
            ddouble[] expected_dist_eta2theta1 = [
                0.0,
                0.03100460281852668,
                0.06151350683481438,
                0.09151963509597505,
                0.1210162934665561,
                0.1499971825855893,
                0.1784564093230829,
                0.2063884976952842,
                0.2337883991982995,
                0.2606515025201175,
                0.2869736425916785,
                0.3127511089384337,
                0.3379806532948084,
                0.3626594964451488,
                0.3867853342560817,
                0.4103563428667705,
                0.4333711830052949,
                0.4558290034013206,
                0.4777294432673675,
                0.499072633823325,
                0.5198591988413792,
                0.5400902541912489,
                0.5597674063685145,
                0.5788927499919156,
                0.5974688642587314,
                0.6154988083507698,
                0.6329861157870484,
                0.6499347877229444,
                0.6663492851994006,
                0.6822345203497044,
                0.697595846575368,
                0.7124390477067171,
                0.7267703261679443,
                0.7405962901705391,
                0.7539239399631992,
                0.7667606531704794,
                0.7791141692565741,
                0.7909925731546836,
                0.8024042781064001,
                0.813358007759404,
                0.8238627775754894,
                0.8339278756044847,
                0.8435628426829995,
                0.8527774521200663,
                0.8615816889346372,
                0.8699857287125214,
                0.8779999161526741,
                0.8856347433747531,
                0.8929008280615385,
                0.8998088915111082,
                0.906369736674608,
                0.9125942262559883,
                0.9184932609502214,
                0.9240777578962338,
                0.9293586294200861,
                0.9343467621428073,
                0.9390529965257333,
                0.9434881069242163,
                0.9476627822181738,
                0.9515876070851257,
                0.9552730439781629,
                0.9587294158676909,
                0.9619668898018375,
                0.9649954613361154,
                0.9678249398783226,
            ];
            ddouble[] expected_dist_eta1theta2 = [
                0.0,
                0.007812420372407813,
                0.01562436174929982,
                0.0234353417798826,
                0.03124487450894708,
                0.03905247044024374,
                0.04685763660194764,
                0.05465987661422889,
                0.0624586907589384,
                0.0702535760514238,
                0.07804402631448248,
                0.08582953225446699,
                0.0936095815395459,
                0.1013836588801328,
                0.1091512461114873,
                0.1169118222784969,
                0.1246648637226423,
                0.1324098441711501,
                0.1401462348283387,
                0.1478735044691543,
                0.1555911195349012,
                0.1632985442311633,
                0.1709952406279164,
                0.1786806687618256,
                0.1863542867407267,
                0.1940155508502812,
                0.2016639156628011,
                0.2092988341482326,
                0.2169197577872889,
                0.2245261366867218,
                0.232117419696716,
                0.239693054530393,
                0.2472524878854099,
                0.2547951655676307,
                0.2623205326168557,
                0.2698280334345834,
                0.2773171119137852,
                0.2847872115706667,
                0.2922377756783886,
                0.2996682474027225,
                0.3070780699396066,
                0.3144666866545747,
                0.3218335412240213,
                0.3291780777782701,
                0.3364997410464066,
                0.3437979765028373,
                0.3510722305155339,
                0.3583219504959199,
                0.3655465850503533,
                0.3727455841331617,
                0.3799183992011775,
                0.3870644833697251,
                0.3941832915700066,
                0.4012742807078302,
                0.4083369098236265,
                0.4153706402536923,
                0.4223749357926031,
                0.429349262856728,
                0.436293090648788,
                0.4432058913233838,
                0.4500871401534311,
                0.4569363156974282,
                0.4637528999674844,
                0.4705363785980349,
                0.4772862410151655,
            ];
            ddouble[] expected_dist_eta2theta2 = [
                0.0,
                0.0155638068327405,
                0.03100460281852668,
                0.04632146831542527,
                0.06151350683481438,
                0.07657984543300134,
                0.09151963509597505,
                0.106332051116975,
                0.1210162934665561,
                0.1355715871548343,
                0.1499971825855893,
                0.1642923559019134,
                0.1784564093230829,
                0.1924886714723425,
                0.2063884976952842,
                0.220155270368515,
                0.2337883991982995,
                0.247287321508872,
                0.2606515025201175,
                0.2738804356143195,
                0.2869736425916785,
                0.2999306739143095,
                0.3127511089384337,
                0.325434556134478,
                0.3379806532948084,
                0.3503890677288244,
                0.3626594964451488,
                0.3747916663206557,
                0.3867853342560817,
                0.3986402873179791,
                0.4103563428667705,
                0.4219333486706761,
                0.4333711830052949,
                0.4446697547386251,
                0.4558290034013206,
                0.4668488992419921,
                0.4777294432673675,
                0.4884706672671376,
                0.499072633823325,
                0.5095354363040256,
                0.5198591988413792,
                0.5300440762936429,
                0.5400902541912489,
                0.5499979486667433,
                0.5597674063685145,
                0.5693989043582292,
                0.5788927499919156,
                0.5882492807846393,
                0.5974688642587314,
                0.6065518977755515,
                0.6154988083507698,
                0.6243100524531779,
                0.6329861157870484,
                0.641527513058074,
                0.6499347877229444,
                0.6582085117226224,
                0.6663492851994006,
                0.6743577361978403,
                0.6822345203497044,
                0.6899803205430125,
                0.697595846575368,
                0.7050818347917167,
                0.7124390477067171,
                0.7196682736119165,
                0.7267703261679443,

            ];
            ddouble[] expected_dist_eta3theta4 = [
                0.0,
                0.01167300363867851,
                0.02325463620357604,
                0.03474508381498687,
                0.04614453743749658,
                0.05745319285642636,
                0.06867125065358981,
                0.079798916182369,
                0.09083639954210254,
                0.1017839155517958,
                0.1126416837231543,
                0.1234099282329386,
                0.1340888778946526,
                0.1446787661295548,
                0.155179830937012,
                0.1655923148641886,
                0.175916464975074,
                0.1861525328188599,
                0.1963007743976637,
                0.2063614501336095,
                0.2163348248352638,
                0.2262211676634338,
                0.2360207520963348,
                0.2457338558941302,
                0.2553607610628452,
                0.2649017538176678,
                0.2743571245456333,
                0.2837271677677052,
                0.2930121821002534,
                0.3022124702159404,
                0.3113283388040137,
                0.3203600985300225,
                0.3293080639949512,
                0.3381725536937884,
                0.3469538899735313,
                0.3556523989906354,
                0.3642684106679125,
                0.3728022586508927,
                0.3812542802636458,
                0.3896248164640825,
                0.3979142117987318,
                0.4061228143570118,
                0.4142509757249944,
                0.4222990509386785,
                0.4302673984367749,
                0.4381563800130144,
                0.4459663607679872,
                0.4536977090605228,
                0.4613507964586163,
                0.4689259976899174,
                0.4764236905917812,
                0.4838442560609006,
                0.4911880780025215,
                0.4984555432792562,
                0.5056470416595017,
                0.5127629657654732,
                0.5198037110208658,
                0.526769675598149,
                0.5336612603655102,
                0.5404788688334546,
                0.5472229071010708,
                0.5538937838019776,
                0.5604919100499575,
                0.5670176993842911,
                0.5734715677148023,
            ];

            foreach ((GompertzDistribution dist, ddouble[] expecteds) in new[]{
                (dist_eta1theta1, expected_dist_eta1theta1), (dist_eta2theta1, expected_dist_eta2theta1),
                (dist_eta1theta2, expected_dist_eta1theta2), (dist_eta2theta2, expected_dist_eta2theta2),
                (dist_eta3theta4, expected_dist_eta3theta4),
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