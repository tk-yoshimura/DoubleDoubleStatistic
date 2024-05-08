
using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.LinearityDistribution {
    [TestClass()]
    public class JohnsonSUDistributionTests {
        readonly JohnsonSUDistribution dist_gamma1delta1mu0sigma1 = new(gamma: 1, delta: 1, mu: 0, sigma: 1);
        readonly JohnsonSUDistribution dist_gamma2delta1mu1sigma1 = new(gamma: 2, delta: 1, mu: 1, sigma: 1);
        readonly JohnsonSUDistribution dist_gamma1delta1mu0sigma2 = new(gamma: 1, delta: 1, mu: 0, sigma: 2);
        readonly JohnsonSUDistribution dist_gamma2delta2mu2sigma2 = new(gamma: 2, delta: 2, mu: 2, sigma: 2);
        readonly JohnsonSUDistribution dist_gamma3delta2mu0sigma4 = new(gamma: 3, delta: 2, mu: 0, sigma: 4);

        JohnsonSUDistribution[] Dists => [
            dist_gamma1delta1mu0sigma1,
            dist_gamma2delta1mu1sigma1,
            dist_gamma1delta1mu0sigma2,
            dist_gamma2delta2mu2sigma2,
            dist_gamma3delta2mu0sigma4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (JohnsonSUDistribution dist in Dists) {
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
            foreach (JohnsonSUDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (!ddouble.IsFinite(dist.Mean)) {
                    continue;
                }

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (JohnsonSUDistribution dist in Dists) {
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
            foreach (JohnsonSUDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (JohnsonSUDistribution dist in Dists) {
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
            foreach (JohnsonSUDistribution dist in Dists) {
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
            foreach (JohnsonSUDistribution dist in Dists) {
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
            foreach (JohnsonSUDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (JohnsonSUDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (JohnsonSUDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -2; x <= 6; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (JohnsonSUDistribution dist in Dists) {
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
            foreach (JohnsonSUDistribution dist in Dists) {
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
            foreach (JohnsonSUDistribution dist in Dists) {
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

            foreach (JohnsonSUDistribution dist in Dists) {

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
            foreach (JohnsonSUDistribution dist in Dists) {
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
                0.05314437641077076,
                0.05484007252460024,
                0.05660226486053547,
                0.05843398302103565,
                0.06033840534172218,
                0.0623188655628945,
                0.06437885957431538,
                0.066522052177893,
                0.06875228379832934,
                0.07107357705425225,
                0.07349014308121059,
                0.07600638747249455,
                0.07862691567321582,
                0.08135653762644864,
                0.08420027142632254,
                0.08716334568036031,
                0.09025120022041165,
                0.09346948472626636,
                0.09682405473610721,
                0.1003209644106239,
                0.1039664552896147,
                0.1077669401274576,
                0.1117289807125072,
                0.1158592583601675,
                0.1201645355142129,
                0.1246516065892589,
                0.1293272358316748,
                0.1341980795585914,
                0.1392705896464032,
                0.144550894572611,
                0.1500446536599219,
                0.1557568794227652,
                0.1616917220707499,
                0.1678522092838868,
                0.1742399333531991,
                0.1808546767053998,
                0.1876939657525031,
                0.1947525420115569,
                0.2020217386605047,
                0.2094887503376029,
                0.2171357843544447,
                0.2249390830071634,
                0.232867809937783,
                0.240882799333195,
                0.2489351762233012,
                0.2569648706016506,
                0.2648990691492326,
                0.2726506777538798,
                0.2801169074311687,
                0.2871781466907388,
                0.2936973443763202,
                0.2995201952120277,
                0.3044764876961689,
                0.3083830257512602,
                0.3110485481706711,
                0.3122810110511727,
                0.3118974303530713,
                0.3097361710927577,
                0.3056711049333472,
                0.299626473242224,
                0.2915906897941372,
                0.2816268701029246,
                0.2698777972646877,
                0.2565635057333895,
                0.2419707245191434,
                0.2264348933809417,
                0.2103169671138663,
                0.1939782915918822,
                0.1777571243349117,
                0.1619498104054871,
                0.1467984408159261,
                0.1324854261197197,
                0.1191342173911201,
                0.1068146518383341,
                0.09555114289842805,
                0.0853320738316758,
                0.07611912192404223,
                0.06785567856843354,
                0.0604739289296672,
                0.0539004596397759,
                0.04806046321517093,
                0.04288071780994945,
                0.03829156434101754,
                0.03422810372567122,
                0.03063081415546914,
                0.0274457554601166,
                0.02462449306726796,
                0.02212384238466871,
                0.01990550760738313,
                0.01793566743534953,
                0.01618454360873093,
                0.01462597580720443,
                0.01323701750170894,
                0.01199756102981537,
                0.01088999583457494,
                0.009898900931795592,
                0.009010770839949748,
                0.008213773108745574,
                0.007497534986742677,
                0.006852956508734092,
                0.006272047241548391,
                0.005747784019720015,
                0.005273987173374141,
                0.004845212961530917,
                0.004456660149865219,
                0.004104088896935682,
                0.003783750327485411,
                0.003492325370322975,
                0.00322687161902706,
                0.002984777135595675,
                0.0027637202606121,
                0.002561634619606065,
                0.002376678625473708,
                0.002207208872646858,
                0.002051756901748002,
                0.001909008885244477,
                0.001777787846530446,
                0.001657038078196801,
                0.001545811471142312,
                0.001443255505647496,
                0.001348602689464534,
                0.001261161257147107,
                0.001180306969925029,
                0.001105475877000011,
                0.001036157917700773,
                9.718912599170589e-4,
                9.122572840010579e-4,
                8.568761331968794e-4,
                8.054027619037526e-4,
            ];
            ddouble[] expected_dist_gamma2delta1mu1sigma1 = [
                0.07451195767556179,
                0.07570358410470811,
                0.07691681602136685,
                0.07815196818295854,
                0.07940934780474468,
                0.08068925298562177,
                0.08199197096270665,
                0.08331777617740335,
                0.08466692813387276,
                0.0860396690288685,
                0.0874362211297418,
                0.08885678387502909,
                0.0903015306693962,
                0.09177060534179576,
                0.09326411823247253,
                0.09478214187089024,
                0.09632470620272422,
                0.09789179331972446,
                0.09948333164146715,
                0.1010991894927382,
                0.1027391680144796,
                0.1044029933398364,
                0.1060903079598196,
                0.1078006611953931,
                0.109533498684357,
                0.1112881507821934,
                0.113063819766001,
                0.1148595657197714,
                0.1166742909674959,
                0.1185067229079655,
                0.1203553950916456,
                0.1222186263657468,
                0.1240944978986849,
                0.1259808278797272,
                0.1278751436740351,
                0.1297746511979697,
                0.1316762012650004,
                0.1335762526396852,
                0.1354708315270695,
                0.1373554872190088,
                0.1392252436193363,
                0.1410745463791588,
                0.142897205395352,
                0.144686332464122,
                0.1464342739432464,
                0.1481325383689983,
                0.1497717191066474,
                0.1513414122994568,
                0.1528301306362231,
                0.1542252138017655,
                0.1555127369334042,
                0.1566774190102223,
                0.1577025338881478,
                0.1585698277070161,
                0.1592594476870628,
                0.1597498889586564,
                0.1600179680897835,
                0.1600388344460405,
                0.1597860334771627,
                0.1592316394770986,
                0.1583464792509255,
                0.1571004722708993,
                0.1554631169704607,
                0.1534041562032656,
                0.1508944566114945,
                0.1479071352494146,
                0.1444189602653784,
                0.1404120381437031,
                0.1358757749221552,
                0.1308090600062746,
                0.1252225668990464,
                0.1191409963514129,
                0.1126050101634906,
                0.1056725315483845,
                0.09841904259769778,
                0.09093651967346912,
                0.083330742038146,
                0.07571690352876552,
                0.06821373931154867,
                0.06093669880932991,
                0.05399096651318806,
                0.0474652600411029,
                0.04142725564199259,
                0.0359212075676222,
                0.03096791557218093,
                0.02656677586657485,
                0.02269934094370579,
                0.01933367771715933,
                0.01642884957328757,
                0.01393900367610896,
                0.01181674916835863,
                0.01001570432493076,
                0.008492235932634698,
                0.00720650212283474,
                0.00612294763308637,
                0.005210402648669009,
                0.004441918072677757,
                0.003794443202234228,
                0.003248424106548123,
                0.00278737665598612,
                0.002397468737751924,
                0.002067131793099192,
                0.001786711759273881,
                0.001548162881460627,
                0.001344783778083444,
                0.001170992835711807,
                0.001022138885646683,
                8.943427396898709e-4,
                7.843652351195557e-4,
                6.894977559214261e-4,
                6.074716286461099e-4,
                5.363832569431217e-4,
                4.74632312620599e-4,
                4.20870718587777e-4,
                3.739605294736561e-4,
                3.329391364439804e-4,
                2.969904956961059e-4,
                2.65421309554323e-4,
                2.376412802193833e-4,
                2.131467143935133e-4,
                1.915068872635315e-4,
                1.723526812357366e-4,
                1.553671023551757e-4,
                1.402773489092932e-4,
                1.268481651751741e-4,
                1.148762610056753e-4,
                1.041856169334276e-4,
                9.462352632196521e-5,
                8.605725213440126e-5,
            ];
            ddouble[] expected_dist_gamma1delta1mu0sigma2 = [
                0.08084586103537496,
                0.08237181144122903,
                0.08392610464194342,
                0.08550881651940488,
                0.08711996667659956,
                0.08875951178860342,
                0.09042733835269991,
                0.09212325479336618,
                0.09384698287625154,
                0.09559814838398158,
                0.09737627100577845,
                0.09918075339263104,
                0.1010108693302524,
                0.1028657509835243,
                0.1047443751688014,
                0.1066455486146043,
                0.1085678921772224,
                0.110509823985953,
                0.1124695415035817,
                0.1144450025017732,
                0.1164339049688915,
                0.118433665990044,
                0.1204413996665975,
                0.1224538941758334,
                0.1244675881116506,
                0.1264785462951877,
                0.1284824353008253,
                0.1304744990091204,
                0.1324495345746163,
                0.1344018692837817,
                0.1363253388769399,
                0.1382132680179055,
                0.1400584537155843,
                0.1418531526316593,
                0.1435890733453694,
                0.1452573747867221,
                0.1468486721881601,
                0.1483530520348183,
                0.1497600976060138,
                0.1510589267841618,
                0.1522382438480845,
                0.1532864069494623,
                0.1541915128756301,
                0.1549415005092451,
                0.1555242740853355,
                0.1559278468998671,
                0.1561405055255863,
                0.1561509938309785,
                0.1559487151765357,
                0.1555239500918541,
                0.1548680855463788,
                0.1539738506640596,
                0.1528355524666736,
                0.1514493040513409,
                0.149813236621112,
                0.1479276861104084,
                0.1457953448970686,
                0.1434213693748195,
                0.1408134350514623,
                0.1379817323731943,
                0.1349388986323438,
                0.1316998840081411,
                0.1282817528666947,
                0.1247034246991929,
                0.1209853622595717,
                0.1171492173140128,
                0.1132174466904708,
                0.1092129128190631,
                0.1051584835569331,
                0.1010766457545623,
                0.09698914579594109,
                0.09291666836649924,
                0.08887856216745584,
                0.08489261843733815,
                0.08097490520274354,
                0.07713965738381844,
                0.07339922040796304,
                0.06976404296704326,
                0.06624271305985985,
                0.06284203050865701,
                0.05956710869556007,
                0.0564214982673398,
                0.05340732591916703,
                0.05052544199578735,
                0.04777557144921402,
                0.0451564635821192,
                0.0426660369158379,
                0.04030151639729638,
                0.03805956096202111,
                0.03593638017683806,
                0.03392783928421677,
                0.03202955245826006,
                0.0302369644648336,
                0.02854542120461444,
                0.02695022981988795,
                0.02544670917680772,
                0.02403023160758546,
                0.02269625682392972,
                0.02144035890497472,
                0.0202582472294144,
                0.01914578217050877,
                0.01809898631050014,
                0.01711405186283561,
                0.01618734492031582,
                0.01531540707773457,
                0.01449495491075815,
                0.0137228777300583,
                0.01299623397187273,
                0.01231224653363398,
                0.01166829731618906,
                0.01106192119233436,
                0.01049079958465888,
                0.009952753803691566,
                0.009445738269689946,
                0.008967833717674765,
                0.008517240465091818,
                0.008092271804365465,
                0.007691347568216886,
                0.007312987903602214,
                0.006955807280159813,
                0.006618508750854469,
                0.006299878475810123,
                0.005998780514907684,
                0.005714151890389594,
                0.005444997917287471,
                0.005190387796823006,
                0.004949450465897796,
                0.004721370694275364,
                0.004505385419974874,
            ];
            ddouble[] expected_dist_gamma2delta2mu2sigma2 = [
                0.03304317322703193,
                0.03445094957335294,
                0.03591820577168132,
                0.03744725195854678,
                0.0390404616736203,
                0.04070027045093105,
                0.04242917393363534,
                0.04422972545468772,
                0.04610453302050303,
                0.04805625562908846,
                0.05008759884816138,
                0.05220130957245148,
                0.05440016987272078,
                0.05668698984203642,
                0.05906459933751652,
                0.06153583850817176,
                0.06410354699162767,
                0.06677055165448746,
                0.06953965274296266,
                0.07241360830224752,
                0.07539511671505912,
                0.0784867972019599,
                0.08169116811868922,
                0.08501062287897442,
                0.08844740332542149,
                0.09200357036640584,
                0.09568097169374598,
                0.09948120639477524,
                0.1034055862737061,
                0.1074550937014853,
                0.1116303358213064,
                0.1159314949493286,
                0.1203582750278033,
                0.1249098440116674,
                0.1295847721008274,
                0.1343809657700012,
                0.1392955975974418,
                0.1443250319545739,
                0.1494647466920955,
                0.1547092510460913,
                0.1600520000919456,
                0.1654853061961099,
                0.1710002480579286,
                0.1765865780975032,
                0.1822326291326541,
                0.1879252214998276,
                0.1936495720113728,
                0.1993892064055151,
                0.205125877235417,
                0.2108394894588088,
                0.2165080363274201,
                0.2221075485318881,
                0.227612059927077,
                0.2329935935364941,
                0.238222171901651,
                0.2432658561882996,
                0.2480908187682619,
                0.2526614542405282,
                0.2569405340113139,
                0.2608894095878371,
                0.2644682696178788,
                0.2676364553853322,
                0.2703528389057849,
                0.2725762669081386,
                0.27426607279045,
                0.2753826570552977,
                0.2758881347234349,
                0.2757470457674368,
                0.2749271216907427,
                0.2734000980186105,
                0.2711425587163114,
                0.2681367944987322,
                0.2643716527867872,
                0.2598433528985367,
                0.2545562361953874,
                0.2485234176534982,
                0.2417673030653505,
                0.2343199351981598,
                0.2262231331535808,
                0.2175283922661898,
                0.2082965174483851,
                0.1985969711062789,
                0.1885069275949651,
                0.1781100393902378,
                0.1674949351852833,
                0.1567534861399607,
                0.1459788924115577,
                0.135263656569164,
                0.1246975221347194,
                0.1143654629537409,
                0.1043458112590405,
                0.0947086084376527,
                0.08551425248225555,
                0.07681250037986624,
                0.06864186338153726,
                0.06102940987753414,
                0.05399096651318806,
                0.0475316853834861,
                0.04164692564878913,
                0.03632338332088191,
                0.03154039427905236,
                0.02727133309536031,
                0.02348503359798534,
                0.02014716532941844,
                0.01722151181875626,
                0.0146711103648372,
                0.01245922733281869,
                0.01055015651788074,
                0.008909839958095814,
                0.007506320078671162,
                0.00631003896825835,
                0.005294004969958317,
                0.004433848880918533,
                0.003707792296899874,
                0.003096549471741304,
                0.002583181941502292,
                0.00215292249628981,
                0.001792982203468513,
                0.001492351346203048,
                0.001241602513821958,
                0.001032701768698387,
                8.588318658728336e-4,
                7.142299228761846e-4,
                5.940407066008621e-4,
                4.94185783374912e-4,
                4.112481214141203e-4,
                3.423712936664346e-4,
                2.851721582686586e-4,
                2.376657528077389e-4,
            ];
            ddouble[] expected_dist_gamma3delta2mu0sigma4 = [
                0.06560806938241914,
                0.06431661935126193,
                0.06300539684991181,
                0.06167571542716004,
                0.06032894661256766,
                0.05896651790521011,
                0.05758991044392066,
                0.05620065635196696,
                0.05480033575086826,
                0.05339057344004827,
                0.05197303524122537,
                0.05054942400885896,
                0.04912147531060228,
                0.04769095278454014,
                0.04625964318300981,
                0.04482935111599251,
                0.04340189351040255,
                0.04197909380506167,
                0.04056277590469992,
                0.03915475791992887,
                0.03775684572374984,
                0.03637082635873874,
                0.03499846133253919,
                0.03364147984263746,
                0.032301571974528,
                0.0309803819202391,
                0.02967950126671242,
                0.02840046240564413,
                0.02714473211803831,
                0.02591370538782333,
                0.02470869949938172,
                0.02353094847368208,
                0.02238159789683236,
                0.0212617001932529,
                0.02017221039326664,
                0.01911398244170731,
                0.01808776609014933,
                0.01709420441058211,
                0.01613383196281705,
                0.01520707364167679,
                0.01431424422314377,
                0.01345554862122237,
                0.01263108285940467,
                0.01184083575244108,
                0.01108469128574265,
                0.01036243167132588,
                0.009673741050911998,
                0.009018209808769335,
                0.008395339449301398,
                0.007804547987393077,
                0.007245175793280453,
                0.006716491828342952,
                0.006217700203848622,
                0.00574794699141072,
                0.005306327211809304,
                0.004891891927939501,
                0.004503655367983929,
                0.004140602006454642,
                0.003801693533465715,
                0.003485875646405734,
                0.003192084602979938,
                0.002919253480260331,
                0.002666318090774451,
                0.002432222513622474,
                0.002215924205969004,
                0.002016398667839375,
                0.001832643640789383,
                0.00166368282854795,
                0.00150856913500151,
                0.001366387421758873,
                0.001236256793886923,
                0.001117332428142686,
                0.001008806963070809,
                9.099114746351923e-4,
                8.19916064580798e-4,
                7.381300914695941e-4,
                6.63902076317029e-4,
                5.966193160040151e-4,
                5.357072382012299e-4,
                4.806285314772899e-4,
                4.308820836386522e-4,
                3.860017602421013e-4,
                3.455550537086762e-4,
                3.091416316298251e-4,
                2.763918107695425e-4,
                2.469649810036877e-4,
                2.205480010675008e-4,
                1.96853585565998e-4,
                1.756187002933705e-4,
                1.566029805529317e-4,
                1.395871849067163e-4,
                1.24371694643926e-4,
                1.107750672629695e-4,
                9.86326504289557e-5,
                8.779526120702326e-5,
                7.812793388628677e-5,
                6.950873839925958e-5,
                6.182767020331486e-5,
                5.498561151698361e-5,
                4.889336298518484e-5,
                4.347074417270451e-5,
                3.864576074211706e-5,
                3.435383574803635e-5,
                3.053710216106716e-5,
                2.714375350925392e-5,
                2.412744937978919e-5,
                2.144677244753199e-5,
                1.906473367854332e-5,
                1.694832238636573e-5,
                1.506809788707426e-5,
                1.339781959811728e-5,
                1.191411254850891e-5,
                1.059616540781312e-5,
                9.425458293250352e-6,
                8.38551777362129e-6,
                7.461696651813303e-6,
                6.640976271317014e-6,
                5.911789253927374e-6,
                5.263860733651614e-6,
                4.688066304267357e-6,
                4.176305043840582e-6,
                3.721386118020534e-6,
                3.316927594563985e-6,
                2.957266224013739e-6,
                2.63737705566769e-6,
                2.352801863979035e-6,
                2.099585458492383e-6,
                1.874219040607654e-6,
                1.673589853199086e-6,
            ];

            foreach ((JohnsonSUDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1delta1mu0sigma1, expected_dist_gamma1delta1mu0sigma1),
                (dist_gamma2delta1mu1sigma1, expected_dist_gamma2delta1mu1sigma1),
                (dist_gamma1delta1mu0sigma2, expected_dist_gamma1delta1mu0sigma2),
                (dist_gamma2delta2mu2sigma2, expected_dist_gamma2delta2mu2sigma2),
                (dist_gamma3delta2mu0sigma4, expected_dist_gamma3delta2mu0sigma4),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 1d / 16, i++) {
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
                0.1368212949397061,
                0.1401954702864504,
                0.1436776892445191,
                0.1472722017454186,
                0.1509834516842389,
                0.1548160864234584,
                0.15877496671626,
                0.1628651770522797,
                0.167092036424804,
                0.1714611095135309,
                0.1759782182708975,
                0.1806494538924251,
                0.1854811891422286,
                0.1904800909934326,
                0.1956531335293476,
                0.2010076110343296,
                0.2065511511827421,
                0.2122917282096114,
                0.2182376759165968,
                0.2243977003307766,
                0.230780891790291,
                0.2373967361786895,
                0.2442551249672205,
                0.2513663636493332,
                0.2587411780620211,
                0.2663907179816289,
                0.27432655725424,
                0.2825606895691407,
                0.2911055188039592,
                0.2999738426572499,
                0.309178828033281,
                0.3187339763489517,
                0.3286530765880972,
                0.3389501435279539,
                0.349639338100691,
                0.3607348663253451,
                0.3722508526503187,
                0.384201182886063,
                0.396599311190717,
                0.409458024817771,
                0.4227891595793896,
                0.4366032582795994,
                0.4509091638179157,
                0.4657135383899347,
                0.4810203004095078,
                0.4968299717154372,
                0.5131389296650191,
                0.5299385623240649,
                0.5472143307135505,
                0.564944750628946,
                0.5831003185945148,
                0.6016424226258102,
                0.6205222988799151,
                0.6396801195020793,
                0.6590443233655048,
                0.6785313265330044,
                0.6980457675087095,
                0.7174814457660572,
                0.7367230911911224,
                0.7556490481007765,
                0.7741348656761261,
                0.7920576609544854,
                0.8093009773257653,
                0.8257597299399906,
                0.8413447460685429,
                0.855986405749447,
                0.8696369806313666,
                0.8822714441303069,
                0.8938867448296337,
                0.9044997448803976,
                0.9141441786000402,
                0.9228670570387448,
                0.9307249319541754,
                0.9377803584175477,
                0.9440987896228247,
                0.9497460289097825,
                0.9547862722097756,
                0.9592807082363665,
                0.9632866044204519,
                0.9668567893929036,
                0.9700394413761244,
                0.9728780999860734,
                0.9754118317479729,
                0.9776754936952433,
                0.9797000527115418,
                0.9815129297700658,
                0.9831383475794161,
                0.9845976674149327,
                0.9859097063502997,
                0.9870910300297695,
                0.9881562188596007,
                0.9891181073353872,
                0.9899879973980107,
                0.9907758474133854,
                0.9914904387439111,
                0.9921395220295429,
                0.9927299453011981,
                0.9932677659638881,
                0.9937583485497249,
                0.9942064499777055,
                0.9946162938850212,
                0.9949916354244199,
                0.9953358177604135,
                0.9956518213474735,
                0.9959423069374667,
                0.9962096531418286,
                0.996455989266009,
                0.9966832240386864,
                0.9968930707750913,
                0.9970870694413244,
                0.9972666060236393,
                0.9974329295521605,
                0.997587167081385,
                0.9977303368891222,
                0.9978633601204149,
                0.9979870710727031,
                0.9981022262923838,
                0.9982095126304016,
                0.9983095543850894,
                0.9984029196437173,
                0.9984901259197421,
                0.99857164517025,
                0.9986479082672776,
                0.9987193089873421,
                0.9987862075754143,
                0.998848933932541,
                0.9989077904702318,
                0.998963054669426,
                0.9990149813772521,
            ];
            ddouble[] expected_dist_gamma2delta1mu1sigma1 = [
                0.3773537077931407,
                0.3820478317756101,
                0.3868171059247931,
                0.3916628904552752,
                0.396586565026067,
                0.4015895282208766,
                0.4066731969247261,
                0.4118390055856859,
                0.4170884053493535,
                0.4224228630524564,
                0.4278438600605791,
                0.4333528909334839,
                0.4389514618998216,
                0.4446410891211731,
                0.4504232967233138,
                0.4562996145703411,
                0.462271575754808,
                0.4683407137742618,
                0.4745085593615497,
                0.4807766369329052,
                0.4871464606141321,
                0.4936195298011349,
                0.5001973242065363,
                0.5068812983391828,
                0.5136728753578644,
                0.5205734402345801,
                0.5275843321560549,
                0.5347068360849612,
                0.5419421733943145,
                0.549291491479788,
                0.5567558522451348,
                0.56433621934549,
                0.572033444061993,
                0.5798482496688927,
                0.5877812141410468,
                0.5958327510355083,
                0.6040030883657357,
                0.6122922452709469,
                0.6207000062663932,
                0.6292258928430818,
                0.6378691321680412,
                0.6466286226190829,
                0.6555028958718221,
                0.6644900752423997,
                0.6735878299781273,
                0.682793325181837,
                0.6921031670562487,
                0.7015133431651028,
                0.7110191574319269,
                0.7206151596400393,
                0.7302950692650905,
                0.7400516935721412,
                0.7498768400532264,
                0.7597612234812554,
                0.769694368127755,
                0.7796645060544707,
                0.7896584728650224,
                0.7996616029190748,
                0.8096576267972688,
                0.8196285747914598,
                0.8295546914109909,
                0.8394143673642953,
                0.8491840972028062,
                0.8588384727784882,
                0.868350224796081,
                0.8776903268898966,
                0.8868281785677599,
                0.8957318846409925,
                0.90436864882102,
                0.9127052972421438,
                0.9207089428437314,
                0.9283477928553945,
                0.9355920883250854,
                0.9424151465595042,
                0.9487944554619797,
                0.9547127456137721,
                0.9601589459497797,
                0.9651289178918969,
                0.969625866872711,
                0.973660353410255,
                0.9772498680518208,
                0.9804179895613535,
                0.9831932024817979,
                0.9856074954042012,
                0.9876948841694196,
                0.9894900004803089,
                0.9910268593754686,
                0.9923378778923654,
                0.9934531731739068,
                0.994400130550006,
                0.9952032060187053,
                0.9958839140973403,
                0.9964609492239136,
                0.9969503933833975,
                0.9973659709915532,
                0.9977193215918378,
                0.9980202698676837,
                0.9982770799562912,
                0.9984966867990711,
                0.9986849013730875,
                0.9988465893860865,
                0.9989858247013392,
                0.9991060196825068,
                0.9992100350566709,
                0.9993002719716275,
                0.9993787488053599,
                0.9994471650628542,
                0.9995069544277084,
                0.9995593287601302,
                0.9996053145700733,
                0.9996457832553338,
                0.9996814761838464,
                0.999713025517746,
                0.9997409715223989,
                0.9997657769739055,
                0.9997878391704607,
                0.9998074999633946,
                0.999825054149808,
                0.9998407565079404,
                0.9998548277065062,
                0.999867459278334,
                0.9998788178151178,
                0.9998890485126373,
                0.9998982781732895,
                0.9999066177543278,
                0.9999141645350499,
                0.9999210039637332,
                0.9999272112348715,
                0.9999328526388316,
            ];
            ddouble[] expected_dist_gamma1delta1mu0sigma2 = [
                0.3286530765880972,
                0.3337534815174116,
                0.3389501435279539,
                0.3442448366924565,
                0.349639338100691,
                0.3551354241372714,
                0.3607348663253451,
                0.3664394266972394,
                0.3722508526503187,
                0.3781708712433955,
                0.384201182886063,
                0.3903434543703028,
                0.396599311190717,
                0.4029703290967909,
                0.409458024817771,
                0.4160638458981413,
                0.4227891595793896,
                0.4296352406619055,
                0.4366032582795994,
                0.4436942615193514,
                0.4509091638179157,
                0.4582487260706723,
                0.4657135383899347,
                0.473304000455733,
                0.4810203004095078,
                0.4888623922514116,
                0.4968299717154372,
                0.504922450613952,
                0.5131389296650191,
                0.5214781698427997,
                0.5299385623240649,
                0.5385180971430764,
                0.5472143307135505,
                0.5560243524306954,
                0.564944750628946,
                0.5739715782423359,
                0.5831003185945148,
                0.5923258518339268,
                0.6016424226258102,
                0.6110436098150125,
                0.6205222988799151,
                0.6300706581048561,
                0.6396801195020793,
                0.6493413656089742,
                0.6590443233655048,
                0.6687781663323151,
                0.6785313265330044,
                0.6882915171846544,
                0.6980457675087095,
                0.707780470679946,
                0.7174814457660572,
                0.7271340142281841,
                0.7367230911911224,
                0.7462332912534184,
                0.7556490481007765,
                0.7649547466268015,
                0.7741348656761261,
                0.7831741289365202,
                0.7920576609544854,
                0.8007711447725423,
                0.8093009773257653,
                0.8176344185266322,
                0.8257597299399906,
                0.8336662991216824,
                0.8413447460685429,
                0.8487870087915613,
                0.855986405749447,
                0.862937673724256,
                0.8696369806313666,
                0.8760819136750235,
                0.8822714441303069,
                0.8882058708010441,
                0.8938867448296337,
                0.8993167789914029,
                0.9044997448803976,
                0.9094403614874687,
                0.9141441786000402,
                0.9186174582406107,
                0.9228670570387448,
                0.9269003120326261,
                0.9307249319541754,
                0.9343488955963302,
                0.9377803584175477,
                0.9410275681264865,
                0.9440987896228247,
                0.9470022393564154,
                0.9497460289097825,
                0.9523381174077926,
                0.9547862722097756,
                0.957098037238057,
                0.9592807082363665,
                0.9613413142250294,
                0.9632866044204519,
                0.965123039907861,
                0.9668567893929036,
                0.9684937284046883,
                0.9700394413761244,
                0.9714992260837474,
                0.9728780999860734,
                0.9741808080550464,
                0.9754118317479729,
                0.9765753988166094,
                0.9776754936952433,
                0.9787158682504349,
                0.9797000527115418,
                0.9806313666333265,
                0.9815129297700658,
                0.9823476727649453,
                0.9831383475794161,
                0.9838875376049869,
                0.9845976674149327,
                0.9852710121259621,
                0.9859097063502997,
                0.9865157527271888,
                0.9870910300297695,
                0.9876373008488626,
                0.9881562188596007,
                0.9886493356802771,
                0.9891181073353872,
                0.9895639003367575,
                0.9899879973980107,
                0.9903916027985034,
                0.9907758474133854,
                0.9911417934266377,
                0.9914904387439111,
                0.991822721121766,
                0.9921395220295429,
                0.9924416702596129,
                0.9927299453011981,
            ];
            ddouble[] expected_dist_gamma2delta2mu2sigma2 = [
                0.05082642293313729,
                0.05293531038989864,
                0.05513403074864479,
                0.05742637335872176,
                0.05981627394054333,
                0.06230781850660116,
                0.06490524717982682,
                0.06761295787774851,
                0.07043550982712357,
                0.07337762686962795,
                0.0764442005147039,
                0.0796402926908296,
                0.08297113814120867,
                0.08644214640419944,
                0.09005890331266794,
                0.0938271719398566,
                0.09775289291227687,
                0.1018421840025578,
                0.1061013389071057,
                0.1105368251048338,
                0.1151552806841269,
                0.1199635100156085,
                0.1249684781382141,
                0.1301773037155504,
                0.1355972504086218,
                0.1412357164997581,
                0.1471002225910915,
                0.1531983971893097,
                0.1595379599767829,
                0.1661267025577175,
                0.1729724664569116,
                0.1800831181382525,
                0.1874665208005816,
                0.1951305027003497,
                0.2030828217439725,
                0.2113311260885106,
                0.2198829104877759,
                0.2287454681228852,
                0.2379258376623846,
                0.2474307453082129,
                0.2572665416009063,
                0.2674391327816617,
                0.2779539065413369,
                0.2888156520284747,
                0.3000284740414152,
                0.3115957013949837,
                0.3235197895317043,
                0.3358022175426369,
                0.3484433798753712,
                0.3614424731380931,
                0.3747973785604202,
                0.3885045408452144,
                0.4025588443417985,
                0.4169534876904925,
                0.4316798583311022,
                0.4467274085330682,
                0.4620835348905821,
                0.4777334635289365,
                0.4936601435840222,
                0.5098441518386846,
                0.5262636117189164,
                0.5428941301584581,
                0.5597087561184704,
                0.5766779647827762,
                0.5937696716190104,
                0.6109492805791614,
                0.628179770684181,
                0.6454218250690903,
                0.6626340062288831,
                0.6797729806730917,
                0.6967937954419233,
                0.7136502079374697,
                0.7302950692650905,
                0.7466807597587175,
                0.7627596735898792,
                0.7784847473620602,
                0.7938100254194574,
                0.8086912523262608,
                0.8230864806987851,
                0.8369566804225136,
                0.8502663334052639,
                0.8629839965672599,
                0.875082814914794,
                0.8865409664437728,
                0.8973420214051452,
                0.9074752002253418,
                0.9169355171398654,
                0.92572380031959,
                0.9338465838141333,
                0.9413158717849808,
                0.9481487809549093,
                0.9543670726036279,
                0.959996590408357,
                0.9650666245877274,
                0.9696092258334896,
                0.9736584941705005,
                0.9772498680518208,
                0.9804194376858141,
                0.9832033039509132,
                0.9856370005403,
                0.9877549925343031,
                0.989590259805055,
                0.9911739688954802,
                0.9925352326201682,
                0.9937009528735569,
                0.9946957391748423,
                0.9955418934079919,
                0.9962594500206556,
                0.9968662605484144,
                0.9973781116021637,
                0.9978088662413532,
                0.998170619792969,
                0.9984738625144066,
                0.9987276429078887,
                0.9989397268722745,
                0.9991167491513204,
                0.9992643546591272,
                0.9993873282103389,
                0.9994897119496702,
                0.9995749103707928,
                0.9996457832553338,
                0.9997047271702011,
                0.9997537463587386,
                0.9997945139708729,
                0.9998284246200725,
                0.999856639248545,
                0.9998801232416353,
                0.9998996786698982,
                0.9999159714621576,
            ];
            ddouble[] expected_dist_gamma3delta2mu0sigma4 = [
                0.8920033823420986,
                0.8960636351696349,
                0.9000425477930732,
                0.9039389251747214,
                0.9077516561952443,
                0.9114797172162497,
                0.9151221755072553,
                0.9186781925169079,
                0.9221470269679232,
                0.9255280377549651,
                0.9288206866245181,
                0.9320245406158417,
                0.9351392742422335,
                0.9381646713921781,
                0.9411006269304706,
                0.9439471479801147,
                0.9467043548667073,
                0.9493724817081621,
                0.9519518766339523,
                0.9544430016196438,
                0.9568464319242732,
                0.9591628551201553,
                0.9613930697069428,
                0.963537983304225,
                0.9655986104196027,
                0.9675760697920368,
                0.9694715813132707,
                0.9712864625333006,
                0.9730221247591393,
                0.9746800687594872,
                0.9762618800913417,
                0.9777692240680003,
                0.9792038403913063,
                0.9805675374743102,
                0.9818621864837009,
                0.9830897151343854,
                0.984252101271385,
                0.9853513662767325,
                0.986389568341251,
                0.9873687956429196,
                0.9882911594749505,
                0.9891587873676692,
                0.9899738162487803,
                0.9907383856865917,
                0.9914546312602373,
                0.992124678099887,
                0.9927506346383509,
                0.9933345866133962,
                0.9938785913575192,
                0.9943846724088771,
                0.9948548144736387,
                0.9952909587662006,
                0.9956949987495942,
                0.99606877629406,
                0.9964140782672334,
                0.9967326335647734,
                0.9970261105856271,
                0.9972961151515464,
                0.9975441888660314,
                0.9977718079036341,
                0.9979803822165881,
                0.9981712551420897,
                0.9983457033902876,
                0.9985049373901994,
                0.9986501019683699,
                0.9987822773331588,
                0.9989024803360947,
                0.9990116659807484,
                0.9991107291490707,
                0.9992005065150619,
                0.9992817786159831,
                0.9993552720520426,
                0.9994216617865402,
                0.9994815735198089,
                0.9995355861118806,
                0.9995842340305975,
                0.9996280098038227,
                0.9996673664564479,
                0.999702719914984,
                0.9997344513646357,
                0.9997629095458448,
                0.9997884129793229,
                0.999811252110543,
                0.9998316913665076,
                0.9998499711193302,
                0.9998663095527587,
                0.9998809044292056,
                0.9998939347561453,
                0.9999055623518748,
                0.9999159333116321,
                0.9999251793759067,
                0.9999334192034902,
                0.9999407595523929,
                0.9999472963722142,
                0.9999531158119045,
                0.9999582951471081,
                0.9999629036314445,
                0.9999670032761722,
                0.9999706495627025,
                0.9999738920923991,
                0.9999767751780154,
                0.9999793383810118,
                0.9999816169988374,
                0.9999836425060972,
                0.9999854429533354,
                0.9999870433269652,
                0.9999884658736696,
                0.9999897303923903,
                0.9999908544968092,
                0.9999918538510274,
                0.9999927423809372,
                0.9999935324636001,
                0.9999942350967481,
                0.9999948600503575,
                0.9999954160020725,
                0.9999959106581041,
                0.9999963508610797,
                0.9999967426861869,
                0.9999970915268243,
                0.9999974021708625,
                0.9999976788685055,
                0.9999979253926474,
                0.9999981450925305,
                0.9999983409414283,
                0.9999985155790031,
                0.9999986713489193,
                0.9999988103322347,
                0.9999989343770351,
                0.9999990451247301,
            ];

            foreach ((JohnsonSUDistribution dist, ddouble[] expecteds) in new[]{
                (dist_gamma1delta1mu0sigma1, expected_dist_gamma1delta1mu0sigma1),
                (dist_gamma2delta1mu1sigma1, expected_dist_gamma2delta1mu1sigma1),
                (dist_gamma1delta1mu0sigma2, expected_dist_gamma1delta1mu0sigma2),
                (dist_gamma2delta2mu2sigma2, expected_dist_gamma2delta2mu2sigma2),
                (dist_gamma3delta2mu0sigma4, expected_dist_gamma3delta2mu0sigma4),
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 1d / 16, i++) {
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