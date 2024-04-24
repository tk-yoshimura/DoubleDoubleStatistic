using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class StudentTDistributionTests {
        readonly StudentTDistribution dist_nu0p5 = new(nu: 0.5);
        readonly StudentTDistribution dist_nu1 = new(nu: 1);
        readonly StudentTDistribution dist_nu1p5 = new(nu: 1.5);
        readonly StudentTDistribution dist_nu2 = new(nu: 2);
        readonly StudentTDistribution dist_nu2p5 = new(nu: 2.5);
        readonly StudentTDistribution dist_nu3 = new(nu: 3);
        readonly StudentTDistribution dist_nu3p5 = new(nu: 3.5);
        readonly StudentTDistribution dist_nu4 = new(nu: 4);
        readonly StudentTDistribution dist_nu4p5 = new(nu: 4.5);
        readonly StudentTDistribution dist_nu5 = new(nu: 5);
        readonly StudentTDistribution dist_nu8 = new(nu: 8);
        readonly StudentTDistribution dist_nu16 = new(nu: 16);
        readonly StudentTDistribution dist_nu32 = new(nu: 32);
        readonly StudentTDistribution dist_nu64 = new(nu: 64);
        readonly StudentTDistribution dist_nu128 = new(nu: 128);
        readonly StudentTDistribution dist_nu129 = new(nu: 129);


        StudentTDistribution[] Dists => [
            dist_nu0p5,
            dist_nu1,
            dist_nu1p5,
            dist_nu2,
            dist_nu2p5,
            dist_nu3,
            dist_nu3p5,
            dist_nu4,
            dist_nu4p5,
            dist_nu5,
            dist_nu8,
            dist_nu16,
            dist_nu32,
            dist_nu64,
            dist_nu128,
            dist_nu129
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Nu={dist.Nu}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (StudentTDistribution dist in Dists) {
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
            foreach (StudentTDistribution dist in Dists) {
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
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (StudentTDistribution dist in Dists) {
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
            foreach (StudentTDistribution dist in Dists) {
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
            foreach (StudentTDistribution dist in Dists) {
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
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (StudentTDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (StudentTDistribution dist in Dists) {
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
            foreach (StudentTDistribution dist in Dists) {
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
        public void PDFExpectedTest() {
            ddouble[] expected_nu1 = [
                "0.01872411095198768656104514863206051318052",
                "0.02402338763651250351228434164113424332596",
                "0.03183098861837906715377675267450287240689",
                "0.04390481188741940297072655541310741021640",
                "0.06366197723675813430755350534900574481378",
                "0.09794150344116636047315923899847037663659",
                "0.1591549430918953357688837633725143620345",
                "0.2546479089470325372302140213960229792551",
                "0.3183098861837906715377675267450287240689",
                "0.2546479089470325372302140213960229792551",
                "0.1591549430918953357688837633725143620345",
                "0.09794150344116636047315923899847037663659",
                "0.06366197723675813430755350534900574481378",
                "0.04390481188741940297072655541310741021640",
                "0.03183098861837906715377675267450287240689",
                "0.02402338763651250351228434164113424332596",
                "0.01872411095198768656104514863206051318052"
            ];
            ddouble[] expected_nu2 = [
                "0.01309457002197310230371934003897868591268",
                "0.01858992781845675517940704169146473883219",
                "0.02741012223434214751334655154273294780105",
                "0.04220064386804796077025242584550177644239",
                "0.06804138174397716939436900207516364977683",
                "0.1141344117818037522441912762899398484470",
                "0.1924500897298752548363829268339858185492",
                "0.2962962962962962962962962962962962962963",
                "0.3535533905932737622004221810524245196424",
                "0.2962962962962962962962962962962962962963",
                "0.1924500897298752548363829268339858185492",
                "0.1141344117818037522441912762899398484470",
                "0.06804138174397716939436900207516364977683",
                "0.04220064386804796077025242584550177644239",
                "0.02741012223434214751334655154273294780105",
                "0.01858992781845675517940704169146473883219",
                "0.01309457002197310230371934003897868591268"
            ];
            ddouble[] expected_nu4 = [
                "0.006708203932499369089227521006193828706322",
                "0.01127321611414344311883033037937787205956",
                "0.01969349809083653687639083077362127963360",
                "0.03567562436955665030413771690855873091766",
                "0.06629126073623883041257915894732959743295",
                "0.1228800000000000000000000000000000000000",
                "0.2146625258399798108552806721982025186023",
                "0.3222618685603870651600694859951242779679",
                "0.3750000000000000000000000000000000000000",
                "0.3222618685603870651600694859951242779679",
                "0.2146625258399798108552806721982025186023",
                "0.1228800000000000000000000000000000000000",
                "0.06629126073623883041257915894732959743295",
                "0.03567562436955665030413771690855873091766",
                "0.01969349809083653687639083077362127963360",
                "0.01127321611414344311883033037937787205956",
                "0.006708203932499369089227521006193828706322"
            ];
            ddouble[] expected_nu5 = [
                "0.005123727051917914253153098487591714825869",
                "0.009244354092520921647506930056678994317208",
                "0.01729257880022296060439170739562203753731",
                "0.03332623888702283072154572385173728114295",
                "0.06509031032621646625301899189792437723233",
                "0.1245173446463551375415496365570363583645",
                "0.2196797973509805736039390976554947731591",
                "0.3279185313227465122017983032058697488556",
                "0.3796066898224944311876067607486949680190",
                "0.3279185313227465122017983032058697488556",
                "0.2196797973509805736039390976554947731591",
                "0.1245173446463551375415496365570363583645",
                "0.06509031032621646625301899189792437723233",
                "0.03332623888702283072154572385173728114295",
                "0.01729257880022296060439170739562203753731",
                "0.009244354092520921647506930056678994317208",
                "0.005123727051917914253153098487591714825869"
            ];
            ddouble[] expected_nu8 = [
                "0.002756305973425001075003373926655934886793",
                "0.005920595490240063168161454672057883856525",
                "0.01300941799263384922791514649345945758809",
                "0.02878134758931450349116644254251926324410",
                "0.06236808463468179554882780036286577452983",
                "0.1267712053742722929704834199883744711569",
                "0.2276075801453030533963318599806939998984",
                "0.3366938979282274901089900531046628639442",
                "0.3866990209613931774067117605260893183589",
                "0.3366938979282274901089900531046628639442",
                "0.2276075801453030533963318599806939998984",
                "0.1267712053742722929704834199883744711569",
                "0.06236808463468179554882780036286577452983",
                "0.02878134758931450349116644254251926324410",
                "0.01300941799263384922791514649345945758809",
                "0.005920595490240063168161454672057883856525",
                "0.002756305973425001075003373926655934886793"
            ];

            foreach ((StudentTDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4),
                (dist_nu5, expected_nu5), (dist_nu8, expected_nu8)
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_nu1 = [
                "0.07797913037736932546051288977313013511652",
                "0.08858553278290474887587605290700607921396",
                "0.1024163823495667258245989237752594740489",
                "0.1211189415908433987235825893719092475981",
                "0.1475836176504332741754010762247405259511",
                "0.1871670418109988161862527475647785523348",
                "0.2500000000000000000000000000000000000000",
                "0.3524163823495667258245989237752594740489",
                "0.5000000000000000000000000000000000000000",
                "0.6475836176504332741754010762247405259511",
                "0.7500000000000000000000000000000000000000",
                "0.8128329581890011838137472524352214476652",
                "0.8524163823495667258245989237752594740489",
                "0.8788810584091566012764174106280907524019",
                "0.8975836176504332741754010762247405259511",
                "0.9114144672170952511241239470929939207860",
                "0.9220208696226306745394871102268698648835"
            ];
            ddouble[] expected_nu2 = [
                "0.02859547920896831706610375859676730714344",
                "0.03641367502723466771353689781909807537232",
                "0.04773298313335456602978189954490636128267",
                "0.06480586011075540455677185846826293043786",
                "0.09175170953613698363378598754901810133901",
                "0.1361965624455005397216403068258167330752",
                "0.2113248654051871177454256097490212721762",
                "0.3333333333333333333333333333333333333333",
                "0.5000000000000000000000000000000000000000",
                "0.6666666666666666666666666666666666666667",
                "0.7886751345948128822545743902509787278238",
                "0.8638034375544994602783596931741832669248",
                "0.9082482904638630163662140124509818986610",
                "0.9351941398892445954432281415317370695621",
                "0.9522670168666454339702181004550936387173",
                "0.9635863249727653322864631021809019246277",
                "0.9714045207910316829338962414032326928566"
            ];
            ddouble[] expected_nu4 = [
                "0.008065044950046266789981792879119228203064",
                "0.01244808173011137605214918808458851502844",
                "0.01997098403585941363797349989298130893104",
                "0.03338327240599407251944874173636399470330",
                "0.05805826175840779724947227368446935044698",
                "0.1040000000000000000000000000000000000000",
                "0.1869504831500294425027156863776213270383",
                "0.3216649815909316371184511307969689868016",
                "0.5000000000000000000000000000000000000000",
                "0.6783350184090683628815488692030310131984",
                "0.8130495168499705574972843136223786729617",
                "0.8960000000000000000000000000000000000000",
                "0.9419417382415922027505277263155306495530",
                "0.9666167275940059274805512582636360052967",
                "0.9800290159641405863620265001070186910690",
                "0.9875519182698886239478508119154114849716",
                "0.9919349550499537332100182071208807717969"
            ];
            ddouble[] expected_nu5 = [
                "0.005161707740415726902199251692191844072254",
                "0.008642215892646677330846717858437413191096",
                "0.01504962394873128692355240594743758371204",
                "0.02724504967118812055775637958030632022379",
                "0.05096973941492917812268055292114371531070",
                "0.09695184012123671606642657958600844110842",
                "0.1816087338245613128000742599987684447992",
                "0.3191494358204645033534055269763769204442",
                "0.5000000000000000000000000000000000000000",
                "0.6808505641795354966465944730236230795558",
                "0.8183912661754386871999257400012315552008",
                "0.9030481598787632839335734204139915588916",
                "0.9490302605850708218773194470788562846893",
                "0.9727549503288118794422436204196936797762",
                "0.9849503760512687130764475940525624162880",
                "0.9913577841073533226691532821415625868089",
                "0.9948382922595842730978007483078081559277"
            ];
            ddouble[] expected_nu8 = [
                "0.001974886401722662905104665366510507883467",
                "0.004039541130205945303011581300234226899652",
                "0.008535840616891325490031898311204082519930",
                "0.01847101885681205240168000127732461034142",
                "0.04025811897863133566864078589658943346580",
                "0.08600164597595564370290795949129554263514",
                "0.1732967535436671239140374942844078646548",
                "0.3152680377784881708086855327987857555905",
                "0.5000000000000000000000000000000000000000",
                "0.6847319622215118291913144672012142444095",
                "0.8267032464563328760859625057155921353452",
                "0.9139983540240443562970920405087044573649",
                "0.9597418810213686643313592141034105665342",
                "0.9815289811431879475983199987226753896586",
                "0.9914641593831086745099681016887959174801",
                "0.9959604588697940546969884186997657731003",
                "0.9980251135982773370948953346334894921165"
            ];

            foreach ((StudentTDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4),
                (dist_nu5, expected_nu5), (dist_nu8, expected_nu8)
            }) {
                for ((ddouble x, int i) = (-4, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }
    }
}