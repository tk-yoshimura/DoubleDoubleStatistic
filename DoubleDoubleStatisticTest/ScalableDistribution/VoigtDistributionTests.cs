using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ScalableDistribution {
    [TestClass()]
    public class VoigtDistributionTests {
        readonly VoigtDistribution dist1 = new() { EnableCDFErrorException = true };
        readonly VoigtDistribution dist2 = new(gamma: 2, sigma: 3) { EnableCDFErrorException = true };
        readonly VoigtDistribution dist3 = new(gamma: 3, sigma: 4) { EnableCDFErrorException = true };
        readonly VoigtDistribution dist4 = new(gamma: 1d / 16, sigma: 3d / 8) { EnableCDFErrorException = true };

        VoigtDistribution[] Dists => [
            dist1,
            dist2,
            dist3,
            dist4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (VoigtDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Gamma={dist.Gamma}");
                Console.WriteLine($"Sigma={dist.Sigma}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);

                for (ddouble t = 0; t <= 1d / 64; t += 1d / 256) {
                    ddouble y = dist.Integrand(t);
                    Console.WriteLine($"{t},{y}");
                }
            }
        }

        [TestMethod()]
        public void ModeTest() {
            foreach (VoigtDistribution dist in Dists) {
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
            foreach (VoigtDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (VoigtDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (VoigtDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (VoigtDistribution dist in Dists) {
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
            Assert.Inconclusive();

            foreach (VoigtDistribution dist in Dists) {
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
            Assert.Inconclusive();

            foreach (VoigtDistribution dist in Dists) {
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
            foreach (VoigtDistribution dist in Dists) {
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

                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(0d, Interval.Lower)), "quantile(0)");
                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(1d, Interval.Lower)), "quantile(1)");

                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(0d, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(0d, Interval.Upper)), "cquantile(0)");
                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(1d, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(1d, Interval.Upper)), "cquantile(1)");

                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Lower)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.Epsilon, Interval.Lower)), "quantile(0+eps)");
                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Lower)), "quantile(1-eps)");

                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.Epsilon, Interval.Upper)) || ddouble.IsPositiveInfinity(dist.Quantile(ddouble.Epsilon, Interval.Upper)), "cquantile(0+eps)");
                //Assert.IsTrue(ddouble.IsFinite(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)) || ddouble.IsNegativeInfinity(dist.Quantile(ddouble.One - ddouble.Epsilon, Interval.Upper)), "cquantile(1-eps)");
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist1 = [
                "0.02281363525870709106014838385691384365070",
                "0.03086120886094135735354018509974065535798",
                "0.04338582232367966337165545402298233253914",
                "0.06268641077529936343580986731424007272669",
                "0.09071519942627542976993106314883366450545",
                "0.12692380186284846359460955172649948351772",
                "0.1657956626891664570736114255348777648542",
                "0.1967698598754764545622833689698698093889",
                "0.2087092805203676891488309954152961385611",
                "0.1967698598754764545622833689698698093889",
                "0.1657956626891664570736114255348777648542",
                "0.12692380186284846359460955172649948351772",
                "0.09071519942627542976993106314883366450545",
                "0.06268641077529936343580986731424007272669",
                "0.04338582232367966337165545402298233253914",
                "0.03086120886094135735354018509974065535798",
                "0.02281363525870709106014838385691384365070",
            ];
            ddouble[] expected_dist2 = [
                "0.05044703043423117923950337640335990536865",
                "0.05655720731611957570940483672002156187293",
                "0.06259989296472363038561271294590173068131",
                "0.06833209582317445882858664007346427669894",
                "0.07349491459470745580309853336469662715798",
                "0.07783338019412799034119106607626315274213",
                "0.08111833862031692803987871040222532646898",
                "0.08316771422941912690099521913621715852549",
                "0.08386432168347958410847372581382211452530",
                "0.08316771422941912690099521913621715852549",
                "0.08111833862031692803987871040222532646898",
                "0.07783338019412799034119106607626315274213",
                "0.07349491459470745580309853336469662715798",
                "0.06833209582317445882858664007346427669894",
                "0.06259989296472363038561271294590173068131",
                "0.05655720731611957570940483672002156187293",
                "0.05044703043423117923950337640335990536865",
            ];

            foreach ((VoigtDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
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
            VoigtDistribution dist_g1s1 = new(gamma: 1, sigma: 1) { EnableCDFErrorException = true };
            VoigtDistribution dist_g2s1 = new(gamma: 2, sigma: 1) { EnableCDFErrorException = true };
            VoigtDistribution dist_g1s2 = new(gamma: 1, sigma: 2) { EnableCDFErrorException = true };
            VoigtDistribution dist_g2s4 = new(gamma: 2, sigma: 4) { EnableCDFErrorException = true };
            VoigtDistribution dist_grcp4srcp2 = new(gamma: 0.25, sigma: 0.5) { EnableCDFErrorException = true };

            ddouble[] expected_g1s1 = [
                "0.4869597013016863855493486350777892009578",
                "0.4739435568556796986247981620292283431457",
                "0.4480793400783183843992284017161498487239",
                "0.3976648579099691697800158222234655216494",
                "0.3064401674999125967887899353900927599618",
                "0.1790409327778761353365974401950742448947",
                "0.08328272536980600233054061457933769058255",
                "0.04021529507115189064236405693848566329169",
                "0.01994654381367084457887805557836011988423",
                "0.009953671384601179528618251422775685390544",
                "0.004974401831469220797942613207441455946383",
                "0.002486897184917028493253508809926723095236",
                "0.001243410641751309304735504342881727173909",
                "0.0006217005775258804205857745601026920493978",
                "0.0003108496958594858217752055979881431940715",
                "0.0001554247738172882012124121730753294983900",
                "0.00007771237764460217079130850317433000279144"
            ];
            ddouble[] expected_g2s1 = [
                "0.4916185091137958384384874964353325139788",
                "0.4832453124814591083014916393316362618196",
                "0.4665567634069671761578017096833965855896",
                "0.4336358240716738541493406479303013578695",
                "0.3712429038426146428180240570777320775097",
                "0.2687854132756434357853757551296878614970",
                "0.1545337038000290456821373244543304943798",
                "0.07912613623757779599037165953559429098244",
                "0.03973582003283893405114894593735683701760",
                "0.01988785759032445722620027711466421296366",
                "0.009946373371235215700774375968323319837866",
                "0.004973490750276544830237497348299145113309",
                "0.002486783336270143203999176459631249376449",
                "0.001243396411810555133436435872694719114419",
                "0.0006216987988189098903714258919294779408190",
                "0.0003108494735222277115802527196912789408301",
                "0.0001554247460251657248481785250370311504974"
            ];
            ddouble[] expected_g1s2 = [
                "0.4912836054889874478914361351257229883410",
                "0.4825729920053764314750036531207066138594",
                "0.4651921542119972575619608996180705862789",
                "0.4307511881902131416803893933472823469520",
                "0.3643598436553283668160205703702042280208",
                "0.2493158045160731102028505942649372905811",
                "0.1118402087461845147872020735721041155661",
                "0.04265177186286753680506962190071474912861",
                "0.02019248339985736858284054540665159359771",
                "0.009983190420698109661683645319413963807314",
                "0.004978056213637081662075552290874402335525",
                "0.002487352892571979595530059451587076168682",
                "0.001243467571288159353707147262003417050295",
                "0.0006217076926591205443280663330422101309758",
                "0.0003108505852180601199715288348130937917367",
                "0.0001554248849860762851814257021676231852828",
                "0.00007771239154066837859900149848833640747037"
            ];
            ddouble[] expected_g2s4 = [
                "0.4956414412776763679432773675653477200330",
                "0.4912836054889874478914361351257229883410",
                "0.4825729920053764314750036531207066138594",
                "0.4651921542119972575619608996180705862789",
                "0.4307511881902131416803893933472823469520",
                "0.3643598436553283668160205703702042280208",
                "0.2493158045160731102028505942649372905811",
                "0.1118402087461845147872020735721041155661",
                "0.04265177186286753680506962190071474912861",
                "0.02019248339985736858284054540665159359771",
                "0.009983190420698109661683645319413963807314",
                "0.004978056213637081662075552290874402335525",
                "0.002487352892571979595530059451587076168682",
                "0.001243467571288159353707147262003417050295",
                "0.0006217076926591205443280663330422101309758",
                "0.0003108505852180601199715288348130937917367",
                "0.0001554248849860762851814257021676231852828"
            ];
            ddouble[] expected_grcp4srcp2 = [
                "0.4651921542119972575619608996180705862789",
                "0.4307511881902131416803893933472823469520",
                "0.3643598436553283668160205703702042280208",
                "0.2493158045160731102028505942649372905811",
                "0.1118402087461845147872020735721041155661",
                "0.04265177186286753680506962190071474912861",
                "0.02019248339985736858284054540665159359771",
                "0.009983190420698109661683645319413963807314",
                "0.004978056213637081662075552290874402335525",
                "0.002487352892571979595530059451587076168682",
                "0.001243467571288159353707147262003417050295",
                "0.0006217076926591205443280663330422101309758",
                "0.0003108505852180601199715288348130937917367",
                "0.0001554248849860762851814257021676231852828",
                "0.00007771239154066837859900149848833640747037",
                "0.00003885618940130357659186950317678167819859",
                "0.00001942809390452344944887601246048244791697"
            ];

            foreach ((VoigtDistribution dist, ddouble[] expecteds) in new[]{
                (dist_g1s1, expected_g1s1), (dist_g2s1, expected_g2s1), (dist_g1s2, expected_g1s2),
                (dist_g2s4, expected_g2s4), (dist_grcp4srcp2, expected_grcp4srcp2)
            }) {
                for ((ddouble x, int i) = (-1d / 16, 0); i < expecteds.Length; x *= 2, i++) {
                    ddouble expected_xm = expecteds[i];
                    ddouble actual_xm = dist.CDF(x);
                    ddouble expected_xp = 1 - expecteds[i];
                    ddouble actual_xp = dist.CDF(-x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected_xm);
                    Console.WriteLine(actual_xm);

                    Console.WriteLine($"{dist} cdf({-x})");
                    Console.WriteLine(expected_xp);
                    Console.WriteLine(actual_xp);

                    Assert.IsTrue(ddouble.Abs(expected_xm - actual_xm) / expected_xm < 1e-30, $"{dist} cdf({x})\n{expected_xm}\n{actual_xm}");
                    Assert.IsTrue(ddouble.Abs(expected_xp - actual_xp) / expected_xp < 1e-30, $"{dist} cdf({-x})\n{expected_xp}\n{actual_xp}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest2() {
            ddouble[] expected_g2x1 = [
                "0.3539727153992217121416004416892697545005",
                "0.3582519319668703389051389987021591696915",
                "0.3712429038426146428180240570777320775097",
                "0.3976648579099691697800158222234655216494",
                "0.4307511881902131416803893933472823469520",
                "0.4587944436679107186848754519074873095959",
                "0.4773836738787166894380940343186936816256",
                "0.4881330274419626437864104964641970678957"
            ];

            ddouble[] expected_g1x4096 = [
                "0.00007771237330208453833412432734577847804988",
                "0.00007771237417058794835015914801789503228459",
                "0.00007771237764460217079130850317433000279144",
                "0.00007771239154066837859900149848833640747037",
                "0.00007771244712508229921910990108440449697049",
                "0.00007771266946512345671778006930215489737940",
                "0.00007771355886345855379591260816406757796795",
                "0.00007771711706770996823904829679178151386965"
            ];

            foreach ((ddouble gamma, ddouble x, ddouble[] expecteds) in new[]{
                (2, -1, expected_g2x1), (1, -4096, expected_g1x4096)
            }) {
                for ((ddouble sigma, int i) = (1d / 4, 0); i < expecteds.Length; sigma *= 2, i++) {
                    VoigtDistribution dist = new(gamma, sigma) { EnableCDFErrorException = true };
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }

            ddouble[] expected_s1x1 = [
                "0.2095260165231521639564897266992738040323",
                "0.2493158045160731102028505942649372905811",
                "0.3064401674999125967887899353900927599618",
                "0.3712429038426146428180240570777320775097",
                "0.4258662262433281673993812035398189748086",
                "0.4609941230973476287681320011995394073356",
                "0.4802077080883765617852977901435431"
            ];

            ddouble[] expected_s2x4096 = [
                "0.00001942809824704404232609240896592444392918",
                "0.00003885619634933730169912129725720105545307",
                "0.00007771239154066837859900149848833640747037",
                "0.0001554247738172882012124121730753294983900",
                "0.0003108494735222277115802527196912789408301",
                "0.0006216983541469381182765169906667901848670",
                "0.001243391965154448460625084278297583421037",
                "0.002486745986496160029397009656499699234470"
            ];

            foreach ((ddouble sigma, ddouble x, ddouble[] expecteds) in new[]{
                (1, -1, expected_s1x1), (2, -4096, expected_s2x4096)
            }) {
                for ((ddouble gamma, int i) = (1d / 4, 0); i < expecteds.Length; gamma *= 2, i++) {
                    VoigtDistribution dist = new(gamma, sigma) { EnableCDFErrorException = true };
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