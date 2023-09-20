using DoubleDouble;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;

namespace DoubleDoubleDistribution.Tests {
    [TestClass()]
    public class VoigtDistributionTests {
        readonly VoigtDistribution dist1 = new() { EnableCDFErrorException = true };
        readonly VoigtDistribution dist2 = new(2, 3) { EnableCDFErrorException = true };
        readonly VoigtDistribution dist3 = new(gamma: 0, sigma: 1);
        readonly VoigtDistribution dist4 = new(gamma: 1, sigma: 0);

        [TestMethod()]
        public void InfoTest() {
            Console.WriteLine(dist1);
            Console.WriteLine($"Support={dist1.Support}");
            Console.WriteLine($"Gamma={dist1.Gamma}");
            Console.WriteLine($"Sigma={dist1.Sigma}");
            Console.WriteLine($"Mean={dist1.Mean}");
            Console.WriteLine($"Median={dist1.Median}");
            Console.WriteLine($"Mode={dist1.Mode}");
            Console.WriteLine($"Variance={dist1.Variance}");
            Console.WriteLine($"Skewness={dist1.Skewness}");
            Console.WriteLine($"Kurtosis={dist1.Kurtosis}");
            Console.WriteLine($"Entropy={dist1.Entropy}");

            Console.WriteLine(dist2);
            Console.WriteLine($"Support={dist2.Support}");
            Console.WriteLine($"Gamma={dist2.Gamma}");
            Console.WriteLine($"Sigma={dist2.Sigma}");
            Console.WriteLine($"Mean={dist2.Mean}");
            Console.WriteLine($"Median={dist2.Median}");
            Console.WriteLine($"Mode={dist2.Mode}");
            Console.WriteLine($"Variance={dist2.Variance}");
            Console.WriteLine($"Skewness={dist2.Skewness}");
            Console.WriteLine($"Kurtosis={dist2.Kurtosis}");
            Console.WriteLine($"Entropy={dist2.Entropy}");
        }

        [TestMethod()]
        public void PDFTest() {
            Console.WriteLine(dist1);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist1.PDF(x)}");
            }
            Console.WriteLine(dist2);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist2.PDF(x)}");
            }
            Console.WriteLine(dist3);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist3.PDF(x)}");
            }
            Console.WriteLine(dist4);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist4.PDF(x)}");
            }
        }

        [TestMethod()]
        public void CDFTest() {
            Console.WriteLine(dist1);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"cdf({x})={dist1.CDF(x)}");
            }
            Console.WriteLine(dist2);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"cdf({x})={dist2.CDF(x)}");
            }
            Console.WriteLine(dist3);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist3.CDF(x)}");
            }
            Console.WriteLine(dist4);
            for (ddouble x = -4; x <= 4; x += 0.125) {
                Console.WriteLine($"pdf({x})={dist4.CDF(x)}");
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            VoigtDistribution dist_g1s1 = new(gamma: 1, sigma: 1) { EnableCDFErrorException = true };
            VoigtDistribution dist_g2s1 = new(gamma: 2, sigma: 1) { EnableCDFErrorException = true };
            VoigtDistribution dist_g1s2 = new(gamma: 1, sigma: 2) { EnableCDFErrorException = true };
            VoigtDistribution dist_g2s4 = new(gamma: 2, sigma: 4) { EnableCDFErrorException = true };
            VoigtDistribution dist_grcp4srcp2 = new(gamma: 0.25, sigma: 0.5) { EnableCDFErrorException = true };

            ddouble[] expected_g1s1 = new ddouble[] {
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
            };
            ddouble[] expected_g2s1 = new ddouble[] {
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
            };
            ddouble[] expected_g1s2 = new ddouble[] {
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
            };
            ddouble[] expected_g2s4 = new ddouble[] {
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
            };
            ddouble[] expected_grcp4srcp2 = new ddouble[] {
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
            };

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

                    Assert.IsTrue(ddouble.Abs(expected_xm - actual_xm) < 1e-31, $"{dist} cdf({x})\n{expected_xm}\n{actual_xm}");
                    Assert.IsTrue(ddouble.Abs(expected_xp - actual_xp) < 1e-31, $"{dist} cdf({-x})\n{expected_xp}\n{actual_xp}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest2() {
            ddouble[] expected_g2x1 = new ddouble[] {
                "0.3539727153992217121416004416892697545005",
                "0.3582519319668703389051389987021591696915",
                "0.3712429038426146428180240570777320775097",
                "0.3976648579099691697800158222234655216494",
                "0.4307511881902131416803893933472823469520",
                "0.4587944436679107186848754519074873095959",
                "0.4773836738787166894380940343186936816256",
                "0.4881330274419626437864104964641970678957"
            };

            ddouble[] expected_g1x4096 = new ddouble[] {
                "0.00007771237330208453833412432734577847804988",
                "0.00007771237417058794835015914801789503228459",
                "0.00007771237764460217079130850317433000279144",
                "0.00007771239154066837859900149848833640747037",
                "0.00007771244712508229921910990108440449697049",
                "0.00007771266946512345671778006930215489737940",
                "0.00007771355886345855379591260816406757796795",
                "0.00007771711706770996823904829679178151386965"
            };

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

                    Assert.IsTrue(ddouble.Abs(expected - actual) < 1e-31, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }

            ddouble[] expected_s1x1 = new ddouble[] {
                "0.2095260165231521639564897266992738040323",
                "0.2493158045160731102028505942649372905811",
                "0.3064401674999125967887899353900927599618",
                "0.3712429038426146428180240570777320775097",
                "0.4258662262433281673993812035398189748086",
                "0.4609941230973476287681320011995394073356",
                "0.4802077080883765617852977901435431"
            };

            ddouble[] expected_s2x4096 = new ddouble[] {
                "0.00001942809824704404232609240896592444392918",
                "0.00003885619634933730169912129725720105545307",
                "0.00007771239154066837859900149848833640747037",
                "0.0001554247738172882012124121730753294983900",
                "0.0003108494735222277115802527196912789408301",
                "0.0006216983541469381182765169906667901848670",
                "0.001243391965154448460625084278297583421037",
                "0.002486745986496160029397009656499699234470"
            };

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

                    Assert.IsTrue(ddouble.Abs(expected - actual) < 1e-31, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void QuantileTest() {
            Assert.Inconclusive();

            //Console.WriteLine(dist1);
            //for (int i = 0; i <= 10; i++) {
            //    ddouble p = (ddouble)i / 10;
            //    ddouble x = dist1.Quantile(p);

            //    Console.WriteLine($"quantile({p})={x}, cdf({x})={dist1.CDF(x)}");
            //}
            //Console.WriteLine(dist2);
            //for (int i = 0; i <= 10; i++) {
            //    ddouble p = (ddouble)i / 10;
            //    ddouble x = dist2.Quantile(p);

            //    Console.WriteLine($"quantile({p})={x}, cdf({x})={dist2.CDF(x)}");
            //}
        }
    }
}