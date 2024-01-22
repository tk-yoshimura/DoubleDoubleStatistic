using DoubleDouble;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;

namespace DoubleDoubleDistribution.Tests {
    [TestClass()]
    public class ChiSquareDistributionTests {
        readonly ChiSquareDistribution dist_nu1 = new(k: 1);
        readonly ChiSquareDistribution dist_nu2 = new(k: 2);
        readonly ChiSquareDistribution dist_nu3 = new(k: 3);
        readonly ChiSquareDistribution dist_nu4 = new(k: 4);
        readonly ChiSquareDistribution dist_nu5 = new(k: 5);
        readonly ChiSquareDistribution dist_nu8 = new(k: 8);
        readonly ChiSquareDistribution dist_nu16 = new(k: 16);
        readonly ChiSquareDistribution dist_nu32 = new(k: 32);
        readonly ChiSquareDistribution dist_nu64 = new(k: 64);
        readonly ChiSquareDistribution dist_nu128 = new(k: 128);


        ChiSquareDistribution[] Dists => new[]{
            dist_nu1,
            dist_nu2,
            dist_nu3,
            dist_nu4,
            dist_nu5,
            dist_nu8,
            dist_nu16,
            dist_nu32,
            dist_nu64,
            dist_nu128,
        };

        [TestMethod()]
        public void InfoTest() {
            foreach (ChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.K}");
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
        public void PDFTest() {
            foreach (ChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);
                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (ChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);

                Console.WriteLine($"cdf(median)={cdf_median}");

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (ChiSquareDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble cdf_median = dist.CDF(dist.Median, Interval.Lower);
                ddouble ccdf_median = dist.CDF(dist.Median, Interval.Upper);

                Console.WriteLine($"ccdf(median)={ccdf_median}");

                Assert.IsTrue(ddouble.Abs(cdf_median + ccdf_median - 1) < 1e-28);

                foreach (ddouble x in new ddouble[] { 0, double.Epsilon, ddouble.Epsilon }) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 0.125; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
                for (ddouble x = 8; x <= ddouble.Ldexp(1, 128); x *= 2) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (ChiSquareDistribution dist in Dists) {
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
            foreach (ChiSquareDistribution dist in Dists) {
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
                ddouble.PositiveInfinity,
                "0.4393912894677223970468619774122289491813",
                "0.2419707245191433497978301929355606548287",
                "0.1538663228054552595710179894167888707926",
                "0.1037768743551486758350670623603343413422",
                "0.07228895706727250797599023913101522913360",
                "0.05139344326792309227040685259379788210172",
                "0.03705618452374812274221024500559180756934",
                "0.02699548325659402597528210020535679086991",
                "0.01982171487060489525516959162023486879074",
                "0.01464498256192648711320029563820125074252",
                "0.01087474033728314171408276418482327795923",
                "0.008108695554940243370931735685399711574086",
                "0.006067311902576735170812033228672466316908",
                "0.004553342921640173367468696285302704084376",
                "0.003425903510139482983574464022795052535111",
                "0.002583373169261506732133617663323184471166"
            ];
            ddouble[] expected_nu2 = [
                0.5,
                "0.3894003915357024341225851334891603236484",
                "0.3032653298563167118018997674955902267210",
                "0.2361832763705073535690232754716339564851",
                "0.1839397205857211607977618850807304337229",
                "0.1432523984300950501624427133239188013966",
                "0.1115650800742149144666402353820062606711",
                "0.08688697172522256334035862933318550800736",
                "0.06766764161830634594699974748624220170382",
                "0.05269961228093216839160884462034904863425",
                "0.04104249931194939758476433723357990391890",
                "0.03196393060335378635121501277897587465432",
                "0.02489353418393197148967120782503088831585",
                "0.01938710391586100494344991763379807163007",
                "0.01509869171115925036989314618180992253583",
                "0.01175887292800455411807559255021646973503",
                "0.009157819444367090146859010636620621105956"
            ];
            ddouble[] expected_nu4 = [
                "0",
                "0.09735009788392560853064628337229008091210",
                "0.1516326649281583559009498837477951133605",
                "0.1771374572778805151767674566037254673638",
                "0.1839397205857211607977618850807304337229",
                "0.1790654980376188127030533916548985017457",
                "0.1673476201113223716999603530730093910066",
                "0.1520522005191394858456276013330746390129",
                "0.1353352832366126918939994949724844034076",
                "0.1185741276320973788811199003957853594271",
                "0.1026062482798734939619108430839497597973",
                "0.08790080915922291246584128514218365529937",
                "0.07468060255179591446901362347509266494755",
                "0.06300808772654826606621223230984373279773",
                "0.05284542098905737629462601163633472887541",
                "0.04409577348001707794278347206331176150638",
                "0.03663127777746836058743604254648248442382"
            ];

            foreach ((ChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else if (ddouble.IsPositiveInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_nu1 = [
                "0",
                "0.5204998778130465376827466538919645287365",
                "0.6826894921370858971704650912640758449558",
                "0.7793286380801532073955805384702048031390",
                "0.8427007929497148693412206350826092592961",
                "0.8861537019933419497213981646064152039057",
                "0.9167354833364495981450806795009063640672",
                "0.9386311708605978269800109631245451217772",
                "0.9544997361036415855994347256669331250564",
                "0.9661051464753107270669762616459478586814",
                "0.9746526813225317360683963041989610076361",
                "0.9809835263276994562359137099868087143623",
                "0.9856941215645703604741521935777493064263",
                "0.9892125507453296434291745225994498093645",
                "0.9918490284064972996870207396720704572571",
                "0.9938301006794558377864011315176458124986",
                "0.9953222650189527341620692563672529286109"
            ];
            ddouble[] expected_nu2 = [
                "0",
                "0.2211992169285951317548297330216793527032",
                "0.3934693402873665763962004650088195465581",
                "0.5276334472589852928619534490567320870298",
                "0.6321205588285576784044762298385391325542",
                "0.7134952031398098996751145733521623972068",
                "0.7768698398515701710667195292359874786578",
                "0.8262260565495548733192827413336289839853",
                "0.8646647167633873081060005050275155965924",
                "0.8946007754381356632167823107593019027315",
                "0.9179150013761012048304713255328401921622",
                "0.9360721387932924272975699744420482506914",
                "0.9502129316321360570206575843499382233683",
                "0.9612257921682779901131001647324038567399",
                "0.9698026165776814992602137076363801549283",
                "0.9764822541439908917638488148995670605299",
                "0.9816843611112658197062819787267587577881"
            ];
            ddouble[] expected_nu4 = [
                "0",
                "0.02649902116074391469353716627709919087903",
                "0.09020401043104986459430069751322931983712",
                "0.1733585327032242625084185358492811523021",
                "0.2642411176571153568089524596770782651084",
                "0.3553642070645722742690077900423653937154",
                "0.4421745996289254276667988230899686966446",
                "0.5221216555112759016280275386674797059595",
                "0.5939941502901619243180015150825467897771",
                "0.6574525201739409054545425099677311838774",
                "0.7127025048163542169066496393649406725677",
                "0.7602705204748466023658874041576809400926",
                "0.8008517265285442280826303373997528934732",
                "0.8352096167151814579806757001127163911444",
                "0.8641117745995667466709616843637106971775",
                "0.8882907071839567358782818707729435375172",
                "0.9084218055563290985314098936337937889404"
            ];

            foreach ((ChiSquareDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }
    }
}