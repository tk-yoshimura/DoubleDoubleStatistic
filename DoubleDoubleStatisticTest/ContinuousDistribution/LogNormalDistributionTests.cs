using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class LogNormalDistributionTests {
        readonly LogNormalDistribution dist1 = new();
        readonly LogNormalDistribution dist2 = new(mu: 1, sigma: 3);

        LogNormalDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (LogNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Sigma={dist.Sigma}");
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
            foreach (LogNormalDistribution dist in Dists) {
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
            foreach (LogNormalDistribution dist in Dists) {
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
            foreach (LogNormalDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (LogNormalDistribution dist in Dists) {
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
            foreach (LogNormalDistribution dist in Dists) {
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
            foreach (LogNormalDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-15, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (LogNormalDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (LogNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (LogNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (LogNormalDistribution dist in Dists) {
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
            foreach (LogNormalDistribution dist in Dists) {
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
            foreach (LogNormalDistribution dist in Dists) {
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
            ddouble[] expected_dist1 = [
                "0",
                "0.6274960771159243653767016268086090314327",
                "0.3989422804014326779399460599343818684759",
                "0.2449736517105099259576674746893011503214",
                "0.1568740192789810913441754067021522578582",
                "0.1048710668896498317667280827745522595137",
                "0.07272825613999471130773276476702907840881",
                "0.05200533188970584454374667096712120183632",
                "0.03815345651188644191736556483467039275977",
                "0.02860595581010994995864455818464238864633",
                "0.02185071483032719963219291008170479355059",
                "0.01696217092022252956815860543455322797117",
                "0.01335453835505393386570121990985099521830",
                "0.01064609197543415135004081142890611226199",
                "0.008581626313996365138005315542859441828836",
                "0.006986618202810428328793408133349640712022",
                "0.005739296497825189868811749135958922862285"
            ];
            ddouble[] expected_dist2 = [
                "0",
                "0.2268043824438617779486782383869347967530",
                "0.1257944092309977213393073414497281190758",
                "0.08692989856582284381499599480553812555520",
                "0.06614347460597307915614152808726912049225",
                "0.05317160078779778619300895891341093313722",
                "0.04430297918053212449051939720253560665085",
                "0.03785988442778966455942473508482959854336",
                "0.03297072052485991499643052762444811191728",
                "0.02913705482374524447458068286306945425861",
                "0.02605298638028731444575405447124664361479",
                "0.02352029332185975215412238994559444945475",
                "0.02140486326054810024095609501886699868397",
                "0.01961271125285002211521057536245095946559",
                "0.01807601997098268390467754340597386231331",
                "0.01674463068814000353671330055109345118506",
                "0.01558065238968434687204464227372714240348"
            ];

            foreach ((LogNormalDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
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
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist1 = [
                "0",
                "0.2441085957855827340050393064789443013828",
                "0.5000000000000000000000000000000000000000",
                "0.6574321694851541102495800791665482749989",
                "0.7558914042144172659949606935210556986172",
                "0.8202427861042146188088320181890927788824",
                "0.8640313923585755422570248647012813433986",
                "0.8948540088998852514357896947445200426377",
                "0.9171714809983015146501108793751879425037",
                "0.9337192802504506853118121644324196919089",
                "0.9462396895483368740690252109505285577360",
                "0.9558792918709237749230494555973365010856",
                "0.9634142480829571043034700301782009184052",
                "0.9693830114995456860171269002048222475960",
                "0.9741672331954078505835383577482155202229",
                "0.9780425942704737106193288901884205235667",
                "0.9812116071859449016891344979915675588882"
            ];
            ddouble[] expected_dist2 = [
                "0",
                "0.2862469592460356079058171470118215438077",
                "0.3694413401817636382727922820695735833283",
                "0.4214527834864947207089855106318382553129",
                "0.4592655190218048085962334216970826218286",
                "0.4888697222195098741161547264345252317990",
                "0.5131111759871297159183908951819066445512",
                "0.5335728858438851177431707674584030710673",
                "0.5512281153063942832333920046757719086016",
                "0.5667185084960893349554465173487764253336",
                "0.5804895293730796755086416074648829655614",
                "0.5928630465166395999334223852935736577021",
                "0.6040791497815755181942754549767157413761",
                "0.6143216474785582817354265214633669633767",
                "0.6237343416814816954282485679648993560892",
                "0.6324318174750614327676950177306392825792",
                "0.6405068263898791577843407099134968863039"
            ];

            foreach ((LogNormalDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
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