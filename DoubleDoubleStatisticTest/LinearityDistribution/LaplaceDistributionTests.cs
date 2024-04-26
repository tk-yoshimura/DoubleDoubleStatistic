using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.LinearityDistribution {
    [TestClass()]
    public class LaplaceDistributionTests {
        readonly LaplaceDistribution dist1 = new();
        readonly LaplaceDistribution dist2 = new(mu: 1, sigma: 3);

        LaplaceDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (LaplaceDistribution dist in Dists) {
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
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (LaplaceDistribution dist in Dists) {
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
            foreach (LaplaceDistribution dist in Dists) {
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
            foreach (LaplaceDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (LaplaceDistribution dist in Dists) {
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
            foreach (LaplaceDistribution dist in Dists) {
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
            foreach (LaplaceDistribution dist in Dists) {
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
            foreach (LaplaceDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (LaplaceDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (LaplaceDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (LaplaceDistribution dist in Dists) {
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
            foreach (LaplaceDistribution dist in Dists) {
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
            foreach (LaplaceDistribution dist in Dists) {
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
                1.677313139512559266e-04,
                2.153712702878437645e-04,
                2.765421850739168135e-04,
                3.550871944212745141e-04,
                4.559409827772581219e-04,
                5.854398103955872146e-04,
                7.517195964887861837e-04,
                9.652270681138546511e-04,
                1.239376088333179245e-03,
                1.591390398254833444e-03,
                2.043385719232033298e-03,
                2.623759199590692311e-03,
                3.368973499542733500e-03,
                4.325847601560317036e-03,
                5.554498269121153041e-03,
                7.132116954499627751e-03,
                9.157819444367089334e-03,
                1.175887292800455348e-02,
                1.509869171115925043e-02,
                1.938710391586100437e-02,
                2.489353418393197223e-02,
                3.196393060335378511e-02,
                4.104249931194939999e-02,
                5.269961228093216626e-02,
                6.766764161830635116e-02,
                8.688697172522256984e-02,
                1.115650800742149090e-01,
                1.432523984300950459e-01,
                1.839397205857211670e-01,
                2.361832763705073446e-01,
                3.032653298563167121e-01,
                3.894003915357024392e-01,
                5.000000000000000000e-01,
                3.894003915357024392e-01,
                3.032653298563167121e-01,
                2.361832763705073446e-01,
                1.839397205857211670e-01,
                1.432523984300950459e-01,
                1.115650800742149090e-01,
                8.688697172522256984e-02,
                6.766764161830635116e-02,
                5.269961228093216626e-02,
                4.104249931194939999e-02,
                3.196393060335378511e-02,
                2.489353418393197223e-02,
                1.938710391586100437e-02,
                1.509869171115925043e-02,
                1.175887292800455348e-02,
                9.157819444367089334e-03,
                7.132116954499627751e-03,
                5.554498269121153041e-03,
                4.325847601560317036e-03,
                3.368973499542733500e-03,
                2.623759199590692311e-03,
                2.043385719232033298e-03,
                1.591390398254833444e-03,
                1.239376088333179245e-03,
                9.652270681138546511e-04,
                7.517195964887861837e-04,
                5.854398103955872146e-04,
                4.559409827772581219e-04,
                3.550871944212745141e-04,
                2.765421850739168135e-04,
                2.153712702878437645e-04,

            ];
            ddouble[] expected_dist2 = [
                8.297844727977324655e-03,
                9.018961037136934838e-03,
                9.802745273738313128e-03,
                1.065464353445126112e-02,
                1.158057520380025712e-02,
                1.258697408479562789e-02,
                1.368083310398313275e-02,
                1.486975290154335209e-02,
                1.616199464406751010e-02,
                1.756653742697738760e-02,
                1.909314066544795518e-02,
                2.075241190735382588e-02,
                2.255588053943545154e-02,
                2.451607789882946642e-02,
                2.664662434661565008e-02,
                2.896232390840752444e-02,
                3.147926713959364048e-02,
                3.421494292998487896e-02,
                3.718836002473830532e-02,
                4.042017910594145058e-02,
                4.393285635262112604e-02,
                4.775079947669835095e-02,
                5.190053731909961138e-02,
                5.641090418445703775e-02,
                6.131324019524039132e-02,
                6.664160905747455732e-02,
                7.243303475117969514e-02,
                7.872775879016911948e-02,
                8.556951983876533163e-02,
                9.300585762834118198e-02,
                1.010884432854389087e-01,
                1.098734383667406278e-01,
                1.194218850956315497e-01,
                1.298001305119008131e-01,
                1.410802874817690122e-01,
                1.533407357715538821e-01,
                1.666666666666666574e-01,
                1.533407357715538821e-01,
                1.410802874817690122e-01,
                1.298001305119008131e-01,
                1.194218850956315497e-01,
                1.098734383667406278e-01,
                1.010884432854389087e-01,
                9.300585762834118198e-02,
                8.556951983876533163e-02,
                7.872775879016911948e-02,
                7.243303475117969514e-02,
                6.664160905747455732e-02,
                6.131324019524039132e-02,
                5.641090418445703775e-02,
                5.190053731909961138e-02,
                4.775079947669835095e-02,
                4.393285635262112604e-02,
                4.042017910594145058e-02,
                3.718836002473830532e-02,
                3.421494292998487896e-02,
                3.147926713959364048e-02,
                2.896232390840752444e-02,
                2.664662434661565008e-02,
                2.451607789882946642e-02,
                2.255588053943545154e-02,
                2.075241190735382588e-02,
                1.909314066544795518e-02,
                1.756653742697738760e-02,

            ];

            foreach ((LaplaceDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (-8, 0); i < expecteds.Length; x += 0.25, i++) {
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
            ddouble[] expected_dist1 = [
                1.677313139512559266e-04,
                2.153712702878437645e-04,
                2.765421850739168135e-04,
                3.550871944212745141e-04,
                4.559409827772581219e-04,
                5.854398103955872146e-04,
                7.517195964887861837e-04,
                9.652270681138546511e-04,
                1.239376088333179245e-03,
                1.591390398254833444e-03,
                2.043385719232033298e-03,
                2.623759199590692311e-03,
                3.368973499542733500e-03,
                4.325847601560317036e-03,
                5.554498269121153041e-03,
                7.132116954499627751e-03,
                9.157819444367089334e-03,
                1.175887292800455348e-02,
                1.509869171115925043e-02,
                1.938710391586100437e-02,
                2.489353418393197223e-02,
                3.196393060335378511e-02,
                4.104249931194939999e-02,
                5.269961228093216626e-02,
                6.766764161830635116e-02,
                8.688697172522256984e-02,
                1.115650800742149090e-01,
                1.432523984300950459e-01,
                1.839397205857211670e-01,
                2.361832763705073446e-01,
                3.032653298563167121e-01,
                3.894003915357024392e-01,
                5.000000000000000000e-01,
                6.105996084642975053e-01,
                6.967346701436832879e-01,
                7.638167236294925999e-01,
                8.160602794142788330e-01,
                8.567476015699049263e-01,
                8.884349199257850493e-01,
                9.131130282747774718e-01,
                9.323323583816935933e-01,
                9.473003877190678823e-01,
                9.589575006880506347e-01,
                9.680360693966462149e-01,
                9.751064658160680798e-01,
                9.806128960841390407e-01,
                9.849013082888407045e-01,
                9.882411270719954066e-01,
                9.908421805556328898e-01,
                9.928678830455003324e-01,
                9.944455017308788669e-01,
                9.956741523984397046e-01,
                9.966310265004573177e-01,
                9.973762408004093194e-01,
                9.979566142807679840e-01,
                9.984086096017451561e-01,
                9.987606239116668672e-01,
                9.990347729318861392e-01,
                9.992482804035112132e-01,
                9.994145601896043951e-01,
                9.995440590172227635e-01,
                9.996449128055787670e-01,
                9.997234578149261086e-01,
                9.997846287297121881e-01,

            ];
            ddouble[] expected_dist2 = [
                2.489353418393197223e-02,
                2.705688311141080452e-02,
                2.940823582121494112e-02,
                3.196393060335378511e-02,
                3.474172561140077137e-02,
                3.776092225438688194e-02,
                4.104249931194939999e-02,
                4.460925870463005455e-02,
                4.848598393220252684e-02,
                5.269961228093216626e-02,
                5.727942199634386555e-02,
                6.225723572206148110e-02,
                6.766764161830635116e-02,
                7.354823369648839926e-02,
                7.993987303984695370e-02,
                8.688697172522256984e-02,
                9.443780141878091450e-02,
                1.026448287899546369e-01,
                1.115650800742149090e-01,
                1.212605373178243517e-01,
                1.317985690578633851e-01,
                1.432523984300950459e-01,
                1.557016119572988411e-01,
                1.692327125533711063e-01,
                1.839397205857211670e-01,
                1.999248271724236858e-01,
                2.172991042535390993e-01,
                2.361832763705073446e-01,
                2.567085595162960088e-01,
                2.790175728850235459e-01,
                3.032653298563167121e-01,
                3.296203151002218834e-01,
                3.582656552868946354e-01,
                3.894003915357024392e-01,
                4.232408624453070645e-01,
                4.600222073146616464e-01,
                5.000000000000000000e-01,
                5.399777926853384091e-01,
                5.767591375546929910e-01,
                6.105996084642975053e-01,
                6.417343447131054202e-01,
                6.703796848997780611e-01,
                6.967346701436832879e-01,
                7.209824271149765096e-01,
                7.432914404837039912e-01,
                7.638167236294925999e-01,
                7.827008957464609562e-01,
                8.000751728275763419e-01,
                8.160602794142788330e-01,
                8.307672874466288659e-01,
                8.442983880427011867e-01,
                8.567476015699049263e-01,
                8.682014309421366427e-01,
                8.787394626821756205e-01,
                8.884349199257850493e-01,
                8.973551712100453770e-01,
                9.055621985812191133e-01,
                9.131130282747774718e-01,
                9.200601269601530463e-01,
                9.264517663035115591e-01,
                9.323323583816935933e-01,
                9.377427642779385675e-01,
                9.427205780036561622e-01,
                9.473003877190678823e-01,

            ];

            foreach ((LaplaceDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (-8, 0); i < expecteds.Length; x += 0.25, i++) {
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