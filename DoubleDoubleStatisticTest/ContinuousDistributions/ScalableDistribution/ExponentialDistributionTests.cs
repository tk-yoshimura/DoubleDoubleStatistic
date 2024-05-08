using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class ExponentialDistributionTests {
        readonly ExponentialDistribution dist1 = new(theta: 1);
        readonly ExponentialDistribution dist2 = new(theta: 2);

        ExponentialDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (ExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
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
            foreach (ExponentialDistribution dist in Dists) {
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
            foreach (ExponentialDistribution dist in Dists) {
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
            foreach (ExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (ExponentialDistribution dist in Dists) {
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
            foreach (ExponentialDistribution dist in Dists) {
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
            foreach (ExponentialDistribution dist in Dists) {
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
            foreach (ExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (ExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (ExponentialDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (ExponentialDistribution dist in Dists) {
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
            foreach (ExponentialDistribution dist in Dists) {
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
            foreach (ExponentialDistribution dist in Dists) {
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

            foreach (ExponentialDistribution dist in Dists) {

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
            foreach (ExponentialDistribution dist in Dists) {
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
            ddouble[] expected_dist1 = [
                1.000000000000000000e+00,
                7.788007830714048785e-01,
                6.065306597126334243e-01,
                4.723665527410146892e-01,
                3.678794411714423340e-01,
                2.865047968601900918e-01,
                2.231301601484298180e-01,
                1.737739434504451397e-01,
                1.353352832366127023e-01,
                1.053992245618643325e-01,
                8.208499862389879997e-02,
                6.392786120670757022e-02,
                4.978706836786394446e-02,
                3.877420783172200874e-02,
                3.019738342231850087e-02,
                2.351774585600910697e-02,
                1.831563888873417867e-02,
                1.426423390899925550e-02,
                1.110899653824230608e-02,
                8.651695203120634073e-03,
                6.737946999085467001e-03,
                5.247518399181384623e-03,
                4.086771438464066597e-03,
                3.182780796509666887e-03,
                2.478752176666358491e-03,
                1.930454136227709302e-03,
                1.503439192977572367e-03,
                1.170879620791174429e-03,
                9.118819655545162437e-04,
                7.101743888425490281e-04,
                5.530843701478336269e-04,
                4.307425405756875290e-04,
                3.354626279025118532e-04,
                2.612585573016675401e-04,
                2.034683690106441685e-04,
                1.584613251157512555e-04,
                1.234098040866795612e-04,
                9.611165206139469499e-05,
                7.485182988770059808e-05,
                5.829466373086881077e-05,
                4.539992976248485417e-05,
                3.535750085040998108e-05,
                2.753644934974715811e-05,
                2.144540831658916372e-05,
                1.670170079024565926e-05,
                1.300729765406762047e-05,
                1.013009359863071055e-05,
                7.889324827200222879e-06,
                6.144212353328209806e-06,
                4.785117392129008755e-06,
                3.726653172078670938e-06,
                2.902320408650404051e-06,
                2.260329406981054240e-06,
                1.760346312156169343e-06,
                1.370959086384084452e-06,
                1.067704010034782666e-06,
                8.315287191035678771e-07,
                6.475952175842209272e-07,
                5.043476625678880272e-07,
                3.927863545481039239e-07,
                3.059023205018257859e-07,
                2.382369667501817973e-07,
                1.855391362615978357e-07,
                1.444980246109244809e-07,
            ];
            ddouble[] expected_dist2 = [
                5.000000000000000000e-01,
                4.412484512922977276e-01,
                3.894003915357024392e-01,
                3.436446393954861178e-01,
                3.032653298563167121e-01,
                2.676307142594951394e-01,
                2.361832763705073446e-01,
                2.084310098392541943e-01,
                1.839397205857211670e-01,
                1.623262336791748695e-01,
                1.432523984300950459e-01,
                1.264197979023732321e-01,
                1.115650800742149090e-01,
                9.845583760209702939e-02,
                8.688697172522256984e-02,
                7.667748342246423487e-02,
                6.766764161830635116e-02,
                5.971648413335980954e-02,
                5.269961228093216626e-02,
                4.650724460533174620e-02,
                4.104249931194939999e-02,
                3.621987851712572815e-02,
                3.196393060335378511e-02,
                2.820806975188867513e-02,
                2.489353418393197223e-02,
                2.196846681170371018e-02,
                1.938710391586100437e-02,
                1.710905915583301601e-02,
                1.509869171115925043e-02,
                1.332454866817774260e-02,
                1.175887292800455348e-02,
                1.037716893684987110e-02,
                9.157819444367089334e-03,
                8.081747294082937058e-03,
                7.132116954499627751e-03,
                6.294071121216999153e-03,
                5.554498269121153041e-03,
                4.901827517910913902e-03,
                4.325847601560317036e-03,
                3.817547109429980830e-03,
                3.368973499542733500e-03,
                2.973108678236047100e-03,
                2.623759199590692311e-03,
                2.315459366766622895e-03,
                2.043385719232033298e-03,
                1.803281568007865246e-03,
                1.591390398254833444e-03,
                1.404397097262756380e-03,
                1.239376088333179245e-03,
                1.093745559091442547e-03,
                9.652270681138546511e-04,
                8.518098979012869894e-04,
                7.517195964887861837e-04,
                6.633902155134957512e-04,
                5.854398103955872146e-04,
                5.166488193238185115e-04,
                4.559409827772581219e-04,
                4.023665050623066231e-04,
                3.550871944212745141e-04,
                3.133633492242287872e-04,
                2.765421850739168135e-04,
                2.440476217617074985e-04,
                2.153712702878437645e-04,
                1.900644789347318410e-04,
            ];

            foreach ((ExponentialDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.25, i++) {
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
                0.000000000000000000e+00,
                2.211992169285951493e-01,
                3.934693402873665757e-01,
                5.276334472589853108e-01,
                6.321205588285576660e-01,
                7.134952031398098526e-01,
                7.768698398515702097e-01,
                8.262260565495548326e-01,
                8.646647167633872977e-01,
                8.946007754381356536e-01,
                9.179150013761011584e-01,
                9.360721387932924298e-01,
                9.502129316321360486e-01,
                9.612257921682779704e-01,
                9.698026165776815199e-01,
                9.764822541439909243e-01,
                9.816843611112657797e-01,
                9.857357660910007757e-01,
                9.888910034617577338e-01,
                9.913483047968794093e-01,
                9.932620530009145243e-01,
                9.947524816008186388e-01,
                9.959132285615359681e-01,
                9.968172192034903123e-01,
                9.975212478233336233e-01,
                9.980695458637722783e-01,
                9.984965608070224263e-01,
                9.988291203792087902e-01,
                9.990881180344455270e-01,
                9.992898256111574229e-01,
                9.994469156298522172e-01,
                9.995692574594242652e-01,
                9.996645373720974836e-01,
                9.997387414426983643e-01,
                9.997965316309893602e-01,
                9.998415386748842337e-01,
                9.998765901959133506e-01,
                9.999038883479386408e-01,
                9.999251481701123545e-01,
                9.999417053362691377e-01,
                9.999546000702375093e-01,
                9.999646424991496252e-01,
                9.999724635506502812e-01,
                9.999785545916833884e-01,
                9.999832982992097330e-01,
                9.999869927023459226e-01,
                9.999898699064013741e-01,
                9.999921106751727962e-01,
                9.999938557876466572e-01,
                9.999952148826078968e-01,
                9.999962733468279463e-01,
                9.999970976795913291e-01,
                9.999977396705930222e-01,
                9.999982396536878371e-01,
                9.999986290409136336e-01,
                9.999989322959900173e-01,
                9.999991684712808970e-01,
                9.999993524047824378e-01,
                9.999994956523374778e-01,
                9.999996072136454472e-01,
                9.999996940976795257e-01,
                9.999997617630332902e-01,
                9.999998144608637229e-01,
                9.999998555019753432e-01,
            ];
            ddouble[] expected_dist2 = [
                0.000000000000000000e+00,
                1.175030974154046004e-01,
                2.211992169285951493e-01,
                3.127107212090278199e-01,
                3.934693402873665757e-01,
                4.647385714810097213e-01,
                5.276334472589853108e-01,
                5.831379803214915558e-01,
                6.321205588285576660e-01,
                6.753475326416502611e-01,
                7.134952031398098526e-01,
                7.471604041952535358e-01,
                7.768698398515702097e-01,
                8.030883247958059412e-01,
                8.262260565495548326e-01,
                8.466450331550715580e-01,
                8.646647167633872977e-01,
                8.805670317332803254e-01,
                8.946007754381356536e-01,
                9.069855107893365354e-01,
                9.179150013761011584e-01,
                9.275602429657485715e-01,
                9.360721387932924298e-01,
                9.435838604962226706e-01,
                9.502129316321360486e-01,
                9.560630663765925519e-01,
                9.612257921682779704e-01,
                9.657818816883340096e-01,
                9.698026165776815199e-01,
                9.733509026636445460e-01,
                9.764822541439909243e-01,
                9.792456621263002647e-01,
                9.816843611112657797e-01,
                9.838365054118340947e-01,
                9.857357660910007757e-01,
                9.874118577575660138e-01,
                9.888910034617577338e-01,
                9.901963449641781878e-01,
                9.913483047968794093e-01,
                9.923649057811400054e-01,
                9.932620530009145243e-01,
                9.940537826435279500e-01,
                9.947524816008186388e-01,
                9.953690812664667265e-01,
                9.959132285615359681e-01,
                9.963934368639842942e-01,
                9.968172192034903123e-01,
                9.971912058054744898e-01,
                9.975212478233336233e-01,
                9.978125088818171617e-01,
                9.980695458637722783e-01,
                9.982963802041974199e-01,
                9.984965608070224263e-01,
                9.986732195689730363e-01,
                9.988291203792087902e-01,
                9.989667023613523211e-01,
                9.990881180344455270e-01,
                9.991952669898753747e-01,
                9.992898256111574229e-01,
                9.993732733015515679e-01,
                9.994469156298522172e-01,
                9.995119047564765546e-01,
                9.995692574594242652e-01,
                9.996198710421305700e-01,
            ];

            foreach ((ExponentialDistribution dist, ddouble[] expecteds) in new[]{
                (dist1, expected_dist1), (dist2, expected_dist2)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.25, i++) {
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