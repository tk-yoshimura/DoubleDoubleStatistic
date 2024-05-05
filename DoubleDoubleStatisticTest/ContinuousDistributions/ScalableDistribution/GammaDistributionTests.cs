using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class GammaDistributionTests {
        readonly GammaDistribution dist_k1theta1 = new(k: 1, theta: 1);
        readonly GammaDistribution dist_k2theta1 = new(k: 2, theta: 1);
        readonly GammaDistribution dist_k1theta2 = new(k: 1, theta: 2);
        readonly GammaDistribution dist_k2theta2 = new(k: 2, theta: 2);
        readonly GammaDistribution dist_k3theta4 = new(k: 3, theta: 4);

        GammaDistribution[] Dists => [
            dist_k1theta1,
            dist_k2theta1,
            dist_k1theta2,
            dist_k2theta2,
            dist_k3theta4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.K}");
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
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (GammaDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            foreach (GammaDistribution dist in Dists) {
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
            ddouble[] expected_dist_k1theta1 = [
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
            ddouble[] expected_dist_k2theta1 = [
                0.000000000000000000e+00,
                1.947001957678512196e-01,
                3.032653298563167121e-01,
                3.542749145557611001e-01,
                3.678794411714423340e-01,
                3.581309960752376287e-01,
                3.346952402226447409e-01,
                3.041044010382789597e-01,
                2.706705664732254046e-01,
                2.371482552641947517e-01,
                2.052124965597469930e-01,
                1.758016183184458181e-01,
                1.493612051035918542e-01,
                1.260161754530965683e-01,
                1.056908419781147634e-01,
                8.819154696003417282e-02,
                7.326255555493672855e-02,
                6.062299411324682807e-02,
                4.999048442209038517e-02,
                4.109555221482301401e-02,
                3.368973499542733674e-02,
                2.754947159570227469e-02,
                2.247724291155236975e-02,
                1.830098957993058287e-02,
                1.487251305999814401e-02,
                1.206533835142318406e-02,
                9.772354754354214967e-03,
                7.903437440340426584e-03,
                6.383173758881614465e-03,
                5.148764319108482812e-03,
                4.148132776108753340e-03,
                3.338254689461577808e-03,
                2.683701023220095259e-03,
                2.155383097738755756e-03,
                1.729481136590475419e-03,
                1.386536594762823689e-03,
                1.110688236780116647e-03,
                8.890327815679011049e-04,
                7.110923839331554242e-04,
                5.683729713759712353e-04,
                4.539992976248486095e-04,
                3.624143837167023230e-04,
                2.891327181723449941e-04,
                2.305381394033333533e-04,
                1.837187086927022654e-04,
                1.463320986082606841e-04,
                1.164960763842531611e-04,
                9.269956671960267092e-05,
                7.373054823993848718e-05,
                5.861768805358037016e-05,
                4.658316465098339689e-05,
                3.700458521029266013e-05,
                2.938428229075369495e-05,
                2.332458863606925184e-05,
                1.850794766618513957e-05,
                1.468093013797827633e-05,
                1.164140206744995240e-05,
                9.228231850575139742e-06,
                7.313041107234380206e-06,
                5.793598729584527954e-06,
                4.588534807527388906e-06,
                3.633113742940268909e-06,
                2.875856612054768743e-06,
                2.275843887622061785e-06,

            ];
            ddouble[] expected_dist_k1theta2 = [
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
            ddouble[] expected_dist_k2theta2 = [
                0.000000000000000000e+00,
                5.515605641153722288e-02,
                9.735009788392560981e-02,
                1.288667397733072872e-01,
                1.516326649281583561e-01,
                1.672691964121844344e-01,
                1.771374572778805501e-01,
                1.823771336093474027e-01,
                1.839397205857211670e-01,
                1.826170128890717004e-01,
                1.790654980376188143e-01,
                1.738272221157631803e-01,
                1.673476201113223705e-01,
                1.599907361034076658e-01,
                1.520522005191394799e-01,
                1.437702814171204213e-01,
                1.353352832366127023e-01,
                1.268975287833896126e-01,
                1.185741276320973758e-01,
                1.104547059376628920e-01,
                1.026062482798734965e-01,
                9.507718110745504159e-02,
                8.790080915922290905e-02,
                8.109820053667993189e-02,
                7.468060255179592710e-02,
                6.865145878657409606e-02,
                6.300808772654828416e-02,
                5.774307465093642122e-02,
                5.284542098905738172e-02,
                4.830148892214432971e-02,
                4.409577348001708641e-02,
                4.021152963029324401e-02,
                3.663127777746836428e-02,
                3.333720758809212426e-02,
                3.031149705662341404e-02,
                2.753656115532437249e-02,
                2.499524221104519259e-02,
                2.267095227033797886e-02,
                2.054777610741150701e-02,
                1.861054215847115709e-02,
                1.684486749771366837e-02,
                1.523718197595974079e-02,
                1.377473579785113734e-02,
                1.244559409637059980e-02,
                1.123862145577618488e-02,
                1.014345882004424187e-02,
                9.150494789965291434e-03,
                8.250832946418696495e-03,
                7.436256529999072003e-03,
                6.699191549435086143e-03,
                6.032669175711592030e-03,
                5.430288099120705696e-03,
                4.886177377177107484e-03,
                4.394960177776910219e-03,
                3.951718720170213292e-03,
                3.551960632851252145e-03,
                3.191586879440807233e-03,
                2.866861348568934195e-03,
                2.574382159554241406e-03,
                2.311054700528687462e-03,
                2.074066388054376670e-03,
                1.860863115933019405e-03,
                1.669127344730788904e-03,
                1.496757771611013834e-03,

            ];
            ddouble[] expected_dist_k3theta4 = [
                0.000000000000000000e+00,
                4.586977845768925384e-04,
                1.723626762860538432e-03,
                3.643194367003711916e-03,
                6.084381117745350613e-03,
                8.930854845540060066e-03,
                1.208125685374755905e-02,
                1.544764540769858786e-02,
                1.895408311601979798e-02,
                2.253535586093982634e-02,
                2.613581193940382394e-02,
                2.970831100316594087e-02,
                3.321327323960260314e-02,
                3.661782002133911335e-02,
                3.989499797704475215e-02,
                4.302307910267567231e-02,
                4.598493014643028481e-02,
                4.876744506579377275e-02,
                5.136103487505143395e-02,
                5.375916968002559387e-02,
                5.595796813675588122e-02,
                5.795582997537603032e-02,
                5.975310760229360363e-02,
                6.135181313549675775e-02,
                6.275535754174588199e-02,
                6.396831883273250419e-02,
                6.499623654200936684e-02,
                6.584542994747306921e-02,
                6.652283772712351551e-02,
                6.703587694032177302e-02,
                6.739231941427518968e-02,
                6.760018378735187239e-02,
                6.766764161830635116e-02,
                6.760293611473962294e-02,
                6.741431216617571609e-02,
                6.710995648799858437e-02,
                6.669794679305476870e-02,
                6.618620900882724034e-02,
                6.558248165048732869e-02,
                6.489428654454804102e-02,
                6.412890517492093878e-02,
                6.329335999352557474e-02,
                6.239440010176736801e-02,
                6.143849076769719586e-02,
                6.043180629696574130e-02,
                5.938022582421804973e-02,
                5.828933163573869064e-02,
                5.716440967432272896e-02,
                5.601045191384693839e-02,
                5.483216032416284025e-02,
                5.363395217701100171e-02,
                5.241996647092745515e-02,
                5.119407127782046873e-02,
                4.995987183625787142e-02,
                4.872071923672761234e-02,
                4.747971956240238517e-02,
                4.623974336542520380e-02,
                4.500343537359270674e-02,
                4.377322433569329208e-02,
                4.255133292578674314e-02,
                4.133978763751600810e-02,
                4.014042860922721612e-02,
                3.895491932934658230e-02,
                3.778475617921163432e-02,

            ];

            foreach ((GammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1theta1, expected_dist_k1theta1), (dist_k2theta1, expected_dist_k2theta1),
                (dist_k1theta2, expected_dist_k1theta2), (dist_k2theta2, expected_dist_k2theta2),
                (dist_k3theta4, expected_dist_k3theta4),
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
            ddouble[] expected_dist_k1theta1 = [
                0.000000000000000000e+00,
                2.211992169285951215e-01,
                3.934693402873665202e-01,
                5.276334472589850888e-01,
                6.321205588285576660e-01,
                7.134952031398098526e-01,
                7.768698398515702097e-01,
                8.262260565495549436e-01,
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
            ddouble[] expected_dist_k2theta1 = [
                0.000000000000000000e+00,
                2.649902116074391231e-02,
                9.020401043104986361e-02,
                1.733585327032242662e-01,
                2.642411176571152764e-01,
                3.553642070645721684e-01,
                4.421745996289251912e-01,
                5.221216555112759838e-01,
                5.939941502901615600e-01,
                6.574525201739409574e-01,
                7.127025048163542209e-01,
                7.602705204748466672e-01,
                8.008517265285441944e-01,
                8.352096167151814576e-01,
                8.641117745995667843e-01,
                8.882907071839567514e-01,
                9.084218055563291205e-01,
                9.251127719777538783e-01,
                9.389005190396673139e-01,
                9.502527525820563259e-01,
                9.595723180054872570e-01,
                9.672030100051163259e-01,
                9.734359856499835706e-01,
                9.785162296235597745e-01,
                9.826487347633354741e-01,
                9.860042075123490735e-01,
                9.887242060526681975e-01,
                9.909256829388684018e-01,
                9.927049442755638831e-01,
                9.941410612920489331e-01,
                9.952987828537434023e-01,
                9.962310027699626813e-01,
                9.969808363488773528e-01,
                9.975833583449595920e-01,
                9.980670504943989219e-01,
                9.984550020801213899e-01,
                9.987659019591331733e-01,
                9.990148555663707075e-01,
                9.992140557861791672e-01,
                9.993733323648931188e-01,
                9.995006007726127129e-01,
                9.996022281154328670e-01,
                9.996833308324779566e-01,
                9.997480164522800639e-01,
                9.997995795905170180e-01,
                9.998406606037376632e-01,
                9.998733738300170648e-01,
                9.998994111084531511e-01,
                9.999201252394067652e-01,
                9.999365971945542997e-01,
                9.999496901821769423e-01,
                9.999600930943810528e-01,
                9.999683553883023102e-01,
                9.999749150650517615e-01,
                9.999801210932474094e-01,
                9.999842513658520060e-01,
                9.999875270692134555e-01,
                9.999901241729318713e-01,
                9.999921826112302403e-01,
                9.999938136149159051e-01,
                9.999951055628719665e-01,
                9.999961286492903278e-01,
                9.999969386042516506e-01,
                9.999975796580877763e-01,

            ];
            ddouble[] expected_dist_k1theta2 = [
                0.000000000000000000e+00,
                1.175030974154046004e-01,
                2.211992169285951215e-01,
                3.127107212090278754e-01,
                3.934693402873665202e-01,
                4.647385714810097768e-01,
                5.276334472589850888e-01,
                5.831379803214914448e-01,
                6.321205588285576660e-01,
                6.753475326416502611e-01,
                7.134952031398098526e-01,
                7.471604041952535358e-01,
                7.768698398515702097e-01,
                8.030883247958059412e-01,
                8.262260565495549436e-01,
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
            ddouble[] expected_dist_k2theta2 = [
                0.000000000000000000e+00,
                7.190984592330174584e-03,
                2.649902116074391231e-02,
                5.497724166241324539e-02,
                9.020401043104986361e-02,
                1.302001786566408525e-01,
                1.733585327032242662e-01,
                2.183837131027967227e-01,
                2.642411176571152764e-01,
                3.101135068635068048e-01,
                3.553642070645721684e-01,
                3.995059599637270087e-01,
                4.421745996289251912e-01,
                4.831068525889904430e-01,
                5.221216555112759838e-01,
                5.591044703208305489e-01,
                5.939941502901615600e-01,
                6.267719741665012112e-01,
                6.574525201739409574e-01,
                6.860760989140107791e-01,
                7.127025048163542209e-01,
                7.374058807508385716e-01,
                7.602705204748466672e-01,
                7.813874594228628068e-01,
                8.008517265285441944e-01,
                8.187601488034443875e-01,
                8.352096167151814576e-01,
                8.502957323864610562e-01,
                8.641117745995667843e-01,
                8.767479248193559282e-01,
                8.882907071839567514e-01,
                8.988226028657138045e-01,
                9.084218055563291205e-01,
                9.171620902356498739e-01,
                9.251127719777538783e-01,
                9.323387354469172550e-01,
                9.389005190396673139e-01,
                9.448544404235021954e-01,
                9.502527525820563259e-01,
                9.551438214641977398e-01,
                9.595723180054872570e-01,
                9.635794186916084580e-01,
                9.672030100051163259e-01,
                9.704778930737255616e-01,
                9.734359856499835706e-01,
                9.761065192238957966e-01,
                9.785162296235597745e-01,
                9.806895399126370760e-01,
                9.826487347633354741e-01,
                9.844141257829469582e-01,
                9.860042075123490735e-01,
                9.874358040059559860e-01,
                9.887242060526681975e-01,
                9.898832992134192210e-01,
                9.909256829388684018e-01,
                9.918627810956498481e-01,
                9.927049442755638831e-01,
                9.934615442927374751e-01,
                9.941410612920489331e-01,
                9.947511639004941175e-01,
                9.952987828537434023e-01,
                9.957901785246104920e-01,
                9.962310027699626813e-01,
                9.966263554989085449e-01,

            ];
            ddouble[] expected_dist_k3theta4 = [
                0.000000000000000000e+00,
                3.882962237440737890e-05,
                2.964775408880204223e-04,
                9.551446927597426933e-04,
                2.161496689762512895e-03,
                4.031067625372406794e-03,
                6.652214247422994463e-03,
                1.008966162911084478e-02,
                1.438767796697068037e-02,
                1.957291291417348825e-02,
                2.565693089902557067e-02,
                3.263846816137334211e-02,
                4.050543974481386755e-02,
                4.923672039268629091e-02,
                5.880372119461777652e-02,
                6.917178190299930196e-02,
                8.030139707139417882e-02,
                9.214929254681498982e-02,
                1.046693673633011384e-01,
                1.178135147244554343e-01,
                1.315323345175487824e-01,
                1.457757486622581755e-01,
                1.604935295545526497e-01,
                1.756357509128337602e-01,
                1.911531694619418298e-01,
                2.069975450943819240e-01,
                2.231219064209531699e-01,
                2.394807679592275418e-01,
                2.560303046027818108e-01,
                2.727284884619547989e-01,
                2.895351926637298456e-01,
                3.064122662400736719e-01,
                3.233235838169365439e-01,
                3.402350734366825780e-01,
                3.571147255017984024e-01,
                3.739325855144102206e-01,
                3.906607330017219382e-01,
                4.072732487595164819e-01,
                4.237461723120614088e-01,
                4.400574512750505818e-01,
                4.561868841166703548e-01,
                4.721160576387694596e-01,
                4.878282803437689608e-01,
                5.033085127119861779e-01,
                5.185432952869835077e-01,
                5.335206753522663270e-01,
                5.482301328799077389e-01,
                5.626625063394955495e-01,
                5.768099188731565796e-01,
                5.906657052685152509e-01,
                6.042243400954002697e-01,
                6.174813673132628633e-01,
                6.304333316038995827e-01,
                6.430777116375850921e-01,
                6.554128554395506345e-01,
                6.674379179873777268e-01,
                6.791528011378658025e-01,
                6.905580959538131225e-01,
                7.016550274765827044e-01,
                7.124454019689148154e-01,
                7.229315566338927468e-01,
                7.331163117999781598e-01,
                7.430029255483273642e-01,
                7.525950507469507667e-01,
            ];

            foreach ((GammaDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1theta1, expected_dist_k1theta1), (dist_k2theta1, expected_dist_k2theta1),
                (dist_k1theta2, expected_dist_k1theta2), (dist_k2theta2, expected_dist_k2theta2),
                (dist_k3theta4, expected_dist_k3theta4),
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