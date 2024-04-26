using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ScalableDistribution {
    [TestClass()]
    public class InverseGaussDistributionTests {
        readonly InverseGaussDistribution dist_mu1lambda1 = new(mu: 1, lambda: 1);
        readonly InverseGaussDistribution dist_mu2lambda1 = new(mu: 2, lambda: 1);
        readonly InverseGaussDistribution dist_mu1lambda2 = new(mu: 1, lambda: 2);
        readonly InverseGaussDistribution dist_mu2lambda2 = new(mu: 2, lambda: 2);
        readonly InverseGaussDistribution dist_mu3lambda4 = new(mu: 3, lambda: 4);

        InverseGaussDistribution[] Dists => [
            dist_mu1lambda1,
            dist_mu2lambda1,
            dist_mu1lambda2,
            dist_mu2lambda2,
            dist_mu3lambda4,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Lambda={dist.Lambda}");
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
            foreach (InverseGaussDistribution dist in Dists) {
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
            foreach (InverseGaussDistribution dist in Dists) {
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
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (InverseGaussDistribution dist in Dists) {
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
            foreach (InverseGaussDistribution dist in Dists) {
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
            foreach (InverseGaussDistribution dist in Dists) {
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
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (InverseGaussDistribution dist in Dists) {
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
            foreach (InverseGaussDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-28);
                    }
                }

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-280"; p /= 10) {
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
            foreach (InverseGaussDistribution dist in Dists) {
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

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-280"; p /= 10) {
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"quantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (InverseGaussDistribution dist in Dists) {
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
            ddouble[] expected_dist_mu1lambda1 = [
                0.000000000000000000e+00,
                1.036140765327133950e+00,
                8.787825789354447581e-01,
                5.891455034058555862e-01,
                3.989422804014327029e-01,
                2.784118295636137663e-01,
                1.997937831333951031e-01,
                1.467424160755297080e-01,
                1.098478223669305948e-01,
                8.352944456601557599e-02,
                6.435281305249991590e-02,
                5.012887535439480152e-02,
                3.941835796981973949e-02,
                3.124902255392518827e-02,
                2.494853617603773333e-02,
                2.004239244292750979e-02,
                1.618969945823646797e-02,
                1.314171093271665480e-02,
                1.071440973308502784e-02,
                8.769998723340765093e-03,
                7.204168934430731953e-03,
                5.937193432662460010e-03,
                4.907602417365810807e-03,
                4.067629680936425858e-03,
                3.379893528659050263e-03,
                2.814950142747340515e-03,
                2.349471486503027645e-03,
                1.964875558634609533e-03,
                1.646287825811404253e-03,
                1.381747531512703651e-03,
                1.161596728020189906e-03,
                9.780068075341530592e-04,
                8.246093114085995554e-04,
                6.962063968934366224e-04,
                5.885425638618871445e-04,
                4.981237883213677662e-04,
                4.220735564369423355e-04,
                3.580177776539568250e-04,
                3.039924142295118560e-04,
                2.583690646675928875e-04,
                2.197948003186267626e-04,
                1.871433645698393445e-04,
                1.594754658697860672e-04,
                1.360063757282683344e-04,
                1.160794151360577703e-04,
                9.914420310432021195e-05,
                8.473876836709907337e-05,
                7.247480414217717410e-05,
                6.202548713590911515e-05,
                5.311539403566850316e-05,
                4.551213794020290081e-05,
                3.901941844570693137e-05,
                3.347123623334795887e-05,
                2.872706894436740117e-05,
                2.466784218516771378e-05,
                2.119255948163845706e-05,
                1.821547931876543155e-05,
                1.566374717883371959e-05,
                1.347540661700952356e-05,
                1.159772659245619569e-05,
                9.985793069201679729e-06,
                8.601321763960161383e-06,
                7.411656209265385145e-06,
                6.388921310276009443e-06,
            ];
            ddouble[] expected_dist_mu1lambda2 = [
                0.000000000000000000e+00,
                4.757211568945174385e-01,
                9.678828980765734613e-01,
                7.991751325335804124e-01,
                5.641895835477562793e-01,
                3.840124634976435636e-01,
                2.599595409636549226e-01,
                1.767149677246939166e-01,
                1.209853622595716827e-01,
                8.347528226774900140e-02,
                5.802965878713834541e-02,
                4.062373919236163256e-02,
                2.862093862528111710e-02,
                2.028162966750505711e-02,
                1.444764748647170763e-02,
                1.034072424244746047e-02,
                7.433143076476834976e-03,
                5.364041724500038036e-03,
                3.884724310179186145e-03,
                2.822567955769658514e-03,
                2.056968850540707040e-03,
                1.503163715533159701e-03,
                1.101254534531053830e-03,
                8.087038568981529379e-04,
                5.951656347464258151e-04,
                4.389007999938188457e-04,
                3.242763861535133776e-04,
                2.400107489341667729e-04,
                1.779358778510661137e-04,
                1.321202347423509936e-04,
                9.824435080151186592e-05,
                7.315453284076124413e-05,
                5.454266844036000962e-05,
                4.071569698115229613e-05,
                3.042909767449812440e-05,
                2.276625396289778854e-05,
                1.705080877219965110e-05,
                1.278282124964620158e-05,
                9.592147156891585598e-06,
                7.204318118713548677e-06,
                5.415514964427321914e-06,
                4.074174221199812821e-06,
                3.067443977082190384e-06,
                2.311192736361016781e-06,
                1.742626773717489291e-06,
                1.314826257300864994e-06,
                9.926953498398672652e-07,
                7.499560840063275041e-07,
                5.669139825773298525e-07,
                4.287951656150410112e-07,
                3.245072660955757244e-07,
                2.457150703583711327e-07,
                1.861501825364718718e-07,
                1.410949216768250763e-07,
                1.069960900045910472e-07,
                8.117565485245966300e-08,
                6.161382814191385523e-08,
                4.678628457585224964e-08,
                3.554190428115991185e-08,
                2.701087761963807523e-08,
                2.053557864940365728e-08,
                1.561852748669927227e-08,
                1.188318920018147882e-08,
                9.044422619773291474e-09,
            ];
            ddouble[] expected_dist_mu2lambda1 = [
                0.000000000000000000e+00,
                6.902185506120922520e-01,
                6.429310691952073720e-01,
                4.733917111740103545e-01,
                3.520653267642995243e-01,
                2.698459980844180039e-01,
                2.126793750717355658e-01,
                1.715592654693821884e-01,
                1.410473958869390698e-01,
                1.177953974098791129e-01,
                9.967158665804999595e-02,
                8.527208945504161086e-02,
                7.364318792573194827e-02,
                6.411883108237549644e-02,
                5.622239416274822654e-02,
                4.960539770317289737e-02,
                4.400816584553744054e-02,
                3.923390170884066830e-02,
                3.513119701275929652e-02,
                3.158197861168046683e-02,
                2.849304086055453114e-02,
                2.578999266895671758e-02,
                2.341286169643149936e-02,
                2.131285563998111185e-02,
                1.944994434463581667e-02,
                1.779103271238377526e-02,
                1.630856456900397350e-02,
                1.497944479587816173e-02,
                1.378419920147079883e-02,
                1.270631387246776503e-02,
                1.173171136329217643e-02,
                1.084833217626189963e-02,
                1.004579795617511519e-02,
                9.315138613801393977e-03,
                8.648569842276843805e-03,
                8.039310638866255937e-03,
                7.481432798298479070e-03,
                6.969736118456758105e-03,
                6.499644407848904556e-03,
                6.067118416986721878e-03,
                5.668582612248956529e-03,
                5.300863329207994475e-03,
                4.961136325659808687e-03,
                4.646882134555253840e-03,
                4.355847917508967604e-03,
                4.086014758503164583e-03,
                3.835569528414331866e-03,
                3.602880604456028344e-03,
                3.386476852525222968e-03,
                3.185029380921655575e-03,
                2.997335655763581734e-03,
                2.822305635384200762e-03,
                2.658949635992997277e-03,
                2.506367686237443403e-03,
                2.363740165835533265e-03,
                2.230319554627222759e-03,
                2.105423144379286234e-03,
                1.988426587409891297e-03,
                1.878758174331154451e-03,
                1.775893748550979407e-03,
                1.679352178125939838e-03,
                1.588691316519330859e-03,
                1.503504393123794060e-03,
                1.423416782327881778e-03,
            ];
            ddouble[] expected_dist_mu2lambda2 = [
                0.000000000000000000e+00,
                2.110999837206014862e-01,
                5.180703826635669751e-01,
                5.159862466781215407e-01,
                4.393912894677223790e-01,
                3.607463244928694190e-01,
                2.945727517029277931e-01,
                2.415407226594145085e-01,
                1.994711402007163514e-01,
                1.660104216298307067e-01,
                1.392059147818068832e-01,
                1.175486869641090421e-01,
                9.989689156669755155e-02,
                8.538885031107104617e-02,
                7.337120803776485400e-02,
                6.334461564671768630e-02,
                5.492391118346529738e-02,
                4.780918347757081438e-02,
                4.176472228300778800e-02,
                3.660361819097512903e-02,
                3.217640652624995795e-02,
                2.836261857166412839e-02,
                2.506443767719740076e-02,
                2.220189199763087165e-02,
                1.970917898490986975e-02,
                1.753183095513924372e-02,
                1.562451127696259413e-02,
                1.394928748449480310e-02,
                1.247426808801886666e-02,
                1.117251895359625685e-02,
                1.002119622146375490e-02,
                9.000848160442194781e-03,
                8.094849729118233986e-03,
                7.288942066714144269e-03,
                6.570855466358327401e-03,
                5.929999159172338538e-03,
                5.357204866542513919e-03,
                4.844513854417346714e-03,
                4.384999361670382546e-03,
                3.972617931570225809e-03,
                3.602084467215365977e-03,
                3.268766844696518667e-03,
                2.968596716331230005e-03,
                2.697993769239011557e-03,
                2.453801208682905403e-03,
                2.233230639130947930e-03,
                2.033814840468212929e-03,
                1.853367198859700053e-03,
                1.689946764329525132e-03,
                1.541828080237454780e-03,
                1.407475071373670257e-03,
                1.285518393551454914e-03,
                1.174735743251513823e-03,
                1.074034704947402288e-03,
                9.824377793173047664e-04,
                8.990692901030755213e-04,
                8.231439129057021263e-04,
                7.539566073144778578e-04,
                6.908737657563518257e-04,
                6.333254193792266706e-04,
                5.807983640100949528e-04,
                5.328300884606009110e-04,
                4.890034037670765296e-04,
                4.489416858252333656e-04,
            ];
            ddouble[] expected_dist_mu3lambda4 = [
                0.000000000000000000e+00,
                7.684330425933131280e-03,
                1.403173887779394224e-01,
                2.740983640955898526e-01,
                3.280201493519872558e-01,
                3.312270708806419051e-01,
                3.111991095141730668e-01,
                2.826268696902350386e-01,
                2.524295107478290445e-01,
                2.236345052995515337e-01,
                1.974145053476377420e-01,
                1.740794934456247611e-01,
                1.535529553205935704e-01,
                1.355998036574632415e-01,
                1.199345836296211304e-01,
                1.062714673486700040e-01,
                9.434580692324828766e-02,
                8.392185813018626650e-02,
                7.479392911046786607e-02,
                6.678457531802944880e-02,
                5.974164090831513008e-02,
                5.353502417299397381e-02,
                4.805357505482682479e-02,
                4.320229473575631224e-02,
                3.889988868927163335e-02,
                3.507666782681964635e-02,
                3.167276749716796402e-02,
                2.863664542015600251e-02,
                2.592381890211338177e-02,
                2.349580457043429910e-02,
                2.131922814951821982e-02,
                1.936507636224974946e-02,
                1.760806735335330303e-02,
                1.602611987552777367e-02,
                1.459990479628144218e-02,
                1.331246528944465997e-02,
                1.214889442044397340e-02,
                1.109606077961539441e-02,
                1.014237442407047753e-02,
                9.277586711711211659e-03,
                8.492618699693962880e-03,
                7.779413675443496053e-03,
                7.130810125888894374e-03,
                6.540432058602263336e-03,
                6.002594090508919789e-03,
                5.512219135051125291e-03,
                5.064766862820432348e-03,
                4.656171396510626163e-03,
                4.282786938993591447e-03,
                3.941340231863900249e-03,
                3.628888907847955484e-03,
                3.342784939675547129e-03,
                3.080642504980089617e-03,
                2.840309685310173689e-03,
                2.619843500501745689e-03,
                2.417487850029597339e-03,
                2.231653992631533068e-03,
                2.060903246216370856e-03,
                1.903931633265186149e-03,
                1.759556233807298186e-03,
                1.626703039590190027e-03,
                1.504396130092544668e-03,
                1.391748014242079464e-03,
                1.287951001672949826e-03,
            ];

            foreach ((InverseGaussDistribution dist, ddouble[] expecteds) in new[]{
                (dist_mu1lambda1, expected_dist_mu1lambda1),
                (dist_mu2lambda1, expected_dist_mu2lambda1),
                (dist_mu1lambda2, expected_dist_mu1lambda2),
                (dist_mu2lambda2, expected_dist_mu2lambda2),
                (dist_mu3lambda4, expected_dist_mu3lambda4),
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
            ddouble[] expected_dist_mu1lambda1 = [
                0.000000000000000000e+00,
                1.126907667166023863e-01,
                3.649755481729599071e-01,
                5.464181447269447212e-01,
                6.681020012231705385e-01,
                7.516606358677878630e-01,
                8.107679929999789259e-01,
                8.536695284085747959e-01,
                8.854754259860063703e-01,
                9.094690723566644097e-01,
                9.278319592945425409e-01,
                9.420561866350805058e-01,
                9.531879207427882417e-01,
                9.619763181346883130e-01,
                9.689676318770072960e-01,
                9.745664175082131564e-01,
                9.790763641788862159e-01,
                9.827281489880055343e-01,
                9.856988170539068284e-01,
                9.881255029119014877e-01,
                9.901152973996735218e-01,
                9.917524417612406662e-01,
                9.931036384906714032e-01,
                9.942220159302010973e-01,
                9.951501178662978120e-01,
                9.959221785838520136e-01,
                9.965658685895064650e-01,
                9.971036443357288448e-01,
                9.975537990175050718e-01,
                9.979312858500274386e-01,
                9.982483668600516946e-01,
                9.985151269243504046e-01,
                9.987398830675029426e-01,
                9.989295118597472101e-01,
                9.990897124196070900e-01,
                9.992252185239535178e-01,
                9.993399703042274984e-01,
                9.994372537067789386e-01,
                9.995198141333510033e-01,
                9.995899493202395902e-01,
                9.996495854627911459e-01,
                9.997003397724986140e-01,
                9.997435720124210334e-01,
                9.997804270520169734e-01,
                9.998118700837438366e-01,
                9.998387158274147923e-01,
                9.998616527962892464e-01,
                9.998812634973759383e-01,
                9.998980412767494164e-01,
                9.999124043905218695e-01,
                9.999247077769934933e-01,
                9.999352529203550111e-01,
                9.999442961271486352e-01,
                9.999520554803543737e-01,
                9.999587166899617241e-01,
                9.999644380212262540e-01,
                9.999693544509085097e-01,
                9.999735811763859061e-01,
                9.999772165815901115e-01,
                9.999803447464313821e-01,
                9.999830375720731857e-01,
                9.999853565825690405e-01,
                9.999873544535290426e-01,
                9.999890763103047586e-01,
            ];
            ddouble[] expected_dist_mu1lambda2 = [
                0.000000000000000000e+00,
                2.805684041471992463e-02,
                2.323571891918430332e-01,
                4.580233401504606450e-01,
                6.276978381552530406e-01,
                7.449252444352248981e-01,
                8.244079562051370713e-01,
                8.782964106873687449e-01,
                9.150466813289288570e-01,
                9.402986873915673005e-01,
                9.577838788517624158e-01,
                9.699797202393830586e-01,
                9.785435738738854639e-01,
                9.845939019468854525e-01,
                9.888921337919750476e-01,
                9.919609659008418179e-01,
                9.941619868900306578e-01,
                9.957471051447951416e-01,
                9.968929489647028186e-01,
                9.977240887556957372e-01,
                9.983288481644420065e-01,
                9.987701547267583413e-01,
                9.990930405480935583e-01,
                9.993298618494389496e-01,
                9.995039532898249979e-01,
                9.996322007290925349e-01,
                9.997268618053186673e-01,
                9.997968602880534794e-01,
                9.998487102798229165e-01,
                9.998871787424528090e-01,
                9.999157621797879480e-01,
                9.999370307703006056e-01,
                9.999528775879920062e-01,
                9.999646996145654665e-01,
                9.999735295673766178e-01,
                9.999801321500559892e-01,
                9.999850744931266799e-01,
                9.999887778194359589e-01,
                9.999915554170647480e-01,
                9.999936406026810332e-01,
                9.999952073513723727e-01,
                9.999963855423228409e-01,
                9.999972722437623363e-01,
                9.999979400789152928e-01,
                9.999984434369370012e-01,
                9.999988230902518982e-01,
                9.999991096316239725e-01,
                9.999993260357957769e-01,
                9.999994895708995690e-01,
                9.999996132262704007e-01,
                9.999997067801430495e-01,
                9.999997775988611615e-01,
                9.999998312356886965e-01,
                9.999998718798849273e-01,
                9.999999026937780311e-01,
                9.999999260659854228e-01,
                9.999999438017931341e-01,
                9.999999572663964553e-01,
                9.999999674927583238e-01,
                9.999999752628775385e-01,
                9.999999811690670759e-01,
                9.999999856601878845e-01,
                9.999999890765552957e-01,
                9.999999916763068208e-01,
            ];
            ddouble[] expected_dist_mu2lambda1 = [
                0.000000000000000000e+00,
                7.328871874451958757e-02,
                2.492117733417387737e-01,
                3.879443205731936128e-01,
                4.901383399453297929e-01,
                5.672289594870942864e-01,
                6.271337295439404214e-01,
                6.748961601790895992e-01,
                7.137917880779035551e-01,
                7.460215894948594162e-01,
                7.731148799422212781e-01,
                7.961665277826092124e-01,
                8.159810287041544541e-01,
                8.331629859307949637e-01,
                8.481757552279055590e-01,
                8.613805606424956851e-01,
                8.730632624933560848e-01,
                8.834530998720541151e-01,
                8.927360851898458183e-01,
                9.010647539384023519e-01,
                9.085653794893191471e-01,
                9.153433921545227614e-01,
                9.214875048026222037e-01,
                9.270728925983647839e-01,
                9.321636713955138154e-01,
                9.368148494603656395e-01,
                9.410738790614736615e-01,
                9.449819007761690592e-01,
                9.485747494626777243e-01,
                9.518837736654777171e-01,
                9.549365077218867626e-01,
                9.577572266417246905e-01,
                9.603674069952539138e-01,
                9.627861119122266276e-01,
                9.650303144065867089e-01,
                9.671151702702615038e-01,
                9.690542494908528370e-01,
                9.708597333716650901e-01,
                9.725425831438917479e-01,
                9.741126847678741596e-01,
                9.755789737546814999e-01,
                9.769495431494203963e-01,
                9.782317372647832965e-01,
                9.794322333077748466e-01,
                9.805571126816183813e-01,
                9.816119234510042313e-01,
                9.826017352185029141e-01,
                9.835311874624994566e-01,
                9.844045322240596585e-01,
                9.852256718951269976e-01,
                9.859981927481494823e-01,
                9.867253947534588843e-01,
                9.874103181521456385e-01,
                9.880557671860935587e-01,
                9.886643313310832237e-01,
                9.892384043316867492e-01,
                9.897802012966043783e-01,
                9.902917740789753953e-01,
                9.907750251370596173e-01,
                9.912317200457309019e-01,
                9.916634988077969171e-01,
                9.920718860957202789e-01,
                9.924583005383927903e-01,
                9.928240631538514105e-01,
            ];
            ddouble[] expected_dist_mu2lambda2 = [
                0.000000000000000000e+00,
                1.206821184832046825e-02,
                1.126907667166023863e-01,
                2.451369679139390745e-01,
                3.649755481729599071e-01,
                4.647908770819937363e-01,
                5.464181447269447212e-01,
                6.131782423371359148e-01,
                6.681020012231705385e-01,
                7.136299705785378311e-01,
                7.516606358677878630e-01,
                7.836608480386715314e-01,
                8.107679929999789259e-01,
                8.338707006492269702e-01,
                8.536695284085747959e-01,
                8.707219118653545253e-01,
                8.854754259860063703e-01,
                8.982925223538784643e-01,
                9.094690723566644097e-01,
                9.192483935189047495e-01,
                9.278319592945425409e-01,
                9.353876530021798530e-01,
                9.420561866350805058e-01,
                9.479561356903500879e-01,
                9.531879207427882417e-01,
                9.578369803897189838e-01,
                9.619763181346883130e-01,
                9.656685606629201724e-01,
                9.689676318770072960e-01,
                9.719201225888594031e-01,
                9.745664175082131564e-01,
                9.769416274379841836e-01,
                9.790763641788862159e-01,
                9.809973876964377393e-01,
                9.827281489880055343e-01,
                9.842892473502977690e-01,
                9.856988170539068284e-01,
                9.869728555331794828e-01,
                9.881255029119014877e-01,
                9.891692808691076033e-01,
                9.901152973996735218e-01,
                9.909734228612473439e-01,
                9.917524417612406662e-01,
                9.924601839778377510e-01,
                9.931036384906714032e-01,
                9.936890521914656382e-01,
                9.942220159302010973e-01,
                9.947075396106255241e-01,
                9.951501178662978120e-01,
                9.955537876137404707e-01,
                9.959221785838520136e-01,
                9.962585577694037520e-01,
                9.965658685895064650e-01,
                9.968467654567696457e-01,
                9.971036443357288448e-01,
                9.973386697989481409e-01,
                9.975537990175050718e-01,
                9.977508030632918912e-01,
                9.979312858500274386e-01,
                9.980967009966866810e-01,
                9.982483668600516946e-01,
                9.983874799513172693e-01,
                9.985151269243504046e-01,
                9.986322952996269375e-01,
            ];
            ddouble[] expected_dist_mu3lambda4 = [
                0.000000000000000000e+00,
                2.285523296768642617e-04,
                1.617263726906272014e-02,
                6.964192865095256146e-02,
                1.463377354671119568e-01,
                2.294579648482432677e-01,
                3.100525901763213077e-01,
                3.843662809374400746e-01,
                4.512407860318277120e-01,
                5.107027262515501098e-01,
                5.632748612063291915e-01,
                6.096513926215343382e-01,
                6.505490513900129557e-01,
                6.866424412737403005e-01,
                7.185396065188452974e-01,
                7.467765120019125291e-01,
                7.718200458880227766e-01,
                7.940744737290883482e-01,
                8.138888977978904204e-01,
                8.315645825709803152e-01,
                8.473616591988195612e-01,
                8.615050436266348255e-01,
                8.741895571750475380e-01,
                8.855843072694069562e-01,
                8.958364119337960396e-01,
                9.050741565148509027e-01,
                9.134096661862046718e-01,
                9.209411689079367580e-01,
                9.277549135840780226e-01,
                9.339267985626997781e-01,
                9.395237569536392108e-01,
                9.446049376910270468e-01,
                9.492227148361314093e-01,
                9.534235522064576918e-01,
                9.572487459029183077e-01,
                9.607350635572972086e-01,
                9.639152960153762884e-01,
                9.668187345991674508e-01,
                9.694715849625229476e-01,
                9.718973267902469848e-01,
                9.741170271272731052e-01,
                9.761496139081219381e-01,
                9.780121152440802712e-01,
                9.797198691803783444e-01,
                9.812867079288165595e-01,
                9.827251199887353517e-01,
                9.840463930713018170e-01,
                9.852607403226666394e-01,
                9.863774119874021418e-01,
                9.874047943538998284e-01,
                9.883504975691045358e-01,
                9.892214336937275032e-01,
                9.900238861847650007e-01,
                9.907635718347104214e-01,
                9.914456960620519066e-01,
                9.920750023320012989e-01,
                9.926558163869667606e-01,
                9.931920858806160268e-01,
                9.936874159354209812e-01,
                9.941451010796097387e-01,
                9.945681539640132129e-01,
                9.949593312111669752e-01,
                9.953211567071619248e-01,
                9.956559426102649102e-01,
            ];

            foreach ((InverseGaussDistribution dist, ddouble[] expecteds) in new[]{
                (dist_mu1lambda1, expected_dist_mu1lambda1),
                (dist_mu2lambda1, expected_dist_mu2lambda1),
                (dist_mu1lambda2, expected_dist_mu1lambda2),
                (dist_mu2lambda2, expected_dist_mu2lambda2),
                (dist_mu3lambda4, expected_dist_mu3lambda4),
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