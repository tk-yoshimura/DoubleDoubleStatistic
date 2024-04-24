using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ScalableDistribution {
    [TestClass()]
    public class HalfNormalDistributionTests {
        readonly HalfNormalDistribution dist1 = new();
        readonly HalfNormalDistribution dist2 = new(sigma: 3);

        HalfNormalDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (HalfNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
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
            foreach (HalfNormalDistribution dist in Dists) {
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
            foreach (HalfNormalDistribution dist in Dists) {
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
            foreach (HalfNormalDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (HalfNormalDistribution dist in Dists) {
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
            foreach (HalfNormalDistribution dist in Dists) {
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
            foreach (HalfNormalDistribution dist in Dists) {
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
            foreach (HalfNormalDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (HalfNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (HalfNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 8; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (HalfNormalDistribution dist in Dists) {
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
            foreach (HalfNormalDistribution dist in Dists) {
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
            foreach (HalfNormalDistribution dist in Dists) {
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
                7.978845608028654057e-01,
                7.733362336056984754e-01,
                7.041306535285990487e-01,
                6.022748643096088594e-01,
                4.839414490382867307e-01,
                3.652981707780438292e-01,
                2.590351913317834875e-01,
                1.725546376530230630e-01,
                1.079819330263761257e-01,
                6.347930367133483576e-02,
                3.505660098713708067e-02,
                1.818712500318210579e-02,
                8.863696823876015055e-03,
                4.058096114599536210e-03,
                1.745365390091520308e-03,
                7.051913647348908134e-04,
                2.676604515297707353e-04,
                9.543727308240990354e-05,
                3.196748221381094922e-05,
                1.005901457718489078e-05,
                2.973439029468595356e-06,
                8.256941977259996767e-07,
                2.153952008508655215e-07,
                5.278486407141147058e-08,
                1.215176569964657217e-08,
                2.628003636311768021e-09,
                5.339113229525703836e-10,
                1.018987591768736839e-10,
                1.826944081672918728e-11,
                3.077075901122550034e-12,
                4.868641066058019208e-13,
                7.236588902225034600e-14,
                1.010454216707378546e-14,
                1.325427491193750452e-15,
                1.633247126333910039e-16,
                1.890620776380570607e-17,
                2.055954714333783305e-18,
                2.100289965994073927e-19,
                2.015587078860001973e-20,
                1.817106862395333075e-21,
                1.538919725341283988e-22,
                1.224356973049424402e-23,
                9.150751181041611095e-25,
                6.424835743242879729e-26,
                4.237638507018707516e-27,
                2.625683521230806162e-28,
                1.528331082317440644e-29,
                8.356975089085292674e-31,
                4.292767471326120941e-32,
                2.071487019241670135e-33,
                9.390390715950292712e-35,
                3.998910684437090640e-36,
                1.599765551401362534e-37,
                6.012119052199490967e-39,
                2.122537627830432506e-40,
                7.039467675990735119e-42,
                2.193213118777942537e-43,
                6.419163620839914003e-45,
                1.764950994918964953e-46,
                4.558725601310310314e-48,
                1.106141909968883277e-49,
                2.521359877877516829e-51,
                5.399026049177175356e-53,
                1.086056959502566845e-54,
            ];
            ddouble[] expected_dist2 = [
                2.659615202676218204e-01,
                2.650396441722279728e-01,
                2.622931440679599491e-01,
                2.577787445352328066e-01,
                2.515888184619954338e-01,
                2.438482442760016022e-01,
                2.347102178428663588e-01,
                2.243512151633021146e-01,
                2.129653370149015013e-01,
                2.007582881032029531e-01,
                1.879412502735350321e-01,
                1.747249020522009533e-01,
                1.613138163460955676e-01,
                1.479014364704906714e-01,
                1.346657903693725677e-01,
                1.217660569260146097e-01,
                1.093400497839957519e-01,
                9.750263618635714169e-02,
                8.634506377726115789e-02,
                7.593512941346125533e-02,
                6.631809252849911462e-02,
                5.751821255100768998e-02,
                4.954077570599540320e-02,
                4.237447101130822208e-02,
                3.599397767545870624e-02,
                3.036263524860371640e-02,
                2.543508233761198820e-02,
                2.115976789044494641e-02,
                1.748125939580632074e-02,
                1.434229334451577083e-02,
                1.168553366237902631e-02,
                9.455022502958199820e-03,
                7.597324015864961172e-03,
                6.062375001060701929e-03,
                4.804066509739481593e-03,
                3.780587479630594105e-03,
                2.954565607958671829e-03,
                2.293042218121663870e-03,
                1.767317302953403486e-03,
                1.352698704866511998e-03,
                1.028185997527403380e-03,
                7.761155091442563424e-04,
                5.817884633638401026e-04,
                4.330997312797091998e-04,
                3.201804344138804519e-04,
                2.350637882449636135e-04,
                1.713802367082042281e-04,
                1.240851505304960969e-04,
                8.922015050992357843e-05,
                6.370744050009290941e-05,
                4.517533992587842738e-05,
                3.181242436080329892e-05,
                2.224724159709228277e-05,
                1.545039693236913760e-05,
                1.065582740460365031e-05,
                7.298250737423571095e-06,
                4.964030580419993453e-06,
                3.353004859061630400e-06,
                2.249147737431210964e-06,
                1.498255394818001837e-06,
                9.911463431561983813e-07,
                6.511391104739706200e-07,
                4.248091347780569753e-07,
                2.752313992419999099e-07,

            ];

            foreach ((HalfNormalDistribution dist, ddouble[] expecteds) in new[]{
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
                1.974126513658474025e-01,
                3.829249225480262364e-01,
                5.467452952462634741e-01,
                6.826894921370858516e-01,
                7.887004526662892978e-01,
                8.663855974622838296e-01,
                9.198816862723657728e-01,
                9.544997361036415828e-01,
                9.755510546899106572e-01,
                9.875806693484476817e-01,
                9.940404735298908889e-01,
                9.973002039367397931e-01,
                9.988459499152184673e-01,
                9.995347418419289198e-01,
                9.998231654295983706e-01,
                9.999366575163337600e-01,
                9.999786229484501909e-01,
                9.999932046537505226e-01,
                9.999979658335149679e-01,
                9.999994266968561529e-01,
                9.999998479007896623e-01,
                9.999999620208750439e-01,
                9.999999910756551813e-01,
                9.999999980268245992e-01,
                9.999999995895472171e-01,
                9.999999999196800271e-01,
                9.999999999852153820e-01,
                9.999999999974402698e-01,
                9.999999999995832223e-01,
                9.999999999999362732e-01,
                9.999999999999908962e-01,
                9.999999999999986677e-01,
                9.999999999999997780e-01,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_dist2 = [
                0.000000000000000000e+00,
                6.641350370524468083e-02,
                1.323676652218073180e-01,
                1.974126513658474025e-01,
                2.611173196364726667e-01,
                3.230777609786206206e-01,
                3.829249225480262364e-01,
                4.403310728005884744e-01,
                4.950149249061541745e-01,
                5.467452952462634741e-01,
                5.953432380727137208e-01,
                6.406826616704293098e-01,
                6.826894921370858516e-01,
                7.213395051007558401e-01,
                7.566549908512376632e-01,
                7.887004526662892978e-01,
                8.175775605482642483e-01,
                8.434195929103651679e-01,
                8.663855974622838296e-01,
                8.866544907804740916e-01,
                9.044192954543706087e-01,
                9.198816862723657728e-01,
                9.332469848303655002e-01,
                9.447197074139792594e-01,
                9.544997361036415828e-01,
                9.627791496202273702e-01,
                9.697397199795283473e-01,
                9.755510546899106572e-01,
                9.803693427427093710e-01,
                9.843366463571023850e-01,
                9.875806693484476817e-01,
                9.902149267955300260e-01,
                9.923392388648204498e-01,
                9.940404735298908889e-01,
                9.953934677366083061e-01,
                9.964620635217777966e-01,
                9.973002039367397931e-01,
                9.979530426714149360e-01,
                9.984580304310599619e-01,
                9.988459499152184673e-01,
                9.991418793336064308e-01,
                9.993660715301604380e-01,
                9.995347418419289198e-01,
                9.996607627506988170e-01,
                9.997542672200696590e-01,
                9.998231654295983706e-01,
                9.998735815362631563e-01,
                9.999102182470136491e-01,
                9.999366575163337600e-01,
                9.999556057278193855e-01,
                9.999690914062353553e-01,
                9.999786229484501909e-01,
                9.999853131523261940e-01,
                9.999899765379487437e-01,
                9.999932046537505226e-01,
                9.999954237834260962e-01,
                9.999969387465268422e-01,
                9.999979658335149679e-01,
                9.999986573430883485e-01,
                9.999991196974553631e-01,
                9.999994266968561529e-01,
                9.999996291319894759e-01,
                9.999997616942932765e-01,
                9.999998479007896623e-01,
            ];

            foreach ((HalfNormalDistribution dist, ddouble[] expecteds) in new[]{
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