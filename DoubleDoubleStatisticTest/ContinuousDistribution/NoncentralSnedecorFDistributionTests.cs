using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class NoncentralSnedecorFDistributionTests {
        readonly NoncentralSnedecorFDistribution dist_n_1_m_1_lambda_1 = new(n: 1, m: 1, lambda: 1);
        readonly NoncentralSnedecorFDistribution dist_n_2_m_1_lambda_2 = new(n: 2, m: 1, lambda: 2);
        readonly NoncentralSnedecorFDistribution dist_n_1_m_2_lambda_3 = new(n: 1, m: 2, lambda: 3);
        readonly NoncentralSnedecorFDistribution dist_n_2_m_2_lambda_4 = new(n: 2, m: 2, lambda: 4);
        readonly NoncentralSnedecorFDistribution dist_n_3_m_4_lambda_1 = new(n: 3, m: 4, lambda: 1);
        readonly NoncentralSnedecorFDistribution dist_n_4_m_2_lambda_2 = new(n: 4, m: 2, lambda: 2);
        readonly NoncentralSnedecorFDistribution dist_n_2_m_4_lambda_3 = new(n: 2, m: 4, lambda: 3);
        readonly NoncentralSnedecorFDistribution dist_n_4_m_4_lambda_4 = new(n: 4, m: 4, lambda: 4);
        readonly NoncentralSnedecorFDistribution dist_n_6_m_8_lambda_1 = new(n: 6, m: 8, lambda: 1);
        readonly NoncentralSnedecorFDistribution dist_n_8_m_6_lambda_2 = new(n: 8, m: 6, lambda: 2);
        readonly NoncentralSnedecorFDistribution dist_n_8_m_8_lambda_3 = new(n: 8, m: 8, lambda: 3);
        readonly NoncentralSnedecorFDistribution dist_n_10_m_10_lambda_4 = new(n: 10, m: 10, lambda: 4);

        NoncentralSnedecorFDistribution[] Dists => [
            dist_n_1_m_1_lambda_1,
            dist_n_2_m_1_lambda_2,
            dist_n_1_m_2_lambda_3,
            dist_n_2_m_2_lambda_4,
            dist_n_3_m_4_lambda_1,
            dist_n_4_m_2_lambda_2,
            dist_n_2_m_4_lambda_3,
            dist_n_4_m_4_lambda_4,
            dist_n_6_m_8_lambda_1,
            dist_n_8_m_6_lambda_2,
            dist_n_8_m_8_lambda_3,
            dist_n_10_m_10_lambda_4
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"D1={dist.N}");
                Console.WriteLine($"D2={dist.M}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Median={dist.Median}");
                //Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Variance={dist.Variance}");
                Console.WriteLine($"Skewness={dist.Skewness}");
                Console.WriteLine($"Kurtosis={dist.Kurtosis}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            Assert.Inconclusive();

            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= 1; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-10"; p /= 10) {
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-25);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 1000; i++) {
                    ddouble p = (ddouble)i / 1000;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-28);
                    }
                }

                for (ddouble p = (ddouble)1 / 1000; p >= "1e-10"; p /= 10) {
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"quantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-25);
                    }
                }
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (NoncentralSnedecorFDistribution dist in Dists) {
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
            ddouble[] expected_dist_n_3_m_4_lambda_1 = [
                0,
                3.313886646841897599e-01,
                4.213709925081267627e-01,
                4.650038330323429947e-01,
                4.848783841413655815e-01,
                4.906306425292091822e-01,
                4.874927551842822004e-01,
                4.786374273126058898e-01,
                4.661126012710236632e-01,
                4.512901330119833410e-01,
                4.351099487300002511e-01,
                4.182244054538874489e-01,
                4.010887903409254274e-01,
                3.840205230788706348e-01,
                3.672469496403262812e-01,
                3.509029409952401646e-01,
                3.350936887115685270e-01,
                3.198809706367322647e-01,
                3.052998781963300878e-01,
                2.913661844732018391e-01,
                2.780817606407026688e-01,
                2.654385776637568917e-01,
                2.534216724518756347e-01,
                2.420113503057152726e-01,
                2.311848209229916129e-01,
                2.209174100396837648e-01,
                2.111834547721680166e-01,
                2.019569597305758968e-01,
                1.932120735251285637e-01,
                1.849234296635859642e-01,
                1.770663849487430219e-01,
                1.696171807410995491e-01,
                1.625530464149926857e-01,
                1.558522591586786188e-01,
                1.494941718555020671e-01,
                1.434592168496340303e-01,
                1.377288926372732814e-01,
                1.322857381122477527e-01,
                1.271132984086154050e-01,
                1.221960847667957623e-01,
                1.175195311384213820e-01,
                1.130699485554775702e-01,
                1.088344792260498156e-01,
                1.048010503143514821e-01,
                1.009583290426108032e-01,
                9.729567862770681363e-02,
                9.380311631335123967e-02,
                9.047127274564024935e-02,
                8.729135372309569652e-02,
                8.425510366338873069e-02,
                8.135477115034372986e-02,
                7.858307662089725198e-02,
                7.593318172571805746e-02,
                7.339866090194657744e-02,
                7.097347439893875176e-02,
                6.865194339678026836e-02,
                6.642872646400732606e-02,
                6.429879767372304045e-02,
                6.225742631149610434e-02,
                6.030015765878803807e-02,
                5.842279527656923238e-02,
                5.662138446160813049e-02,
                5.489219688792345425e-02,
                5.323171597959585472e-02,
            ];
            ddouble[] expected_dist_n_2_m_4_lambda_3 = [
                0,
                2.324907367292141924e-01,
                2.398318992480677203e-01,
                2.453977988771046947e-01,
                2.494112193329733607e-01,
                2.520804605863485870e-01,
                2.535754691839961494e-01,
                2.540595454103219963e-01,
                2.536743870781843491e-01,
                2.525452673566783157e-01,
                2.507826427477349851e-01,
                2.484950950230532529e-01,
                2.457493173724928093e-01,
                2.426280146593240150e-01,
                2.391965325365674466e-01,
                2.355118998792482266e-01,
                2.316237979298196581e-01,
                2.275754275690017037e-01,
                2.234042821888659602e-01,
                2.191428342856321654e-01,
                2.148191444631408487e-01,
                2.104573997220582082e-01,
                2.060889463289388990e-01,
                2.017119190275695473e-01,
                1.973507008503488647e-01,
                1.930182024573179056e-01,
                1.887253383077436009e-01,
                1.844812789003613396e-01,
                1.802936722708481099e-01,
                1.761688389938231569e-01,
                1.721119430336814637e-01,
                1.681271420965468888e-01,
                1.642177195240734910e-01,
                1.603862000307265667e-01,
                1.566344513953543682e-01,
                1.529637734615230116e-01,
                1.493749763582796763e-01,
                1.458684491306705411e-01,
                1.424442196688857631e-01,
                1.391020074341386970e-01,
                1.358412695322275265e-01,
                1.326612410243455109e-01,
                1.295609702911537653e-01,
                1.265393497755507024e-01,
                1.235951429770509730e-01,
                1.207270079128796997e-01,
                1.179335176702822480e-01,
                1.152131781020915152e-01,
                1.125644433816463813e-01,
                1.099857294850625639e-01,
                1.074826740860013036e-01,
                1.050393373711044770e-01,
                1.026611418097900102e-01,
                1.003464515800245177e-01,
                9.809363688204063292e-02,
                9.590107906548571426e-02,
                9.376717538950306297e-02,
                9.169034280026576444e-02,
                8.966902124291074860e-02,
                8.770167618452061331e-02,
                8.578680092269362623e-02,
                8.392291830222609050e-02,
                8.210858199186876760e-02,
                8.034237776455288582e-02,
            ];
            ddouble[] expected_dist_n_6_m_8_lambda_1 = [
                0,
                4.584450195367118647e-02,
                1.414836177749366486e-01,
                2.477615859196400372e-01,
                3.456368738049414024e-01,
                4.270774918322880476e-01,
                4.898937993199313223e-01,
                5.348307534177448863e-01,
                5.639543899055993181e-01,
                5.797788718467367097e-01,
                5.848223649719541672e-01,
                5.813635774894881925e-01,
                5.713955793251040971e-01,
                5.565864334457215623e-01,
                5.383034372582334859e-01,
                5.176470201276317518e-01,
                4.954892051094172412e-01,
                4.725111854550489432e-01,
                4.492374495637662668e-01,
                4.260654830751108846e-01,
                4.032909272970452430e-01,
                3.811284959590043875e-01,
                3.597291329086870770e-01,
                3.391939445124281249e-01,
                3.195854238635664757e-01,
                3.009364339404130284e-01,
                2.832573567833340800e-01,
                2.665417530675617819e-01,
                2.507708210162939988e-01,
                2.359168902427155323e-01,
                2.219461451458126511e-01,
                2.088207336034541839e-01,
                1.965003885801008543e-01,
                1.849436633588607826e-01,
                1.741088626106001847e-01,
                1.639547337659030823e-01,
                1.544409710091398935e-01,
                1.455285727230970494e-01,
                1.371800855118343598e-01,
                1.293597607157637164e-01,
                1.220336435098223621e-01,
                1.151696111853506149e-01,
                1.087373735803032065e-01,
                1.027084447297044040e-01,
                9.705609472498986923e-02,
                9.175528701277356480e-02,
                8.678260609851216889e-02,
                8.211617951447980346e-02,
                7.773559651258477032e-02,
                7.362182565251451649e-02,
                6.975985773899351372e-02,
                6.612766593194763232e-02,
                6.271265846347073847e-02,
                5.950052340275568258e-02,
                5.647794220908774337e-02,
                5.363252237738923983e-02,
                5.095273322153203582e-02,
                4.842784591535530953e-02,
                4.604787684492350069e-02,
                4.380353470501140850e-02,
                4.168617128291751683e-02,
                3.968773518714385773e-02,
                3.780072917042298641e-02,
                3.601817003956453433e-02,
            ];
            ddouble[] expected_dist_n_8_m_6_lambda_2 = [
                0,
                1.111858220853592302e-02,
                5.931536218352903567e-02,
                1.361717137984483350e-01,
                2.236295244903906965e-01,
                3.077905712467031707e-01,
                3.807099255324213893e-01,
                4.390489884018716205e-01,
                4.823475479536972754e-01,
                5.117454409050981035e-01,
                5.290981487235043579e-01,
                5.364674329159241362e-01,
                5.358421927670331542e-01,
                5.290047088095384265e-01,
                5.174815595559989490e-01,
                5.025407635680934071e-01,
                4.852119493810313489e-01,
                4.663161443305141041e-01,
                4.464977999704688449e-01,
                4.262692878141682296e-01,
                4.059829946723158556e-01,
                3.859357440898213976e-01,
                3.663354699354192667e-01,
                3.473302496913299175e-01,
                3.290213957361787234e-01,
                3.114739922115139503e-01,
                2.947253131088678479e-01,
                2.787915041956023909e-01,
                2.636728591576265046e-01,
                2.493579638884746152e-01,
                2.358269360874798160e-01,
                2.230539438805023422e-01,
                2.110091536389191169e-01,
                1.996602252467283012e-01,
                1.889734520692243069e-01,
                1.789146205230518205e-01,
                1.694496508025666870e-01,
                1.605450662189866762e-01,
                1.521683297450726613e-01,
                1.442880769614296455e-01,
                1.368742699431324539e-01,
                1.298982904740508815e-01,
                1.233329864291920902e-01,
                1.171526842619341341e-01,
                1.113331754343249574e-01,
                1.058516841442869033e-01,
                1.006868218494938322e-01,
                9.581853271795058580e-02,
                9.122803346778284028e-02,
                8.689774949055538433e-02,
                8.281124978248843860e-02,
                7.895318171885490344e-02,
                7.530920674847507712e-02,
                7.186593798946816491e-02,
                6.861087957360512135e-02,
                6.553236911455151414e-02,
                6.261952260472325449e-02,
                5.986218238052120472e-02,
                5.725086785479670953e-02,
                5.477672936210486315e-02,
                5.243150480727098300e-02,
                5.020747895767562613e-02,
                4.809744561379280386e-02,
                4.609467215976880555e-02,
            ];
            ddouble[] expected_dist_n_8_m_8_lambda_3 = [
                0,
                5.593424595878982927e-03,
                3.301213636650268679e-02,
                8.275933746789386480e-02,
                1.467836813825486308e-01,
                2.161352968876655845e-01,
                2.836968671004503895e-01,
                3.447539238766658132e-01,
                3.966923722054787604e-01,
                4.384827433912690808e-01,
                4.701274497337370462e-01,
                4.923012566276571134e-01,
                5.060209340067289840e-01,
                5.124480462751629384e-01,
                5.127586581193943616e-01,
                5.080659101180939885e-01,
                4.993787096202461639e-01,
                4.875840326181573592e-01,
                4.734439537433520750e-01,
                4.576199745517927631e-01,
                4.406109386956291174e-01,
                4.228704605355181889e-01,
                4.047519222200035016e-01,
                3.865370894473584840e-01,
                3.684473666151765747e-01,
                3.506538316197515548e-01,
                3.332859953672495656e-01,
                3.164393282698019227e-01,
                3.001816493812536901e-01,
                2.845584977720361231e-01,
                2.695976120764553841e-01,
                2.553126417254025782e-01,
                2.417062022939952814e-01,
                2.287723778603145408e-01,
                2.164987598696210547e-01,
                2.048680985994310866e-01,
                1.938596354347787709e-01,
                1.834501699535251162e-01,
                1.736149101483341617e-01,
                1.643281443017285959e-01,
                1.555637677813637154e-01,
                1.472956916409584593e-01,
                1.394981552563234750e-01,
                1.321459616426823303e-01,
                1.252146507924045693e-01,
                1.186806226613490306e-01,
                1.125212209865389251e-01,
                1.067147853556216081e-01,
                1.012406788583608414e-01,
                9.607929618848930209e-02,
                9.121205741258185673e-02,
                8.662138997750366831e-02,
                8.229070293241935730e-02,
                7.820435455468954800e-02,
                7.434761622460150932e-02,
                7.070663325670079047e-02,
                6.726838411297286768e-02,
                6.402063922605893964e-02,
                6.095191939781140905e-02,
                5.805145484316653198e-02,
                5.530914529705554389e-02,
                5.271552115243548914e-02,
                5.026170597499656639e-02,
                4.793938078034853251e-02,
            ];
            ddouble[] expected_dist_n_10_m_10_lambda_4 = [
                0,
                8.957744028843447068e-04,
                9.920214956654974103e-03,
                3.500091702580672426e-02,
                7.770067153860422959e-02,
                1.343591198806818876e-01,
                1.990259305859661543e-01,
                2.656563859484896706e-01,
                3.292907751365620195e-01,
                3.864301583868245404e-01,
                4.349987224873743674e-01,
                4.740414920113497077e-01,
                5.034893777100187817e-01,
                5.238411178817847791e-01,
                5.359439786720848753e-01,
                5.408222524354994176e-01,
                5.395564686674209076e-01,
                5.332036732471484575e-01,
                5.227488134425251065e-01,
                5.091050413214004422e-01,
                4.930003870154408552e-01,
                4.751249138651436055e-01,
                4.560411871773472625e-01,
                4.362158787234171053e-01,
                4.160301939587495723e-01,
                3.957907242724267305e-01,
                3.757399835921381714e-01,
                3.560662015039683381e-01,
                3.369121584800849334e-01,
                3.183829842122454234e-01,
                3.005529287690433815e-01,
                2.834711634952524761e-01,
                2.671666984442599979e-01,
                2.516525121454282754e-01,
                2.369289927528494821e-01,
                2.229867839348420311e-01,
                2.098091231084175667e-01,
                1.973737496274718672e-01,
                1.856544508099799706e-01,
                1.746223072543462151e-01,
                1.642466860657387784e-01,
                1.544960285371210773e-01,
                1.453384660177192700e-01,
                1.367422976558674286e-01,
                1.286867402791602721e-01,
                1.211204908951879350e-01,
                1.140247220643430953e-01,
                1.073712095955858675e-01,
                1.011329590461240574e-01,
                9.528425909816061468e-02,
                8.980070548697537125e-02,
                8.465920233458046340e-02,
                7.983794609173200385e-02,
                7.531639644453580473e-02,
                7.107523722071906747e-02,
                6.709633129237113636e-02,
                6.336267013573593232e-02,
                5.985832109839383380e-02,
                5.656837280265668982e-02,
                5.347888060031057850e-02,
                5.057681195241858080e-02,
                4.784999353690810170e-02,
                4.528705909306873084e-02,
                4.287739979458304179e-02,
            ];

            foreach ((NoncentralSnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_n_3_m_4_lambda_1,   expected_dist_n_3_m_4_lambda_1),
                 (dist_n_2_m_4_lambda_3,   expected_dist_n_2_m_4_lambda_3),
                 (dist_n_6_m_8_lambda_1,   expected_dist_n_6_m_8_lambda_1),
                 (dist_n_8_m_6_lambda_2,   expected_dist_n_8_m_6_lambda_2),
                 (dist_n_8_m_8_lambda_3,   expected_dist_n_8_m_8_lambda_3),
                 (dist_n_10_m_10_lambda_4, expected_dist_n_10_m_10_lambda_4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (ddouble.IsInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-4, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist_n_3_m_4_lambda_1 = [
                0.000000000000000000e+00,
                1.442311879319116072e-02,
                3.829767468787535767e-02,
                6.615898973970207875e-02,
                9.593551496802137091e-02,
                1.264774521838201438e-01,
                1.570804703341167063e-01,
                1.872956800815069600e-01,
                2.168341551828511238e-01,
                2.455121800160625645e-01,
                2.732173315650455181e-01,
                2.998863630142523107e-01,
                3.254902666288236324e-01,
                3.500239318780347109e-01,
                3.735007513958670811e-01,
                3.959403660679001180e-01,
                4.173747850074125409e-01,
                4.378395232645354573e-01,
                4.573730784617968670e-01,
                4.760155050505990859e-01,
                4.938073847667145411e-01,
                5.107890890690212071e-01,
                5.270002576676572348e-01,
                5.424794373848282891e-01,
                5.572638400866062325e-01,
                5.713891889682961578e-01,
                5.848896302157750027e-01,
                5.977976927908235094e-01,
                6.101442833524863829e-01,
                6.219587065203555643e-01,
                6.332687030896926306e-01,
                6.441005006253452958e-01,
                6.544788722394714320e-01,
                6.644272004060706172e-01,
                6.739675434638975338e-01,
                6.831207030686665682e-01,
                6.919062913206666110e-01,
                7.003427966489740308e-01,
                7.084476478041655811e-01,
                7.162372755173710415e-01,
                7.237271715396209615e-01,
                7.309319448932044549e-01,
                7.378653752549494405e-01,
                7.445404634567929225e-01,
                7.509694791368767275e-01,
                7.571640056088571047e-01,
                7.631349820411695850e-01,
                7.688927430541750985e-01,
                7.744470558532574600e-01,
                7.798071550215506242e-01,
                7.849817750981613607e-01,
                7.899791810674124948e-01,
                7.948071968824519384e-01,
                7.994732321430710886e-01,
                8.039843070431532679e-01,
                8.083470756981809835e-01,
                8.125678479578379854e-01,
                8.166526098032022007e-01,
                8.206070424224188864e-01,
                8.244365400531996624e-01,
                8.281462266750662105e-01,
                8.317409716290184996e-01,
                8.352254042372712073e-01,
                8.386039274908966279e-01,
            ];
            ddouble[] expected_dist_n_2_m_4_lambda_3 = [
                0.000000000000000000e+00,
                1.424936564468136380e-02,
                2.901933267629702717e-02,
                4.419140317194911183e-02,
                5.966172994109367883e-02,
                7.534016866943313295e-02,
                9.114757977413118217e-02,
                1.070160537799991346e-01,
                1.228869252990042460e-01,
                1.387098632516091112e-01,
                1.544418949215388415e-01,
                1.700480994895767206e-01,
                1.854953783568628056e-01,
                2.007589472373087214e-01,
                2.158174238845152471e-01,
                2.306532449293314802e-01,
                2.452521769183528821e-01,
                2.596028847021524255e-01,
                2.736965510325476414e-01,
                2.875265417303659032e-01,
                3.010881112957001782e-01,
                3.143781443312623214e-01,
                3.273982244615369908e-01,
                3.401419564031794018e-01,
                3.526125452770699775e-01,
                3.648113937203374402e-01,
                3.767406475053265424e-01,
                3.884030788285185443e-01,
                3.998019843855886712e-01,
                4.109410964414870193e-01,
                4.218245053174225490e-01,
                4.324565919052579255e-01,
                4.428419689865960018e-01,
                4.529854302810081279e-01,
                4.628919062776204463e-01,
                4.725664260185210019e-01,
                4.820140841030262635e-01,
                4.912400122702677563e-01,
                5.002493549953358487e-01,
                5.090472486025470200e-01,
                5.176388034594652243e-01,
                5.260290888680864141e-01,
                5.342231203159621966e-01,
                5.422258487907971602e-01,
                5.500421518978644286e-01,
                5.576768265510657629e-01,
                5.651345830361383937e-01,
                5.724200402688427447e-01,
                5.795377220923687567e-01,
                5.864920544770231281e-01,
                5.932948204489835176e-01,
                5.999357897682024987e-01,
                6.064260947285483594e-01,
                6.127697555182536338e-01,
                6.189706902438206759e-01,
                6.250327154706335486e-01,
                6.309595470750526403e-01,
                6.367548013595744072e-01,
                6.424219963886480533e-01,
                6.479645535080147400e-01,
                6.533857990150974127e-01,
                6.586889659520760443e-01,
                6.638771959968954306e-01,
                6.689535414306485572e-01,
            ];
            ddouble[] expected_dist_n_6_m_8_lambda_1 = [
                0.000000000000000000e+00,
                1.022094267160775791e-03,
                6.743820934077797077e-03,
                1.891578437511768490e-02,
                3.753163945960842629e-02,
                6.177412532720727434e-02,
                9.052638365852673508e-02,
                1.226375047052903206e-01,
                1.570505057100305968e-01,
                1.928547589263958362e-01,
                2.292992994516888372e-01,
                2.657813813629573141e-01,
                3.018343626667774338e-01,
                3.371051815400520346e-01,
                3.713354806070857106e-01,
                4.043438234376707796e-01,
                4.360102112186732271e-01,
                4.662629949262995788e-01,
                4.950680364616455997e-01,
                5.224198666923031142e-01,
                5.483345571241990957e-01,
                5.728440292136580680e-01,
                5.959915507610658025e-01,
                6.178282010956480663e-01,
                6.384101198263167909e-01,
                6.577963847886926230e-01,
                6.760473921672942588e-01,
                6.932236352549587810e-01,
                7.093847980551607479e-01,
                7.245890962914341982e-01,
                7.388928118002007217e-01,
                7.523499771946213777e-01,
                7.650121765120356798e-01,
                7.769284346629129479e-01,
                7.881451741996181592e-01,
                7.987062224833020307e-01,
                8.086528559660759452e-01,
                8.180238712031686621e-01,
                8.268556745124685525e-01,
                8.351823840252956543e-01,
                8.430359393180723027e-01,
                8.504462149566622831e-01,
                8.574411351851762841e-01,
                8.640467876983944917e-01,
                8.702875349910786928e-01,
                8.761861222099016633e-01,
                8.817637807696308005e-01,
                8.870403272546658613e-01,
                8.920342573260146013e-01,
                8.967628345050380201e-01,
                9.012486945226921842e-01,
                9.054940094767189862e-01,
                9.095191766990750049e-01,
                9.133373176048368292e-01,
                9.169606905197739133e-01,
                9.204007506431957086e-01,
                9.236682058977605525e-01,
                9.267730689002749989e-01,
                9.297247052865442729e-01,
                9.325318786189582543e-01,
                9.352027920986434317e-01,
                9.377451272954393513e-01,
                9.401660800992697187e-01,
                9.424723940860902927e-01,
            ];
            ddouble[] expected_dist_n_8_m_6_lambda_2 = [
                0.000000000000000000e+00,
                1.893797413861743665e-04,
                2.197020981229057257e-03,
                8.204960169997053618e-03,
                1.943565513909807294e-02,
                3.608563741067687819e-02,
                5.767207541018615213e-02,
                8.336914955628632606e-02,
                1.122391659412702997e-01,
                1.433725246083777105e-01,
                1.759562804653261869e-01,
                2.093018674484079522e-01,
                2.428483071936321791e-01,
                2.761529010002604090e-01,
                3.088764591055468967e-01,
                3.407670612410094169e-01,
                3.716444491744592682e-01,
                4.013860317377899722e-01,
                4.299148488243671462e-01,
                4.571945707758097166e-01,
                4.832018677523677863e-01,
                5.079474776268259051e-01,
                5.314531890894335087e-01,
                5.537518405912966735e-01,
                5.748840029179815847e-01,
                5.948953977787440062e-01,
                6.138349065931212323e-01,
                6.317530492931653363e-01,
                6.487008352515971632e-01,
                6.647289072943225463e-01,
                6.798869153949097521e-01,
                6.942230694553861436e-01,
                7.077838309662761729e-01,
                7.206137117073098119e-01,
                7.327551543563853853e-01,
                7.442484752270853710e-01,
                7.551318536141340410e-01,
                7.654413556066914648e-01,
                7.752109829068754010e-01,
                7.844727393077485234e-01,
                7.932567091556629313e-01,
                8.015911434381798317e-01,
                8.095025501740050755e-01,
                8.170157865939194908e-01,
                8.241541512379781942e-01,
                8.309394745911408142e-01,
                8.373922072662710603e-01,
                8.435315050432726780e-01,
                8.493753103044237029e-01,
                8.549404295832707668e-01,
                8.602426070792903445e-01,
                8.652965940920739341e-01,
                8.701162144042146052e-01,
                8.747144256971594034e-01,
                8.791033771234902927e-01,
                8.832944631860358831e-01,
                8.872983740916752238e-01,
                8.911251427579218731e-01,
                8.947841886550997881e-01,
                8.982843586675368641e-01,
                9.016339651547549128e-01,
                9.048408213890040130e-01,
                9.079122745392914640e-01,
                9.108552363648194650e-01,
            ];
            ddouble[] expected_dist_n_8_m_8_lambda_3 = [
                0.000000000000000000e+00,
                9.307563546801616552e-05,
                1.170573852200218869e-03,
                4.689930948724297180e-03,
                1.181318380661470872e-02,
                2.314758512111799588e-02,
                3.879116960222520755e-02,
                5.847270782984054921e-02,
                8.169435025228209690e-02,
                1.078479566882942775e-01,
                1.362935476011780667e-01,
                1.664163432990637126e-01,
                1.976549824890467266e-01,
                2.295170258544419195e-01,
                2.615836078099524342e-01,
                2.935076957483480076e-01,
                3.250087146734440546e-01,
                3.558653827838441086e-01,
                3.859079393740711406e-01,
                4.150158300342952633e-01,
                4.430904081998912059e-01,
                4.700770171840058520e-01,
                4.959413891596590740e-01,
                5.206690431386526985e-01,
                5.442611666386658609e-01,
                5.667311633089406886e-01,
                5.881017878606443761e-01,
                6.084027895395731633e-01,
                6.276689898818590407e-01,
                6.459387273301130383e-01,
                6.632526090401308361e-01,
                6.796525180268978561e-01,
                6.951808312023995962e-01,
                7.098798106010691944e-01,
                7.237911360721870357e-01,
                7.369555529304284702e-01,
                7.494126125340467848e-01,
                7.612004875684202920e-01,
                7.723558470252093189e-01,
                7.829137785600099875e-01,
                7.929077481566838559e-01,
                8.023695888907439455e-01,
                8.113295121264800791e-01,
                8.198161357544648808e-01,
                8.278565251222474997e-01,
                8.354762431692946967e-01,
                8.426994069797284492e-01,
                8.495487485400129435e-01,
                8.560456779560885465e-01,
                8.622103477642707459e-01,
                8.680617172781209812e-01,
                8.736176161624082948e-01,
                8.788948066258288616e-01,
                8.839090437851843118e-01,
                8.886751338824603330e-01,
                8.932069901386691102e-01,
                8.975176861093130398e-01,
                9.016195064699163986e-01,
                9.055239952095527345e-01,
                9.092420012483627412e-01,
                9.127837215239148483e-01,
                9.161587416127348060e-01,
                9.193760739688543460e-01,
                9.224441938720476797e-01,
            ];
            ddouble[] expected_dist_n_10_m_10_lambda_4 = [
                0.000000000000000000e+00,
                1.192939885636593661e-05,
                2.818148653197692578e-04,
                1.591563952548513423e-03,
                5.027703135932504522e-03,
                1.159646698862771066e-02,
                2.198958971079717148e-02,
                3.651492522908599686e-02,
                5.513316904737839558e-02,
                7.753996684306500853e-02,
                1.032586796091825021e-01,
                1.317166323199703837e-01,
                1.623135768785827771e-01,
                1.944630642634488926e-01,
                2.276217462600976238e-01,
                2.613055137107053172e-01,
                2.950965454810965949e-01,
                3.286441447802977422e-01,
                3.616616145037225283e-01,
                3.939277185758850996e-01,
                4.252543704615989184e-01,
                4.555159478780765880e-01,
                4.846198794655687814e-01,
                5.125057060349013094e-01,
                5.391394082506553165e-01,
                5.645084041939915087e-01,
                5.886172182557033450e-01,
                6.114837869226339429e-01,
                6.331363474028894878e-01,
                6.536108461943811454e-01,
                6.729488027469421318e-01,
                6.911955656047025220e-01,
                7.083989030152026967e-01,
                7.246078757465599374e-01,
                7.398719459888111194e-01,
                7.542402822535971829e-01,
                7.677612258520769561e-01,
                7.804818896820590046e-01,
                7.924478646351935662e-01,
                8.037030129403776035e-01,
                8.142893312182994947e-01,
                8.242468689792731285e-01,
                8.336136908045982397e-01,
                8.424258725648193780e-01,
                8.507327890280742544e-01,
                8.585367391488236688e-01,
                8.658826512735237069e-01,
                8.727990422228226741e-01,
                8.793127028687612867e-01,
                8.854487768323215935e-01,
                8.912308415615541879e-01,
                8.966809901698352681e-01,
                9.018199127860404563e-01,
                9.066669764698782608e-01,
                9.112403029887057304e-01,
                9.155568439472916431e-01,
                9.196324529177993812e-01,
                9.234819543408973042e-01,
                9.271192090663119556e-01,
                9.305571764771813914e-01,
                9.338079732012851686e-01,
                9.368829284569499816e-01,
                9.397926361148378760e-01,
                9.425470035812155700e-01,
            ];

            foreach ((NoncentralSnedecorFDistribution dist, ddouble[] expecteds) in new[]{
                 (dist_n_3_m_4_lambda_1,   expected_dist_n_3_m_4_lambda_1),
                 (dist_n_2_m_4_lambda_3,   expected_dist_n_2_m_4_lambda_3),
                 (dist_n_6_m_8_lambda_1,   expected_dist_n_6_m_8_lambda_1),
                 (dist_n_8_m_6_lambda_2,   expected_dist_n_8_m_6_lambda_2),
                 (dist_n_8_m_8_lambda_3,   expected_dist_n_8_m_8_lambda_3),
                 (dist_n_10_m_10_lambda_4, expected_dist_n_10_m_10_lambda_4),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 16, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-4, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}