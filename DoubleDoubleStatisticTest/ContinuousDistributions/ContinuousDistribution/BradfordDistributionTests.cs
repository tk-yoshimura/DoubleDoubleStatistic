using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class BradfordDistributionTests {
        readonly BradfordDistribution dist_c1 = new(c: 1);
        readonly BradfordDistribution dist_c2 = new(c: 2);

        BradfordDistribution[] Dists => [
            dist_c1,
            dist_c2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"C={dist.C}");
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
            foreach (BradfordDistribution dist in Dists) {
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
            foreach (BradfordDistribution dist in Dists) {
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
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Variance)) {
                    continue;
                }

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Skewness)) {
                    continue;
                }

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);

                if (ddouble.IsNaN(dist.Kurtosis)) {
                    continue;
                }

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void EntropyTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (BradfordDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-28);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (BradfordDistribution dist in Dists) {
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
            foreach (BradfordDistribution dist in Dists) {
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

            foreach (BradfordDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 5) * 0.1, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (BradfordDistribution dist in Dists) {
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
            ddouble[] expected_dist_c1 = [
                1.442695040888963387e+00,
                1.437081441508072599e+00,
                1.431511358401452050e+00,
                1.425984287519593163e+00,
                1.420499732567594764e+00,
                1.415057204856607687e+00,
                1.409656223158681776e+00,
                1.404296313564922682e+00,
                1.398977009346873634e+00,
                1.393697850821036477e+00,
                1.388458385216446045e+00,
                1.383258166545223489e+00,
                1.378096755476024837e+00,
                1.372973719210314636e+00,
                1.367888631361387608e+00,
                1.362841071836068840e+00,
                1.357830626719024325e+00,
                1.352856888159614135e+00,
                1.347919454261221395e+00,
                1.343017928972998654e+00,
                1.338151921983966153e+00,
                1.333321048619403149e+00,
                1.328524929739477223e+00,
                1.323763191640052517e+00,
                1.319035465955623598e+00,
                1.314341389564322560e+00,
                1.309680604494945610e+00,
                1.305052757835952759e+00,
                1.300457501646389513e+00,
                1.295894492868683079e+00,
                1.291363393243267987e+00,
                1.286863869224998735e+00,
                1.282395591901300813e+00,
                1.277958236912023038e+00,
                1.273551484370947007e+00,
                1.269175018788916365e+00,
                1.264828528998543256e+00,
                1.260511708080459536e+00,
                1.256224253291070125e+00,
                1.251965865991778504e+00,
                1.247736251579644184e+00,
                1.243535119419443280e+00,
                1.239362182777096200e+00,
                1.235217158754430278e+00,
                1.231099768225248914e+00,
                1.227009735772673249e+00,
                1.222946789627730624e+00,
                1.218910661609157176e+00,
                1.214901087064390151e+00,
                1.210917804811720178e+00,
                1.206960557083577301e+00,
                1.203029089470927238e+00,
                1.199123150868748988e+00,
                1.195242493422571828e+00,
                1.191386872476047243e+00,
                1.187556046519532504e+00,
                1.183749777139662340e+00,
                1.179967828969886945e+00,
                1.176209969641957453e+00,
                1.172475969738332235e+00,
                1.168765602745489351e+00,
                1.165078645008121860e+00,
                1.161414875684197101e+00,
                1.157774076700860855e+00,
                1.154156032711170843e+00,
                1.150560531051634472e+00,
                1.146987361700542385e+00,
                1.143436317237073174e+00,
                1.139907192801156155e+00,
                1.136399786054075767e+00,
                1.132913897139799619e+00,
                1.129449328647017303e+00,
                1.126005885571873977e+00,
                1.122583375281381857e+00,
                1.119181607477498952e+00,
                1.115800394161856834e+00,
                1.112439549601128563e+00,
                1.109098890293017003e+00,
                1.105778234932858206e+00,
                1.102477404380819870e+00,
                1.099196221629686443e+00,
                1.095934511773218656e+00,
                1.092692101975072827e+00,
                1.089468821438273283e+00,
                1.086264501375219460e+00,
                1.083078974978224807e+00,
                1.079912077390569047e+00,
                1.076763645678060266e+00,
                1.073633518801089037e+00,
                1.070521537587172833e+00,
                1.067427544703972941e+00,
                1.064351384632779896e+00,
                1.061292903642455876e+00,
                1.058251949763824085e+00,
                1.055228372764498879e+00,
                1.052222024124144228e+00,
                1.049232757010155170e+00,
                1.046260426253752485e+00,
                1.043304888326481938e+00,
                1.040366001317111655e+00,
                1.037443624908917617e+00,
                1.034537620357351972e+00,
                1.031647850468085714e+00,
                1.028774179575416836e+00,
                1.025916473521040650e+00,
                1.023074599633170711e+00,
                1.020248426706007239e+00,
                1.017437824979544603e+00,
                1.014642666119710546e+00,
                1.011862823198834693e+00,
                1.009098170676433481e+00,
                1.006348584380312383e+00,
                1.003613941487974559e+00,
                1.000894120508332374e+00,
                9.981890012637153475e-01,
                9.954984648721688645e-01,
                9.928223937300394430e-01,
                9.901606714948382182e-01,
                9.875131830683813172e-01,
                9.848798145801990200e-01,
                9.822604533712091524e-01,
                9.796549879776514924e-01,
                9.770633081152767518e-01,
                9.744853046637853566e-01,
                9.719208696515122536e-01,
                9.693698962403534125e-01,
                9.668322787109284722e-01,
                9.643079124479755349e-01,
                9.617966939259755543e-01,
                9.592985206949991017e-01,
                9.568132913667737327e-01,
                9.543409056009680702e-01,
                9.518812640916871626e-01,
                9.494342685541763416e-01,
                9.469998217117298056e-01,
                9.445778272827997091e-01,
                9.421681899683025385e-01,
                9.397708154391213320e-01,
                9.373856103237934168e-01,
                9.350124821963915478e-01,
                9.326513395645824600e-01,
                9.303020918578706056e-01,
                9.279646494160167514e-01,
                9.256389234776306596e-01,
                9.233248261689366299e-01,
                9.210222704927047976e-01,
                9.187311703173498545e-01,
                9.164514403661901953e-01,
                9.141829962068679372e-01,
                9.119257542409250350e-01,
                9.096796316935336080e-01,
                9.074445466033774776e-01,
                9.052204178126829204e-01,
                9.030071649573951920e-01,
                9.008047084574991592e-01,
                8.986129695074809298e-01,
                8.964318700669287043e-01,
                8.942613328512704074e-01,
                8.921012813226439908e-01,
                8.899516396809028507e-01,
                8.878123328547468107e-01,
                8.856832864929846671e-01,
                8.835644269559201902e-01,
                8.814556813068608054e-01,
                8.793569773037492876e-01,
                8.772682433909136934e-01,
                8.751894086909352000e-01,
                8.731204029966302960e-01,
                8.710611567631476593e-01,
                8.690116011001756124e-01,
                8.669716677642597125e-01,
                8.649412891512286983e-01,
                8.629203982887256874e-01,
                8.609089288288452879e-01,
                8.589068150408712965e-01,
                8.569139918041174253e-01,
                8.549303946008671717e-01,
                8.529559595094102775e-01,
                8.509906231971765544e-01,
                8.490343229139646342e-01,
                8.470869964852629774e-01,
                8.451485823056628099e-01,
                8.432190193323622074e-01,
                8.412982470787576661e-01,
                8.393862056081241585e-01,
                8.374828355273801206e-01,
                8.355880779809381487e-01,
                8.337018746446379636e-01,
                8.318241677197627526e-01,
                8.299548999271341154e-01,
                8.280940145012883891e-01,
                8.262414551847306887e-01,
                8.243971662222647767e-01,
                8.225610923554000919e-01,
                8.207331788168324982e-01,
                8.189133713249991953e-01,
                8.171016160787049065e-01,
                8.152978597518204529e-01,
                8.135020494880499387e-01,
                8.117141328957684365e-01,
                8.099340580429268410e-01,
                8.081617734520232244e-01,
                8.063972280951411387e-01,
                8.046403713890516451e-01,
                8.028911531903796250e-01,
                8.011495237908343725e-01,
                7.994154339124992514e-01,
                7.976888347031850790e-01,
                7.959696777318419070e-01,
                7.942579149840315322e-01,
                7.925534988574562956e-01,
                7.908563821575474995e-01,
                7.891665180931083379e-01,
                7.874838602720141134e-01,
                7.858083626969672997e-01,
                7.841399797613050060e-01,
                7.824786662448615093e-01,
                7.808243773098829665e-01,
                7.791770684969928640e-01,
                7.775366957212097585e-01,
                7.759032152680139793e-01,
                7.742765837894647341e-01,
                7.726567583003652873e-01,
                7.710436961744772111e-01,
                7.694373551407804879e-01,
                7.678376932797810195e-01,
                7.662446690198644328e-01,
                7.646582411336949603e-01,
                7.630783687346583966e-01,
                7.615050112733497967e-01,
                7.599381285341042513e-01,
                7.583776806315701835e-01,
                7.568236280073251665e-01,
                7.552759314265330426e-01,
                7.537345519746422307e-01,
                7.521994510541234469e-01,
                7.506705903812493919e-01,
                7.491479319829099648e-01,
                7.476314381934709985e-01,
                7.461210716516659680e-01,
                7.446167952975295545e-01,
                7.431185723693655154e-01,
                7.416263664007521905e-01,
                7.401401412175844241e-01,
                7.386598609351493483e-01,
                7.371854899552388041e-01,
                7.357169929632960681e-01,
                7.342543349255956642e-01,
                7.327974810864575916e-01,
                7.313463969654943053e-01,
                7.299010483548905581e-01,
                7.284614013167153290e-01,
                7.270274221802650594e-01,
                7.255990775394393077e-01,
                7.241763342501463807e-01,
                7.227591594277389397e-01,
                7.213475204444816935e-01,
            ];
            ddouble[] expected_dist_c2 = [
                1.820478453253674633e+00,
                1.806366217181940792e+00,
                1.792471092434387403e+00,
                1.778788106995956797e+00,
                1.765312439518714749e+00,
                1.752039413657671796e+00,
                1.738964492660226568e+00,
                1.726083274196076767e+00,
                1.713391485415223237e+00,
                1.700884978222411359e+00,
                1.688559724757031644e+00,
                1.676411813068132028e+00,
                1.664437442974788128e+00,
                1.652632922102626578e+00,
                1.640994662087819522e+00,
                1.629519174940352144e+00,
                1.618203069558821872e+00,
                1.607043048389450668e+00,
                1.596035904222399493e+00,
                1.585178517118845942e+00,
                1.574467851462637569e+00,
                1.563900953130673566e+00,
                1.553474946776469023e+00,
                1.543187033221657956e+00,
                1.533034486950462849e+00,
                1.523014653702420729e+00,
                1.513124948158898420e+00,
                1.503362851719163507e+00,
                1.493725910361989317e+00,
                1.484211732588983201e+00,
                1.474817987446014778e+00,
                1.465542402619310458e+00,
                1.456382762602939751e+00,
                1.447336906934598266e+00,
                1.438402728496730454e+00,
                1.429578171880186321e+00,
                1.420861231807746039e+00,
                1.412249951614971799e+00,
                1.403742421785965977e+00,
                1.395336778541738765e+00,
                1.387031202478990144e+00,
                1.378823917257220977e+00,
                1.370713188332178589e+00,
                1.362697321733744582e+00,
                1.354774662886455427e+00,
                1.346943595470926791e+00,
                1.339202540324542223e+00,
                1.331549954379830591e+00,
                1.323984329639036117e+00,
                1.316504192183448385e+00,
                1.309108101216125641e+00,
                1.301794648136705979e+00,
                1.294562455647057631e+00,
                1.287410176886576441e+00,
                1.280336494595990970e+00,
                1.273340120308581236e+00,
                1.266419793567773677e+00,
                1.259574281170110055e+00,
                1.252802376432636366e+00,
                1.246102898483798516e+00,
                1.239474691576969878e+00,
                1.232916624425768992e+00,
                1.226427589560370235e+00,
                1.220006502704033435e+00,
                1.213652302169116348e+00,
                1.207363948271867216e+00,
                1.201140422765311078e+00,
                1.194980728289591454e+00,
                1.188883887839134346e+00,
                1.182848944246042366e+00,
                1.176874959679143240e+00,
                1.170961015158142438e+00,
                1.165106210082351712e+00,
                1.159309661773484379e+00,
                1.153570505032031512e+00,
                1.147887891706750541e+00,
                1.142260990276815491e+00,
                1.136688985446196787e+00,
                1.131171077749856035e+00,
                1.125706483171354355e+00,
                1.120294432771492099e+00,
                1.114934172327609385e+00,
                1.109624961983192160e+00,
                1.104366075907442513e+00,
                1.099156801964482844e+00,
                1.093996441391879682e+00,
                1.088884308488179187e+00,
                1.083819730309164342e+00,
                1.078802046372547840e+00,
                1.073830608370831108e+00,
                1.068904779892065893e+00,
                1.064023936148266403e+00,
                1.059187463711228849e+00,
                1.054394760255522145e+00,
                1.049645234308425046e+00,
                1.044938305006593460e+00,
                1.040273401859242552e+00,
                1.035649964517646016e+00,
                1.031067442550753777e+00,
                1.026525295226741630e+00,
                1.022022991300308492e+00,
                1.017560008805547422e+00,
                1.013135834854218986e+00,
                1.008749965439265539e+00,
                1.004401905243406778e+00,
                1.000091167452662466e+00,
                9.958172735746597670e-01,
                9.915797532615759913e-01,
                9.873781441375862888e-01,
                9.832119916306765184e-01,
                9.790808488086989447e-01,
                9.749842762195410906e-01,
                9.709218417352931674e-01,
                9.668931204002919610e-01,
                9.628976942829353680e-01,
                9.589351523311537839e-01,
                9.550050902314358714e-01,
                9.511071102713076764e-01,
                9.472408212051641740e-01,
                9.434058381233617618e-01,
                9.396017823244772194e-01,
                9.358282811906438736e-01,
                9.320849680658814806e-01,
                9.283714821373321158e-01,
                9.246874683193267996e-01,
                9.210325771401989270e-01,
                9.174064646317729999e-01,
                9.138087922214523928e-01,
                9.102392266268373167e-01,
                9.066974397528029250e-01,
                9.031831085909703960e-01,
                8.996959151215071504e-01,
                8.962355462171937015e-01,
                8.928016935496947415e-01,
                8.893940534979783985e-01,
                8.860123270588226019e-01,
                8.826562197593573744e-01,
                8.793254415715862082e-01,
                8.760197068288358979e-01,
                8.727387341440837609e-01,
                8.694822463301132842e-01,
                8.662499703214511237e-01,
                8.630416370980383833e-01,
                8.598569816105916530e-01,
                8.566957427076116183e-01,
                8.535576630639939433e-01,
                8.504424891112056795e-01,
                8.473499709689831239e-01,
                8.442798623785158219e-01,
                8.412319206370770797e-01,
                8.382059065340660142e-01,
                8.352015842884241703e-01,
                8.322187214873940642e-01,
                8.292570890265847927e-01,
                8.263164610513132891e-01,
                8.233966148991884726e-01,
                8.204973310439097611e-01,
                8.176183930402468603e-01,
                8.147595874701760721e-01,
                8.119207038901405937e-01,
                8.091015347794109358e-01,
                8.063018754895168172e-01,
                8.035215241947253340e-01,
                8.007602818435407555e-01,
                7.980179521111997465e-01,
                7.952943413531411432e-01,
                7.925892585594229711e-01,
                7.899025153100689423e-01,
                7.872339257313187844e-01,
                7.845833064527620859e-01,
                7.819504765653367828e-01,
                7.793352575801684168e-01,
                7.767374733882345117e-01,
                7.741569502208317521e-01,
                7.715935166108289778e-01,
                7.690470033546876749e-01,
                7.665172434752314246e-01,
                7.640040721851486527e-01,
                7.615073268512103644e-01,
                7.590268469591868739e-01,
                7.565624740794492098e-01,
                7.541140518332374310e-01,
                7.516814258595817533e-01,
                7.492644437828628323e-01,
                7.468629551809946587e-01,
                7.444768115542182985e-01,
                7.421058662944916007e-01,
                7.397499746554614397e-01,
                7.374089937230073888e-01,
                7.350827823863418375e-01,
                7.327712013096552290e-01,
                7.304741129042957581e-01,
                7.281913813014698755e-01,
                7.259228723254528282e-01,
                7.236684534672991331e-01,
                7.214279938590414387e-01,
                7.192013642483652269e-01,
                7.169884369737549168e-01,
                7.147890859400931607e-01,
                7.126031865947105581e-01,
                7.104306159038730195e-01,
                7.082712523296971208e-01,
                7.061249758074858995e-01,
                7.039916677234754339e-01,
                7.018712108929829885e-01,
                6.997634895389499565e-01,
                6.976683892708693824e-01,
                6.955857970640906274e-01,
                6.935156012394950720e-01,
                6.914576914435321964e-01,
                6.894119586286104884e-01,
                6.873782950338358511e-01,
                6.853565941660892946e-01,
                6.833467507814380282e-01,
                6.813486608668722910e-01,
                6.793622216223624832e-01,
                6.773873314432277137e-01,
                6.754238899028126575e-01,
                6.734717977354633955e-01,
                6.715309568197992407e-01,
                6.696012701622711116e-01,
                6.676826418810039021e-01,
                6.657749771899152957e-01,
                6.638781823831064743e-01,
                6.619921648195180586e-01,
                6.601168329078480612e-01,
                6.582520960917241926e-01,
                6.563978648351277423e-01,
                6.545540506080628207e-01,
                6.527205658724659632e-01,
                6.508973240683529893e-01,
                6.490842396001960113e-01,
                6.472812278235288153e-01,
                6.454882050317738429e-01,
                6.437050884432882203e-01,
                6.419317961886236157e-01,
                6.401682472979954852e-01,
                6.384143616889599304e-01,
                6.366700601542906179e-01,
                6.349352643500554283e-01,
                6.332098967838868386e-01,
                6.314938808034427087e-01,
                6.297871405850550275e-01,
                6.280896011225616249e-01,
                6.264011882163181832e-01,
                6.247218284623869966e-01,
                6.230514492418992578e-01,
                6.213899787105875427e-01,
                6.197373457884849390e-01,
                6.180934801497887099e-01,
                6.164583122128844961e-01,
                6.148317731305286138e-01,
                6.132137947801851174e-01,
                6.116043097545154072e-01,
                6.100032513520167177e-01,
                6.084105535678077104e-01,
                6.068261510845581741e-01,
            ];

            foreach ((BradfordDistribution dist, ddouble[] expecteds) in new[]{
                (dist_c1, expected_dist_c1),
                (dist_c2, expected_dist_c2),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 256, i++) {
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
            ddouble[] expected_dist_c1 = [
                0.000000000000000000e+00,
                5.624549193878106666e-03,
                1.122725542325412121e-02,
                1.680828768655389174e-02,
                2.236781302845450986e-02,
                2.790599656988448571e-02,
                3.342300153745027952e-02,
                3.891898929230235005e-02,
                4.439411935845343632e-02,
                4.984854945056153147e-02,
                5.528243550118960847e-02,
                6.069593168755394591e-02,
                6.608919045777243706e-02,
                7.146236255662415104e-02,
                7.681559705083089440e-02,
                8.214904135387156303e-02,
                8.746284125033940149e-02,
                9.275714091985244603e-02,
                9.803208296052672022e-02,
                1.032878084120219531e-01,
                1.085244567781690622e-01,
                1.137421660491883257e-01,
                1.189410727235074294e-01,
                1.241213118291875778e-01,
                1.292830169449664668e-01,
                1.344263202209260988e-01,
                1.395513523987935711e-01,
                1.446582428318823332e-01,
                1.497471195046820580e-01,
                1.548181090521040248e-01,
                1.598713367783894113e-01,
                1.649069266756877927e-01,
                1.699250014423123734e-01,
                1.749256825006788030e-01,
                1.799090900149344918e-01,
                1.848753429082838629e-01,
                1.898245588800172301e-01,
                1.947568544222478548e-01,
                1.996723448363643960e-01,
                2.045711442492036403e-01,
                2.094533656289497836e-01,
                2.143191208007658166e-01,
                2.191685204621615646e-01,
                2.240016741981050441e-01,
                2.288186904958808554e-01,
                2.336196767597020496e-01,
                2.384047393250789404e-01,
                2.431739834729509386e-01,
                2.479275134435855177e-01,
                2.526654324502486393e-01,
                2.573878426926517471e-01,
                2.620948453701794079e-01,
                2.667865406949013751e-01,
                2.714630279043745986e-01,
                2.761244052742375943e-01,
                2.807707701306025871e-01,
                2.854022188622483136e-01,
                2.900188469326183127e-01,
                2.946207488916270378e-01,
                2.992080183872788401e-01,
                3.037807481771029328e-01,
                3.083390301394072219e-01,
                3.128829552843553907e-01,
                3.174126137648694002e-01,
                3.219280948873623482e-01,
                3.264294871223031302e-01,
                3.309168781146169525e-01,
                3.353903546939249192e-01,
                3.398500028846247467e-01,
                3.442959079158168789e-01,
                3.487281542310775584e-01,
                3.531468254980825727e-01,
                3.575520046180836742e-01,
                3.619437737352415030e-01,
                3.663222142458157915e-01,
                3.706874068072176787e-01,
                3.750394313469247454e-01,
                3.793783670712621570e-01,
                3.837042924740522443e-01,
                3.880172853451347437e-01,
                3.923174227787603052e-01,
                3.966047811818584345e-01,
                4.008794362821843094e-01,
                4.051414631363438601e-01,
                4.093909361377017775e-01,
                4.136279290241725026e-01,
                4.178525148858978633e-01,
                4.220647661728123556e-01,
                4.262647547020979588e-01,
                4.304525516655314243e-01,
                4.346282276367246511e-01,
                4.387918525782609214e-01,
                4.429434958487283303e-01,
                4.470832262096522958e-01,
                4.512111118323288150e-01,
                4.553272203045607425e-01,
                4.594316186372972566e-01,
                4.635243732711803455e-01,
                4.676055500829974787e-01,
                4.716752143920444618e-01,
                4.757334309663977523e-01,
                4.797802640290997345e-01,
                4.838157772642563970e-01,
                4.878400338230514111e-01,
                4.918530963296747216e-01,
                4.958550268871710354e-01,
                4.998458870832053758e-01,
                5.038257379957506998e-01,
                5.077946401986963565e-01,
                5.117526537673796616e-01,
                5.156998382840425332e-01,
                5.196362528432127981e-01,
                5.235619560570129449e-01,
                5.274770060603960475e-01,
                5.313814605163120763e-01,
                5.352753766208033781e-01,
                5.391588111080313217e-01,
                5.430318202552377738e-01,
                5.468944598876366303e-01,
                5.507467853832432869e-01,
                5.545888516776373844e-01,
                5.584207132686643815e-01,
                5.622424242210726231e-01,
                5.660540381710917890e-01,
                5.698556083309478382e-01,
                5.736471874933219972e-01,
                5.774288280357486869e-01,
                5.812005819249571603e-01,
                5.849625007211561867e-01,
                5.887146355822636679e-01,
                5.924570372680804109e-01,
                5.961897561444102767e-01,
                5.999128421871277039e-01,
                6.036263449861919428e-01,
                6.073303137496106618e-01,
                6.110247973073522632e-01,
                6.147098441152082371e-01,
                6.183855022586065564e-01,
                6.220518194563762204e-01,
                6.257088430644652810e-01,
                6.293566200796096854e-01,
                6.329951971429578217e-01,
                6.366246205436488781e-01,
                6.402449362223457952e-01,
                6.438561897747246965e-01,
                6.474584264549203549e-01,
                6.510516911789285821e-01,
                6.546360285279674285e-01,
                6.582114827517947520e-01,
                6.617780977719871505e-01,
                6.653359171851762621e-01,
                6.688849842662470957e-01,
                6.724253419714956159e-01,
                6.759570329417488033e-01,
                6.794800995054460779e-01,
                6.829945836816828653e-01,
                6.865005271832184119e-01,
                6.899979714194454106e-01,
                6.934869574993252073e-01,
                6.969675262342871491e-01,
                7.004397181410921824e-01,
                7.039035734446636994e-01,
                7.073591320808827465e-01,
                7.108064336993515919e-01,
                7.142455176661226535e-01,
                7.176764230663961186e-01,
                7.210991887071851458e-01,
                7.245138531199497578e-01,
                7.279204545631992040e-01,
                7.313190310250641257e-01,
                7.347096202258381892e-01,
                7.380922596204904096e-01,
                7.414669864011469436e-01,
                7.448338374995455702e-01,
                7.481928495894603071e-01,
                7.515440590890981598e-01,
                7.548875021634685600e-01,
                7.582232147267249367e-01,
                7.615512324444793091e-01,
                7.648715907360906785e-01,
                7.681843247769263305e-01,
                7.714894695005983793e-01,
                7.747870596011733335e-01,
                7.780771295353582362e-01,
                7.813597135246597158e-01,
                7.846348455575206104e-01,
                7.879025593914316117e-01,
                7.911628885550183732e-01,
                7.944158663501059703e-01,
                7.976615258537601560e-01,
                8.008998999203047475e-01,
                8.041310211833178068e-01,
                8.073549220576040630e-01,
                8.105716347411469069e-01,
                8.137811912170371809e-01,
                8.169836232553809863e-01,
                8.201789624151878400e-01,
                8.233672400462351826e-01,
                8.265484872909150127e-01,
                8.297227350860586492e-01,
                8.328900141647417321e-01,
                8.360503550580697940e-01,
                8.392037880969439589e-01,
                8.423503434138079893e-01,
                8.454900509443752377e-01,
                8.486229404293380574e-01,
                8.517490414160575618e-01,
                8.548683832602363974e-01,
                8.579809951275720881e-01,
                8.610869059953937255e-01,
                8.641861446542802305e-01,
                8.672787397096619610e-01,
                8.703647195834045558e-01,
                8.734441125153765695e-01,
                8.765169465649997882e-01,
                8.795832496127832245e-01,
                8.826430493618413475e-01,
                8.856963733393953264e-01,
                8.887432488982591750e-01,
                8.917837032183102419e-01,
                8.948177633079436033e-01,
                8.978454560055116884e-01,
                9.008668079807485851e-01,
                9.038818457361802450e-01,
                9.068905956085184794e-01,
                9.098930837700419660e-01,
                9.128893362299616010e-01,
                9.158793788357731946e-01,
                9.188632372745945132e-01,
                9.218409370744899967e-01,
                9.248125036057808224e-01,
                9.277779620823420892e-01,
                9.307373375628863466e-01,
                9.336906549522338006e-01,
                9.366379390025706408e-01,
                9.395792143146930453e-01,
                9.425145053392399719e-01,
                9.454438363779116283e-01,
                9.483672315846776169e-01,
                9.512847149669719782e-01,
                9.541963103868752460e-01,
                9.571020415622862876e-01,
                9.600019320680810431e-01,
                9.628960053372604966e-01,
                9.657842846620869892e-01,
                9.686667931952085420e-01,
                9.715435539507719653e-01,
                9.744145898055270871e-01,
                9.772799234999165474e-01,
                9.801395776391571557e-01,
                9.829935746943102570e-01,
                9.858419370033404405e-01,
                9.886846867721658105e-01,
                9.915218460756953789e-01,
                9.943534368588580197e-01,
                9.971794809376214319e-01,
                1.000000000000000000e+00,
            ];
            ddouble[] expected_dist_c2 = [
                0.000000000000000000e+00,
                7.083609497477530784e-03,
                1.411251876197523736e-02,
                2.108756612364194594e-02,
                2.800957078685089099e-02,
                3.487933340759554607e-02,
                4.169763664925946883e-02,
                4.846524571772510559e-02,
                5.518290887673613782e-02,
                6.185135794438179069e-02,
                6.847130877152649009e-02,
                7.504346170296732899e-02,
                8.156850202206224199e-02,
                8.804710037953442092e-02,
                9.447991320712462171e-02,
                1.008675831167287068e-01,
                1.072107392856276875e-01,
                1.135099978283877487e-01,
                1.197659621559798576e-01,
                1.259792233226424418e-01,
                1.321503603609856525e-01,
                1.382799406058125014e-01,
                1.443685200071097019e-01,
                1.504166434326395796e-01,
                1.564248449605455948e-01,
                1.623936481623638184e-01,
                1.683235663768156976e-01,
                1.742151029747399382e-01,
                1.800687516155054790e-01,
                1.858849964952324629e-01,
                1.916643125871325126e-01,
                1.974071658742676538e-01,
                2.031140135750122977e-01,
                2.087853043614912862e-01,
                2.144214785712553750e-01,
                2.200229684124429264e-01,
                2.255901981626681208e-01,
                2.311235843618631713e-01,
                2.366235359992947263e-01,
                2.420904546949641800e-01,
                2.475247348755924837e-01,
                2.529267639453833261e-01,
                2.582969224517484563e-01,
                2.636355842461732268e-01,
                2.689431166403923879e-01,
                2.742198805580379761e-01,
                2.794662306819180042e-01,
                2.846825155970745258e-01,
                2.898690779297657394e-01,
                2.950262544825113520e-01,
                3.001543763653332642e-01,
                3.052537691233204176e-01,
                3.103247528606399297e-01,
                3.153676423611132540e-01,
                3.203827472054703307e-01,
                3.253703718853913629e-01,
                3.303308159144413558e-01,
                3.352643739359978947e-01,
                3.401713358282701383e-01,
                3.450519868065018980e-01,
                3.499066075224492867e-01,
                3.547354741612201434e-01,
                3.595388585355578925e-01,
                3.643170281776503816e-01,
                3.690702464285425255e-01,
                3.737987725252256976e-01,
                3.785028616854770878e-01,
                3.831827651905177490e-01,
                3.878387304655573353e-01,
                3.924710011582881464e-01,
                3.970798172153934269e-01,
                4.016654149571269317e-01,
                4.062280271500245954e-01,
                4.107678830778019874e-01,
                4.152852086104933949e-01,
                4.197802262718828281e-01,
                4.242531553052786841e-01,
                4.287042117376804184e-01,
                4.331336084423835775e-01,
                4.375415552000689878e-01,
                4.419282587584203448e-01,
                4.462939228903113342e-01,
                4.506387484506047536e-01,
                4.549629334316021567e-01,
                4.592666730171825473e-01,
                4.635501596356671472e-01,
                4.678135830114467653e-01,
                4.720571302154046855e-01,
                4.762809857141702685e-01,
                4.804853314182350177e-01,
                4.846703467289625311e-01,
                4.888362085845223692e-01,
                4.929830915047780371e-01,
                4.971111676351565034e-01,
                5.012206067895281780e-01,
                5.053115764921227715e-01,
                5.093842420185072939e-01,
                5.134387664356522274e-01,
                5.174753106411080772e-01,
                5.214940334013178358e-01,
                5.254950913890881203e-01,
                5.294786392202384118e-01,
                5.334448294894535980e-01,
                5.373938128053582508e-01,
                5.413257378248328422e-01,
                5.452407512865928840e-01,
                5.491389980440479768e-01,
                5.530206210974615288e-01,
                5.568857616254262455e-01,
                5.607345590156750381e-01,
                5.645671508952435635e-01,
                5.683836731599990477e-01,
                5.721842600035548232e-01,
                5.759690439455809052e-01,
                5.797381558595314788e-01,
                5.834917249997978450e-01,
                5.872298790283061454e-01,
                5.909527440405696330e-01,
                5.946604445912107018e-01,
                5.983531037189658841e-01,
                6.020308429711850318e-01,
                6.056937824278373350e-01,
                6.093420407250368376e-01,
                6.129757350780984382e-01,
                6.165949813041350369e-01,
                6.201998938442071507e-01,
                6.237905857850361002e-01,
                6.273671688802909818e-01,
                6.309297535714574190e-01,
                6.344784490083014150e-01,
                6.380133630689348578e-01,
                6.415346023794930019e-01,
                6.450422723334326980e-01,
                6.485364771104605852e-01,
                6.520173196950993511e-01,
                6.554849018948999406e-01,
                6.589393243583082649e-01,
                6.623806865921947340e-01,
                6.658090869790529442e-01,
                6.692246227938757341e-01,
                6.726273902207168254e-01,
                6.760174843689423785e-01,
                6.793949992891824552e-01,
                6.827600279889870505e-01,
                6.861126624481935776e-01,
                6.894529936340128007e-01,
                6.927811115158390987e-01,
                6.960971050797902793e-01,
                6.994010623429838258e-01,
                7.026930703675555723e-01,
                7.059732152744246925e-01,
                7.092415822568127748e-01,
                7.124982555935196471e-01,
                7.157433186619641674e-01,
                7.189768539509917566e-01,
                7.221989430734570004e-01,
                7.254096667785820962e-01,
                7.286091049641003625e-01,
                7.317973366881861397e-01,
                7.349744401811755257e-01,
                7.381404928570850510e-01,
                7.412955713249296252e-01,
                7.444397513998451954e-01,
                7.475731081140196688e-01,
                7.506957157274373182e-01,
                7.538076477384384466e-01,
                7.569089768940998608e-01,
                7.599997752004384877e-01,
                7.630801139324430160e-01,
                7.661500636439360079e-01,
                7.692096941772699203e-01,
                7.722590746728617006e-01,
                7.752982735785670654e-01,
                7.783273586588997928e-01,
                7.813463970040969153e-01,
                7.843554550390360314e-01,
                7.873545985320029583e-01,
                7.903438926033183876e-01,
                7.933234017338212096e-01,
                7.962931897732138387e-01,
                7.992533199482730888e-01,
                8.022038548709260475e-01,
                8.051448565461972739e-01,
                8.080763863800267677e-01,
                8.109985051869628148e-01,
                8.139112731977324877e-01,
                8.168147500666897987e-01,
                8.197089948791472791e-01,
                8.225940661585898761e-01,
                8.254700218737757167e-01,
                8.283369194457250728e-01,
                8.311948157545985261e-01,
                8.340437671464696612e-01,
                8.368838294399892908e-01,
                8.397150579329487607e-01,
                8.425375074087391170e-01,
                8.453512321427127940e-01,
                8.481562859084448291e-01,
                8.509527219839003731e-01,
                8.537405931575050566e-01,
                8.565199517341255397e-01,
                8.592908495409566916e-01,
                8.620533379333206181e-01,
                8.648074678003778715e-01,
                8.675532895707521730e-01,
                8.702908532180707590e-01,
                8.730202082664215713e-01,
                8.757414037957290676e-01,
                8.784544884470498749e-01,
                8.811595104277891721e-01,
                8.838565175168406896e-01,
                8.865455570696504362e-01,
                8.892266760232058198e-01,
                8.918999209009507156e-01,
                8.945653378176307013e-01,
                8.972229724840646847e-01,
                8.998728702118498068e-01,
                9.025150759179961790e-01,
                9.051496341294953396e-01,
                9.077765889878228744e-01,
                9.103959842533754232e-01,
                9.130078633098452912e-01,
                9.156122691685320003e-01,
                9.182092444725906688e-01,
                9.207988315012231029e-01,
                9.233810721738073823e-01,
                9.259560080539688265e-01,
                9.285236803535943384e-01,
                9.310841299367906831e-01,
                9.336373973237859225e-01,
                9.361835226947778921e-01,
                9.387225458937266120e-01,
                9.412545064320972932e-01,
                9.437794434925486087e-01,
                9.462973959325706730e-01,
                9.488084022880739488e-01,
                9.513125007769276387e-01,
                9.538097293024495604e-01,
                9.563001254568486154e-01,
                9.587837265246211826e-01,
                9.612605694858986638e-01,
                9.637306910197531717e-01,
                9.661941275074552582e-01,
                9.686509150356898967e-01,
                9.711010893997275017e-01,
                9.735446861065545354e-01,
                9.759817403779591505e-01,
                9.784122871535794186e-01,
                9.808363610939065946e-01,
                9.832539965832531870e-01,
                9.856652277326775069e-01,
                9.880700883828725800e-01,
                9.904686121070152005e-01,
                9.928608322135786812e-01,
                9.952467817491076341e-01,
                9.976264935009581114e-01,
                1.000000000000000000e+00,
            ];

            foreach ((BradfordDistribution dist, ddouble[] expecteds) in new[]{
                (dist_c1, expected_dist_c1),
                (dist_c2, expected_dist_c2),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 256, i++) {
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