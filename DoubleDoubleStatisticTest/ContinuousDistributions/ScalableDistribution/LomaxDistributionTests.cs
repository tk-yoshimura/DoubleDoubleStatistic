﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class LomaxDistributionTests {
        readonly LomaxDistribution dist_alpha2theta1 = new(alpha: 2, theta: 1);
        readonly LomaxDistribution dist_alpha3theta2 = new(alpha: 3, theta: 2);
        readonly LomaxDistribution dist_alpha4theta2 = new(alpha: 4, theta: 2);
        readonly LomaxDistribution dist_alpha4theta3 = new(alpha: 4, theta: 3);
        readonly LomaxDistribution dist_alpha5theta6 = new(alpha: 5, theta: 6);

        LomaxDistribution[] Dists => [
            dist_alpha2theta1,
            dist_alpha3theta2,
            dist_alpha4theta2,
            dist_alpha4theta3,
            dist_alpha5theta6,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (LomaxDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"Alpha={dist.Alpha}");
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
            foreach (LomaxDistribution dist in Dists) {
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
            foreach (LomaxDistribution dist in Dists) {
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
            foreach (LomaxDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (LomaxDistribution dist in Dists) {
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
            foreach (LomaxDistribution dist in Dists) {
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
            foreach (LomaxDistribution dist in Dists) {
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
            foreach (LomaxDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (LomaxDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (LomaxDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (LomaxDistribution dist in Dists) {
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
            foreach (LomaxDistribution dist in Dists) {
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
            foreach (LomaxDistribution dist in Dists) {
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

            foreach (LomaxDistribution dist in Dists) {

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
        public void FitTest() {
            Random random = new(1234);

            foreach (LomaxDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                (LomaxDistribution? dist_fit, ddouble error) = LomaxDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (LomaxDistribution dist in Dists) {
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
            ddouble[] expected_dist_alpha2theta1 = [
                2.000000000000000000e+00,
                1.024000000000000021e+00,
                5.925925925925925597e-01,
                3.731778425655976617e-01,
                2.500000000000000000e-01,
                1.755829903978051998e-01,
                1.280000000000000027e-01,
                9.616829451540195861e-02,
                7.407407407407406996e-02,
                5.826126536185707860e-02,
                4.664723032069970771e-02,
                3.792592592592592543e-02,
                3.125000000000000000e-02,
                2.605332790555668632e-02,
                2.194787379972564997e-02,
                1.866161247995334546e-02,
                1.600000000000000033e-02,
                1.382140157650361793e-02,
                1.202103681442524483e-02,
                1.052025971891181085e-02,
                9.259259259259258745e-03,
                8.191999999999999629e-03,
                7.282658170232134826e-03,
                6.503073718437229686e-03,
                5.830903790087463463e-03,
                5.248267661650743779e-03,
                4.740740740740740679e-03,
                4.296599644187841863e-03,
                3.906250000000000000e-03,
                3.561788685755627824e-03,
                3.256665988194585790e-03,
                2.985422740524781415e-03,
                2.743484224965706247e-03,
                2.526997413776084153e-03,
                2.332701559994168183e-03,
                2.157824643031743764e-03,
                2.000000000000000042e-03,
                1.857198821839497289e-03,
                1.727675197062952242e-03,
                1.609921139019205921e-03,
                1.502629601803155603e-03,
                1.404663923182441739e-03,
                1.315032464863976356e-03,
                1.232867476378066455e-03,
                1.157407407407407343e-03,
                1.087982048296203116e-03,
                1.023999999999999954e-03,
                9.649380705761735231e-04,
                9.103322712790168532e-04,
                8.597701458250771008e-04,
                8.128842148046537108e-04,
                7.693463561232156229e-04,
                7.288629737609329329e-04,
                6.911708325908646989e-04,
                6.560334577063429724e-04,
                6.232380136235934898e-04,
                5.925925925925925849e-04,
                5.639238526572708807e-04,
                5.370749555234802328e-04,
                5.119037620927265660e-04,
                4.882812500000000000e-04,
                4.660901228948566223e-04,
                4.452235857194534780e-04,
                4.255842640218377684e-04,
            ];
            ddouble[] expected_dist_alpha3theta2 = [
                1.500000000000000000e+00,
                9.364426154549611026e-01,
                6.143999999999999462e-01,
                4.196434669762994507e-01,
                2.962962962962962798e-01,
                2.151185182591645872e-01,
                1.599333610995418709e-01,
                1.213629629629629642e-01,
                9.375000000000000000e-02,
                7.356233761568946317e-02,
                5.852766346593506891e-02,
                4.714512626514529708e-02,
                3.839999999999999664e-02,
                3.159177503200826659e-02,
                2.622771668601871567e-02,
                2.195532463077247420e-02,
                1.851851851851851749e-02,
                1.572863999999999859e-02,
                1.344490739119778670e-02,
                1.156101994388840891e-02,
                9.995835068721366928e-03,
                8.686787853766749640e-03,
                7.585185185185185260e-03,
                6.652799449065045605e-03,
                5.859375000000000000e-03,
                5.180783542917276677e-03,
                4.597646100980591448e-03,
                4.094294044148271902e-03,
                3.657978966620941807e-03,
                3.278266915168974189e-03,
                2.946570391571581068e-03,
                2.655784176039069182e-03,
                2.399999999999999790e-03,
                2.174281547519411450e-03,
                1.974485939500516662e-03,
                1.797121271463299587e-03,
                1.639232292876169729e-03,
                1.498308184727937869e-03,
                1.372207789423279638e-03,
                1.259098699279727600e-03,
                1.157407407407407343e-03,
                1.065778333024852134e-03,
                9.830399999999999121e-04,
                9.081770076011045625e-04,
                8.403067119498616689e-04,
                7.786597547095037812e-04,
                7.225637464930255568e-04,
                6.714295471620790713e-04,
                6.247396917950854330e-04,
                5.820385958659913140e-04,
                5.429242408604218525e-04,
                5.070410958293641510e-04,
                4.740740740740740787e-04,
                4.437433594680164121e-04,
                4.157999655665653503e-04,
                3.900219139754107067e-04,
                3.662109375000000000e-04,
                3.441896292146633444e-04,
                3.237989714323297923e-04,
                3.048961891499733291e-04,
                2.873528813112869655e-04,
                2.710533905033639065e-04,
                2.558933777592669939e-04,
                2.417785741919237729e-04,
            ];
            ddouble[] expected_dist_alpha4theta2 = [
                2.000000000000000000e+00,
                1.109857914613287200e+00,
                6.553600000000000536e-01,
                4.069269982800479135e-01,
                2.633744855967078413e-01,
                1.765075021613658202e-01,
                1.218539894091747455e-01,
                8.630255144032922265e-02,
                6.250000000000000000e-02,
                4.615676085690319230e-02,
                3.468305983166522499e-02,
                2.646743930674823678e-02,
                2.048000000000000168e-02,
                1.604661588927404001e-02,
                1.271646869625149730e-02,
                1.018217953890897427e-02,
                8.230452674897120041e-03,
                6.710886400000000244e-03,
                5.515859442542681880e-03,
                4.567316521042334362e-03,
                3.807937169036710798e-03,
                3.195140360006160692e-03,
                2.696954732510288208e-03,
                2.289135294301951166e-03,
                1.953125000000000000e-03,
                1.674596700740937916e-03,
                1.442398776778224759e-03,
                1.247784851549949524e-03,
                1.083845619739538281e-03,
                9.450859575261907376e-04,
                8.271074783358823993e-04,
                7.263683216517111477e-04,
                6.400000000000000524e-04,
                5.656667440700908480e-04,
                5.014567465398137502e-04,
                4.457975247040742986e-04,
                3.973896467578592905e-04,
                3.551545326762519522e-04,
                3.181931105909054460e-04,
                2.857528962904346231e-04,
                2.572016460905350013e-04,
                2.320061677333011129e-04,
                2.097152000000000076e-04,
                1.899455179296427891e-04,
                1.723706075794588087e-04,
                1.567113971742397564e-04,
                1.427286412825729488e-04,
                1.302166394496153301e-04,
                1.189980365323972124e-04,
                1.089195033199515859e-04,
                9.984813625019252163e-05,
                9.166844670361385436e-05,
                8.427983539094650650e-05,
                7.759446722937991727e-05,
                7.153547794693597395e-05,
                6.603545633446107740e-05,
                6.103515625000000000e-05,
                5.648240069163706185e-05,
                5.233114689815430987e-05,
                4.854068682984650238e-05,
                4.507496177431952373e-05,
                4.190197341114803664e-05,
                3.899327661093592262e-05,
                3.632354166263644039e-05,

            ];
            ddouble[] expected_dist_alpha4theta3 = [
                1.333333333333333259e+00,
                8.935692296919147681e-01,
                6.168858213839468752e-01,
                4.369066666666667209e-01,
                3.164062500000001110e-01,
                2.336686018380723651e-01,
                1.755829903978052275e-01,
                1.339914114904128817e-01,
                1.036800000000000499e-01,
                8.123599293944983035e-02,
                6.437712277477318501e-02,
                5.154728391572669516e-02,
                4.166666666666666435e-02,
                3.397386240000002666e-02,
                2.792403842787230528e-02,
                2.312203988777681782e-02,
                1.927768191824836067e-02,
                1.617539807253117953e-02,
                1.365333333333333503e-02,
                1.158874742740363499e-02,
                9.887695312499994796e-03,
                8.477645797500998198e-03,
                7.302143807439767480e-03,
                6.316910810971615615e-03,
                5.486968449931413361e-03,
                4.784497659976339606e-03,
                4.187231609075405155e-03,
                3.677239628361788065e-03,
                3.239999999999999391e-03,
                2.863687891854835386e-03,
                2.538624779357807199e-03,
                2.256849968814375954e-03,
                2.011785086711662899e-03,
                1.797969821673525400e-03,
                1.610852622366458356e-03,
                1.446624037470325324e-03,
                1.302083333333333261e-03,
                1.174531224149835988e-03,
                1.061683200000000833e-03,
                9.615991845188165063e-04,
                8.726262008710095400e-04,
                7.933514481945892717e-04,
                7.225637464930255568e-04,
                6.592217372136771310e-04,
                6.024275599452612708e-04,
                5.514049855572549690e-04,
                5.054811897665993604e-04,
                4.640715114370454197e-04,
                4.266666666666667196e-04,
                3.928219903487359735e-04,
                3.621483571063632682e-04,
                3.343044976932091849e-04,
                3.089904785156251084e-04,
                2.859421535014125570e-04,
                2.649264311719061937e-04,
                2.457372270760979920e-04,
                2.281919939824925440e-04,
                2.121287403939369550e-04,
                1.974034628428631506e-04,
                1.838879296670968956e-04,
                1.714677640603566675e-04,
                1.600407825375088514e-04,
                1.495155518742606127e-04,
                1.398101333333333384e-04,
            ];

            foreach ((LomaxDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha2theta1, expected_dist_alpha2theta1),
                (dist_alpha3theta2, expected_dist_alpha3theta2),
                (dist_alpha4theta2, expected_dist_alpha4theta2),
                (dist_alpha4theta3, expected_dist_alpha4theta3),
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
            ddouble[] expected_dist_alpha2theta1 = [
                0.000000000000000000e+00,
                3.599999999999999867e-01,
                5.555555555555555802e-01,
                6.734693877551020114e-01,
                7.500000000000000000e-01,
                8.024691358024691468e-01,
                8.400000000000000799e-01,
                8.677685950413223104e-01,
                8.888888888888889506e-01,
                9.053254437869822091e-01,
                9.183673469387755306e-01,
                9.288888888888888751e-01,
                9.375000000000000000e-01,
                9.446366782006920815e-01,
                9.506172839506172867e-01,
                9.556786703601107824e-01,
                9.599999999999999645e-01,
                9.637188208616780383e-01,
                9.669421487603305776e-01,
                9.697542533081285887e-01,
                9.722222222222222099e-01,
                9.744000000000000439e-01,
                9.763313609467455523e-01,
                9.780521262002743743e-01,
                9.795918367346938549e-01,
                9.809750297265160679e-01,
                9.822222222222222188e-01,
                9.833506763787721594e-01,
                9.843750000000000000e-01,
                9.853076216712580715e-01,
                9.861591695501730204e-01,
                9.869387755102040360e-01,
                9.876543209876543772e-01,
                9.883126369612855733e-01,
                9.889196675900276956e-01,
                9.894806048652202701e-01,
                9.899999999999999911e-01,
                9.904818560380725945e-01,
                9.909297052154194541e-01,
                9.913466738777717557e-01,
                9.917355371900826722e-01,
                9.920987654320987392e-01,
                9.924385633270321749e-01,
                9.927569035762788818e-01,
                9.930555555555555802e-01,
                9.933361099541857531e-01,
                9.936000000000000387e-01,
                9.938485198000769349e-01,
                9.940828402366863603e-01,
                9.943040227839088363e-01,
                9.945130315500685381e-01,
                9.947107438016529191e-01,
                9.948979591836735192e-01,
                9.950754078177901363e-01,
                9.952437574316289615e-01,
                9.954036196495259503e-01,
                9.955555555555555269e-01,
                9.957000806234883594e-01,
                9.958376690946930676e-01,
                9.959687578735197944e-01,
                9.960937500000000000e-01,
                9.962130177514793017e-01,
                9.963269054178145456e-01,
                9.964357317888170584e-01,
            ];
            ddouble[] expected_dist_alpha3theta2 = [
                0.000000000000000000e+00,
                2.976680384087791453e-01,
                4.879999999999999893e-01,
                6.153268219383921656e-01,
                7.037037037037037202e-01,
                7.669549385525716856e-01,
                8.134110787172010859e-01,
                8.482962962962963260e-01,
                8.750000000000000000e-01,
                8.957866883777733102e-01,
                9.122085048010973862e-01,
                9.253535500801866043e-01,
                9.360000000000000542e-01,
                9.447143936939854658e-01,
                9.519158527422990623e-01,
                9.579189611243527080e-01,
                9.629629629629630205e-01,
                9.672319999999999807e-01,
                9.708693673190714746e-01,
                9.739877051262511021e-01,
                9.766763848396501357e-01,
                9.790069293533970596e-01,
                9.810370370370370408e-01,
                9.828136014232485840e-01,
                9.843750000000000000e-01,
                9.857528452569774835e-01,
                9.869733360472217054e-01,
                9.880583090379009281e-01,
                9.890260631001371872e-01,
                9.898920103448957120e-01,
                9.906691937600233810e-01,
                9.913687014278730336e-01,
                9.919999999999999929e-01,
                9.925712047126420412e-01,
                9.930892992117481555e-01,
                9.935603154439232032e-01,
                9.939894815927873273e-01,
                9.943813443072702096e-01,
                9.947398701405441024e-01,
                9.950685300944877021e-01,
                9.953703703703703498e-01,
                9.956480718068151381e-01,
                9.959040000000000115e-01,
                9.961402477176952708e-01,
                9.963586709148839482e-01,
                9.965609194166996643e-01,
                9.967484631407813600e-01,
                9.969226145755071267e-01,
                9.970845481049562808e-01,
                9.972353166696364957e-01,
                9.973758661691746186e-01,
                9.975070479455055983e-01,
                9.976296296296296440e-01,
                9.977443045893709117e-01,
                9.978517001779060314e-01,
                9.979523849516290790e-01,
                9.980468750000000000e-01,
                9.981356395084205468e-01,
                9.982191056571222409e-01,
                9.982976629439126448e-01,
                9.983716670059027409e-01,
                9.984414430046056887e-01,
                9.985072886297375883e-01,
                9.985694767693644724e-01,
            ];
            ddouble[] expected_dist_alpha4theta2 = [
                0.000000000000000000e+00,
                3.757049230300259501e-01,
                5.904000000000000359e-01,
                7.202376886824670699e-01,
                8.024691358024691468e-01,
                8.565876544938902937e-01,
                8.933777592669720491e-01,
                9.190913580246913295e-01,
                9.375000000000000000e-01,
                9.509584415895403486e-01,
                9.609815576893766531e-01,
                9.685699158232364825e-01,
                9.744000000000000439e-01,
                9.789388166453277806e-01,
                9.825148555426541641e-01,
                9.853631169128183043e-01,
                9.876543209876543772e-01,
                9.895142399999999894e-01,
                9.910367284058680948e-01,
                9.922926533707410179e-01,
                9.933361099541857531e-01,
                9.942088080974887943e-01,
                9.949432098765431665e-01,
                9.955648003672900037e-01,
                9.960937500000000000e-01,
                9.965461443047217704e-01,
                9.969349025993462510e-01,
                9.972704706372345074e-01,
                9.975613473555859922e-01,
                9.978144887232206495e-01,
                9.980356197389522732e-01,
                9.982294772159739926e-01,
                9.983999999999999542e-01,
                9.985504789683203875e-01,
                9.986836760403330349e-01,
                9.988019191523578311e-01,
                9.989071784714158575e-01,
                9.990011278768480274e-01,
                9.990851948070511579e-01,
                9.991606008671468642e-01,
                9.992283950617284471e-01,
                9.992894811113167686e-01,
                9.993446400000000063e-01,
                9.993945486615992690e-01,
                9.994397955253667698e-01,
                9.994808934968603120e-01,
                9.995182908356713414e-01,
                9.995523803018919740e-01,
                9.995835068721365957e-01,
                9.996119742694227206e-01,
                9.996380505060930011e-01,
                9.996619726027804465e-01,
                9.996839506172839895e-01,
                9.997041710936880321e-01,
                9.997228000229556599e-01,
                9.997399853906830947e-01,
                9.997558593750000000e-01,
                9.997705402471902536e-01,
                9.997841340190450898e-01,
                9.997967358739000554e-01,
                9.998084314124591199e-01,
                9.998192977396643810e-01,
                9.998294044148271498e-01,
                9.998388142838720016e-01,
            ];
            ddouble[] expected_dist_alpha4theta3 = [
                0.000000000000000000e+00,
                2.739750008753195076e-01,
                4.602249062890462206e-01,
                5.904000000000000359e-01,
                6.835937500000000000e-01,
                7.517271105470481052e-01,
                8.024691358024691468e-01,
                8.408851988551346857e-01,
                8.703999999999999515e-01,
                8.933777592669720491e-01,
                9.114814561846869001e-01,
                9.259007793711429013e-01,
                9.375000000000000000e-01,
                9.469158399999999531e-01,
                9.546234375547074658e-01,
                9.609815576893766531e-01,
                9.662640566430653610e-01,
                9.706820909935371944e-01,
                9.744000000000000439e-01,
                9.775468018594054564e-01,
                9.802246093750000000e-01,
                9.825148555426541641e-01,
                9.844829444091904858e-01,
                9.861817576009995756e-01,
                9.876543209876543772e-01,
                9.889358491613047253e-01,
                9.900553249284459456e-01,
                9.910367284058680948e-01,
                9.919000000000000039e-01,
                9.926617997771219892e-01,
                9.933361099541857531e-01,
                9.939347157088114226e-01,
                9.944675910115429174e-01,
                9.949432098765431665e-01,
                9.953687987106963897e-01,
                9.957505418899309513e-01,
                9.960937500000000000e-01,
                9.964029981260411617e-01,
                9.966822400000000526e-01,
                9.969349025993462510e-01,
                9.971639648471691819e-01,
                9.973720233278554614e-01,
                9.975613473555859922e-01,
                9.977339252783279377e-01,
                9.978915035401916267e-01,
                9.980356197389522732e-01,
                9.981676306870961302e-01,
                9.982887363015758853e-01,
                9.983999999999999542e-01,
                9.985023661617954405e-01,
                9.985966751162128618e-01,
                9.986836760403330349e-01,
                9.987640380859375000e-01,
                9.988383600014004715e-01,
                9.989071784714158575e-01,
                9.989709753616188292e-01,
                9.990301840255744192e-01,
                9.990851948070511579e-01,
                9.991363598500624388e-01,
                9.991839973121022789e-01,
                9.992283950617284471e-01,
                9.992698139296726234e-01,
                9.993084905725815315e-01,
                9.993446400000000063e-01,
            ];

            foreach ((LomaxDistribution dist, ddouble[] expecteds) in new[]{
                (dist_alpha2theta1, expected_dist_alpha2theta1),
                (dist_alpha3theta2, expected_dist_alpha3theta2),
                (dist_alpha4theta2, expected_dist_alpha4theta2),
                (dist_alpha4theta3, expected_dist_alpha4theta3),
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