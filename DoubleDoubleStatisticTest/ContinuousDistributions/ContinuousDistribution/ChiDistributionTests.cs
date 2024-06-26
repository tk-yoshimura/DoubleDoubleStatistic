﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class ChiDistributionTests {
        readonly ChiDistribution dist_nu1 = new(nu: 1);
        readonly ChiDistribution dist_nu2 = new(nu: 2);
        readonly ChiDistribution dist_nu3 = new(nu: 3);
        readonly ChiDistribution dist_nu4 = new(nu: 4);
        readonly ChiDistribution dist_nu5 = new(nu: 5);
        readonly ChiDistribution dist_nu8 = new(nu: 8);
        readonly ChiDistribution dist_nu16 = new(nu: 16);
        readonly ChiDistribution dist_nu32 = new(nu: 32);
        readonly ChiDistribution dist_nu64 = new(nu: 64);
        readonly ChiDistribution dist_nu128 = new(nu: 128);


        ChiDistribution[] Dists => [
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
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (ChiDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.Nu}");
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
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-28, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
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
            foreach (ChiDistribution dist in Dists) {
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

            foreach (ChiDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 1) * 0.02, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (ChiDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                (ChiDistribution? dist_fit, ddouble error) = ChiDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (ChiDistribution dist in Dists) {
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
            ddouble[] expected_nu1 = [
                7.978845608028654057e-01,
                7.041306535285990487e-01,
                4.839414490382867862e-01,
                2.590351913317834875e-01,
                1.079819330263761257e-01,
                3.505660098713708067e-02,
                8.863696823876011585e-03,
                1.745365390091519657e-03,
                2.676604515297708437e-04,
                3.196748221381096955e-05,
                2.973439029468597050e-06,
                2.153952008508656539e-07,
                1.215176569964657879e-08,
                5.339113229525706938e-10,
                1.826944081672919698e-11,
                4.868641066058022237e-13,
                1.010454216707379020e-14,
                1.633247126333911025e-16,
                2.055954714333784075e-18,
                2.015587078860003177e-20,
                1.538919725341284693e-22,
                9.150751181041616605e-25,
                4.237638507018709669e-27,
                1.528331082317430835e-29,
                4.292767471326093025e-32,
                9.390390715950230704e-35,
                1.599765551401351885e-37,
                2.122537627830418232e-40,
                2.193213118777928199e-43,
                1.764950994918953285e-46,
                1.106141909968876060e-49,
                5.399026049177139191e-53,
                2.052326145583793601e-56,
                6.075803397579806523e-60,
                1.400836426863707071e-63,
                2.515344765756370985e-67,
                3.517499085190184561e-71,
                3.830864983343845943e-75,
                3.249272073547194703e-79,
                2.146355668136516533e-83,
                1.104189672431945563e-87,
                4.423968760421112081e-92,
                1.380405884025434955e-96,
                3.354501045708761354e-101,
                6.348563105650482805e-106,
                9.357273634500031866e-111,
                1.074112073004141771e-115,
                9.602321571813973107e-121,
                6.685428883589061166e-126,
                3.625005587103129117e-131,
                1.530785947283912051e-136,
                5.034387410386920289e-142,
                1.289451994279598679e-147,
                2.572113348142795001e-153,
                3.995778518336646233e-159,
                4.834361157239206160e-165,
                4.555154957473421706e-171,
                3.342670730260233052e-177,
                1.910338908389809412e-183,
                8.502636707055287153e-190,
                2.947292269757159835e-196,
                7.956446854564414275e-203,
                1.672790321171329046e-209,
                2.738988616159594583e-216,
                3.492732513517630450e-223,
                3.468700527212947009e-230,
                2.682839334698930853e-237,
                1.616028082774524761e-244,
                7.581052800185902807e-252,
                2.769724156675932308e-259,
                7.880792554272220875e-267,
                1.746343420124140599e-274,
                3.013809435240855385e-282,
            ];
            ddouble[] expected_nu2 = [
                0.000000000000000000e+00,
                4.412484512922977276e-01,
                6.065306597126334243e-01,
                4.869787010375246084e-01,
                2.706705664732254046e-01,
                1.098423340585185509e-01,
                3.332698961472693039e-02,
                7.656218913640098697e-03,
                1.341850511610047630e-03,
                1.802938382682798938e-04,
                1.863326586039334114e-05,
                1.484768176849658036e-06,
                9.137987846827572860e-08,
                4.349530959340313648e-09,
                1.602814391951884374e-10,
                4.576452508203986385e-12,
                1.013133243927532528e-13,
                1.739925956189575219e-15,
                2.319081398239488409e-17,
                2.399855591879407124e-19,
                1.928749847963913166e-21,
                1.204220411336873109e-23,
                5.842201474647015780e-26,
                2.202800794762216041e-28,
                6.456223392025375090e-31,
                1.471138679902095907e-33,
                2.606511416550167373e-36,
                3.591278661524400795e-39,
                3.848299011074287901e-42,
                3.207455148721440257e-45,
                2.079514940461756944e-48,
                1.048834729651103570e-51,
                4.115534996227836169e-55,
                1.256456898466489112e-58,
                2.984669766353212839e-62,
                5.516905021503770788e-66,
                7.935356396633784607e-70,
                8.882362897277790443e-74,
                7.737481388946301390e-78,
                5.245612910036339328e-82,
                2.767793053473453870e-86,
                1.136647631047969725e-90,
                3.633172640333462948e-95,
                9.039123706086775562e-100,
                1.750483656229302195e-104,
                2.638710749890953298e-109,
                3.096259646161858263e-114,
                2.828160463595815623e-119,
                2.010946208116654937e-124,
                1.113101333790180226e-129,
                4.796389172336852996e-135,
                1.608965573111042324e-140,
                4.201829875932749026e-146,
                8.542714958313600422e-152,
                1.352150740785444300e-157,
                1.666217625395605313e-163,
                1.598531229641923210e-169,
                1.193983697548320211e-175,
                6.943338806751297535e-182,
                3.143660062875926378e-188,
                1.108164920546163607e-194,
                3.041437834315525263e-201,
                6.499248450694511811e-208,
                1.081336143692345248e-214,
                1.400797131856989401e-221,
                1.412895707882557016e-228,
                1.109605353887011686e-235,
                6.785059322174537774e-243,
                3.230489821070749636e-250,
                1.197610382499015242e-257,
                3.456988052532051726e-265,
                7.769944984525434308e-273,
                1.359809990050363609e-280,
            ];
            ddouble[] expected_nu4 = [
                0.000000000000000000e+00,
                5.515605641153721594e-02,
                3.032653298563167121e-01,
                5.478510386672151844e-01,
                5.413411329464505872e-01,
                3.432572939328704664e-01,
                1.499714532662711208e-01,
                4.689434084604560127e-02,
                1.073480409288037410e-02,
                1.825475112466332485e-03,
                2.329158232549168083e-04,
                2.245711867485107612e-05,
                1.644837812428964226e-06,
                9.188384151606385894e-08,
                3.926895260282124562e-09,
                1.287127267932374345e-10,
                3.242026380568111361e-12,
                6.285482516734839779e-14,
                9.392279662869889490e-16,
                1.082934835835587529e-17,
                9.643749239819569213e-20,
                6.638265017494517876e-22,
                3.534531892161426696e-24,
                1.456602025536517466e-26,
                4.648480842258273848e-29,
                1.149327093673502543e-31,
                2.202502146984885647e-34,
                3.272552680314092122e-37,
                3.771333030852823720e-40,
                3.371837225093399001e-43,
                2.339454308019474091e-46,
                1.259912718993387743e-49,
                5.267884795171623806e-53,
                1.710351953037521555e-56,
                4.312847812380408732e-60,
                8.447760814177801473e-64,
                1.285527736254663483e-67,
                1.519994350796657164e-71,
                1.396615390704789194e-75,
                9.973221545206620472e-80,
                5.535586106946969318e-84,
                2.388380834739571891e-88,
                8.011145671935149646e-93,
                2.089167466569321894e-97,
                4.236170448074875103e-102,
                6.679236585661412932e-107,
                8.189606764097983396e-112,
                7.809258080103969493e-117,
                5.791525079375992955e-122,
                3.340695378037690311e-127,
                1.498871616355280434e-132,
                5.231149319577251667e-138,
                1.420218498065250292e-143,
                2.999560789737902634e-149,
                4.928589450162865634e-155,
                6.300385396027246236e-161,
                6.266242420196487392e-167,
                4.849066291668150192e-173,
                2.919673968239001439e-179,
                1.367885084858897422e-185,
                4.986742142457891749e-192,
                1.414648772686048817e-198,
                3.122888880558635133e-205,
                5.364778942893664183e-212,
                7.172081315107906507e-219,
                7.461855457255276766e-226,
                6.041801151914676278e-233,
                3.807266412155075870e-240,
                1.867223116578800829e-247,
                7.127278788846886562e-255,
                2.117405182175867496e-262,
                4.896036583374261906e-270,
                8.811568735526823467e-278,
            ];
            ddouble[] expected_nu8 = [
                0.000000000000000000e+00,
                1.436355635717114366e-04,
                1.263605541067986532e-02,
                1.155623284688657021e-01,
                3.608940886309671692e-01,
                5.586870018438649943e-01,
                5.061536547736653624e-01,
                2.932117509670717492e-01,
                1.145045769907239996e-01,
                3.118995367940525024e-02,
                6.065516230596796829e-03,
                8.562361315585801302e-04,
                8.882124187116408390e-05,
                6.834099993594543231e-06,
                3.928531466640572118e-07,
                1.696896300496779470e-08,
                5.533058356169564660e-10,
                1.367108815833882008e-11,
                2.567614452837063344e-13,
                3.675238300545026021e-15,
                4.018228849924824886e-17,
                3.362021924133733300e-19,
                2.156211726380652772e-21,
                1.061502519344182240e-23,
                4.016287447711129574e-26,
                1.169155979078690417e-28,
                2.621069325834791461e-31,
                4.529085075465645908e-34,
                6.036647071385073368e-37,
                6.210511469795047373e-40,
                4.934786430978571357e-43,
                3.030093370201803520e-46,
                1.438483741401549196e-49,
                5.282141402339099045e-53,
                1.500889008907612116e-56,
                3.301281529108153248e-60,
                5.622898318377975703e-64,
                7.418526386675514507e-68,
                7.583679763834975598e-72,
                6.008453027282422214e-76,
                3.690390737964626064e-80,
                1.757550368738146594e-84,
                6.491731755931902199e-89,
                1.860012454134003672e-93,
                4.134784768684331643e-98,
                7.132563678145138110e-103,
                9.549115610299732999e-108,
                9.923616738953141397e-113,
                8.006204269729450203e-118,
                5.015219806251861234e-123,
                2.439569688078302569e-128,
                9.216087658321099218e-134,
                2.704190701549370523e-139,
                6.163535786399103440e-145,
                1.091356044160029174e-150,
                1.501366252982308513e-156,
                1.604826458761773770e-162,
                1.332988245414479893e-168,
                8.604291349708530342e-175,
                4.316447234310194312e-181,
                1.683025473079548548e-187,
                5.100781764962547650e-194,
                1.201688942442672005e-200,
                2.200811288047543687e-207,
                3.133530140446007226e-214,
                3.468718322178737744e-221,
                2.985457859949818466e-228,
                1.997934535171277521e-235,
                1.039682279465202693e-242,
                4.207158612175187036e-250,
                1.323929646459425336e-257,
                3.240013537006024823e-265,
                6.166688263871047712e-273,
                9.128450313405448755e-281,
                1.050988510522451479e-288,
            ];

            foreach ((ChiDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4), (dist_nu8, expected_nu8)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
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
            ddouble[] expected_nu1 = [
                0.000000000000000000e+00,
                3.829249225480261254e-01,
                6.826894921370858516e-01,
                8.663855974622841627e-01,
                9.544997361036414718e-01,
                9.875806693484476817e-01,
                9.973002039367397931e-01,
                9.995347418419289198e-01,
                9.999366575163337600e-01,
                9.999932046537505226e-01,
                9.999994266968562640e-01,
                9.999999620208750439e-01,
                9.999999980268247102e-01,
                9.999999999196800271e-01,
                9.999999999974403808e-01,
                9.999999999999361622e-01,
                9.999999999999987788e-01,
            ];
            ddouble[] expected_nu2 = [
                0.000000000000000000e+00,
                1.175030974154046004e-01,
                3.934693402873665202e-01,
                6.753475326416502611e-01,
                8.646647167633872977e-01,
                9.560630663765925519e-01,
                9.888910034617577338e-01,
                9.978125088818171617e-01,
                9.996645373720974836e-01,
                9.999599347026070228e-01,
                9.999962733468279463e-01,
                9.999997300421497037e-01,
                9.999999847700202782e-01,
                9.999999993308413826e-01,
                9.999999999771026493e-01,
                9.999999999993898214e-01,
                9.999999999999873435e-01,
                9.999999999999997780e-01,
            ];
            ddouble[] expected_nu4 = [
                0.000000000000000000e+00,
                7.190984592330174584e-03,
                9.020401043104986361e-02,
                3.101135068635068048e-01,
                5.939941502901615600e-01,
                8.187601488034443875e-01,
                9.389005190396673139e-01,
                9.844141257829469582e-01,
                9.969808363488773528e-01,
                9.995542735665033929e-01,
                9.999496901821769423e-01,
                9.999956469296633621e-01,
                9.999997106303848415e-01,
                9.999999851948657703e-01,
                9.999999994161176131e-01,
                9.999999999822281049e-01,
                9.999999999995821121e-01,
                9.999999999999924505e-01,
                9.999999999999998890e-01,
            ];
            ddouble[] expected_nu8 = [
                0.000000000000000000e+00,
                9.206413744597217042e-06,
                1.751622556290823940e-03,
                2.762781505072401181e-02,
                1.428765395014529593e-01,
                3.807495393578545495e-01,
                6.577040441654089520e-01,
                8.596067915520650971e-01,
                9.576198880083159892e-01,
                9.905695132473332221e-01,
                9.984454421569889382e-01,
                9.998090890948383436e-01,
                9.999824398333543307e-01,
                9.999987844836243367e-01,
                9.999999364221727749e-01,
                9.999999974783633228e-01,
                9.999999999239348458e-01,
                9.999999999982505106e-01,
                9.999999999999692468e-01,
                9.999999999999995559e-01,
            ];

            foreach ((ChiDistribution dist, ddouble[] expecteds) in new[]{
                (dist_nu1, expected_nu1), (dist_nu2, expected_nu2), (dist_nu4, expected_nu4), (dist_nu8, expected_nu8)
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.5, i++) {
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