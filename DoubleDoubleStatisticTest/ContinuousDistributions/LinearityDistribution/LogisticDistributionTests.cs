﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.LinearityDistribution {
    [TestClass()]
    public class LogisticDistributionTests {
        readonly LogisticDistribution dist1 = new();
        readonly LogisticDistribution dist2 = new(mu: 1, sigma: 3);

        LogisticDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (LogisticDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (LogisticDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
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
            foreach (LogisticDistribution dist in Dists) {
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

            foreach (LogisticDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 2) * 0.02, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (LogisticDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                (LogisticDistribution? dist_fit, ddouble error) = LogisticDistribution.Fit(xs, (0.1, 0.9));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (LogisticDistribution dist in Dists) {
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
                3.352376707564740750e-04,
                4.303717019243360141e-04,
                5.524730727021605293e-04,
                7.091667670252287115e-04,
                9.102211801218266455e-04,
                1.168142510795444405e-03,
                1.498928708569049854e-03,
                1.923022356864607452e-03,
                2.466509291360047176e-03,
                3.162616926066188182e-03,
                4.053571694869767365e-03,
                5.192875981301849249e-03,
                6.648056670790152001e-03,
                8.503912157689539134e-03,
                1.086622972222523083e-02,
                1.386584143763396754e-02,
                1.766270621329111071e-02,
                2.244941038204345887e-02,
                2.845302387973556307e-02,
                3.593359082532811971e-02,
                4.517665973091213039e-02,
                5.647624464487405876e-02,
                7.010371654510814954e-02,
                8.625794444256298932e-02,
                1.049935854035065202e-01,
                1.261292251866551750e-01,
                1.491464520703328356e-01,
                1.731047869924969840e-01,
                1.966119332414818510e-01,
                2.178949937618140098e-01,
                2.350037122015944946e-01,
                2.461340827375982920e-01,
                2.500000000000000000e-01,
                2.461340827375982920e-01,
                2.350037122015944946e-01,
                2.178949937618140098e-01,
                1.966119332414818510e-01,
                1.731047869924969840e-01,
                1.491464520703328356e-01,
                1.261292251866551750e-01,
                1.049935854035065202e-01,
                8.625794444256298932e-02,
                7.010371654510814954e-02,
                5.647624464487405876e-02,
                4.517665973091213039e-02,
                3.593359082532811971e-02,
                2.845302387973556307e-02,
                2.244941038204345887e-02,
                1.766270621329111071e-02,
                1.386584143763396754e-02,
                1.086622972222523083e-02,
                8.503912157689539134e-03,
                6.648056670790152001e-03,
                5.192875981301849249e-03,
                4.053571694869767365e-03,
                3.162616926066188182e-03,
                2.466509291360047176e-03,
                1.923022356864607452e-03,
                1.498928708569049854e-03,
                1.168142510795444405e-03,
                9.102211801218266455e-04,
                7.091667670252287115e-04,
                5.524730727021605293e-04,
                4.303717019243360141e-04,
            ];
            ddouble[] expected_dist2 = [
                1.505888657697070955e-02,
                1.623347648684589709e-02,
                1.748784661785123679e-02,
                1.882541488162468510e-02,
                2.024939206379837431e-02,
                2.176270436813676490e-02,
                2.336790551503604985e-02,
                2.506707816846890161e-02,
                2.686172475891542358e-02,
                2.875264814752099760e-02,
                3.073982304635384730e-02,
                3.282225967567142105e-02,
                3.499786180116883777e-02,
                3.726328204519072096e-02,
                3.961377819003224221e-02,
                4.204307506221838936e-02,
                4.454323746508398463e-02,
                4.710456046055013640e-02,
                4.971548402344427853e-02,
                5.236253962428102232e-02,
                5.503033655096791626e-02,
                5.770159566416566133e-02,
                6.035723770630962992e-02,
                6.297653217531155867e-02,
                6.553731108049394571e-02,
                6.801624960756402771e-02,
                7.038921286768033692e-02,
                7.263166458727134123e-02,
                7.471912996707620602e-02,
                7.662770121792261691e-02,
                7.833457073386483616e-02,
                7.981857378808095149e-02,
                8.106072033492357776e-02,
                8.204469424586609272e-02,
                8.275729830837390277e-02,
                8.318882469180895189e-02,
                8.333333333333332871e-02,
                8.318882469180895189e-02,
                8.275729830837390277e-02,
                8.204469424586609272e-02,
                8.106072033492357776e-02,
                7.981857378808095149e-02,
                7.833457073386483616e-02,
                7.662770121792261691e-02,
                7.471912996707620602e-02,
                7.263166458727134123e-02,
                7.038921286768033692e-02,
                6.801624960756402771e-02,
                6.553731108049394571e-02,
                6.297653217531155867e-02,
                6.035723770630962992e-02,
                5.770159566416566133e-02,
                5.503033655096791626e-02,
                5.236253962428102232e-02,
                4.971548402344427853e-02,
                4.710456046055013640e-02,
                4.454323746508398463e-02,
                4.204307506221838936e-02,
                3.961377819003224221e-02,
                3.726328204519072096e-02,
                3.499786180116883777e-02,
                3.282225967567142105e-02,
                3.073982304635384730e-02,
                2.875264814752099760e-02,
            ];

            foreach ((LogisticDistribution dist, ddouble[] expecteds) in new[]{
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
                3.353501304664781085e-04,
                4.305570813246148800e-04,
                5.527786369235995539e-04,
                7.096703991005881085e-04,
                9.110511944006453909e-04,
                1.169510265055514770e-03,
                1.501182256736991686e-03,
                1.926734663327475672e-03,
                2.472623156634774343e-03,
                3.172682842485189270e-03,
                4.070137715896127682e-03,
                5.220125693558397284e-03,
                6.692850924284855438e-03,
                8.577485413711984086e-03,
                1.098694263059317965e-02,
                1.406362704324547533e-02,
                1.798620996209155873e-02,
                2.297736991002561485e-02,
                2.931223075135631906e-02,
                3.732688734412945714e-02,
                4.742587317756678106e-02,
                6.008665017400762615e-02,
                7.585818002124354587e-02,
                9.534946489910949008e-02,
                1.192029220221175467e-01,
                1.480471980316894753e-01,
                1.824255238063563489e-01,
                2.227001388253088410e-01,
                2.689414213699951040e-01,
                3.208213008246070252e-01,
                3.775406687981454068e-01,
                4.378234991142019306e-01,
                5.000000000000000000e-01,
                5.621765008857980694e-01,
                6.224593312018545932e-01,
                6.791786991753929748e-01,
                7.310585786300048960e-01,
                7.772998611746910758e-01,
                8.175744761936436511e-01,
                8.519528019683105802e-01,
                8.807970779778823145e-01,
                9.046505351008905516e-01,
                9.241418199787565513e-01,
                9.399133498259923947e-01,
                9.525741268224333647e-01,
                9.626731126558706331e-01,
                9.706877692486436393e-01,
                9.770226300899743643e-01,
                9.820137900379084517e-01,
                9.859363729567544032e-01,
                9.890130573694068117e-01,
                9.914225145862880506e-01,
                9.933071490757152677e-01,
                9.947798743064416582e-01,
                9.959298622841039617e-01,
                9.968273171575148250e-01,
                9.975273768433653432e-01,
                9.980732653366725105e-01,
                9.984988177432629897e-01,
                9.988304897349444822e-01,
                9.990889488055993972e-01,
                9.992903296008994740e-01,
                9.994472213630763990e-01,
                9.995694429186754437e-01,
            ];
            ddouble[] expected_dist2 = [
                4.742587317756678106e-02,
                5.133579311531624723e-02,
                5.554926015761174618e-02,
                6.008665017400762615e-02,
                6.496916912866407268e-02,
                7.021879183055279583e-02,
                7.585818002124354587e-02,
                8.191057715532519545e-02,
                8.839967720705840803e-02,
                9.534946489910949008e-02,
                1.027840249172517761e-01,
                1.107273179723682194e-01,
                1.192029220221175467e-01,
                1.282337375925192702e-01,
                1.378416569649357215e-01,
                1.480471980316894753e-01,
                1.588691048809151574e-01,
                1.703239186438458286e-01,
                1.824255238063563489e-01,
                1.951846770138401799e-01,
                2.086085273260449569e-01,
                2.227001388253088410e-01,
                2.374580283439025052e-01,
                2.528757327293304491e-01,
                2.689414213699951040e-01,
                2.856375705089441164e-01,
                3.029407160345927164e-01,
                3.208213008246070252e-01,
                3.392436312341828297e-01,
                3.581659549112690133e-01,
                3.775406687981454068e-01,
                3.973146620215083358e-01,
                4.174297935376852786e-01,
                4.378234991142019306e-01,
                4.584295167832001527e-01,
                4.791787146272570297e-01,
                5.000000000000000000e-01,
                5.208212853727429703e-01,
                5.415704832167999028e-01,
                5.621765008857980694e-01,
                5.825702064623147214e-01,
                6.026853379784916642e-01,
                6.224593312018545932e-01,
                6.418340450887310977e-01,
                6.607563687658172258e-01,
                6.791786991753929748e-01,
                6.970592839654073947e-01,
                7.143624294910558836e-01,
                7.310585786300048960e-01,
                7.471242672706694954e-01,
                7.625419716560974948e-01,
                7.772998611746910758e-01,
                7.913914726739550431e-01,
                8.048153229861597646e-01,
                8.175744761936436511e-01,
                8.296760813561542269e-01,
                8.411308951190848981e-01,
                8.519528019683105802e-01,
                8.621583430350643340e-01,
                8.717662624074806743e-01,
                8.807970779778823145e-01,
                8.892726820276318778e-01,
                8.972159750827483071e-01,
                9.046505351008905516e-01,
            ];

            foreach ((LogisticDistribution dist, ddouble[] expecteds) in new[]{
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