﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ContinuousDistribution {
    [TestClass()]
    public class ArcsineDistributionTests {
        readonly ArcsineDistribution dist = new();

        ArcsineDistribution[] Dists => [
            dist
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
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
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Mean;
                ddouble expected = IntegrationStatistics.Mean(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void MedianTest() {
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Variance;
                ddouble expected = IntegrationStatistics.Variance(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void SkewnessTest() {
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Skewness;
                ddouble expected = IntegrationStatistics.Skewness(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void KurtosisTest() {
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Kurtosis;
                ddouble expected = IntegrationStatistics.Kurtosis(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (ArcsineDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (ArcsineDistribution dist in Dists) {
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
            foreach (ArcsineDistribution dist in Dists) {
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
            foreach (ArcsineDistribution dist in Dists) {
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

            foreach (ArcsineDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 95; i++) {
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
            foreach (ArcsineDistribution dist in Dists) {
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
            ddouble[] expected_dist = [
                double.PositiveInfinity,
                5.102934600213477445e+00,
                3.615415673813619168e+00,
                2.957802724672829786e+00,
                2.566609672215115712e+00,
                2.300213932746645984e+00,
                2.103993835729919493e+00,
                1.951827349096080244e+00,
                1.829444583876738362e+00,
                1.728304900124160381e+00,
                1.642943161581896794e+00,
                1.569678493153294108e+00,
                1.505929207181997853e+00,
                1.449823983289395768e+00,
                1.399968912117335984e+00,
                1.355301571191928955e+00,
                1.314996147331199516e+00,
                1.278399762992876720e+00,
                1.244988562149981659e+00,
                1.214336695973995628e+00,
                1.186093957176622915e+00,
                1.159969350395513965e+00,
                1.135718822517168602e+00,
                1.113135963022804376e+00,
                1.092044860609388213e+00,
                1.072294549181830314e+00,
                1.053754641593402486e+00,
                1.036311862223869618e+00,
                1.019867267641508146e+00,
                1.004333999624527296e+00,
                9.896354541107929004e-01,
                9.757037780675937855e-01,
                9.624786270806683364e-01,
                9.499061318644860252e-01,
                9.379380334163418542e-01,
                9.265309552372392732e-01,
                9.156457876727548406e-01,
                9.052471645226233266e-01,
                8.953030160151287387e-01,
                8.857841853232177876e-01,
                8.766640982207998256e-01,
                8.679184773937601571e-01,
                8.595250944458958653e-01,
                8.514635538620765054e-01,
                8.437151041754252789e-01,
                8.362624723826609374e-01,
                8.290897183008352211e-01,
                8.221821060896948863e-01,
                8.155259906002850778e-01,
                8.091087165706638551e-01,
                8.029185289882438958e-01,
                7.969444931868687743e-01,
                7.911764234544589325e-01,
                7.856048191012904303e-01,
                7.802208070856204714e-01,
                7.750160904172447296e-01,
                7.699829016645347579e-01,
                7.651139609797362739e-01,
                7.604024381333962523e-01,
                7.558419181138543719e-01,
                7.514263699035427235e-01,
                7.471501180918249663e-01,
                7.430078170254862391e-01,
                7.389944272337537479e-01,
                7.351051938957228193e-01,
                7.313356271449891199e-01,
                7.276814840297235465e-01,
                7.241387519668823769e-01,
                7.207036335471256328e-01,
                7.173725325626933991e-01,
                7.141420411442430671e-01,
                7.110089279047611122e-01,
                7.079701269993309287e-01,
                7.050227280189617085e-01,
                7.021639666450184558e-01,
                6.993912159981809928e-01,
                6.967019786224182276e-01,
                6.940938790502931832e-01,
                6.915646569011110190e-01,
                6.891121604680512380e-01,
                6.867343407545661860e-01,
                6.844292459240308713e-01,
                6.821950161299463744e-01,
                6.800298786969765086e-01,
                6.779321436257732536e-01,
                6.759001993969523392e-01,
                6.739325090517510608e-01,
                6.720276065288542933e-01,
                6.701840932386439187e-01,
                6.684006348577258283e-01,
                6.666759583280322810e-01,
                6.650088490461163460e-01,
                6.633981482294382115e-01,
                6.618427504475351553e-01,
                6.603416013069427493e-01,
                6.588936952796391466e-01,
                6.574980736655997582e-01,
                6.561538226807984708e-01,
                6.548600716626700180e-01,
                6.536159913856818271e-01,
                6.524207924802281067e-01,
                6.512737239485872909e-01,
                6.501740717721702367e-01,
                6.491211576047225540e-01,
                6.481143375465601153e-01,
                6.471530009952893847e-01,
                6.462365695688120359e-01,
                6.453644960967365174e-01,
                6.445362636766140962e-01,
                6.437513847916931375e-01,
                6.430094004871429458e-01,
                6.423098796019330869e-01,
                6.416524180537789279e-01,
                6.410366381747685249e-01,
                6.404621880954799451e-01,
                6.399287411755799626e-01,
                6.394359954790636014e-01,
                6.389836732924559781e-01,
                6.385715206844464475e-01,
                6.381993071055698241e-01,
                6.378668250266846806e-01,
                6.375738896151271762e-01,
                6.373203384475422117e-01,
                6.371060312585125063e-01,
                6.369308497242192058e-01,
                6.367946972804769956e-01,
                6.366974989745943780e-01,
                6.366392013506133596e-01,
                6.366197723675813824e-01,
                6.366392013506133596e-01,
                6.366974989745943780e-01,
                6.367946972804769956e-01,
                6.369308497242192058e-01,
                6.371060312585125063e-01,
                6.373203384475422117e-01,
                6.375738896151271762e-01,
                6.378668250266846806e-01,
                6.381993071055698241e-01,
                6.385715206844464475e-01,
                6.389836732924559781e-01,
                6.394359954790636014e-01,
                6.399287411755799626e-01,
                6.404621880954799451e-01,
                6.410366381747685249e-01,
                6.416524180537789279e-01,
                6.423098796019330869e-01,
                6.430094004871429458e-01,
                6.437513847916931375e-01,
                6.445362636766140962e-01,
                6.453644960967365174e-01,
                6.462365695688120359e-01,
                6.471530009952893847e-01,
                6.481143375465601153e-01,
                6.491211576047225540e-01,
                6.501740717721702367e-01,
                6.512737239485872909e-01,
                6.524207924802281067e-01,
                6.536159913856818271e-01,
                6.548600716626700180e-01,
                6.561538226807984708e-01,
                6.574980736655997582e-01,
                6.588936952796391466e-01,
                6.603416013069427493e-01,
                6.618427504475351553e-01,
                6.633981482294382115e-01,
                6.650088490461163460e-01,
                6.666759583280322810e-01,
                6.684006348577258283e-01,
                6.701840932386439187e-01,
                6.720276065288542933e-01,
                6.739325090517510608e-01,
                6.759001993969523392e-01,
                6.779321436257732536e-01,
                6.800298786969765086e-01,
                6.821950161299463744e-01,
                6.844292459240308713e-01,
                6.867343407545661860e-01,
                6.891121604680512380e-01,
                6.915646569011110190e-01,
                6.940938790502931832e-01,
                6.967019786224182276e-01,
                6.993912159981809928e-01,
                7.021639666450184558e-01,
                7.050227280189617085e-01,
                7.079701269993309287e-01,
                7.110089279047611122e-01,
                7.141420411442430671e-01,
                7.173725325626933991e-01,
                7.207036335471256328e-01,
                7.241387519668823769e-01,
                7.276814840297235465e-01,
                7.313356271449891199e-01,
                7.351051938957228193e-01,
                7.389944272337537479e-01,
                7.430078170254862391e-01,
                7.471501180918249663e-01,
                7.514263699035427235e-01,
                7.558419181138543719e-01,
                7.604024381333962523e-01,
                7.651139609797362739e-01,
                7.699829016645347579e-01,
                7.750160904172447296e-01,
                7.802208070856204714e-01,
                7.856048191012904303e-01,
                7.911764234544589325e-01,
                7.969444931868687743e-01,
                8.029185289882438958e-01,
                8.091087165706638551e-01,
                8.155259906002850778e-01,
                8.221821060896948863e-01,
                8.290897183008352211e-01,
                8.362624723826609374e-01,
                8.437151041754252789e-01,
                8.514635538620765054e-01,
                8.595250944458958653e-01,
                8.679184773937601571e-01,
                8.766640982207998256e-01,
                8.857841853232177876e-01,
                8.953030160151287387e-01,
                9.052471645226233266e-01,
                9.156457876727548406e-01,
                9.265309552372392732e-01,
                9.379380334163418542e-01,
                9.499061318644860252e-01,
                9.624786270806683364e-01,
                9.757037780675937855e-01,
                9.896354541107929004e-01,
                1.004333999624527296e+00,
                1.019867267641508146e+00,
                1.036311862223869618e+00,
                1.053754641593402486e+00,
                1.072294549181830314e+00,
                1.092044860609388213e+00,
                1.113135963022804376e+00,
                1.135718822517168602e+00,
                1.159969350395513965e+00,
                1.186093957176622915e+00,
                1.214336695973995628e+00,
                1.244988562149981659e+00,
                1.278399762992876720e+00,
                1.314996147331199516e+00,
                1.355301571191928955e+00,
                1.399968912117335984e+00,
                1.449823983289395768e+00,
                1.505929207181997853e+00,
                1.569678493153294108e+00,
                1.642943161581896794e+00,
                1.728304900124160381e+00,
                1.829444583876738362e+00,
                1.951827349096080244e+00,
                2.103993835729919493e+00,
                2.300213932746645984e+00,
                2.566609672215115712e+00,
                2.957802724672829786e+00,
                3.615415673813619168e+00,
                5.102934600213477445e+00,
                double.PositiveInfinity,
            ];

            foreach ((ArcsineDistribution dist, ddouble[] expecteds) in new[]{
                (dist, expected_dist)
            }) {
                for ((ddouble x, int i) = (1d / 256, 1); i < expecteds.Length - 1; x += 1d / 256, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-14, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist = [
                0.000000000000000000e+00,
                3.981468553857747672e-02,
                5.634329647599925495e-02,
                6.905142851401453730e-02,
                7.978617534953647006e-02,
                8.926251030067318404e-02,
                9.784688372416877611e-02,
                1.057568520356759606e-01,
                1.131340822573211613e-01,
                1.200769208317815917e-01,
                1.266569338873447526e-01,
                1.329281193826306329e-01,
                1.389324068593618122e-01,
                1.447031242248240090e-01,
                1.502672808130794402e-01,
                1.556471260933865175e-01,
                1.608612465103324840e-01,
                1.659253574377640983e-01,
                1.708528878228005587e-01,
                1.756554202101504591e-01,
                1.803430275957909334e-01,
                1.849245352126235309e-01,
                1.894077267298194811e-01,
                1.937995086410816414e-01,
                1.981060427563264281e-01,
                2.023328540475890225e-01,
                2.064849192292999847e-01,
                2.105667401180692511e-01,
                2.145824048502942560e-01,
                2.185356393262696817e-01,
                2.224298507221364463e-01,
                2.262681645146422948e-01,
                2.300534561626159102e-01,
                2.337883783581126973e-01,
                2.374753845814266962e-01,
                2.411167495546006345e-01,
                2.447145870782462540e-01,
                2.482708656494191091e-01,
                2.517874221887800301e-01,
                2.552659741494138212e-01,
                2.587081302345012879e-01,
                2.621153999142954993e-01,
                2.654892019027996763e-01,
                2.688308717298289463e-01,
                2.721416685237193489e-01,
                2.754227811029890671e-01,
                2.786753334611161703e-01,
                2.819003897167521866e-01,
                2.850989585917253488e-01,
                2.882719974707778632e-01,
                2.914204160898503093e-01,
                2.945450798936618542e-01,
                2.976468130981602678e-01,
                3.007264014889837744e-01,
                3.037845949832692627e-01,
                3.068221099788626671e-01,
                3.098396315121523248e-01,
                3.128378152432922388e-01,
                3.158172892854451685e-01,
                3.187786558928216718e-01,
                3.217224930206649680e-01,
                3.246493557689114628e-01,
                3.275597777200117688e-01,
                3.304542721803001015e-01,
                3.333333333333333703e-01,
                3.361974373127702420e-01,
                3.390470432016039481e-01,
                3.418825939638941436e-01,
                3.447045173145492636e-01,
                3.475132265321806968e-01,
                3.503091212195793003e-01,
                3.530925880159433428e-01,
                3.558640012646110407e-01,
                3.586237236397116779e-01,
                3.613721067348472649e-01,
                3.641094916166458528e-01,
                3.668362093457795958e-01,
                3.695525814678227716e-01,
                3.722589204761234116e-01,
                3.749555302486831110e-01,
                3.776427064608769424e-01,
                3.803207369756953504e-01,
                3.829899022130595077e-01,
                3.856504754996356032e-01,
                3.883027234004666206e-01,
                3.909469060336364654e-01,
                3.935832773690903785e-01,
                3.962120855126529673e-01,
                3.988335729762055304e-01,
                4.014479769349181826e-01,
                4.040555294723641166e-01,
                4.066564578142878861e-01,
                4.092509845517444678e-01,
                4.118393278542758495e-01,
                4.144217016737474224e-01,
                4.169983159394248795e-01,
                4.195693767448336886e-01,
                4.221350865269065111e-01,
                4.246956442378945207e-01,
                4.272512455104859885e-01,
                4.298020828165493579e-01,
                4.323483456198907726e-01,
                4.348902205233956542e-01,
                4.374278914108987171e-01,
                4.399615395841091625e-01,
                4.424913438948993583e-01,
                4.450174808732452192e-01,
                4.475401248510954555e-01,
                4.500594480824288701e-01,
                4.525756208597452868e-01,
                4.550888116272260864e-01,
                4.575991870907865633e-01,
                4.601069123252318205e-01,
                4.626121508787189862e-01,
                4.651150648747207628e-01,
                4.676158151116733275e-01,
                4.701145611604884977e-01,
                4.726114614601001351e-01,
                4.751066734112099343e-01,
                4.776003534683925245e-01,
                4.800926572307113172e-01,
                4.825837395309975353e-01,
                4.850737545239355852e-01,
                4.875628557730972146e-01,
                4.900511963370613455e-01,
                4.925389288547581401e-01,
                4.950262056301659719e-01,
                4.975131787164974173e-01,
                5.000000000000001110e-01,
                5.024868212835026382e-01,
                5.049737943698340281e-01,
                5.074610711452419709e-01,
                5.099488036629385990e-01,
                5.124371442269028964e-01,
                5.149262454760644703e-01,
                5.174162604690024647e-01,
                5.199073427692887384e-01,
                5.223996465316075311e-01,
                5.248933265887900657e-01,
                5.273885385398999759e-01,
                5.298854388395115578e-01,
                5.323841848883267280e-01,
                5.348849351252792372e-01,
                5.373878491212811248e-01,
                5.398930876747682905e-01,
                5.424008129092134922e-01,
                5.449111883727739691e-01,
                5.474243791402547687e-01,
                5.499405519175711854e-01,
                5.524598751489044890e-01,
                5.549825191267548918e-01,
                5.575086561051006973e-01,
                5.600384604158907820e-01,
                5.625721085891013384e-01,
                5.651097794766043458e-01,
                5.676516543801092274e-01,
                5.701979171834506976e-01,
                5.727487544895140115e-01,
                5.753043557621055903e-01,
                5.778649134730935444e-01,
                5.804306232551662559e-01,
                5.830016840605750650e-01,
                5.855782983262527441e-01,
                5.881606721457242060e-01,
                5.907490154482555322e-01,
                5.933435421857121694e-01,
                5.959444705276359944e-01,
                5.985520230650819284e-01,
                6.011664270237945251e-01,
                6.037879144873471438e-01,
                6.064167226309097325e-01,
                6.090530939663635346e-01,
                6.116972765995333239e-01,
                6.143495245003643968e-01,
                6.170100977869406034e-01,
                6.196792630243046496e-01,
                6.223572935391231686e-01,
                6.250444697513168890e-01,
                6.277410795238768104e-01,
                6.304474185321773394e-01,
                6.331637906542205707e-01,
                6.358905083833541472e-01,
                6.386278932651527906e-01,
                6.413762763602883776e-01,
                6.441359987353889593e-01,
                6.469074119840568793e-01,
                6.496908787804207552e-01,
                6.524867734678193587e-01,
                6.552954826854506809e-01,
                6.581174060361059119e-01,
                6.609529567983960519e-01,
                6.638025626872298135e-01,
                6.666666666666666297e-01,
                6.695457278197000095e-01,
                6.724402222799882312e-01,
                6.753506442310884816e-01,
                6.782775069793350875e-01,
                6.812213441071783837e-01,
                6.841827107145548315e-01,
                6.871621847567077612e-01,
                6.901603684878477862e-01,
                6.931778900211373884e-01,
                6.962154050167307373e-01,
                6.992735985110162256e-01,
                7.023531869018397877e-01,
                7.054549201063382569e-01,
                7.085795839101497462e-01,
                7.117280025292223034e-01,
                7.149010414082747067e-01,
                7.180996102832478689e-01,
                7.213246665388839407e-01,
                7.245772188970109884e-01,
                7.278583314762806511e-01,
                7.311691282701711092e-01,
                7.345107980972003237e-01,
                7.378846000857044451e-01,
                7.412918697654988787e-01,
                7.447340258505861232e-01,
                7.482125778112199699e-01,
                7.517291343505809742e-01,
                7.552854129217538848e-01,
                7.588832504453993932e-01,
                7.625246154185733038e-01,
                7.662116216418872749e-01,
                7.699465438373841453e-01,
                7.737318354853577329e-01,
                7.775701492778636092e-01,
                7.814643606737303738e-01,
                7.854175951497057717e-01,
                7.894332598819309155e-01,
                7.935150807707002096e-01,
                7.976671459524110608e-01,
                8.018939572436737384e-01,
                8.062004913589183586e-01,
                8.105922732701806854e-01,
                8.150754647873764691e-01,
                8.196569724042089833e-01,
                8.243445797898495409e-01,
                8.291471121771993857e-01,
                8.340746425622360682e-01,
                8.391387534896677103e-01,
                8.443528739066136213e-01,
                8.497327191869205043e-01,
                8.552968757751762130e-01,
                8.610675931406380768e-01,
                8.670718806173695059e-01,
                8.733430661126554417e-01,
                8.799230791682186581e-01,
                8.868659177426788665e-01,
                8.942431479643239145e-01,
                9.021531162758312794e-01,
                9.107374896993267743e-01,
                9.202138246504637520e-01,
                9.309485714859854211e-01,
                9.436567035240009949e-01,
                9.601853144614221902e-01,
                1.000000000000000000e+00,
            ];

            foreach ((ArcsineDistribution dist, ddouble[] expecteds) in new[]{
                (dist, expected_dist),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 256, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-14, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                    else {
                        Assert.AreEqual(0, actual);
                    }
                }
            }
        }
    }
}