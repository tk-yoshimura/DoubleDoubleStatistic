using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.InternalUtils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistribution {
    [TestClass()]
    public class IrwinHallDistributionTests {
        readonly IrwinHallDistribution dist_n1 = new(n: 1);
        readonly IrwinHallDistribution dist_n2 = new(n: 2);
        readonly IrwinHallDistribution dist_n3 = new(n: 3);
        readonly IrwinHallDistribution dist_n4 = new(n: 4);
        readonly IrwinHallDistribution dist_n5 = new(n: 5);
        readonly IrwinHallDistribution dist_n8 = new(n: 8);
        readonly IrwinHallDistribution dist_n16 = new(n: 16);
        readonly IrwinHallDistribution dist_n32 = new(n: 32);
        readonly IrwinHallDistribution dist_n64 = new(n: 64);


        IrwinHallDistribution[] Dists => [
            dist_n1,
            dist_n2,
            dist_n3,
            dist_n4,
            dist_n5,
            dist_n8,
            dist_n16,
            dist_n32,
            dist_n64,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"N={dist.N}");
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
            foreach (IrwinHallDistribution dist in Dists) {
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
            foreach (IrwinHallDistribution dist in Dists) {
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
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (IrwinHallDistribution dist in Dists) {
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
            foreach (IrwinHallDistribution dist in Dists) {
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
            foreach (IrwinHallDistribution dist in Dists) {
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
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = 0; x <= dist.N; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);

                for (ddouble x = 0; x <= dist.N; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);

                for (ddouble x = 0; x <= dist.N; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"ccdf({x})={ccdf}");

                    Assert.IsTrue(ddouble.Abs(cdf + ccdf - 1) < 1e-30);
                }
            }
        }

        [TestMethod()]
        public void QuantileLowerTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 100; i++) {
                    ddouble p = (ddouble)i / 100;
                    ddouble x = dist.Quantile(p, Interval.Lower);
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-30);
                    }

                    Assert.IsTrue(ddouble.IsFinite(x));
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (IrwinHallDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int i = 0; i <= 100; i++) {
                    ddouble p = (ddouble)i / 100;
                    ddouble x = dist.Quantile(p, Interval.Upper);
                    ddouble ccdf = dist.CDF(x, Interval.Upper);

                    Console.WriteLine($"cquantile({p})={x}, ccdf({x})={ccdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - ccdf) < 1e-30);
                    }

                    Assert.IsTrue(ddouble.IsFinite(x));
                }
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (IrwinHallDistribution dist in Dists) {
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
            ddouble[] expected_n1 = [
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ];
            ddouble[] expected_n2 = [
                0,
                0.125,
                0.25,
                0.375,
                0.5,
                0.625,
                0.75,
                0.875,
                1,
                0.875,
                0.75,
                0.625,
                0.5,
                0.375,
                0.25,
                0.125,
                0,
            ];
            ddouble[] expected_n3 = [
                0,
                0.0078125,
                0.03125,
                0.0703125,
                0.125,
                0.1953125,
                0.28125,
                0.3828125,
                0.5,
                0.609375,
                0.6875,
                0.734375,
                0.75,
                0.734375,
                0.6875,
                0.609375,
                0.5,
                0.3828125,
                0.28125,
                0.1953125,
                0.125,
                0.0703125,
                0.03125,
                0.0078125,
                0,
            ];
            ddouble[] expected_n4 = [
                "0",
                "3.2552083333333333333333333333333333e-4",
                "0.0026041666666666666666666666666666667",
                "0.0087890625",
                "0.020833333333333333333333333333333333",
                "0.040690104166666666666666666666666667",
                "0.0703125",
                "0.11165364583333333333333333333333333",
                "0.16666666666666666666666666666666667",
                "0.23600260416666666666666666666666667",
                "0.31510416666666666666666666666666667",
                "0.39811197916666666666666666666666667",
                "0.47916666666666666666666666666666667",
                "0.55240885416666666666666666666666667",
                "0.61197916666666666666666666666666667",
                "0.65201822916666666666666666666666667",
                "0.66666666666666666666666666666666667",
                "0.65201822916666666666666666666666667",
                "0.61197916666666666666666666666666667",
                "0.55240885416666666666666666666666667",
                "0.47916666666666666666666666666666667",
                "0.39811197916666666666666666666666667",
                "0.31510416666666666666666666666666667",
                "0.23600260416666666666666666666666667",
                "0.16666666666666666666666666666666667",
                "0.11165364583333333333333333333333333",
                "0.0703125",
                "0.040690104166666666666666666666666667",
                "0.020833333333333333333333333333333333",
                "0.0087890625",
                "0.0026041666666666666666666666666666667",
                "3.2552083333333333333333333333333333e-4",
                "0",
            ];
            ddouble[] expected_n5 = [
                "0",
                "1.0172526041666666666666666666666667e-5",
                "1.6276041666666666666666666666666667e-4",
                "8.23974609375e-4",
                "0.0026041666666666666666666666666666667",
                "0.0063578287760416666666666666666666667",
                "0.01318359375",
                "0.024424235026041666666666666666666667",
                "0.041666666666666666666666666666666667",
                "0.066691080729166666666666666666666667",
                "0.10091145833333333333333333333333333",
                "0.14481608072916666666666666666666667",
                "0.19791666666666666666666666666666667",
                "0.25874837239583333333333333333333333",
                "0.32486979166666666666666666666666667",
                "0.39286295572916666666666666666666667",
                "0.45833333333333333333333333333333333",
                "0.51601155598958333333333333333333333",
                "0.56087239583333333333333333333333333",
                "0.58925374348958333333333333333333333",
                "0.59895833333333333333333333333333333",
                "0.58925374348958333333333333333333333",
                "0.56087239583333333333333333333333333",
                "0.51601155598958333333333333333333333",
                "0.45833333333333333333333333333333333",
                "0.39286295572916666666666666666666667",
                "0.32486979166666666666666666666666667",
                "0.25874837239583333333333333333333333",
                "0.19791666666666666666666666666666667",
                "0.14481608072916666666666666666666667",
                "0.10091145833333333333333333333333333",
                "0.066691080729166666666666666666666667",
                "0.041666666666666666666666666666666667",
                "0.024424235026041666666666666666666667",
                "0.01318359375",
                "0.0063578287760416666666666666666666667",
                "0.0026041666666666666666666666666666667",
                "8.23974609375e-4",
                "1.6276041666666666666666666666666667e-4",
                "1.0172526041666666666666666666666667e-5",
                "0",
            ];
            ddouble[] expected_n8 = [
                "0",
                "9.4610547262524801587301587301587302e-11",
                "1.2110150049603174603174603174603175e-8",
                "2.0691326686314174107142857142857143e-7",
                "1.5500992063492063492063492063492063e-6",
                "7.3914490048847501240079365079365079e-6",
                "2.6484898158482142857142857142857143e-5",
                "7.7915853924221462673611111111111111e-5",
                "1.9841269841269841269841269841269841e-4",
                "4.525185577453128875248015873015873e-4",
                "9.4600859142485119047619047619047619e-4",
                "0.0018420366067734975663442460317460317",
                "0.0033776661706349206349206349206349206",
                "0.0058775399412427629743303571428571429",
                "0.0097613501170324900793650793650793651",
                "0.015541772142289176819816468253968254",
                "0.02380952380952380952380952380952381",
                "0.035202214547566005161830357142857143",
                "0.050353967575799851190476190476190476",
                "0.069825952677499680292038690476190476",
                "0.09402436755952380952380952380952381",
                "0.12311556083815438406808035714285714",
                "0.15694830758231026785714285714285714",
                "0.19499325099445524669828869047619048",
                "0.23630952380952380952380952380952381",
                "0.27954855769399612668960813492063492",
                "0.32300445314437624007936507936507937",
                "0.36471397223926725841703869047619048",
                "0.40259641617063492063492063492063492",
                "0.43461733933479066879030257936507937",
                "0.45895941598074776785714285714285714",
                "0.47418377011541336301773313492063492",
                "0.47936507936507936507936507936507937",
                "0.47418377011541336301773313492063492",
                "0.45895941598074776785714285714285714",
                "0.43461733933479066879030257936507937",
                "0.40259641617063492063492063492063492",
                "0.36471397223926725841703869047619048",
                "0.32300445314437624007936507936507937",
                "0.27954855769399612668960813492063492",
                "0.23630952380952380952380952380952381",
                "0.19499325099445524669828869047619048",
                "0.15694830758231026785714285714285714",
                "0.12311556083815438406808035714285714",
                "0.09402436755952380952380952380952381",
                "0.069825952677499680292038690476190476",
                "0.050353967575799851190476190476190476",
                "0.035202214547566005161830357142857143",
                "0.02380952380952380952380952380952381",
                "0.015541772142289176819816468253968254",
                "0.0097613501170324900793650793650793651",
                "0.0058775399412427629743303571428571429",
                "0.0033776661706349206349206349206349206",
                "0.0018420366067734975663442460317460317",
                "9.4600859142485119047619047619047619e-4",
                "4.525185577453128875248015873015873e-4",
                "1.9841269841269841269841269841269841e-4",
                "7.7915853924221462673611111111111111e-5",
                "2.6484898158482142857142857142857143e-5",
                "7.3914490048847501240079365079365079e-6",
                "1.5500992063492063492063492063492063e-6",
                "2.0691326686314174107142857142857143e-7",
                "1.2110150049603174603174603174603175e-8",
                "9.4610547262524801587301587301587302e-11",
                "0",
            ];

            foreach ((IrwinHallDistribution dist, ddouble[] expecteds) in new[]{
                (dist_n1, expected_n1), (dist_n2, expected_n2), (dist_n3, expected_n3),
                (dist_n4, expected_n4), (dist_n5, expected_n5), (dist_n8, expected_n8),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else if (ddouble.IsPositiveInfinity(expected)) {
                        Assert.IsTrue(ddouble.IsPositiveInfinity(actual));
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_n1 = [
                0,
                0.125,
                0.25,
                0.375,
                0.5,
                0.625,
                0.75,
                0.875,
                1,
            ];
            ddouble[] expected_n2 = [
                0,
                0.0078125,
                0.03125,
                0.0703125,
                0.125,
                0.1953125,
                0.28125,
                0.3828125,
                0.5,
                0.6171875,
                0.71875,
                0.8046875,
                0.875,
                0.9296875,
                0.96875,
                0.9921875,
                1,
            ];
            ddouble[] expected_n3 = [
                "0",
                "3.2552083333333333333333333333333333e-4",
                "0.0026041666666666666666666666666666667",
                "0.0087890625",
                "0.020833333333333333333333333333333333",
                "0.040690104166666666666666666666666667",
                "0.0703125",
                "0.11165364583333333333333333333333333",
                "0.16666666666666666666666666666666667",
                "0.236328125",
                "0.31770833333333333333333333333333333",
                "0.40690104166666666666666666666666667",
                "0.5",
                "0.59309895833333333333333333333333333",
                "0.68229166666666666666666666666666667",
                "0.763671875",
                "0.83333333333333333333333333333333333",
                "0.88834635416666666666666666666666667",
                "0.9296875",
                "0.95930989583333333333333333333333333",
                "0.97916666666666666666666666666666667",
                "0.9912109375",
                "0.99739583333333333333333333333333333",
                "0.99967447916666666666666666666666667",
                "1",
            ];
            ddouble[] expected_n4 = [
                "0",
                "1.0172526041666666666666666666666667e-5",
                "1.6276041666666666666666666666666667e-4",
                "8.23974609375e-4",
                "0.0026041666666666666666666666666666667",
                "0.0063578287760416666666666666666666667",
                "0.01318359375",
                "0.024424235026041666666666666666666667",
                "0.041666666666666666666666666666666667",
                "0.066701253255208333333333333333333333",
                "0.10107421875",
                "0.14564005533854166666666666666666667",
                "0.20052083333333333333333333333333333",
                "0.265106201171875",
                "0.33805338541666666666666666666666667",
                "0.41728719075520833333333333333333333",
                "0.5",
                "0.58271280924479166666666666666666667",
                "0.66194661458333333333333333333333333",
                "0.734893798828125",
                "0.79947916666666666666666666666666667",
                "0.85435994466145833333333333333333333",
                "0.89892578125",
                "0.93329874674479166666666666666666667",
                "0.95833333333333333333333333333333333",
                "0.97557576497395833333333333333333333",
                "0.98681640625",
                "0.99364217122395833333333333333333333",
                "0.99739583333333333333333333333333333",
                "0.999176025390625",
                "0.99983723958333333333333333333333333",
                "0.99998982747395833333333333333333333",
                "1",
            ];
            ddouble[] expected_n5 = [
                "0",
                "2.5431315104166666666666666666666667e-7",
                "8.1380208333333333333333333333333333e-6",
                "6.1798095703125e-5",
                "2.6041666666666666666666666666666667e-4",
                "7.9472859700520833333333333333333333e-4",
                "0.0019775390625",
                "0.0042742411295572916666666666666666667",
                "0.0083333333333333333333333333333333333",
                "0.015015665690104166666666666666666667",
                "0.025390625",
                "0.040648396809895833333333333333333333",
                "0.061979166666666666666666666666666667",
                "0.0904510498046875",
                "0.12688802083333333333333333333333333",
                "0.17174784342447916666666666666666667",
                "0.225",
                "0.28600616455078125",
                "0.353466796875",
                "0.42553558349609375",
                "0.5",
                "0.57446441650390625",
                "0.646533203125",
                "0.71399383544921875",
                "0.775",
                "0.82825215657552083333333333333333333",
                "0.87311197916666666666666666666666667",
                "0.9095489501953125",
                "0.93802083333333333333333333333333333",
                "0.95935160319010416666666666666666667",
                "0.974609375",
                "0.98498433430989583333333333333333333",
                "0.99166666666666666666666666666666667",
                "0.99572575887044270833333333333333333",
                "0.9980224609375",
                "0.99920527140299479166666666666666667",
                "0.99973958333333333333333333333333333",
                "0.999938201904296875",
                "0.99999186197916666666666666666666667",
                "0.99999974568684895833333333333333333",
                "1",
            ];
            ddouble[] expected_n8 = [
                "0",
                "1.4782898009769500248015873015873016e-12",
                "3.7844218905009920634920634920634921e-10",
                "9.6990593842097691127232142857142857e-9",
                "9.6881200396825396825396825396825397e-8",
                "5.7745695350662110343812003968253968e-7",
                "2.4829592023577008928571428571428571e-6",
                "8.5220465229617224799262152777777778e-6",
                "2.4801587301587301587301587301587302e-5",
                "6.3635516793481887332976810515873016e-5",
                "1.4782595256018260168650793650793651e-4",
                "3.1680695505605803595648871527777778e-4",
                "6.3486250620039682539682539682539683e-4",
                "0.0012012667495698209792848617311507937",
                "0.0021617802362593393477182539682539683",
                "0.0037205186997732472798180958581349206",
                "0.0061507936507936507936507936507936508",
                "0.0098031068915530802711607917906746032",
                "0.015108074082268608940972222222222222",
                "0.022571823221173078294784303695436508",
                "0.03276183113219246031746031746031746",
                "0.046282224790267054996793232266865079",
                "0.063738815746610126798115079365079365",
                "0.085695387572345752564687577504960317",
                "0.11262400793650793650793650793650794",
                "0.14485338854913910230000813802083333",
                "0.18252054736727759951636904761904762",
                "0.22553190323598091564481220548115079",
                "0.27353951590401785714285714285714286",
                "0.32593647411447905358814057849702381",
                "0.38187336883847675626240079365079365",
                "0.44029570264032199269249325706845238",
                "0.5",
                "0.55970429735967800730750674293154762",
                "0.61812663116152324373759920634920635",
                "0.67406352588552094641185942150297619",
                "0.72646048409598214285714285714285714",
                "0.77446809676401908435518779451884921",
                "0.81747945263272240048363095238095238",
                "0.85514661145086089769999186197916667",
                "0.88737599206349206349206349206349206",
                "0.91430461242765424743531242249503968",
                "0.93626118425338987320188492063492064",
                "0.95371777520973294500320676773313492",
                "0.96723816886780753968253968253968254",
                "0.97742817677882692170521569630456349",
                "0.98489192591773139105902777777777778",
                "0.9901968931084469197288392082093254",
                "0.99384920634920634920634920634920635",
                "0.99627948130022675272018190414186508",
                "0.99783821976374066065228174603174603",
                "0.99879873325043017902071513826884921",
                "0.9993651374937996031746031746031746",
                "0.99968319304494394196404351128472222",
                "0.99985217404743981739831349206349206",
                "0.99993636448320651811266702318948413",
                "0.99997519841269841269841269841269841",
                "0.99999147795347703827752007378472222",
                "0.99999751704079764229910714285714286",
                "0.99999942254304649337889656187996032",
                "0.9999999031187996031746031746031746",
                "0.99999999030094061579023088727678571",
                "0.99999999962155781094990079365079365",
                "0.99999999999852171019902304997519841",
                "1",
            ];

            foreach ((IrwinHallDistribution dist, ddouble[] expecteds) in new[]{
                (dist_n1, expected_n1), (dist_n2, expected_n2), (dist_n3, expected_n3),
                (dist_n4, expected_n4), (dist_n5, expected_n5), (dist_n8, expected_n8),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 0.125, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0d) {
                        Assert.IsTrue(actual == 0d);
                    }
                    else {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                    }
                }
            }
        }
    }
}