﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class WignerSemicircleDistributionTests {
        readonly WignerSemicircleDistribution dist_r1 = new(r: 1);
        readonly WignerSemicircleDistribution dist_r2 = new(r: 2);
        readonly WignerSemicircleDistribution dist_r3 = new(r: 3);

        WignerSemicircleDistribution[] Dists => [
            dist_r1,
            dist_r2,
            dist_r3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (WignerSemicircleDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"R={dist.R}");
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
            foreach (WignerSemicircleDistribution dist in Dists) {
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
            foreach (WignerSemicircleDistribution dist in Dists) {
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
            foreach (WignerSemicircleDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (WignerSemicircleDistribution dist in Dists) {
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
            foreach (WignerSemicircleDistribution dist in Dists) {
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
            foreach (WignerSemicircleDistribution dist in Dists) {
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
            foreach (WignerSemicircleDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (WignerSemicircleDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (WignerSemicircleDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -4; x <= 4; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (WignerSemicircleDistribution dist in Dists) {
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
            foreach (WignerSemicircleDistribution dist in Dists) {
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

            for (long h = 1000; h <= 1000000000000; h *= 100) {
                for (int i = 1000; i >= 1; i--) {
                    ddouble p = (ddouble)i / h;
                    ddouble x = dist_r1.Quantile(p, Interval.Lower);
                    ddouble cdf = dist_r1.CDF(x, Interval.Lower);

                    Console.WriteLine($"quantile({p})={x}, cdf({x})={cdf}");

                    if (ddouble.IsFinite(x)) {
                        Assert.IsTrue(ddouble.Abs(p - cdf) < 1e-20);
                    }
                }
            }
        }

        [TestMethod()]
        public void QuantileUpperTest() {
            foreach (WignerSemicircleDistribution dist in Dists) {
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

            foreach (WignerSemicircleDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 100000).ToArray();

                double max_error = 0d;

                for (int i = 5; i <= 90; i++) {
                    double p = (double)i / 100;
                    double expected = (double)dist.Quantile(p, Interval.Lower);
                    double actual = xs.Quantile((double)p);

                    max_error = double.Max(max_error, double.Abs(expected - actual));

                    Assert.AreEqual(expected, actual, (double.Abs(expected) + 2) * 0.01, $"{p}\n{expected}\n{actual}");
                }

                Console.WriteLine(max_error);
            }
        }

        [TestMethod()]
        public void FitTest() {
            Random random = new(1234);

            foreach (WignerSemicircleDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                (WignerSemicircleDistribution? dist_fit, ddouble error) = WignerSemicircleDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (WignerSemicircleDistribution dist in Dists) {
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
            ddouble[] expected_dist_r1 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                "3.082022220307499027824845471717387725006e-1",
                "4.210843993477923898408659806294586443604e-1",
                "4.969611505220486728408461856211720634528e-1",
                "5.513288954217920495113264983129694413974e-1",
                "5.901623240859552611248590547677828824232e-1",
                "6.164044440614998055649690943434775450011e-1",
                "6.316265990216885847612989709441879665405e-1",
                "6.366197723675813430755350534900574481378e-1",
                "6.316265990216885847612989709441879665405e-1",
                "6.164044440614998055649690943434775450011e-1",
                "5.901623240859552611248590547677828824232e-1",
                "5.513288954217920495113264983129694413974e-1",
                "4.969611505220486728408461856211720634528e-1",
                "4.210843993477923898408659806294586443604e-1",
                "3.082022220307499027824845471717387725006e-1",
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0
            ];
            ddouble[] expected_dist_r2 = [
                0,
                "1.107671525394118791634956154030184012551e-1",
                "1.541011110153749513912422735858693862503e-1",
                "1.855623102988608703050018388047598334955e-1",
                "2.105421996738961949204329903147293221802e-1",
                "2.311516665230624270868634103788040793754e-1",
                "2.484805752610243364204230928105860317264e-1",
                "2.631777495923702436505412378934116527252e-1",
                "2.756644477108960247556632491564847206987e-1",
                "2.862301099391825992289047673359533925688e-1",
                "2.950811620429776305624295273838914412116e-1",
                "3.023682139416538000385001287382432112869e-1",
                "3.082022220307499027824845471717387725006e-1",
                "3.126645339335944481668424192252557711944e-1",
                "3.158132995108442923806494854720939832703e-1",
                "3.176875788707120502615896026407171984208e-1",
                "3.183098861837906715377675267450287240689e-1",
                "3.176875788707120502615896026407171984208e-1",
                "3.158132995108442923806494854720939832703e-1",
                "3.126645339335944481668424192252557711944e-1",
                "3.082022220307499027824845471717387725006e-1",
                "3.023682139416538000385001287382432112869e-1",
                "2.950811620429776305624295273838914412116e-1",
                "2.862301099391825992289047673359533925688e-1",
                "2.756644477108960247556632491564847206987e-1",
                "2.631777495923702436505412378934116527252e-1",
                "2.484805752610243364204230928105860317264e-1",
                "2.311516665230624270868634103788040793754e-1",
                "2.105421996738961949204329903147293221802e-1",
                "1.855623102988608703050018388047598334955e-1",
                "1.541011110153749513912422735858693862503e-1",
                "1.107671525394118791634956154030184012551e-1",
                0,
            ];
            ddouble[] expected_dist_r3 = [
                "1.581694540927060129948361503119562366746e-1",
                "1.656537168406828909469487285403906878176e-1",
                "1.723611665862753261631557921018833139202e-1",
                "1.783794543649308776651916945091818799576e-1",
                "1.837762984739306831704421661043231471325e-1",
                "1.886050534071877745962358775661181717854e-1",
                "1.9290838397896848394527143831919935134e-1",
                "1.967207746953184203749530182559276274744e-1",
                "2.00070292479356904345599820223720129481e-1",
                "2.029798532306216237824329701183667140626e-1",
                "2.054681480204999351883230314478258483337e-1",
                "2.075503289159608709977375074098210270395e-1",
                "2.09238520268073753928918433455679268566e-1",
                "2.105421996738961949204329903147293221802e-1",
                "2.114684786718165590524885594241174585158e-1",
                "2.120223036586537772065961239525580685695e-1",
                "2.122065907891937810251783511633524827126e-1",
                "2.120223036586537772065961239525580685695e-1",
                "2.114684786718165590524885594241174585158e-1",
                "2.105421996738961949204329903147293221802e-1",
                "2.09238520268073753928918433455679268566e-1",
                "2.075503289159608709977375074098210270395e-1",
                "2.054681480204999351883230314478258483337e-1",
                "2.029798532306216237824329701183667140626e-1",
                "2.00070292479356904345599820223720129481e-1",
                "1.967207746953184203749530182559276274744e-1",
                "1.9290838397896848394527143831919935134e-1",
                "1.886050534071877745962358775661181717854e-1",
                "1.837762984739306831704421661043231471325e-1",
                "1.783794543649308776651916945091818799576e-1",
                "1.723611665862753261631557921018833139202e-1",
                "1.656537168406828909469487285403906878176e-1",
                "1.581694540927060129948361503119562366746e-1"
            ];

            foreach ((WignerSemicircleDistribution dist, ddouble[] expecteds) in new[]{
                (dist_r1, expected_dist_r1), (dist_r2, expected_dist_r2),
                (dist_r3, expected_dist_r3),
            }) {
                for ((ddouble x, int i) = (-2, 0); i < expecteds.Length; x += 1 / 8d, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF(x);

                    Console.WriteLine($"{dist} pdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0 && actual == 0) {
                        continue;
                    }

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} pdf({x})\n{expected}\n{actual}");
                }
            }
        }

        [TestMethod()]
        public void CDFExpectedTest() {
            ddouble[] expected_dist_r1 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                "2.602277437187940507456339726933279819916e-2",
                "7.214680640719373902345582496909590137486e-2",
                "1.297985990535851231209091871360261893886e-1",
                "1.95501109477885320955501708755090972984e-1",
                "2.669872706947602631985512872164145853293e-1",
                "3.425188212371462805334286698785810512908e-1",
                "4.20630249886376233879979949726274038111e-1",
                "5.0e-1",
                "5.79369750113623766120020050273725961889e-1",
                "6.574811787628537194665713301214189487092e-1",
                "7.330127293052397368014487127835854146707e-1",
                "8.04498890522114679044498291244909027016e-1",
                "8.702014009464148768790908128639738106114e-1",
                "9.278531935928062609765441750309040986251e-1",
                "9.739772256281205949254366027306672018009e-1",
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1
            ];
            ddouble[] expected_dist_r2 = [
                0,
                "9.289876751622511909728689064552356662445e-3",
                "2.602277437187940507456339726933279819916e-2",
                "4.733666563850198931170312973072765302603e-2",
                "7.214680640719373902345582496909590137486e-2",
                "9.979135949989585006493082537911942867707e-2",
                "1.297985990535851231209091871360261893886e-1",
                "1.618021473664440823052288505695272663123e-1",
                "1.95501109477885320955501708755090972984e-1",
                "2.306383281662186890187387604703506770973e-1",
                "2.669872706947602631985512872164145853293e-1",
                "3.043435061194386807450351511628247405206e-1",
                "3.425188212371462805334286698785810512908e-1",
                "3.813369394715602548543125461624700856347e-1",
                "4.20630249886376233879979949726274038111e-1",
                "4.602371835512917584438689851100398648338e-1",
                "5.0e-1",
                "5.397628164487082415561310148899601351662e-1",
                "5.79369750113623766120020050273725961889e-1",
                "6.186630605284397451456874538375299143653e-1",
                "6.574811787628537194665713301214189487092e-1",
                "6.956564938805613192549648488371752594794e-1",
                "7.330127293052397368014487127835854146707e-1",
                "7.693616718337813109812612395296493229027e-1",
                "8.04498890522114679044498291244909027016e-1",
                "8.381978526335559176947711494304727336877e-1",
                "8.702014009464148768790908128639738106114e-1",
                "9.002086405001041499350691746208805713229e-1",
                "9.278531935928062609765441750309040986251e-1",
                "9.52663334361498010688296870269272346974e-1",
                "9.739772256281205949254366027306672018009e-1",
                "9.907101232483774880902713109354476433376e-1",
                1
            ];
            ddouble[] expected_dist_r3 = [
                "1.095510187085239977503976456292756548197e-1",
                "1.297985990535851231209091871360261893886e-1",
                "1.509321197222754563669861768829263348968e-1",
                "1.728602027961400797765469233340056643353e-1",
                "1.95501109477885320955501708755090972984e-1",
                "2.187806216937912044889727247807869203246e-1",
                "2.42630491390010943429939025877592311537e-1",
                "2.669872706947602631985512872164145853293e-1",
                "2.917914057909288179980556491105708465761e-1",
                "3.169865181004688999016442609070936774172e-1",
                "3.425188212371462805334286698785810512908e-1",
                "3.683366378772854315907435761042828594529e-1",
                "3.943899908948678868264452369655296780527e-1",
                "4.20630249886376233879979949726274038111e-1",
                "4.470098187926942698829332164051651251385e-1",
                "4.734818534476187126098955715794373601367e-1",
                "5.0e-1",
                "5.265181465523812873901044284205626398633e-1",
                "5.529901812073057301170667835948348748615e-1",
                "5.79369750113623766120020050273725961889e-1",
                "6.056100091051321131735547630344703219473e-1",
                "6.316633621227145684092564238957171405471e-1",
                "6.574811787628537194665713301214189487092e-1",
                "6.830134818995311000983557390929063225828e-1",
                "7.082085942090711820019443508894291534239e-1",
                "7.330127293052397368014487127835854146707e-1",
                "7.57369508609989056570060974122407688463e-1",
                "7.812193783062087955110272752192130796754e-1",
                "8.04498890522114679044498291244909027016e-1",
                "8.271397972038599202234530766659943356647e-1",
                "8.490678802777245436330138231170736651032e-1",
                "8.702014009464148768790908128639738106114e-1",
                "8.904489812914760022496023543707243451803e-1"
            ];

            foreach ((WignerSemicircleDistribution dist, ddouble[] expecteds) in new[]{
                (dist_r1, expected_dist_r1), (dist_r2, expected_dist_r2),
                (dist_r3, expected_dist_r3),
            }) {
                for ((ddouble x, int i) = (-2, 0); i < expecteds.Length; x += 1 / 8d, i++) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.CDF(x);

                    Console.WriteLine($"{dist} cdf({x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected == 0 && actual == 0) {
                        continue;
                    }

                    Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-30, $"{dist} cdf({x})\n{expected}\n{actual}");
                }
            }
        }
    }
}