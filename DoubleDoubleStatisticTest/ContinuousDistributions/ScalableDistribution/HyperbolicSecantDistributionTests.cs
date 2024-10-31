﻿using DoubleDouble;
using DoubleDoubleStatistic;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ContinuousDistributions.ScalableDistribution {
    [TestClass()]
    public class HyperbolicSecantDistributionTests {
        readonly HyperbolicSecantDistribution dist = new();
        readonly HyperbolicSecantDistribution dist_pi2 = new(sigma: ddouble.Pi / 2);

        HyperbolicSecantDistribution[] Dists => [
            dist,
            dist_pi2
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (HyperbolicSecantDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"S={dist.Sigma}");
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (HyperbolicSecantDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (HyperbolicSecantDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
            foreach (HyperbolicSecantDistribution dist in Dists) {
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

            foreach (HyperbolicSecantDistribution dist in Dists) {

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

            foreach (HyperbolicSecantDistribution dist in Dists) {

                Console.WriteLine(dist);

                double[] xs = dist.Sample(random, 10000).ToArray();

                (HyperbolicSecantDistribution? dist_fit, ddouble error) = HyperbolicSecantDistribution.Fit(xs, (0.05, 0.95));

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);
                Console.WriteLine(error);

                Assert.IsTrue(error < 1e-2);
            }
        }

        [TestMethod()]
        public void IrregularValueTest() {
            foreach (HyperbolicSecantDistribution dist in Dists) {
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
                2.062820820908705111e-01,
                2.068960160876329446e-01,
                2.075104395435024562e-01,
                2.081253379978671436e-01,
                2.087406968372008398e-01,
                2.093565012946668469e-01,
                2.099727364497356585e-01,
                2.105893872278191414e-01,
                2.112064383999196215e-01,
                2.118238745822957902e-01,
                2.124416802361444034e-01,
                2.130598396672990225e-01,
                2.136783370259454917e-01,
                2.142971563063550688e-01,
                2.149162813466349298e-01,
                2.155356958284967983e-01,
                2.161553832770439498e-01,
                2.167753270605763949e-01,
                2.173955103904156028e-01,
                2.180159163207482376e-01,
                2.186365277484892955e-01,
                2.192573274131654482e-01,
                2.198782978968185087e-01,
                2.204994216239296867e-01,
                2.211206808613641883e-01,
                2.217420577183377983e-01,
                2.223635341464045845e-01,
                2.229850919394665565e-01,
                2.236067127338053351e-01,
                2.242283780081372746e-01,
                2.248500690836899851e-01,
                2.254717671243036681e-01,
                2.260934531365547118e-01,
                2.267151079699039062e-01,
                2.273367123168685000e-01,
                2.279582467132189871e-01,
                2.285796915382005967e-01,
                2.292010270147801232e-01,
                2.298222332099179588e-01,
                2.304432900348660773e-01,
                2.310641772454924414e-01,
                2.316848744426314721e-01,
                2.323053610724613027e-01,
                2.329256164269083718e-01,
                2.335456196440792720e-01,
                2.341653497087201885e-01,
                2.347847854527046751e-01,
                2.354039055555492421e-01,
                2.360226885449581424e-01,
                2.366411127973965800e-01,
                2.372591565386934775e-01,
                2.378767978446734710e-01,
                2.384940146418190465e-01,
                2.391107847079625426e-01,
                2.397270856730082944e-01,
                2.403428950196859759e-01,
                2.409581900843343893e-01,
                2.415729480577167021e-01,
                2.421871459858670750e-01,
                2.428007607709689875e-01,
                2.434137691722655650e-01,
                2.440261478070019918e-01,
                2.446378731514007587e-01,
                2.452489215416694401e-01,
                2.458592691750410286e-01,
                2.464688921108484365e-01,
                2.470777662716312773e-01,
                2.476858674442770081e-01,
                2.482931712811959624e-01,
                2.488996533015297719e-01,
                2.495052888923948720e-01,
                2.501100533101598966e-01,
                2.507139216817579341e-01,
                2.513168690060332833e-01,
                2.519188701551239862e-01,
                2.525198998758785285e-01,
                2.531199327913086217e-01,
                2.537189434020770129e-01,
                2.543169060880208221e-01,
                2.549137951097114052e-01,
                2.555095846100488011e-01,
                2.561042486158933706e-01,
                2.566977610397329079e-01,
                2.572900956813858886e-01,
                2.578812262297411340e-01,
                2.584711262645342789e-01,
                2.590597692581596556e-01,
                2.596471285775195925e-01,
                2.602331774859097391e-01,
                2.608178891449413617e-01,
                2.614012366164998880e-01,
                2.619831928647403108e-01,
                2.625637307581188962e-01,
                2.631428230714620842e-01,
                2.637204424880718601e-01,
                2.642965616018670971e-01,
                2.648711529195624248e-01,
                2.654441888628830681e-01,
                2.660156417708162691e-01,
                2.665854839018990119e-01,
                2.671536874365424974e-01,
                2.677202244793920882e-01,
                2.682850670617242250e-01,
                2.688481871438788695e-01,
                2.694095566177276413e-01,
                2.699691473091783700e-01,
                2.705269309807145084e-01,
                2.710828793339705167e-01,
                2.716369640123420526e-01,
                2.721891566036314103e-01,
                2.727394286427279324e-01,
                2.732877516143226715e-01,
                2.738340969556578020e-01,
                2.743784360593099492e-01,
                2.749207402760072694e-01,
                2.754609809174803914e-01,
                2.759991292593467205e-01,
                2.765351565440272719e-01,
                2.770690339836967553e-01,
                2.776007327632654120e-01,
                2.781302240433937145e-01,
                2.786574789635375970e-01,
                2.791824686450260495e-01,
                2.797051641941686873e-01,
                2.802255367053945734e-01,
                2.807435572644206290e-01,
                2.812591969514495749e-01,
                2.817724268443976277e-01,
                2.822832180221505616e-01,
                2.827915415678474709e-01,
                2.832973685721938417e-01,
                2.838006701367999374e-01,
                2.843014173775477160e-01,
                2.847995814279826177e-01,
                2.852951334427314412e-01,
                2.857880446009461450e-01,
                2.862782861097705189e-01,
                2.867658292078324456e-01,
                2.872506451687583673e-01,
                2.877327053047109562e-01,
                2.882119809699489887e-01,
                2.886884435644079816e-01,
                2.891620645373025877e-01,
                2.896328153907480885e-01,
                2.901006676834021469e-01,
                2.905655930341247695e-01,
                2.910275631256562523e-01,
                2.914865497083127810e-01,
                2.919425246036983501e-01,
                2.923954597084321150e-01,
                2.928453269978915086e-01,
                2.932920985299687922e-01,
                2.937357464488415948e-01,
                2.941762429887558872e-01,
                2.946135604778207795e-01,
                2.950476713418143659e-01,
                2.954785481079997278e-01,
                2.959061634089500403e-01,
                2.963304899863825614e-01,
                2.967515006949993372e-01,
                2.971691685063361787e-01,
                2.975834665126153022e-01,
                2.979943679306050752e-01,
                2.984018461054814830e-01,
                2.988058745146949247e-01,
                2.992064267718372861e-01,
                2.996034766305113450e-01,
                2.999969979881997317e-01,
                3.003869648901336120e-01,
                3.007733515331595942e-01,
                3.011561322696031939e-01,
                3.015352816111304679e-01,
                3.019107742326028188e-01,
                3.022825849759279149e-01,
                3.026506888539028939e-01,
                3.030150610540505163e-01,
                3.033756769424464372e-01,
                3.037325120675368195e-01,
                3.040855421639453993e-01,
                3.044347431562686168e-01,
                3.047800911628582021e-01,
                3.051215624995897158e-01,
                3.054591336836168236e-01,
                3.057927814371087516e-01,
                3.061224826909723085e-01,
                3.064482145885540909e-01,
                3.067699544893254804e-01,
                3.070876799725464901e-01,
                3.074013688409089617e-01,
                3.077109991241574471e-01,
                3.080165490826872188e-01,
                3.083179972111177447e-01,
                3.086153222418410724e-01,
                3.089085031485437338e-01,
                3.091975191497017827e-01,
                3.094823497120466893e-01,
                3.097629745540024238e-01,
                3.100393736490921759e-01,
                3.103115272293129334e-01,
                3.105794157884781415e-01,
                3.108430200855268333e-01,
                3.111023211477980666e-01,
                3.113573002742694995e-01,
                3.116079390387603842e-01,
                3.118542192930962020e-01,
                3.120961231702350513e-01,
                3.123336330873546784e-01,
                3.125667317488992070e-01,
                3.127954021495836789e-01,
                3.130196275773574621e-01,
                3.132393916163233039e-01,
                3.134546781496128087e-01,
                3.136654713622166191e-01,
                3.138717557437683547e-01,
                3.140735160912818125e-01,
                3.142707375118403146e-01,
                3.144634054252376520e-01,
                3.146515055665686789e-01,
                3.148350239887703905e-01,
                3.150139470651115436e-01,
                3.151882614916296510e-01,
                3.153579542895162402e-01,
                3.155230128074471563e-01,
                3.156834247238599067e-01,
                3.158391780491746625e-01,
                3.159902611279600815e-01,
                3.161366626410421210e-01,
                3.162783716075558971e-01,
                3.164153773869392561e-01,
                3.165476696808680046e-01,
                3.166752385351314092e-01,
                3.167980743414479106e-01,
                3.169161678392207193e-01,
                3.170295101172314611e-01,
                3.171380926152728708e-01,
                3.172419071257185363e-01,
                3.173409457950305823e-01,
                3.174352011252032368e-01,
                3.175246659751431721e-01,
                3.176093335619851743e-01,
                3.176891974623433645e-01,
                3.177642516134973061e-01,
                3.178344903145121081e-01,
                3.178999082272928600e-01,
                3.179605003775729521e-01,
                3.180162621558353941e-01,
                3.180671893181672982e-01,
                3.181132779870470273e-01,
                3.181545246520639525e-01,
                3.181909261705703762e-01,
                3.182224797682652317e-01,
                3.182491830397100596e-01,
                3.182710339487757611e-01,
                3.182880308290217952e-01,
                3.183001723840059305e-01,
                3.183074576875257189e-01,
                3.183098861837906912e-01,
                3.183074576875257189e-01,
                3.183001723840059305e-01,
                3.182880308290217952e-01,
                3.182710339487757611e-01,
                3.182491830397100596e-01,
                3.182224797682652317e-01,
                3.181909261705703762e-01,
                3.181545246520639525e-01,
                3.181132779870470273e-01,
                3.180671893181672982e-01,
                3.180162621558353941e-01,
                3.179605003775729521e-01,
                3.178999082272928600e-01,
                3.178344903145121081e-01,
                3.177642516134973061e-01,
                3.176891974623433645e-01,
                3.176093335619851743e-01,
                3.175246659751431721e-01,
                3.174352011252032368e-01,
                3.173409457950305823e-01,
                3.172419071257185363e-01,
                3.171380926152728708e-01,
                3.170295101172314611e-01,
                3.169161678392207193e-01,
                3.167980743414479106e-01,
                3.166752385351314092e-01,
                3.165476696808680046e-01,
                3.164153773869392561e-01,
                3.162783716075558971e-01,
                3.161366626410421210e-01,
                3.159902611279600815e-01,
                3.158391780491746625e-01,
                3.156834247238599067e-01,
                3.155230128074471563e-01,
                3.153579542895162402e-01,
                3.151882614916296510e-01,
                3.150139470651115436e-01,
                3.148350239887703905e-01,
                3.146515055665686789e-01,
                3.144634054252376520e-01,
                3.142707375118403146e-01,
                3.140735160912818125e-01,
                3.138717557437683547e-01,
                3.136654713622166191e-01,
                3.134546781496128087e-01,
                3.132393916163233039e-01,
                3.130196275773574621e-01,
                3.127954021495836789e-01,
                3.125667317488992070e-01,
                3.123336330873546784e-01,
                3.120961231702350513e-01,
                3.118542192930962020e-01,
                3.116079390387603842e-01,
                3.113573002742694995e-01,
                3.111023211477980666e-01,
                3.108430200855268333e-01,
                3.105794157884781415e-01,
                3.103115272293129334e-01,
                3.100393736490921759e-01,
                3.097629745540024238e-01,
                3.094823497120466893e-01,
                3.091975191497017827e-01,
                3.089085031485437338e-01,
                3.086153222418410724e-01,
                3.083179972111177447e-01,
                3.080165490826872188e-01,
                3.077109991241574471e-01,
                3.074013688409089617e-01,
                3.070876799725464901e-01,
                3.067699544893254804e-01,
                3.064482145885540909e-01,
                3.061224826909723085e-01,
                3.057927814371087516e-01,
                3.054591336836168236e-01,
                3.051215624995897158e-01,
                3.047800911628582021e-01,
                3.044347431562686168e-01,
                3.040855421639453993e-01,
                3.037325120675368195e-01,
                3.033756769424464372e-01,
                3.030150610540505163e-01,
                3.026506888539028939e-01,
                3.022825849759279149e-01,
                3.019107742326028188e-01,
                3.015352816111304679e-01,
                3.011561322696031939e-01,
                3.007733515331595942e-01,
                3.003869648901336120e-01,
                2.999969979881997317e-01,
                2.996034766305113450e-01,
                2.992064267718372861e-01,
                2.988058745146949247e-01,
                2.984018461054814830e-01,
                2.979943679306050752e-01,
                2.975834665126153022e-01,
                2.971691685063361787e-01,
                2.967515006949993372e-01,
                2.963304899863825614e-01,
                2.959061634089500403e-01,
                2.954785481079997278e-01,
                2.950476713418143659e-01,
                2.946135604778207795e-01,
                2.941762429887558872e-01,
                2.937357464488415948e-01,
                2.932920985299687922e-01,
                2.928453269978915086e-01,
                2.923954597084321150e-01,
                2.919425246036983501e-01,
                2.914865497083127810e-01,
                2.910275631256562523e-01,
                2.905655930341247695e-01,
                2.901006676834021469e-01,
                2.896328153907480885e-01,
                2.891620645373025877e-01,
                2.886884435644079816e-01,
                2.882119809699489887e-01,
                2.877327053047109562e-01,
                2.872506451687583673e-01,
                2.867658292078324456e-01,
                2.862782861097705189e-01,
                2.857880446009461450e-01,
                2.852951334427314412e-01,
                2.847995814279826177e-01,
                2.843014173775477160e-01,
                2.838006701367999374e-01,
                2.832973685721938417e-01,
                2.827915415678474709e-01,
                2.822832180221505616e-01,
                2.817724268443976277e-01,
                2.812591969514495749e-01,
                2.807435572644206290e-01,
                2.802255367053945734e-01,
                2.797051641941686873e-01,
                2.791824686450260495e-01,
                2.786574789635375970e-01,
                2.781302240433937145e-01,
                2.776007327632654120e-01,
                2.770690339836967553e-01,
                2.765351565440272719e-01,
                2.759991292593467205e-01,
                2.754609809174803914e-01,
                2.749207402760072694e-01,
                2.743784360593099492e-01,
                2.738340969556578020e-01,
                2.732877516143226715e-01,
                2.727394286427279324e-01,
                2.721891566036314103e-01,
                2.716369640123420526e-01,
                2.710828793339705167e-01,
                2.705269309807145084e-01,
                2.699691473091783700e-01,
                2.694095566177276413e-01,
                2.688481871438788695e-01,
                2.682850670617242250e-01,
                2.677202244793920882e-01,
                2.671536874365424974e-01,
                2.665854839018990119e-01,
                2.660156417708162691e-01,
                2.654441888628830681e-01,
                2.648711529195624248e-01,
                2.642965616018670971e-01,
                2.637204424880718601e-01,
                2.631428230714620842e-01,
                2.625637307581188962e-01,
                2.619831928647403108e-01,
                2.614012366164998880e-01,
                2.608178891449413617e-01,
                2.602331774859097391e-01,
                2.596471285775195925e-01,
                2.590597692581596556e-01,
                2.584711262645342789e-01,
                2.578812262297411340e-01,
                2.572900956813858886e-01,
                2.566977610397329079e-01,
                2.561042486158933706e-01,
                2.555095846100488011e-01,
                2.549137951097114052e-01,
                2.543169060880208221e-01,
                2.537189434020770129e-01,
                2.531199327913086217e-01,
                2.525198998758785285e-01,
                2.519188701551239862e-01,
                2.513168690060332833e-01,
                2.507139216817579341e-01,
                2.501100533101598966e-01,
                2.495052888923948720e-01,
                2.488996533015297719e-01,
                2.482931712811959624e-01,
                2.476858674442770081e-01,
                2.470777662716312773e-01,
                2.464688921108484365e-01,
                2.458592691750410286e-01,
                2.452489215416694401e-01,
                2.446378731514007587e-01,
                2.440261478070019918e-01,
                2.434137691722655650e-01,
                2.428007607709689875e-01,
                2.421871459858670750e-01,
                2.415729480577167021e-01,
                2.409581900843343893e-01,
                2.403428950196859759e-01,
                2.397270856730082944e-01,
                2.391107847079625426e-01,
                2.384940146418190465e-01,
                2.378767978446734710e-01,
                2.372591565386934775e-01,
                2.366411127973965800e-01,
                2.360226885449581424e-01,
                2.354039055555492421e-01,
                2.347847854527046751e-01,
                2.341653497087201885e-01,
                2.335456196440792720e-01,
                2.329256164269083718e-01,
                2.323053610724613027e-01,
                2.316848744426314721e-01,
                2.310641772454924414e-01,
                2.304432900348660773e-01,
                2.298222332099179588e-01,
                2.292010270147801232e-01,
                2.285796915382005967e-01,
                2.279582467132189871e-01,
                2.273367123168685000e-01,
                2.267151079699039062e-01,
                2.260934531365547118e-01,
                2.254717671243036681e-01,
                2.248500690836899851e-01,
                2.242283780081372746e-01,
                2.236067127338053351e-01,
                2.229850919394665565e-01,
                2.223635341464045845e-01,
                2.217420577183377983e-01,
                2.211206808613641883e-01,
                2.204994216239296867e-01,
                2.198782978968185087e-01,
                2.192573274131654482e-01,
                2.186365277484892955e-01,
                2.180159163207482376e-01,
                2.173955103904156028e-01,
                2.167753270605763949e-01,
                2.161553832770439498e-01,
                2.155356958284967983e-01,
                2.149162813466349298e-01,
                2.142971563063550688e-01,
                2.136783370259454917e-01,
                2.130598396672990225e-01,
                2.124416802361444034e-01,
                2.118238745822957902e-01,
                2.112064383999196215e-01,
                2.105893872278191414e-01,
                2.099727364497356585e-01,
                2.093565012946668469e-01,
                2.087406968372008398e-01,
                2.081253379978671436e-01,
                2.075104395435024562e-01,
                2.068960160876329446e-01,
                2.062820820908705111e-01,

            ];

            foreach ((HyperbolicSecantDistribution dist, ddouble[] expecteds) in new[]{
                (dist_pi2, expected_dist)
            }) {
                for ((ddouble x, int i) = (-1, 0); i < expecteds.Length; x += 1d / 256, i++) {
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
            ddouble[] expected_dist = [
                2.244170143285850183e-01,
                2.252240026399228989e-01,
                2.260333900915928174e-01,
                2.268451785673507581e-01,
                2.276593698941673793e-01,
                2.284759658416299910e-01,
                2.292949681213428681e-01,
                2.301163783863260370e-01,
                2.309401982304132572e-01,
                2.317664291876476435e-01,
                2.325950727316768440e-01,
                2.334261302751468037e-01,
                2.342596031690941105e-01,
                2.350954927023377050e-01,
                2.359338001008691732e-01,
                2.367745265272421795e-01,
                2.376176730799614556e-01,
                2.384632407928704290e-01,
                2.393112306345382967e-01,
                2.401616435076467382e-01,
                2.410144802483757120e-01,
                2.418697416257888366e-01,
                2.427274283412184652e-01,
                2.435875410276503450e-01,
                2.444500802491079705e-01,
                2.453150465000370206e-01,
                2.461824402046891846e-01,
                2.470522617165066548e-01,
                2.479245113175060633e-01,
                2.487991892176630027e-01,
                2.496762955542967677e-01,
                2.505558303914552365e-01,
                2.514377937193002510e-01,
                2.523221854534938302e-01,
                2.532090054345844665e-01,
                2.540982534273946714e-01,
                2.549899291204089091e-01,
                2.558840321251629191e-01,
                2.567805619756334257e-01,
                2.576795181276296809e-01,
                2.585808999581855616e-01,
                2.594847067649534433e-01,
                2.603909377655991286e-01,
                2.612995920971988295e-01,
                2.622106688156372600e-01,
                2.631241668950075052e-01,
                2.640400852270132215e-01,
                2.649584226203722803e-01,
                2.658791778002227435e-01,
                2.668023494075308921e-01,
                2.677279359985016982e-01,
                2.686559360439916277e-01,
                2.695863479289237752e-01,
                2.705191699517062176e-01,
                2.714544003236524228e-01,
                2.723920371684048769e-01,
                2.733320785213619875e-01,
                2.742745223291074286e-01,
                2.752193664488435387e-01,
                2.761666086478273829e-01,
                2.771162466028107008e-01,
                2.780682778994832960e-01,
                2.790227000319204786e-01,
                2.799795104020339487e-01,
                2.809387063190268985e-01,
                2.819002849988532788e-01,
                2.828642435636813390e-01,
                2.838305790413612195e-01,
                2.847992883648979290e-01,
                2.857703683719276522e-01,
                2.867438158042001084e-01,
                2.877196273070651289e-01,
                2.886977994289645078e-01,
                2.896783286209290154e-01,
                2.906612112360806854e-01,
                2.916464435291406532e-01,
                2.926340216559427110e-01,
                2.936239416729521379e-01,
                2.946161995367911346e-01,
                2.956107911037695879e-01,
                2.966077121294223296e-01,
                2.976069582680528902e-01,
                2.986085250722831930e-01,
                2.996124079926100747e-01,
                3.006186023769689686e-01,
                3.016271034703034148e-01,
                3.026379064141427877e-01,
                3.036510062461863502e-01,
                3.046663978998948030e-01,
                3.056840762040893278e-01,
                3.067040358825581792e-01,
                3.077262715536711601e-01,
                3.087507777300017575e-01,
                3.097775488179572712e-01,
                3.108065791174171588e-01,
                3.118378628213798187e-01,
                3.128713940156175877e-01,
                3.139071666783404546e-01,
                3.149451746798683782e-01,
                3.159854117823125974e-01,
                3.170278716392661011e-01,
                3.180725477955026470e-01,
                3.191194336866856607e-01,
                3.201685226390865169e-01,
                3.212198078693122016e-01,
                3.222732824840424115e-01,
                3.233289394797772554e-01,
                3.243867717425944486e-01,
                3.254467720479162773e-01,
                3.265089330602880535e-01,
                3.275732473331653960e-01,
                3.286397073087131693e-01,
                3.297083053176147471e-01,
                3.307790335788923786e-01,
                3.318518841997374347e-01,
                3.329268491753533121e-01,
                3.340039203888082731e-01,
                3.350830896109000534e-01,
                3.361643485000322951e-01,
                3.372476886021017473e-01,
                3.383331013503980134e-01,
                3.394205780655144555e-01,
                3.405101099552718114e-01,
                3.416016881146528017e-01,
                3.426953035257501168e-01,
                3.437909470577261706e-01,
                3.448886094667847901e-01,
                3.459882813961565606e-01,
                3.470899533760961053e-01,
                3.481936158238917445e-01,
                3.492992590438894762e-01,
                3.504068732275283926e-01,
                3.515164484533899736e-01,
                3.526279746872605370e-01,
                3.537414417822069002e-01,
                3.548568394786653646e-01,
                3.559741574045444668e-01,
                3.570933850753408301e-01,
                3.582145118942691164e-01,
                3.593375271524061332e-01,
                3.604624200288473190e-01,
                3.615891795908791062e-01,
                3.627177947941638836e-01,
                3.638482544829396148e-01,
                3.649805473902337782e-01,
                3.661146621380911292e-01,
                3.672505872378162284e-01,
                3.683883110902300140e-01,
                3.695278219859413071e-01,
                3.706691081056323056e-01,
                3.718121575203587881e-01,
                3.729569581918654175e-01,
                3.741034979729155330e-01,
                3.752517646076349855e-01,
                3.764017457318718507e-01,
                3.775534288735704624e-01,
                3.787068014531603244e-01,
                3.798618507839599534e-01,
                3.810185640725957668e-01,
                3.821769284194360128e-01,
                3.833369308190393565e-01,
                3.844985581606193414e-01,
                3.856617972285226181e-01,
                3.868266347027236041e-01,
                3.879930571593333544e-01,
                3.891610510711239201e-01,
                3.903306028080674728e-01,
                3.915016986378910269e-01,
                3.926743247266460957e-01,
                3.938484671392930014e-01,
                3.950241118403011731e-01,
                3.962012446942634325e-01,
                3.973798514665266013e-01,
                3.985599178238360407e-01,
                3.997414293349957903e-01,
                4.009243714715439166e-01,
                4.021087296084417950e-01,
                4.032944890247797676e-01,
                4.044816349044967341e-01,
                4.056701523371145646e-01,
                4.068600263184879995e-01,
                4.080512417515684276e-01,
                4.092437834471832625e-01,
                4.104376361248291416e-01,
                4.116327844134801128e-01,
                4.128292128524103100e-01,
                4.140269058920311718e-01,
                4.152258478947422615e-01,
                4.164260231357979070e-01,
                4.176274158041858864e-01,
                4.188300100035223794e-01,
                4.200337897529594078e-01,
                4.212387389881069866e-01,
                4.224448415619684316e-01,
                4.236520812458900997e-01,
                4.248604417305243963e-01,
                4.260699066268059387e-01,
                4.272804594669417089e-01,
                4.284920837054139175e-01,
                4.297047627199964137e-01,
                4.309184798127835281e-01,
                4.321332182112325726e-01,
                4.333489610692182192e-01,
                4.345656914681001459e-01,
                4.357833924178028395e-01,
                4.370020468579073336e-01,
                4.382216376587557694e-01,
                4.394421476225673362e-01,
                4.406635594845666470e-01,
                4.418858559141229936e-01,
                4.431090195159017586e-01,
                4.443330328310268751e-01,
                4.455578783382543873e-01,
                4.467835384551571143e-01,
                4.480099955393198052e-01,
                4.492372318895455074e-01,
                4.504652297470715938e-01,
                4.516939712967964482e-01,
                4.529234386685164204e-01,
                4.541536139381715520e-01,
                4.553844791291023486e-01,
                4.566160162133150457e-01,
                4.578482071127562003e-01,
                4.590810337005964969e-01,
                4.603144778025229922e-01,
                4.615485211980407954e-01,
                4.627831456217818551e-01,
                4.640183327648225164e-01,
                4.652540642760094602e-01,
                4.664903217632921373e-01,
                4.677270867950641398e-01,
                4.689643409015100128e-01,
                4.702020655759606704e-01,
                4.714402422762545286e-01,
                4.726788524261060220e-01,
                4.739178774164788388e-01,
                4.751572986069670934e-01,
                4.763970973271808851e-01,
                4.776372548781377403e-01,
                4.788777525336597174e-01,
                4.801185715417750077e-01,
                4.813596931261243661e-01,
                4.826010984873723153e-01,
                4.838427688046226804e-01,
                4.850846852368373430e-01,
                4.863268289242597131e-01,
                4.875691809898410334e-01,
                4.888117225406702349e-01,
                4.900544346694063469e-01,
                4.912972984557139600e-01,
                4.925402949677010200e-01,
                4.937834052633585102e-01,
                4.950266103920025751e-01,
                4.962698913957172553e-01,
                4.975132293108000536e-01,
                4.987566051692071123e-01,
                5.000000000000000000e-01,
                5.012433948307929432e-01,
                5.024867706892000019e-01,
                5.037301086042828002e-01,
                5.049733896079975359e-01,
                5.062165947366414898e-01,
                5.074597050322989800e-01,
                5.087027015442860955e-01,
                5.099455653305937641e-01,
                5.111882774593298207e-01,
                5.124308190101590776e-01,
                5.136731710757402869e-01,
                5.149153147631627681e-01,
                5.161572311953773751e-01,
                5.173989015126276847e-01,
                5.186403068738756339e-01,
                5.198814284582250478e-01,
                5.211222474663403936e-01,
                5.223627451218623152e-01,
                5.236029026728191704e-01,
                5.248427013930329066e-01,
                5.260821225835212722e-01,
                5.273211475738940335e-01,
                5.285597577237455269e-01,
                5.297979344240394406e-01,
                5.310356590984900427e-01,
                5.322729132049358602e-01,
                5.335096782367079182e-01,
                5.347459357239906508e-01,
                5.359816672351775946e-01,
                5.372168543782182004e-01,
                5.384514788019592046e-01,
                5.396855221974771188e-01,
                5.409189662994037251e-01,
                5.421517928872439107e-01,
                5.433839837866850653e-01,
                5.446155208708977069e-01,
                5.458463860618285590e-01,
                5.470765613314836351e-01,
                5.483060287032035518e-01,
                5.495347702529285172e-01,
                5.507627681104545481e-01,
                5.519900044606802503e-01,
                5.532164615448429412e-01,
                5.544421216617456682e-01,
                5.556669671689732359e-01,
                5.568909804840983524e-01,
                5.581141440858771174e-01,
                5.593364405154334085e-01,
                5.605578523774327193e-01,
                5.617783623412443417e-01,
                5.629979531420926664e-01,
                5.642166075821971605e-01,
                5.654343085318999096e-01,
                5.666510389307818363e-01,
                5.678667817887675939e-01,
                5.690815201872165829e-01,
                5.702952372800036418e-01,
                5.715079162945860825e-01,
                5.727195405330582911e-01,
                5.739300933731941168e-01,
                5.751395582694757147e-01,
                5.763479187541099558e-01,
                5.775551584380317349e-01,
                5.787612610118930689e-01,
                5.799662102470406477e-01,
                5.811699899964776206e-01,
                5.823725841958141691e-01,
                5.835739768642022041e-01,
                5.847741521052577385e-01,
                5.859730941079689392e-01,
                5.871707871475897456e-01,
                5.883672155865199427e-01,
                5.895623638751709139e-01,
                5.907562165528167375e-01,
                5.919487582484316279e-01,
                5.931399736815120560e-01,
                5.943298476628854354e-01,
                5.955183650955032659e-01,
                5.967055109752202879e-01,
                5.978912703915583160e-01,
                5.990756285284561944e-01,
                6.002585706650042097e-01,
                6.014400821761639593e-01,
                6.026201485334734542e-01,
                6.037987553057366785e-01,
                6.049758881596989379e-01,
                6.061515328607069986e-01,
                6.073256752733540154e-01,
                6.084983013621090286e-01,
                6.096693971919325827e-01,
                6.108389489288761354e-01,
                6.120069428406667011e-01,
                6.131733652972765070e-01,
                6.143382027714774374e-01,
                6.155014418393807141e-01,
                6.166630691809606990e-01,
                6.178230715805640427e-01,
                6.189814359274042888e-01,
                6.201381492160401576e-01,
                6.212931985468397311e-01,
                6.224465711264296486e-01,
                6.235982542681282048e-01,
                6.247482353923651255e-01,
                6.258965020270845780e-01,
                6.270430418081346380e-01,
                6.281878424796413229e-01,
                6.293308918943677499e-01,
                6.304721780140587484e-01,
                6.316116889097700415e-01,
                6.327494127621838826e-01,
                6.338853378619089263e-01,
                6.350194526097663328e-01,
                6.361517455170604407e-01,
                6.372822052058361164e-01,
                6.384108204091210048e-01,
                6.395375799711526810e-01,
                6.406624728475940334e-01,
                6.417854881057308836e-01,
                6.429066149246593920e-01,
                6.440258425954556998e-01,
                6.451431605213346909e-01,
                6.462585582177932109e-01,
                6.473720253127395186e-01,
                6.484835515466100819e-01,
                6.495931267724716074e-01,
                6.507007409561106348e-01,
                6.518063841761082555e-01,
                6.529100466239040612e-01,
                6.540117186038434394e-01,
                6.551113905332152099e-01,
                6.562090529422738294e-01,
                6.573046964742499387e-01,
                6.583983118853473648e-01,
                6.594898900447282442e-01,
                6.605794219344856000e-01,
                6.616668986496020421e-01,
                6.627523113978983638e-01,
                6.638356514999678160e-01,
                6.649169103890999466e-01,
                6.659960796111917825e-01,
                6.670731508246466879e-01,
                6.681481158002626763e-01,
                6.692209664211077325e-01,
                6.702916946823851418e-01,
                6.713602926912869417e-01,
                6.724267526668348260e-01,
                6.734910669397120575e-01,
                6.745532279520838337e-01,
                6.756132282574056624e-01,
                6.766710605202227446e-01,
                6.777267175159576995e-01,
                6.787801921306878539e-01,
                6.798314773609135386e-01,
                6.808805663133143948e-01,
                6.819274522044973530e-01,
                6.829721283607338433e-01,
                6.840145882176874581e-01,
                6.850548253201317328e-01,
                6.860928333216596009e-01,
                6.871286059843824123e-01,
                6.881621371786201813e-01,
                6.891934208825829522e-01,
                6.902224511820428399e-01,
                6.912492222699982980e-01,
                6.922737284463288399e-01,
                6.932959641174418763e-01,
                6.943159237959106722e-01,
                6.953336021001051970e-01,
                6.963489937538136498e-01,
                6.973620935858572123e-01,
                6.983728965296966962e-01,
                6.993813976230310869e-01,
                7.003875920073900918e-01,
                7.013914749277169181e-01,
                7.023930417319471653e-01,
                7.033922878705777260e-01,
                7.043892088962305786e-01,
                7.053838004632089209e-01,
                7.063760583270479732e-01,
                7.073659783440573445e-01,
                7.083535564708594023e-01,
                7.093387887639193146e-01,
                7.103216713790710957e-01,
                7.113022005710355478e-01,
                7.122803726929348711e-01,
                7.132561841957999471e-01,
                7.142296316280725144e-01,
                7.152007116351021265e-01,
                7.161694209586387805e-01,
                7.171357564363186610e-01,
                7.180997150011467767e-01,
                7.190612936809732680e-01,
                7.200204895979660513e-01,
                7.209772999680794658e-01,
                7.219317221005167040e-01,
                7.228837533971893548e-01,
                7.238333913521727281e-01,
                7.247806335511565168e-01,
                7.257254776708925714e-01,
                7.266679214786381236e-01,
                7.276079628315952341e-01,
                7.285455996763477993e-01,
                7.294808300482937824e-01,
                7.304136520710763358e-01,
                7.313440639560084833e-01,
                7.322720640014984683e-01,
                7.331976505924692189e-01,
                7.341208221997773675e-01,
                7.350415773796278307e-01,
                7.359599147729867230e-01,
                7.368758331049924948e-01,
                7.377893311843628510e-01,
                7.387004079028012260e-01,
                7.396090622344009269e-01,
                7.405152932350466122e-01,
                7.414191000418145494e-01,
                7.423204818723704301e-01,
                7.432194380243666298e-01,
                7.441159678748371364e-01,
                7.450100708795912574e-01,
                7.459017465726054397e-01,
                7.467909945654155335e-01,
                7.476778145465062808e-01,
                7.485622062806998045e-01,
                7.494441696085448745e-01,
                7.503237044457031768e-01,
                7.512008107823371361e-01,
                7.520754886824940755e-01,
                7.529477382834933730e-01,
                7.538175597953108431e-01,
                7.546849534999630071e-01,
                7.555499197508921405e-01,
                7.564124589723497660e-01,
                7.572725716587815903e-01,
                7.581302583742112189e-01,
                7.589855197516243157e-01,
                7.598383564923533173e-01,
                7.606887693654617033e-01,
                7.615367592071295988e-01,
                7.623823269200385999e-01,
                7.632254734727578205e-01,
                7.640661998991308268e-01,
                7.649045072976623505e-01,
                7.657403968309060005e-01,
                7.665738697248533073e-01,
                7.674049272683232115e-01,
                7.682335708123525508e-01,
                7.690598017695868815e-01,
                7.698836216136739630e-01,
                7.707050318786572429e-01,
                7.715240341583700090e-01,
                7.723406301058326484e-01,
                7.731548214326493529e-01,
                7.739666099084071549e-01,
                7.747759973600771843e-01,
                7.755829856714150372e-01,
            ];

            foreach ((HyperbolicSecantDistribution dist, ddouble[] expecteds) in new[]{
                (dist_pi2, expected_dist),
            }) {
                for ((ddouble x, int i) = (-1, 0); i < expecteds.Length; x += 1d / 256, i++) {
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