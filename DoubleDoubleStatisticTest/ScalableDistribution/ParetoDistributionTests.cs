using DoubleDouble;
using DoubleDoubleStatistic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.ScalableDistribution {
    [TestClass()]
    public class ParetoDistributionTests {
        readonly ParetoDistribution dist_k1alpha2 = new(k: 1, alpha: 2);
        readonly ParetoDistribution dist_k2alpha3 = new(k: 2, alpha: 3);
        readonly ParetoDistribution dist_k2alpha4 = new(k: 2, alpha: 4);
        readonly ParetoDistribution dist_k4alpha3 = new(k: 4, alpha: 3);
        readonly ParetoDistribution dist_k5alpha6 = new(k: 5, alpha: 6);

        ParetoDistribution[] Dists => [
            dist_k1alpha2,
            dist_k2alpha3,
            dist_k2alpha4,
            dist_k4alpha3,
            dist_k5alpha6,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Support={dist.Support}");
                Console.WriteLine($"K={dist.K}");
                Console.WriteLine($"Alpha={dist.Alpha}");
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
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);

                Assert.IsTrue(ddouble.Abs(dist.CDF(dist.Median) - 0.5) < 1e-20, $"{dist}\n{dist.Median}");
            }
        }

        [TestMethod()]
        public void VarianceTest() {
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);

                ddouble actual = dist.Entropy;
                ddouble expected = IntegrationStatistics.Entropy(dist, eps: 1e-28, discontinue_eval_points: 65536);
                Assert.IsTrue(ddouble.Abs(actual - expected) < 1e-20, $"{dist}\n{expected}\n{actual}");
            }
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble pdf = dist.PDF(x);

                    Console.WriteLine($"pdf({x})={pdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFLowerTest() {
            foreach (ParetoDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (ddouble x = -1; x <= 5; x += 0.125) {
                    ddouble cdf = dist.CDF(x, Interval.Lower);

                    Console.WriteLine($"cdf({x})={cdf}");
                }
            }
        }

        [TestMethod()]
        public void CDFUpperTest() {
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
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
            foreach (ParetoDistribution dist in Dists) {
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
            ddouble[] expected_dist_k1alpha2 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                2.0,
                1.404663923182442,
                1.024,
                0.7693463561232157,
                0.5925925925925926,
                0.4660901228948566,
                0.3731778425655977,
                0.3034074074074074,
                0.25,
                0.2084266232444535,
                0.1755829903978052,
                0.1492928998396268,
                0.128,
                0.1105712126120289,
                0.09616829451540196,
                0.08416207775129449,
                0.07407407407407407,
                0.065536,
                0.05826126536185708,
                0.05202458974749784,
                0.04664723032069971,
                0.04198614129320595,
                0.03792592592592593,
                0.03437279715350273,
                0.03125,
                0.02849430948604502,
                0.02605332790555669,
                0.02388338192419825,
                0.02194787379972565,
                0.02021597931020867,
                0.01866161247995335,
                0.01726259714425395,
                0.016,
                0.01485759057471598,
                0.01382140157650362,
                0.01287936911215365,
                0.01202103681442524,
                0.01123731138545953,
                0.01052025971891181,
                0.009862939811024532,
                0.009259259259259259,
                0.008703856386369625,
                0.008192,
                0.007719504564609388,
                0.007282658170232135,
                0.006878161166600617,
                0.00650307371843723,
                0.006154770848985725,
                0.005830903790087463,
                0.005529366660726918,
                0.005248267661650744,
                0.004985904108988748,
                0.004740740740740741,
                0.004511390821258167,
                0.004296599644187842,
                0.004095230096741813,
                0.00390625,
            ];
            ddouble[] expected_dist_k2alpha3 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1.5,
                1.176997401851031,
                0.9364426154549612,
                0.7543220202423248,
                0.6144000000000001,
                0.5054684005121323,
                0.4196434669762994,
                0.3512851940923596,
                0.2962962962962963,
                0.25165824,
                0.2151185182591646,
                0.1849763191022145,
                0.1599333610995419,
                0.138988605660268,
                0.121362962962963,
                0.1064447911850407,
                0.09375,
                0.08289253668667643,
                0.07356233761568946,
                0.06550870470637235,
                0.05852766346593508,
                0.05245227064270359,
                0.0471451262651453,
                0.04249254681662511,
                0.0384,
                0.03478850476031059,
                0.03159177503200827,
                0.02875394034341279,
                0.02622771668601871,
                0.023972930955647,
                0.02195532463077247,
                0.02014557918847564,
                0.01851851851851852,
                0.01705245332839763,
                0.01572864,
                0.01453083212161767,
                0.01344490739119779,
                0.01245855607535206,
                0.01156101994388841,
                0.01074287275459327,
                0.009995835068721367,
                0.009312617533855861,
                0.00868678785376675,
                0.008112657533269826,
                0.007585185185185185,
                0.007099893751488263,
                0.006652799449065045,
                0.006240350623606571,
                0.005859375,
            ];
            ddouble[] expected_dist_k2alpha4 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                2.0,
                1.477016347420902,
                1.109857914613287,
                0.8469580578159436,
                0.65536,
                0.5134917084567693,
                0.4069269982800479,
                0.3258297452450872,
                0.2633744855967078,
                0.2147483648,
                0.1765075021613658,
                0.1461541286733547,
                0.1218539894091747,
                0.1022444915201971,
                0.08630255144032922,
                0.07325232941766244,
                0.0625,
                0.05358709442371001,
                0.04615676085690319,
                0.03992911524959838,
                0.03468305983166522,
                0.0302427506408381,
                0.02646743930674824,
                0.02324378629285476,
                0.02048,
                0.01810133581024291,
                0.01604661588927404,
                0.01426552079053038,
                0.0127164686962515,
                0.01136494504564006,
                0.01018217953890897,
                0.009144092681293908,
                0.00823045267489712,
                0.007424197367465636,
                0.0067108864,
                0.006078256573748569,
                0.005515859442542682,
                0.005014764709575672,
                0.004567316521042334,
                0.004166932462387691,
                0.003807937169036711,
                0.003485424106238451,
                0.003195140360006161,
                0.002933390294515643,
                0.002696954732510288,
                0.002483022951340157,
                0.002289135294301951,
                0.002113134602702754,
                0.001953125,
            ];
            ddouble[] expected_dist_k4alpha3 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.75,
                0.6631402934934114,
                0.5884987009255157,
                0.5240696376509788,
                0.4682213077274806,
                0.4196181651416288,
                0.3771610101211624,
                0.3399403745330009,
                0.3072,
                0.2783080380824847,
                0.2527342002560661,
                0.2300315227473023,
                0.2098217334881497,
                0.191783447645176,
                0.1756425970461798,
                0.1611646335078051,
                0.1481481481481481,
                0.136419626627181,
                0.12582912,
                0.1162466569729414,
                0.1075592591295823,
                0.09966844860281648,
                0.09248815955110727,
                0.08594298203674612,
                0.07996668054977094,
                0.07450094027084689,
                0.069494302830134,
                0.06490126026615861,
                0.06068148148148148,
                0.0567991500119061,
                0.05322239559252036,
                0.04992280498885257,
                0.046875,
            ];

            foreach ((ParetoDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1alpha2, expected_dist_k1alpha2),
                (dist_k2alpha3, expected_dist_k2alpha3),
                (dist_k2alpha4, expected_dist_k2alpha4),
                (dist_k4alpha3, expected_dist_k4alpha3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 8, i++) {
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
            ddouble[] expected_dist_k1alpha2 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.2098765432098766,
                0.3599999999999999,
                0.4710743801652892,
                0.5555555555555556,
                0.621301775147929,
                0.6734693877551021,
                0.7155555555555555,
                0.75,
                0.7785467128027681,
                0.8024691358024691,
                0.8227146814404432,
                0.84,
                0.854875283446712,
                0.8677685950413223,
                0.8790170132325141,
                0.8888888888888888,
                0.8976,
                0.9053254437869822,
                0.9122085048010974,
                0.9183673469387755,
                0.9239001189060643,
                0.9288888888888889,
                0.9334027055150884,
                0.9375,
                0.9412304866850322,
                0.9446366782006921,
                0.9477551020408164,
                0.9506172839506173,
                0.9532505478451424,
                0.9556786703601108,
                0.957922419460881,
                0.96,
                0.9619274241522903,
                0.963718820861678,
                0.965386695511087,
                0.9669421487603306,
                0.9683950617283951,
                0.9697542533081286,
                0.9710276143051154,
                0.9722222222222222,
                0.973344439816743,
                0.9744,
                0.9753940792003075,
                0.9763313609467456,
                0.9772160911356355,
                0.9780521262002744,
                0.9788429752066116,
                0.9795918367346939,
                0.9803016312711603,
                0.9809750297265161,
                0.981614478598104,
                0.9822222222222222,
                0.9828003224939532,
                0.9833506763787722,
                0.9838750314940791,
                0.984375,
            ];
            ddouble[] expected_dist_k2alpha3 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.166293507022186,
                0.2976680384087793,
                0.4028284006414931,
                0.4879999999999999,
                0.5577151495518844,
                0.6153268219383922,
                0.6633516889948221,
                0.7037037037037037,
                0.737856,
                0.7669549385525717,
                0.7919016410100087,
                0.8134110787172012,
                0.8320554348271763,
                0.8482962962962963,
                0.8625088113859891,
                0.875,
                0.8860227620558199,
                0.8957866883777732,
                0.904466472303207,
                0.9122085048010974,
                0.9191360827591653,
                0.9253535500801866,
                0.9309496114229843,
                0.9359999999999999,
                0.9405696377011361,
                0.9447143936939856,
                0.9484825235513854,
                0.9519158527422991,
                0.9550507544581619,
                0.9579189611243528,
                0.9605482407559018,
                0.962962962962963,
                0.9651845744545215,
                0.967232,
                0.9691219817415625,
                0.9708693673190715,
                0.9724873553335975,
                0.9739877051262511,
                0.9753809166040571,
                0.9766763848396501,
                0.9778825333570923,
                0.9790069293533971,
                0.980056383564045,
                0.981037037037037,
                0.9819544367149673,
                0.9828136014232486,
                0.9836190796130327,
                0.984375,
            ];
            ddouble[] expected_dist_k2alpha4 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.2153350654326457,
                0.375704923030026,
                0.4971186531717836,
                0.5903999999999999,
                0.6630210663252452,
                0.7202376886824671,
                0.7658098706050936,
                0.8024691358024691,
                0.83222784,
                0.8565876544938902,
                0.876682453931857,
                0.8933777592669722,
                0.9073409295598214,
                0.9190913580246913,
                0.9290368058766395,
                0.9375,
                0.944738308875549,
                0.9509584415895403,
                0.9563275301957518,
                0.9609815576893767,
                0.965031819571531,
                0.9685699158232365,
                0.9716716354555832,
                0.9744,
                0.9768076634931263,
                0.9789388166453278,
                0.9808307064377249,
                0.9825148555426542,
                0.9840180460295687,
                0.9853631169128183,
                0.9865696138743496,
                0.9876543209876544,
                0.9886316977810683,
                0.98951424,
                0.9903127785855882,
                0.9910367284058681,
                0.9916942959497653,
                0.992292653370741,
                0.9928380848302711,
                0.9933361099541858,
                0.9937915883107628,
                0.9942088080974888,
                0.9945915616444868,
                0.9949432098765432,
                0.9952667374990078,
                0.99556480036729,
                0.995839766250929,
                0.99609375,
            ];
            ddouble[] expected_dist_k4alpha3 = [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.0,
                0.08818209644655917,
                0.166293507022186,
                0.2357317784256561,
                0.2976680384087793,
                0.3530886620733223,
                0.4028284006414931,
                0.4475968913838737,
                0.4879999999999999,
                0.5245571016090886,
                0.5577151495518844,
                0.5878601884110833,
                0.6153268219383922,
                0.6404060356652949,
                0.6633516889948221,
                0.6843859260472149,
                0.7037037037037037,
                0.7214765956361721,
                0.737856,
                0.7529758539324996,
                0.7669549385525717,
                0.7798988426687803,
                0.7919016410100087,
                0.8030473328324568,
                0.8134110787172012,
                0.8230602668567386,
                0.8320554348271763,
                0.8404510685123601,
                0.8482962962962963,
                0.8556354937197386,
                0.8625088113859891,
                0.868952636904262,
                0.875,
            ];

            foreach ((ParetoDistribution dist, ddouble[] expecteds) in new[]{
                (dist_k1alpha2, expected_dist_k1alpha2),
                (dist_k2alpha3, expected_dist_k2alpha3),
                (dist_k2alpha4, expected_dist_k2alpha4),
                (dist_k4alpha3, expected_dist_k4alpha3),
            }) {
                for ((ddouble x, int i) = (0, 0); i < expecteds.Length; x += 1d / 8, i++) {
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