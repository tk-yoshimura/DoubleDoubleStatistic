﻿using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.MultiVariateDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.MultiVariateDistributions {
    [TestClass()]
    public class MultiVariateNormalDistribution2DTests {
        readonly MultiVariateNormalDistribution dist1 = new(
            mu: new Vector(0.5, -0.25),
            cov: new Matrix(new double[,] { { 2.0, 0.25 }, { 0.25, 0.5 } })
        );
        readonly MultiVariateNormalDistribution dist2 = new(
            mu: new Vector(0.25, 1),
            cov: new Matrix(new double[,] { { 0.25, -0.125 }, { -0.125, 0.5 } })
        );

        MultiVariateNormalDistribution[] Dists => [
            dist1,
            dist2,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (MultiVariateNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Mu={dist.Mu}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Covariance={dist.Covariance}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            Assert.AreEqual(new Vector(0.5, -0.25), dist1.Mean);
            Assert.AreEqual(new Vector(0.25, 1), dist2.Mean);
        }

        [TestMethod()]
        public void ModeTest() {
            Assert.AreEqual(new Vector(0.5, -0.25), dist1.Mode);
            Assert.AreEqual(new Vector(0.25, 1), dist2.Mode);
        }

        [TestMethod()]
        public void CovarianceTest() {
            Assert.AreEqual(new Matrix(new double[,] { { 2.0, 0.25 }, { 0.25, 0.5 } }), dist1.Covariance);
            Assert.AreEqual(new Matrix(new double[,] { { 0.25, -0.125 }, { -0.125, 0.5 } }), dist2.Covariance);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(2.80560780584056, (double)dist1.Entropy, 1e-10);
            Assert.AreEqual(1.731390599257166, (double)dist2.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (MultiVariateNormalDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (int y = -5; y <= 5; y++) {
                    for (int x = -5; x <= 5; x++) {
                        ddouble pdf = dist.PDF((x, y));

                        Console.WriteLine($"pdf({x}, {y})={pdf}");
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateAndFitTest() {
            Random random = new(1234);

            foreach (MultiVariateNormalDistribution dist in Dists) {
                double[][] samples = dist.Sample(random, count: 100000).ToArray();

                MultiVariateNormalDistribution? dist_fit = MultiVariateNormalDistribution.Fit(samples);

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);

                Assert.IsTrue((dist.Mu - dist_fit.Mu).Norm < 1e-2, $"{dist},mu");
                Assert.IsTrue((dist.Covariance - dist_fit.Covariance).Norm < 1e-1, $"{dist},cov");
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist1 = [
                1.932368868751906118e-12,
                7.836142172494426923e-12,
                1.864192935532458780e-11,
                2.601690823586798583e-11,
                2.130084287271311014e-11,
                1.023090774944493150e-11,
                2.882755603357201385e-12,
                4.765162962494956094e-13,
                4.620872296677146432e-14,
                2.628732845802659957e-15,
                8.772941091654915375e-17,
                3.861251541367038296e-09,
                2.044335786812987226e-08,
                6.349691769213479416e-08,
                1.156989274936872861e-07,
                1.236751080910180002e-07,
                7.755531040384265924e-08,
                2.853100425124338618e-08,
                6.157421061930475963e-09,
                7.795732413035028348e-10,
                5.790169411175720977e-11,
                2.522910789379786371e-12,
                9.138423117258438155e-07,
                6.316947401974237739e-06,
                2.561648489504544525e-05,
                6.094079091895388835e-05,
                8.504972500106548099e-05,
                6.963282539919765310e-05,
                3.344501516909502718e-05,
                9.423778147965979155e-06,
                1.557740050705668172e-06,
                1.510571181381264090e-07,
                8.593373340517006557e-09,
                2.561648489504544525e-05,
                2.311891219858352498e-04,
                1.224028506130906985e-03,
                3.801823448377126841e-03,
                6.927373981053436096e-03,
                7.404940948484006226e-03,
                4.643557645885823537e-03,
                1.708269391815856342e-03,
                3.686703012621685106e-04,
                4.667627872718187009e-05,
                3.466814238797386414e-06,
                8.504972500106548099e-05,
                1.002149780613473522e-03,
                6.927373981053442167e-03,
                2.809188673828856564e-02,
                6.682969202258061403e-02,
                9.326834855139018443e-02,
                7.636166524781937137e-02,
                3.667691261857235208e-02,
                1.033442759473303905e-02,
                1.708269391815859378e-03,
                1.656542445669143514e-04,
                3.344501516909514238e-05,
                5.145208532017140051e-04,
                4.643557645885831343e-03,
                2.458526975492655170e-02,
                7.636166524781937137e-02,
                1.391400258771782750e-01,
                1.487322148347996476e-01,
                9.326834855139018443e-02,
                3.431150794406873789e-02,
                7.404940948484006226e-03,
                9.375181198118088350e-04,
                1.557740050705668172e-06,
                3.128804530517681519e-05,
                3.686703012621685106e-04,
                2.548438468935527654e-03,
                1.033442759473303905e-02,
                2.458526975492655170e-02,
                3.431150794406872401e-02,
                2.809188673828855176e-02,
                1.349268211801423097e-02,
                3.801823448377130311e-03,
                6.284371892314980119e-04,
                8.593373340516976779e-09,
                2.253507395929717379e-07,
                3.466814238797386414e-06,
                3.128804530517681519e-05,
                1.656542445669140532e-04,
                5.145208532017135714e-04,
                9.375181198118079677e-04,
                1.002149780613474390e-03,
                6.284371892314974698e-04,
                2.311891219858352498e-04,
                4.989409964224299080e-05,
                5.614841220699958368e-12,
                1.922403944755600747e-10,
                3.861251541367038296e-09,
                4.549752965134142354e-08,
                3.145022921782907672e-07,
                1.275369684813981546e-06,
                3.034063323873648770e-06,
                4.234376473296074603e-06,
                3.466814238797386414e-06,
                1.665129256787986475e-06,
                4.691822869363660751e-07,
                4.345266168427482110e-16,
                1.942385446673746955e-14,
                5.093669035869256532e-13,
                7.836142172494426923e-12,
                7.072128889024692668e-11,
                3.744331603858293833e-10,
                1.162986614996817474e-09,
                2.119099775788222663e-09,
                2.265188619320252510e-09,
                1.420475059260476293e-09,
                5.225635709987135560e-10,
                3.982909093715531798e-21,
                2.324509295456534104e-19,
                7.958632601635135438e-18,
                1.598534089782276428e-16,
                1.883569391151342164e-15,
                1.302019902033303638e-14,
                5.279951063556552199e-14,
                1.256083319568730166e-13,
                1.753005487733492973e-13,
                1.435239503121878400e-13,
                6.893531416829272554e-14,
            ];
            ddouble[] expected_dist2 = [
                6.593128850798786826e-62,
                1.688516589896243915e-49,
                4.472732216802843033e-39,
                1.225446230799714353e-30,
                3.472717215971876553e-24,
                1.017884208578177567e-19,
                3.085894898759650030e-17,
                9.676486192865536625e-17,
                3.138398077083700517e-18,
                1.052815266342688138e-21,
                3.653003729392641424e-27,
                7.665853458408099812e-54,
                6.260913892260377238e-42,
                5.288941285836651235e-32,
                4.621187157579207321e-24,
                4.176304601620554900e-18,
                3.903773150036825124e-14,
                3.774250972915355571e-12,
                3.774250972915341839e-12,
                3.903773150036839008e-14,
                4.176304601620584944e-18,
                4.621187157579207321e-24,
                9.064761063726091629e-47,
                2.361005660341507045e-35,
                6.360503449043119894e-26,
                1.772310268203460932e-18,
                5.107893244255588111e-13,
                1.522641516342018295e-09,
                4.694694314641711782e-08,
                1.497168801570806037e-08,
                4.938419909593730144e-11,
                1.684841526689950193e-15,
                5.945438600301428816e-22,
                1.090132049016721114e-40,
                9.054890945602335061e-30,
                7.779311064816736172e-21,
                6.912781789537080919e-14,
                6.353577847813537094e-09,
                6.040010032726765018e-06,
                5.938964940282081471e-05,
                6.040010032726765018e-06,
                6.353577847813537094e-09,
                6.912781789537155388e-14,
                7.779311064816736172e-21,
                1.333302682557747861e-35,
                3.531801759431985529e-25,
                9.676486192865467600e-17,
                2.742160312504970183e-10,
                8.037515023253895255e-06,
                2.436713960186979045e-03,
                7.640839939554988462e-03,
                2.478172023981614386e-04,
                8.313340995593604809e-08,
                2.884519880312162921e-13,
                1.035202412168882581e-20,
                1.658461127365211465e-31,
                1.400994116581004021e-21,
                1.224111910019793965e-13,
                1.106266426437545561e-07,
                1.034075237385296935e-03,
                9.997659497023531072e-02,
                9.997659497023531072e-02,
                1.034075237385297802e-03,
                1.106266426437553502e-07,
                1.224111910019807091e-13,
                1.400994116580994053e-21,
                2.098015723733930010e-28,
                5.652013661426866926e-19,
                1.574894491988696213e-11,
                4.538930389541075205e-06,
                1.353036079732891639e-02,
                4.171757253990641123e-01,
                1.330400743861456236e-01,
                4.388334511332665582e-04,
                1.497168801570816625e-08,
                5.283182449517212733e-15,
                1.928298423487930940e-23,
                2.699224948608472917e-26,
                2.318979945234754344e-17,
                2.060671208819724720e-10,
                1.893974863173119105e-05,
                1.800501614886903731e-02,
                1.770380497349466342e-01,
                1.800501614886902343e-02,
                1.893974863173122493e-05,
                2.060671208819739196e-10,
                2.318979945234770676e-17,
                2.699224948608549255e-26,
                3.531801759432010784e-25,
                9.676486192865536625e-17,
                2.742160312504989829e-10,
                8.037515023253880009e-06,
                2.436713960186974708e-03,
                7.640839939554988462e-03,
                2.478172023981619265e-04,
                8.313340995593604809e-08,
                2.884519880312152824e-13,
                1.035202412168882581e-20,
                3.842649837879750712e-30,
                4.699811680242149883e-25,
                4.106437981820046760e-17,
                3.711110425730611567e-11,
                3.468935965821857818e-07,
                3.353841127746021834e-05,
                3.353841127746015735e-05,
                3.468935965821857818e-07,
                3.711110425730611567e-11,
                4.106437981820075726e-17,
                4.699811680242249984e-25,
                5.563517280601041988e-35,
                6.360503449043119894e-26,
                1.772310268203460932e-18,
                5.107893244255569935e-13,
                1.522641516342018295e-09,
                4.694694314641711782e-08,
                1.497168801570806037e-08,
                4.938419909593747592e-11,
                1.684841526689950193e-15,
                5.945438600301428816e-22,
                2.170013999980846474e-30,
                8.192094812896955560e-41,
            ];

            foreach ((MultiVariateNormalDistribution dist, ddouble[] expecteds) in new[]{
                (dist1,  expected_dist1),
                (dist2,  expected_dist2),
            }) {
                for (int y = -5, i = 0; y <= 5; y++) {
                    for (int x = -5; x <= 5; x++, i++) {
                        ddouble expected = expecteds[i];
                        ddouble actual = dist.PDF((x, y));

                        Console.WriteLine($"{dist} pdf({x},{y})");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        if (expected > 0) {
                            Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x},{y})\n{expected}\n{actual}");
                        }
                        else {
                            Assert.AreEqual(0, actual);
                        }
                    }
                }
            }
        }
    }
}