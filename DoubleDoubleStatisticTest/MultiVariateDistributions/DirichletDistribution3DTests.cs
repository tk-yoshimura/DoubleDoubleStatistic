using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.MultiVariateDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.MultiVariateDistributions {
    [TestClass()]
    public class DirichletDistribution3DTests {
        readonly DirichletDistribution dist1 = new(0.5, 0.25, 0.75);
        readonly DirichletDistribution dist2 = new(2, 3, 4);
        readonly DirichletDistribution dist3 = new(1, 4, 2);

        DirichletDistribution[] Dists => [
            dist1,
            dist2,
            dist3,
        ];

        [TestMethod()]
        public void InfoTest() {
            foreach (DirichletDistribution dist in Dists) {
                Console.WriteLine(dist);
                Console.WriteLine($"Alpha={dist.Alpha}");
                Console.WriteLine($"Mean={dist.Mean}");
                Console.WriteLine($"Mode={dist.Mode}");
                Console.WriteLine($"Covariance={dist.Covariance}");
                Console.WriteLine($"Entropy={dist.Entropy}");
                Console.WriteLine(dist.Formula);
            }
        }

        [TestMethod()]
        public void MeanTest() {
            Assert.IsTrue((new Vector(0.33333333, 0.16666667, 0.5) - dist1.Mean).Norm < 1e-5);
            Assert.IsTrue((new Vector(0.22222222, 0.33333333, 0.44444444) - dist2.Mean).Norm < 1e-5);
        }

        [TestMethod()]
        public void ModeTest() {
            Assert.IsTrue((new Vector(0.1666666, 0.3333333, 0.5) - dist2.Mode).Norm < 1e-5);
        }

        [TestMethod()]
        public void CovarianceTest() {
            Assert.IsTrue((new Matrix(new double[,]
                {{ 0.08888889,-0.02222222,-0.06666667}, {-0.02222222,0.05555556,-0.03333333}, {-0.06666667,-0.03333333, 0.1       }})
                - dist1.Covariance).Norm < 1e-5);
            Assert.IsTrue((new Matrix(new double[,]
                {{ 0.01728395, -0.00740741, -0.00987654}, {-0.00740741,  0.02222222, -0.01481481}, {-0.00987654, -0.01481481,  0.02469136}})
                - dist2.Covariance).Norm < 1e-5);
            Assert.IsTrue((new Matrix(new double[,]
                {{ 0.01530612,-0.01020408,-0.00510204}, {-0.01020408, 0.03061224, -0.02040816}, {-0.00510204, -0.02040816,  0.0255102 }})
                - dist3.Covariance).Norm < 1e-5);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(-2.294094687268076, (double)dist1.Entropy, 1e-10);
            Assert.AreEqual(-1.312553395814394, (double)dist2.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (DirichletDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (double y = 0.125; y <= 0.875; y += 0.125) {
                    for (double x = 0.125; x <= 0.875; x += 0.125) {
                        double z = 1 - x - y;

                        if (z <= 0 || z >= 1) {
                            continue;
                        }

                        ddouble pdf = dist.PDF((x, y, z));

                        Console.WriteLine($"pdf({x}, {y}, {z})={pdf}");
                    }
                }
            }
        }

        [TestMethod()]
        public void RandomGenerateAndFitTest() {
            Random random = new(1234);

            foreach (DirichletDistribution dist in Dists) {
                double[][] samples = dist.Sample(random, count: 100000).ToArray();

                DirichletDistribution? dist_fit = DirichletDistribution.Fit(samples);

                Assert.IsNotNull(dist_fit);

                Console.WriteLine(dist_fit);

                Assert.IsTrue((new Vector(dist.Alpha) - new Vector(dist_fit.Alpha)).Norm < 1e-1, $"{dist},alpha");
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist1 = [
                1.627055254324391464e+00,
                1.204155615684303404e+00,
                1.039595734978234765e+00,
                9.674528424725640230e-01,
                9.576297140341928360e-01,
                1.039595734978234765e+00,
                1.012570140636560367e+00,
                7.570727628509403839e-01,
                6.642425260681785737e-01,
                6.366197723675813824e-01,
                6.771464646593120529e-01,
                7.899219380876583152e-01,
                6.002108774380706668e-01,
                5.423517514414637475e-01,
                5.585591590298033537e-01,
                6.840924653905505748e-01,
                5.353312844635349510e-01,
                5.197978674891174933e-01,
                6.404055870177099985e-01,
                5.385147624316430903e-01,
                6.642425260681785737e-01,
            ];
            ddouble[] expected_dist2 = [
                2.768554687500006661e+00,
                3.204345703125002665e+00,
                2.460937500000003109e+00,
                1.384277343750001998e+00,
                5.126953125000003331e-01,
                7.690429687500019429e-02,
                6.408691406250006217e+00,
                6.562500000000013323e+00,
                4.152832031250003553e+00,
                1.640625000000003109e+00,
                2.563476562500006106e-01,
                7.382812500000004441e+00,
                6.229248046875007105e+00,
                2.768554687500004441e+00,
                4.614257812500009437e-01,
                5.537109375000008882e+00,
                3.281250000000006217e+00,
                6.152343750000006661e-01,
                2.563476562500004441e+00,
                6.408691406250009992e-01,
                4.614257812500009437e-01,
            ];
            ddouble[] expected_dist3 = [
                1.757812500000003331e-01,
                1.464843750000002220e-01,
                1.171875000000001943e-01,
                8.789062500000015266e-02,
                5.859375000000009021e-02,
                2.929687500000004163e-02,
                1.171875000000001110e+00,
                9.375000000000008882e-01,
                7.031250000000007772e-01,
                4.687500000000004441e-01,
                2.343750000000003886e-01,
                3.164062500000002220e+00,
                2.373046875000002220e+00,
                1.582031250000001110e+00,
                7.910156250000012212e-01,
                5.625000000000005329e+00,
                3.750000000000004441e+00,
                1.875000000000001998e+00,
                7.324218750000007105e+00,
                3.662109375000003553e+00,
                6.328125000000005329e+00,
            ];

            foreach ((DirichletDistribution dist, ddouble[] expecteds) in new[]{
                (dist1,  expected_dist1),
                (dist2,  expected_dist2),
                (dist3,  expected_dist3),
            }) {
                int i = 0;
                for (double y = 0.125; y <= 0.875; y += 0.125) {
                    for (double x = 0.125; x <= 0.875; x += 0.125) {
                        double z = 1 - x - y;

                        if (z <= 0 || z >= 1) {
                            continue;
                        }

                        ddouble expected = expecteds[i];
                        ddouble actual = dist.PDF((x, y, z));

                        Console.WriteLine($"{dist} pdf({x},{y},{z})");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        if (expected > 0) {
                            Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x},{y},{z})\n{expected}\n{actual}");
                        }
                        else {
                            Assert.AreEqual(0, actual);
                        }

                        i++;
                    }
                }
            }
        }
    }
}
