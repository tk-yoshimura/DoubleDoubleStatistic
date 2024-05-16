using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.MultiVariateDistributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DoubleDoubleStatisticTest.MultiVariateDistributions {
    [TestClass()]
    public class DirichletDistribution2DTests {
        readonly DirichletDistribution dist1 = new(0.5, 0.25);
        readonly DirichletDistribution dist2 = new(2, 3);
        readonly DirichletDistribution dist3 = new(1, 4);

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
            Assert.IsTrue((new Vector(0.66666667, 0.33333333) - dist1.Mean).Norm < 1e-5);
            Assert.IsTrue((new Vector(0.4, 0.6) - dist2.Mean).Norm < 1e-5);
        }

        [TestMethod()]
        public void ModeTest() {
            Assert.IsTrue((new Vector(0.33333333, 0.66666667) - dist2.Mode).Norm < 1e-5);
        }

        [TestMethod()]
        public void CovarianceTest() {
            Assert.IsTrue((new Matrix(new double[,] { { 0.12698413, -0.12698413 }, { -0.12698413, 0.12698413 } }) - dist1.Covariance).Norm < 1e-5);
            Assert.IsTrue((new Matrix(new double[,] { { 0.04, -0.04 }, { -0.04, 0.04 } }) - dist2.Covariance).Norm < 1e-5);
            Assert.IsTrue((new Matrix(new double[,] { { 0.02666667, -0.02666667 }, { -0.02666667, 0.02666667 } }) - dist3.Covariance).Norm < 1e-5);
        }

        [TestMethod()]
        public void EntropyTest() {
            Assert.AreEqual(-1.1379125471183391, (double)dist1.Entropy, 1e-10);
            Assert.AreEqual(-0.2349066497880008, (double)dist2.Entropy, 1e-10);
        }

        [TestMethod()]
        public void PDFTest() {
            foreach (DirichletDistribution dist in Dists) {
                Console.WriteLine(dist);
                for (double x = 0; x <= 1; x += 0.125) {
                    ddouble pdf = dist.PDF((x, 1 - x));

                    Console.WriteLine($"pdf({x}, {1 - x})={pdf}");
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

                Assert.IsTrue((new Vector(dist.Alpha) - new Vector(dist_fit.Alpha)).Norm < 1e-2, $"{dist},alpha");
            }
        }

        [TestMethod()]
        public void PDFExpectedTest() {
            ddouble[] expected_dist1 = [
                5.961653599113061475e-01,
                4.732186281615705647e-01,
                4.429983791482233202e-01,
                4.535396688970745882e-01,
                5.033433377019234101e-01,
                6.227907389684713690e-01,
                9.697086013278551908e-01,
            ];
            ddouble[] expected_dist2 = [
                1.148437500000000444e+00,
                1.687500000000000222e+00,
                1.757812500000000000e+00,
                1.500000000000000444e+00,
                1.054687500000000000e+00,
                5.625000000000001110e-01,
                1.640625000000000555e-01,
            ];
            ddouble[] expected_dist3 = [
                2.679687500000000444e+00,
                1.687500000000000444e+00,
                9.765625000000002220e-01,
                5.000000000000002220e-01,
                2.109375000000000555e-01,
                6.250000000000002776e-02,
                7.812500000000008674e-03,
            ];

            foreach ((DirichletDistribution dist, ddouble[] expecteds) in new[]{
                (dist1,  expected_dist1),
                (dist2,  expected_dist2),
                (dist3,  expected_dist3),
            }) {
                int i = 0;
                for (double x = 0.125; x <= 0.875; x += 0.125) {
                    ddouble expected = expecteds[i];
                    ddouble actual = dist.PDF((x, 1 - x));

                    Console.WriteLine($"{dist} pdf({x},{1 - x})");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (expected > 0) {
                        Assert.IsTrue(ddouble.Abs(expected - actual) / expected < 1e-10, $"{dist} pdf({x},{1 - x})\n{expected}\n{actual}");
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
