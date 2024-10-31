using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.MultiVariateDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class MultiVariateNormalDistribution :
        MultiVariateDistribution,
        IFittableMultiVariateDistribution<MultiVariateNormalDistribution> {

        public int Dim { get; }
        public Vector Mu => mu.Copy();
        public override Matrix Covariance => cov_matrix.Copy();

        private readonly Vector mu;
        private readonly Matrix cov_matrix, cov_matrix_inv, random_gen_weight;
        private readonly ddouble pdf_norm;

        public MultiVariateNormalDistribution(Vector mu, Matrix cov) {
            int dim = mu.Dim;

            ParamAssert.ValidateLocation(nameof(mu), mu.Dim >= 1 && Vector.IsFinite(mu));
            ParamAssert.ValidateScale(nameof(cov), Matrix.IsSquare(cov) && cov.Shape == (dim, dim) && Matrix.IsSymmetric(cov));

            ddouble det = cov.Det;

            ParamAssert.ValidateScale($"{nameof(det)}({nameof(cov)})", ParamAssert.IsFinitePositive(det));

            Dim = dim;
            this.mu = mu;
            this.cov_matrix = cov;
            this.cov_matrix_inv = cov.Inverse;
            this.pdf_norm = 1d / Sqrt(Pow(2d * Pi, dim) * det);

            (Matrix u, Vector s, _) = Matrix.SVD(cov);
            random_gen_weight = u * Matrix.FromDiagonals(Vector.Func(x => IsPositive(x) ? Sqrt(x) : -Sqrt(-x), s));
        }

        public override ddouble PDF(Vector x) {
            if (x.Dim != Dim) {
                throw new ArgumentException("Invalid argument dim.", nameof(x));
            }

            Vector u = x - Mu;

            ddouble pdf = Exp(Vector.Dot(u * cov_matrix_inv, u) * -0.5d) * pdf_norm;

            return pdf;
        }

        public override double[] Sample(Random random) {
            double[] z = new double[Dim], v = new double[Dim];

            for (int i = 0; i < Dim; i++) {
                z[i] = random.NextGaussian();
            }

            for (int i = 0; i < Dim; i++) {
                double x = (double)mu[i];

                for (int j = 0; j < Dim; j++) {
                    x += (double)random_gen_weight[i, j] * z[j];
                }

                v[i] = x;
            }

            return v;
        }

        public override Vector Mean => Mu;

        public override Vector Mode => Mu;

        public override ddouble Entropy => -Log(pdf_norm) + Dim * 0.5d;

        public static MultiVariateNormalDistribution? Fit(IEnumerable<double[]> samples)
            => Fit(samples.Select(v => new Vector(v)));

        public static MultiVariateNormalDistribution? Fit(IEnumerable<ddouble[]> samples)
            => Fit(samples.Select(v => new Vector(v)));

        public static MultiVariateNormalDistribution? Fit(IEnumerable<Vector> samples) {
            if (!samples.Any()) {
                return null;
            }

            int dim = samples.First().Dim;

            if (dim < 1 || samples.Any(v => v.Dim != dim)) {
                throw new ArgumentException("Invalid argument: mismatch dim", nameof(samples));
            }

            (Vector mu, Matrix cov) = samples.Covariance();

            int n = samples.Count();

            cov *= (ddouble)n / (n - 1);

            try {
                return new MultiVariateNormalDistribution(mu, cov);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string ToString() {
            return $"{typeof(MultiVariateNormalDistribution).Name}[mu={Mu},sigma={cov_matrix}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-(u^T . sigma . u) / 2) / sqrt((2 * pi)^dim * det(sigma)), u = x - mu";
    }
}
