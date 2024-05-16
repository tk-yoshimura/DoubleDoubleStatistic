using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
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
        private readonly ReadOnlyCollection<int> random_gen_index;
        private readonly ddouble pdf_norm;

        public MultiVariateNormalDistribution(Vector mu, Matrix covariance) {
            ParamAssert.ValidateLocation(nameof(mu), mu.Dim >= 1 && Vector.IsFinite(mu));
            ParamAssert.ValidateScale(nameof(covariance), Matrix.IsSquare(covariance) && Matrix.IsSymmetric(covariance));

            int dim = mu.Dim;
            ddouble det = covariance.Det;

            ParamAssert.ValidateScale(nameof(covariance), covariance.Shape == (dim, dim) && ParamAssert.IsFinitePositive(det));

            Dim = dim;
            this.mu = mu;
            this.cov_matrix = covariance;
            this.cov_matrix_inv = covariance.Inverse;
            this.pdf_norm = 1d / Sqrt(Pow(2d * PI, dim) * det);
            (Matrix p, Matrix lower, Matrix upper) = Matrix.LU(covariance);

            random_gen_weight = lower * Matrix.FromDiagonals(Vector.Func(Sqrt, upper.Diagonals));

            int[] generation_index = new int[dim];
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    if (p[i, j] > 0) {
                        generation_index[i] = j;
                        break;
                    }
                }
            }

            this.random_gen_index = new ReadOnlyCollection<int>(generation_index);
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
                double x = (double)mu[i];

                for (int j = 0; j <= i; j++) {
                    x += (double)random_gen_weight[i, j] * z[j];
                }

                v[random_gen_index[i]] = x;
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

            Vector mu = Vector.Zero(dim);
            Matrix cov = Matrix.Zero(dim);

            for (int i = 0; i < dim; i++) {
                mu[i] = samples.Select(v => v[i]).Mean();
            }

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j <= i; j++) {
                    cov[i, j] = cov[j, i] = samples.Select(v => v[i] * v[j]).Mean() - mu[i] * mu[j];
                }
            }

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
