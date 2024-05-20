using Algebra;

namespace DoubleDoubleStatistic.SampleStatistic {
    public static partial class EnumerableExtension {

        public static (Vector mu, Matrix cov) Covariance(this IEnumerable<double[]> xs) {
            return Covariance(xs.Select(x => new Vector(x)));
        }

        public static (Vector mu, Matrix cov) Covariance(this IEnumerable<Vector> xs) {
            if (!xs.Any()) {
                throw new ArgumentException("Invalid vectors", nameof(xs));
            }

            int dim = xs.First().Dim;

            if (dim < 1 || xs.Any(x => x.Dim != dim)) {
                throw new ArgumentException("Invalid vector dim", nameof(xs));
            }

            Vector mu = Vector.Zero(dim);
            Matrix cov = Matrix.Zero(dim);

            for (int i = 0; i < dim; i++) {
                mu[i] = xs.Select(x => x[i]).Mean();
            }
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j <= i; j++) {
                    cov[i, j] = cov[j, i] = xs.Select(x => x[i] * x[j]).Mean() - mu[i] * mu[j];
                }
            }

            return (mu, cov);
        }
    }
}
