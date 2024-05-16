using Algebra;
using DoubleDouble;
using DoubleDoubleStatistic.ContinuousDistributions;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.MultiVariateDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class DirichletDistribution :
        MultiVariateDistribution,
        IFittableMultiVariateDistribution<DirichletDistribution> {

        public int Dim { get; }

        public ReadOnlyCollection<ddouble> Alpha { get; }

        private readonly ddouble pdf_norm, alpha_sum;
        private readonly ReadOnlyCollection<GammaDistribution> random_gen_gamma_dists;

        public DirichletDistribution(params ddouble[] alpha)
            : this(new ReadOnlyCollection<ddouble>(alpha)) { }

        public DirichletDistribution(IEnumerable<ddouble> alpha)
            : this(new ReadOnlyCollection<ddouble>(alpha.ToArray())) { }

        public DirichletDistribution(ReadOnlyCollection<ddouble> alpha) {
            ddouble alpha_sum = alpha.Sum();

            ParamAssert.ValidateLocation(nameof(alpha), alpha.Count >= 2 && alpha.All(ParamAssert.IsFinitePositive) && IsFinite(alpha_sum));

            Dim = alpha.Count;
            Alpha = alpha;

            pdf_norm = alpha.Select(LogGamma).Sum() - LogGamma(alpha_sum);
            this.alpha_sum = alpha_sum;
            random_gen_gamma_dists = new ReadOnlyCollection<GammaDistribution>(alpha.Select(a => new GammaDistribution(a)).ToArray());
        }

        public override ddouble PDF(Vector x) {
            if (x.Dim != Dim) {
                throw new ArgumentException("Invalid argument dim.", nameof(x));
            }

            ddouble u = 0d;

            for (int i = 0; i < Dim; i++) {
                if (!(x[i] >= 0d && x[i] <= 1d)) {
                    return NaN;
                }

                if (x[i] == 0d) {
                    if (Alpha[i] < 1d) {
                        return PositiveInfinity;
                    }
                    if (Alpha[i] > 1d) {
                        return 0d;
                    }

                    continue;
                }

                u += Log(x[i]) * ((double)Alpha[i] - 1d);
            }

            ddouble pdf = Exp(u - pdf_norm);

            return pdf;
        }

        public override double[] Sample(Random random) {
            double[] r = new double[Dim];
            double rs = 0;

            for (int i = 0; i < Dim; i++) {
                r[i] = random_gen_gamma_dists[i].Sample(random);
                rs += r[i];
            }

            rs = Math.Max(rs, double.Epsilon);

            double[] v = new double[Dim];
            double vs = 1d;

            for (int i = 0; i < Dim - 1; i++) {
                v[i] = r[i] / rs;
                vs -= v[i];
            }

            v[Dim - 1] = vs;

            return v;
        }

        public override Vector Mean =>
            Alpha.Select(a => a / alpha_sum).ToArray();

        public override Vector Mode => Alpha.All(a => a > 1d)
            ? Alpha.Select(a => (a - 1d) / (alpha_sum - Dim)).ToArray()
            : Vector.Invalid(Dim);

        public override Matrix Covariance {
            get {
                Vector mean = Mean;
                Matrix cov = Matrix.Zero(Dim);

                for (int i = 0; i < Dim; i++) {
                    for (int j = 0; j < Dim; j++) {
                        cov[i, j] = (i == j) ? (mean[i] * (1d - mean[j])) : (-mean[i] * mean[j]);
                        cov[i, j] /= alpha_sum + 1d;
                    }
                }

                return cov;
            }
        }

        public override ddouble Entropy =>
            pdf_norm + (alpha_sum - Dim) * Digamma(alpha_sum) - Alpha.Select(a => (a - 1d) * Digamma(a)).Sum();

        public static DirichletDistribution? Fit(IEnumerable<double[]> samples)
            => Fit(samples.Select(v => new Vector(v)));

        public static DirichletDistribution? Fit(IEnumerable<ddouble[]> samples)
            => Fit(samples.Select(v => new Vector(v)));

        public static DirichletDistribution? Fit(IEnumerable<Vector> samples) {
            if (!samples.Any()) {
                return null;
            }

            int dim = samples.First().Dim;

            if (dim < 1 || samples.Any(v => v.Dim != dim)) {
                throw new ArgumentException("Invalid argument: mismatch dim", nameof(samples));
            }

            ddouble[] alpha = new ddouble[dim];

            for (int i = 0; i < dim; i++) {
                (ddouble mu, ddouble var) = samples.Select(v => v[i]).MeanVariance();
                alpha[i] = mu * (mu * (1d - mu) - var) / var;
            }

            try {
                return new DirichletDistribution(alpha);
            }
            catch (ArgumentOutOfRangeException) {
                return null;
            }
        }

        public override string ToString() {
            return $"{typeof(DirichletDistribution).Name}[alpha=[{string.Join(',', Alpha)}]]";
        }

        public override string Formula => "p(x; alpha) := prod(x_i * (alpha_i - 1)) / (prod(gamma(alpha_i)) / gamma(sum(alpha)))";
    }
}
