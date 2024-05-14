using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.SampleStatistic;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class StudentTDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<StudentTDistribution> {

        public ddouble Nu { get; }

        private readonly ddouble pdf_norm, nu_inv, nu_half, p;
        private readonly bool is_integer_nu;
        private readonly int n;
        private readonly double zero_thr;

        private readonly ChiSquareDistribution randam_gen_chisq_dist;

        public StudentTDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu > 0d);

            Nu = nu;

            ddouble c = Sqrt(nu * PI);

            pdf_norm = nu < 70d
                ? Gamma((nu + 1d) / 2d) / (Gamma(nu / 2d) * c)
                : Exp(LogGamma((nu + 1d) / 2d) - LogGamma(nu / 2d)) / c;
            nu_inv = 1d / nu;
            nu_half = nu / 2d;
            p = -(nu + 1d) / 2d;
            is_integer_nu = nu <= 1024d && IsInteger(nu);
            n = is_integer_nu ? (int)nu : 0;

            const int zero_thr_log = 710;
            zero_thr = nu < 0.5
                ? double.PositiveInfinity
                : double.Exp((double)(((nu + 1d) * Log(nu) + 2d * zero_thr_log) / (2d * nu + 2d)));

            randam_gen_chisq_dist = new(nu);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (Abs(x) >= zero_thr) {
                return 0d;
            }

            ddouble u = 1d + x * x * nu_inv;
            ddouble v = is_integer_nu ? Pow(Sqrt(u), -(n + 1)) : Pow(u, p);
            ddouble pdf = pdf_norm * v;

            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                return CDF(-x, Interval.Lower);
            }

            if (Abs(x) >= zero_thr) {
                return IsNegative(x) ? 0d : 1d;
            }

            ddouble v = Sqrt(x * x + Nu);
            ddouble t = (x + v) / (2d * v);

            if (IsNaN(t)) {
                return IsNegative(x) ? 0d : 1d;
            }

            ddouble cdf = IncompleteBetaRegularized(t, nu_half, nu_half);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (Abs(p - 0.5d) < 1e-31d) {
                return 0d;
            }

            ddouble x;

            if (p == 0d) {
                x = NegativeInfinity;
            }
            else if (p == 1d) {
                x = PositiveInfinity;
            }
            else if (Nu == 1d) {
                x = Tan(PI * (p - 0.5d));
            }
            else if (Nu == 2d) {
                x = 2d * (p - 0.5d) * Sqrt(1d / (2d * p * (1d - p)));
            }
            else if (Nu == 4d) {
                ddouble a = 4d * p * (1d - p);
                ddouble q = Cos(Acos(Sqrt(a)) / 3d) / Sqrt(a);

                x = Sign(p - 0.5d) * 2d * Sqrt(q - 1d);
            }
            else {
                ddouble t = InverseIncompleteBeta(p, nu_half, nu_half);
                ddouble u = Sqrt(Nu / (t * (1d - t)));

                x = u * (2d * t - 1d) * 0.5d;
            }

            return interval == Interval.Lower ? x : -x;
        }

        public override double Sample(Random random) {
            double c = randam_gen_chisq_dist.Sample(random), z = random.NextGaussian();

            double r = z / double.Max(double.Sqrt(c * (double)nu_inv), double.Epsilon);

            return r;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => Nu > 1d ? 0d : NaN;

        public override ddouble Median => 0d;

        public override ddouble Mode => 0d;

        public override ddouble Variance => Nu > 2d
            ? Nu / (Nu - 2d)
            : Nu > 1d
            ? PositiveInfinity
            : NaN;

        public override ddouble Skewness => Nu > 3d ? 0d : NaN;

        public override ddouble Kurtosis => Nu > 4d
            ? 6d / (Nu - 4d)
            : Nu > 2d
            ? PositiveInfinity
            : NaN;

        public override ddouble Entropy =>
            (Nu + 1d) * 0.5d * (Digamma((Nu + 1d) * 0.5d) - Digamma(Nu * 0.5d)) + Log(Sqrt(Nu) * Beta(Nu * 0.5d, 0.5d));

        public static (StudentTDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (StudentTDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble nu = t / (1d - t);

                    try {
                        StudentTDistribution dist = new(nu);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-4d, 100d / 101d), iter: 32
            );

            try {
                ddouble nu = t / (1d - t);
                StudentTDistribution dist = new(nu);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(StudentTDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := (1 + x^2 / nu)^(- (nu + 1) / 2) * gamma((nu + 1) / 2) / (gamma(nu / 2) * sqrt(pi * nu))";
    }
}
