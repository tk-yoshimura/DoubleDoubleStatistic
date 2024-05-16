using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Optimizer;
using DoubleDoubleStatistic.SampleStatistic;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class InverseChiSquareDistribution : ContinuousDistribution,
        IFittableContinuousDistribution<InverseChiSquareDistribution> {

        public ddouble Nu { get; }

        private readonly ddouble c, pdf_lognorm;

        private readonly ChiSquareDistribution randam_gen_chisq_dist;

        public InverseChiSquareDistribution(ddouble nu) {
            ParamAssert.ValidateShape(nameof(nu), ParamAssert.IsFinitePositive(nu));

            Nu = nu;

            c = nu * 0.5d + 1d;
            pdf_lognorm = nu * 0.5d + LogGamma(nu * 0.5d) * LbE;

            randam_gen_chisq_dist = new(nu);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble pdf = Pow2(-c * Log2(x) - LbE / (2d * x) - pdf_lognorm);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = 1d / (2d * x);

            if (interval == Interval.Lower) {
                if (x <= 0d) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(Nu * 0.5d, u);

                if (IsNaN(cdf)) {
                    return x < Mean ? 0d : 1d;
                }

                return cdf;
            }
            else {
                if (x <= 0d) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(Nu * 0.5d, u);

                if (IsNaN(cdf)) {
                    return x < Mean ? 1d : 0d;
                }

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = 0.5d / InverseUpperIncompleteGamma(Nu * 0.5d, p);

                return x;
            }
            else {
                ddouble x = 0.5d / InverseLowerIncompleteGamma(Nu * 0.5d, p);

                return x;
            }
        }

        public override double Sample(Random random) {
            return 1d / randam_gen_chisq_dist.Sample(random);
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Nu > 2d
            ? 1d / (Nu - 2d)
            : NaN;

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => 1d / (Nu + 2d);

        public override ddouble Variance => Nu > 4d
            ? 2d / (Square(Nu - 2d) * (Nu - 4d))
            : NaN;

        public override ddouble Skewness => Nu > 6d
            ? 4d * Sqrt(2d * (Nu - 4d)) / (Nu - 6d)
            : NaN;

        public override ddouble Kurtosis => Nu > 8d
            ? 12d * (5d * Nu - 22d) / ((Nu - 6d) * (Nu - 8d))
            : NaN;

        public override ddouble Entropy {
            get {
                ddouble nu_half = Nu * 0.5d;

                return nu_half + Log(nu_half * Gamma(nu_half)) - (1d + nu_half) * Digamma(nu_half) - Log(Nu);
            }
        }

        public static (InverseChiSquareDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (InverseChiSquareDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            ddouble[] qs = EnumerableUtil.Linspace(fitting_quantile_range.min, fitting_quantile_range.max, quantile_partitions + 1, end_point: true).ToArray();
            ddouble[] ys = samples.Quantile(qs).ToArray();

            ddouble t = GridMinimizeSearch1D.Search(
                t => {
                    ddouble nu = t / (1d - t);

                    try {
                        InverseChiSquareDistribution dist = new(nu);
                        return EvalFitness.MeanSquaredError(dist, qs, ys);
                    }
                    catch (ArgumentOutOfRangeException) {
                        return NaN;
                    }

                }, (1e-10d, 1000d / 1001d), iter: 32
            );

            try {
                ddouble nu = t / (1d - t);
                InverseChiSquareDistribution dist = new(nu);
                ddouble error = EvalFitness.MeanSquaredError(dist, qs, ys);

                return (dist, error);
            }
            catch (ArgumentOutOfRangeException) {
                return (null, NaN);
            }
        }

        public override string ToString() {
            return $"{typeof(InverseChiSquareDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := x^(-nu / 2 - 1) * exp(-1 / (2 * x)) / (2^(nu / 2) * gamma(nu / 2))";
    }
}
