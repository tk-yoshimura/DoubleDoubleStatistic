using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class MaxwellDistribution : ScalableDistribution<MaxwellDistribution>,
        IMultiplyOperators<MaxwellDistribution, ddouble, MaxwellDistribution>,
        IDivisionOperators<MaxwellDistribution, ddouble, MaxwellDistribution>,
        IFittableContinuousDistribution<MaxwellDistribution> {

        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv, sigma_sq, pdf_norm;

        private readonly GammaDistribution randam_gen_gamma_dist;

        public MaxwellDistribution() : this(sigma: 1d) { }

        public MaxwellDistribution(ddouble sigma) {
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Sigma = sigma;
            sigma_inv = 1d / sigma;
            sigma_sq = sigma * sigma;
            pdf_norm = Sqrt(2d * RcpPi) * sigma_inv;

            randam_gen_gamma_dist = new(kappa: 1.5d, 2d * sigma_sq);
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble u2 = u * u;

            ddouble pdf = pdf_norm * u2 * Exp(u2 * -0.5d);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * sigma_inv;
            ddouble u2 = u * u;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (u <= 0d) {
                    return 0d;
                }
                if (IsPositiveInfinity(u2)) {
                    return 1d;
                }

                ddouble cdf = LowerIncompleteGammaRegularized(1.5d, u2 * 0.5d);

                if (IsNaN(cdf)) {
                    return 1d;
                }

                return cdf;
            }
            else {
                if (u <= 0d) {
                    return 1d;
                }
                if (IsPositiveInfinity(u2)) {
                    return 0d;
                }

                ddouble cdf = UpperIncompleteGammaRegularized(1.5d, u2 * 0.5d);

                if (IsNaN(cdf)) {
                    return 0d;
                }

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Sqrt(InverseLowerIncompleteGamma(1.5d, p) * 2d) * Sigma;

                return x;
            }
            else {
                ddouble x = Sqrt(InverseUpperIncompleteGamma(1.5d, p) * 2d) * Sigma;

                return x;
            }
        }

        public override double Sample(Random random) {
            return double.Sqrt(randam_gen_gamma_dist.Sample(random));
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean =>
            2d * Sigma * Sqrt(2d * RcpPi);

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => Sqrt2 * Sigma;

        public override ddouble Variance =>
            sigma_sq * (3d - 8d * RcpPi);

        public override ddouble Skewness =>
            2 * Sqrt2 * (16d - Pi * 5d) / ExMath.Pow3d2(3d * Pi - 8d);

        public override ddouble Kurtosis =>
            (-384d + Pi * (160d + Pi * -12d)) / Square(3d * Pi - 8d);

        public override ddouble Entropy =>
            -0.5d + Log(Sigma * Sqrt2 * Sqrt(Pi)) + EulerGamma;

        public static MaxwellDistribution operator *(MaxwellDistribution dist, ddouble k) {
            return new(dist.Sigma * k);
        }

        public static MaxwellDistribution operator /(MaxwellDistribution dist, ddouble k) {
            return new(dist.Sigma / k);
        }

        public static (MaxwellDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (MaxwellDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileScaleFitter<MaxwellDistribution>.Fit(new MaxwellDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(MaxwellDistribution).Name}[sigma={Sigma}]";
        }

        public override string Formula => "p(x; sigma) := sqrt(2 / pi) * u^2 * Exp(-u^2 / 2) / sigma, u = x / sigma";
    }
}
