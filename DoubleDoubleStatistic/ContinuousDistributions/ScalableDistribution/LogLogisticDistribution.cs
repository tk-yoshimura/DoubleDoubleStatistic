using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LogLogisticDistribution : ScalableDistribution<LogLogisticDistribution>,
        IMultiplyOperators<LogLogisticDistribution, ddouble, LogLogisticDistribution>,
        IDivisionOperators<LogLogisticDistribution, ddouble, LogLogisticDistribution> {

        public ddouble Gamma { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv, gamma_inv, c;

        public LogLogisticDistribution(ddouble gamma) : this(gamma, sigma: 1d) { }
        public LogLogisticDistribution(ddouble gamma, ddouble sigma) {
            ValidateScale(sigma);
            ValidateShape(gamma, IsFinite);

            Sigma = sigma;
            Gamma = gamma;

            sigma_inv = 1d / sigma;
            gamma_inv = 1d / gamma;
            c = gamma / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * sigma_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(u)) {
                return 0d;
            }

            ddouble v = Pow(u, Gamma);

            if (v == 0d) {
                return Gamma < 1d ? PositiveInfinity : Gamma == 1d ? c : 0d;
            }

            ddouble pdf = c * v / (u * Square(1d + v));
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * sigma_inv;
            ddouble v = Pow(u, -Gamma), vp1 = v + 1d;

            if (interval == Interval.Lower) {
                if (u <= 0d || IsPositiveInfinity(vp1)) {
                    return 0d;
                }
                if (IsPositiveInfinity(u)) {
                    return 1d;
                }

                ddouble cdf = 1d / vp1;

                return cdf;
            }
            else {
                if (u <= 0d || IsPositiveInfinity(vp1)) {
                    return 1d;
                }
                if (IsPositiveInfinity(u)) {
                    return 0d;
                }

                ddouble cdf = v / vp1;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p == 0d) {
                    return 0d;
                }
                if (p == 1d) {
                    return PositiveInfinity;
                }

                ddouble x = Sigma * Pow(p / (1d - p), gamma_inv);

                return x;
            }
            else {
                if (p == 0d) {
                    return PositiveInfinity;
                }
                if (p == 1d) {
                    return 0d;
                }

                ddouble x = Sigma * Pow((1d - p) / p, gamma_inv);

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Gamma > 1d
            ? Sigma * PI / (Gamma * SinPI(gamma_inv))
            : NaN;

        public override ddouble Median => Sigma;

        public override ddouble Mode => Gamma > 1d
            ? Sigma * Pow((Gamma - 1d) / (Gamma + 1d), gamma_inv)
            : 0d;

        public override ddouble Variance => Gamma > 2d
            ? Square(Sigma) * PI * (2d / (Gamma * SinPI(2d * gamma_inv)) - PI / Square(Gamma * SinPI(gamma_inv)))
            : NaN;

        public override ddouble Skewness {
            get {
                if (Gamma <= 3d) {
                    return NaN;
                }

                ddouble csc1 = 1d / SinPI(gamma_inv), csc2 = 1d / SinPI(2d * gamma_inv), csc3 = 1d / SinPI(3d * gamma_inv);

                return (3d * Square(Gamma) * csc3 + 2d * Square(PI) * Cube(csc1) - 6d * PI * Gamma * csc1 * csc2) /
                    (Sqrt(PI) * ExMath.Pow3d2(2d * Gamma * csc2 - PI * Square(csc1)));
            }
        }

        public override ddouble Kurtosis {
            get {
                if (Gamma <= 4d) {
                    return NaN;
                }

                ddouble csc1 = 1d / SinPI(gamma_inv), csc2 = 1d / SinPI(2d * gamma_inv);
                ddouble csc3 = 1d / SinPI(3d * gamma_inv), csc4 = 1d / SinPI(4d * gamma_inv);

                return (4d * Cube(Gamma) * csc4 - 3d * PI * csc1 * (4d * Square(Gamma) * csc3 + Square(PI) * Cube(csc1) - 4d * PI * Gamma * csc1 * csc2)) /
                    (PI * Square(2d * Gamma * csc2 - PI * Square(csc1))) - 3d;
            }
        }

        public override ddouble Entropy => Log(Sigma / Gamma) + 2d;

        public static LogLogisticDistribution operator *(LogLogisticDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Sigma * k);
        }

        public static LogLogisticDistribution operator /(LogLogisticDistribution dist, ddouble k) {
            return new(dist.Gamma, dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(LogLogisticDistribution).Name}[gamma={Gamma},sigma={Sigma}]";
        }

        public override string Formula => "p(x; gamma, sigma) := (u^(gamma - 1) * gamma) / (u^gamma + 1)^2 / sigma, u = x / sigma";
    }
}
