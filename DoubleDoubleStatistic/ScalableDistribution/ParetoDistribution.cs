﻿using DoubleDouble;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    [DebuggerDisplay("{ToString(),nq}")]
    public class ParetoDistribution : ScalableDistribution<ParetoDistribution>,
        IMultiplyOperators<ParetoDistribution, ddouble, ParetoDistribution>,
        IDivisionOperators<ParetoDistribution, ddouble, ParetoDistribution> {

        public ddouble K { get; }
        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm, alpha_inv;

        public ParetoDistribution(ddouble alpha) : this(k: 1d, alpha: alpha) { }

        public ParetoDistribution(ddouble k, ddouble alpha) {
            K = k;
            Alpha = alpha;

            pdf_norm = alpha / k;
            alpha_inv = 1d / alpha;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = K / x;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsNegative(x) || u > 1d) {
                return 0d;
            }

            ddouble pdf = pdf_norm * Pow(u, Alpha + 1d);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = K / x;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(x) || u > 1d) {
                    return 0d;
                }

                ddouble cdf = 1d - Pow(u, Alpha);

                return cdf;
            }
            else {
                if (IsNegative(x) || u > 1d) {
                    return 1d;
                }

                ddouble cdf = Pow(u, Alpha);

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                return Quantile(1d - p, Interval.Upper);
            }

            ddouble x = K / Pow(p, alpha_inv);

            return x;
        }

        public override (ddouble min, ddouble max) Support => (K, PositiveInfinity);

        public override ddouble Mean => (Alpha > 1d) ? (Alpha * K) / (Alpha - 1d) : PositiveInfinity;

        public override ddouble Median => K * Pow(2d, alpha_inv);

        public override ddouble Mode => K;

        public override ddouble Variance => (Alpha > 2d)
            ? (Alpha * K * K) / (Square(Alpha - 1d) * (Alpha - 2d))
            : PositiveInfinity;

        public override ddouble Skewness => (Alpha > 3d)
            ? 2d * (Alpha + 1d) / (Alpha - 3d) * Sqrt((Alpha - 2d) * alpha_inv)
            : NaN;

        public override ddouble Kurtosis => (Alpha > 4d)
            ? 6d * (-2d + Alpha * (-6d + Alpha * (1d + Alpha))) / (Alpha * (Alpha - 3d) * (Alpha - 4d))
            : NaN;

        public override ddouble Entropy => 1d + Log(K * alpha_inv) + alpha_inv;

        public static ParetoDistribution operator *(ParetoDistribution dist, ddouble k) {
            return new(dist.K * k, dist.Alpha);
        }

        public static ParetoDistribution operator /(ParetoDistribution dist, ddouble k) {
            return new(dist.K / k, dist.Alpha);
        }

        public override string ToString() {
            return $"{typeof(ParetoDistribution).Name}[k={K},alpha={Alpha}]";
        }

        public override string Formula => "p(x; k, alpha) := u^(- alpha - 1) * (alpha / k), u = x / K";
    }
}
