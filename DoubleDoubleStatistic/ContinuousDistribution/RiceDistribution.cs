using DoubleDouble;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class RiceDistribution : ContinuousDistribution {

        public ddouble Nu { get; }

        private QuantileBuilder quantile_lower_builder = null, quantile_upper_builder = null;

        public RiceDistribution(ddouble nu) {
            ValidateShape(nu, nu => nu >= 0d);

            Nu = nu;
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d || IsPositiveInfinity(x)) {
                return 0d;
            }

            ddouble pdf = x * Exp(-Square(x - Nu) * 0.5d) * BesselI(0, x * Nu, scale: true);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (x <= 0d) {
                return interval == Interval.Lower ? 0d : 1d;
            }

            if (x < Mode) {
                ddouble cdf_assume = PDF(x) * x;

                ddouble value = 0d;

                if (cdf_assume > 0d) {
                    ddouble eps = interval == Interval.Lower ? cdf_assume * 1e-28 : 1e-28;

                    (value, ddouble error, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(PDF, 0, x, eps, 2048);

                    Debug.WriteLine($"evals: {eval_points}");
                    Debug.WriteLine($"relative_error: {error / value}");
                }

                ddouble cdf = interval == Interval.Lower ? value : (1d - value);

                return cdf;
            }
            else {
                ddouble cdf_assume = PDF(x) / (x * x);
                ddouble value = 0d;

                if (cdf_assume > 0d) {
                    ddouble eps = interval == Interval.Lower ? 1e-28 : cdf_assume * 1e-28;

                    (value, ddouble error, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(PDF, x, PositiveInfinity, eps, 2048);

                    Debug.WriteLine($"evals: {eval_points}");
                    Debug.WriteLine($"relative_error: {error / value}");
                }

                ddouble cdf = interval == Interval.Lower ? (1d - value) : value;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p > 0.5) {
                return (interval == Interval.Lower) ? Quantile(1d - p, Interval.Upper) : Quantile(1d - p, Interval.Lower);
            }

            if (interval == Interval.Lower) {
                this.quantile_lower_builder ??= new QuantileBuilder(0d, 2d * Nu + 1d, x => CDF(x, Interval.Lower));

                (ddouble x, ddouble x0, ddouble x1) = quantile_lower_builder.Estimate(p);

                for (int i = 0; i < 8; i++) {
                    ddouble y = CDF(x, Interval.Lower), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x - dx, x0, x1);

                    if (Abs(dx / x) < 1e-32) {
                        break;
                    }
                }

                return x;
            }
            else {
                if (p <= 0d) {
                    return PositiveInfinity;
                }

                this.quantile_upper_builder ??= new QuantileBuilder(2d * Nu + 32d, 0d, x => CDF(x, Interval.Upper));

                (ddouble x, ddouble x0, ddouble x1) = quantile_upper_builder.Estimate(p);

                for (int i = 0; i < 8; i++) {
                    ddouble y = CDF(x, Interval.Upper), dx = (y - p) / PDF(x);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x = Clamp(x + dx, x1, x0);

                    if (Abs(dx / x) < 1e-32) {
                        break;
                    }
                }

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Median => Quantile(0.5d);

        private ddouble? mode = null;
        public override ddouble Mode {
            get {
                if (mode is not null) {
                    return mode.Value;
                }

                ddouble x = Nu + 1d / (Nu + 1d);

                for (int i = 0; i < 32; i++) {
                    ddouble b0 = BesselI(0, x * Nu, scale: true), b1 = BesselI(1, x * Nu, scale: true);

                    ddouble dx = (x * Nu * b1 - (x * x - 1d) * b0) / (x * (x * x + Nu * Nu - 3d) * b0 + Nu * (1d - 2d * x * x) * b1);

                    if (!IsFinite(dx)) {
                        break;
                    }

                    x -= dx;

                    if (Abs(dx / x) < 1e-30 || Abs(dx) < Epsilon) {
                        break;
                    }
                }

                mode ??= x;

                return x;
            }
        }

        private ddouble? mean = null;
        public override ddouble Mean => mean ??=
            IntegrationStatistics.Mean(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? variance = null;
        public override ddouble Variance => variance ??=
            IntegrationStatistics.Variance(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? skewness = null;
        public override ddouble Skewness => skewness ??=
            IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => kurtosis ??=
            IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 2048);

        private ddouble? entropy = null;
        public override ddouble Entropy => entropy ??=
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384);

        public override string ToString() {
            return $"{typeof(RiceDistribution).Name}[nu={Nu}]";
        }

        public override string Formula => "p(x; nu) := x * exp(-(x^2 + v^2) / 2) * bessel_i(0, x * nu)";
    }
}
