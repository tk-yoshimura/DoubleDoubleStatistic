using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class GompertzDistribution : ScalableDistribution<GompertzDistribution>,
        IMultiplyOperators<GompertzDistribution, ddouble, GompertzDistribution>,
        IDivisionOperators<GompertzDistribution, ddouble, GompertzDistribution> {

        public ddouble Eta { get; }
        public ddouble Theta { get; }

        private readonly ddouble theta_inv;

        public GompertzDistribution(ddouble eta) : this(eta, theta: 1d) { }

        public GompertzDistribution(ddouble eta, ddouble theta) {
            ValidateShape(eta, eta => eta > 0d);
            ValidateScale(theta);

            Eta = eta;
            Theta = theta;

            theta_inv = 1d / theta;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (IsNegative(u) || IsPositiveInfinity(u)) {
                return 0d;
            }

            ddouble pdf = Eta * Exp(Eta * (-Expm1(u)) + u) * theta_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = x * theta_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(u)) {
                    return 0d;
                }

                ddouble cdf = -Expm1(Eta * (-Expm1(u)));

                return cdf;
            }
            else {
                if (IsNegative(u)) {
                    return 1d;
                }

                ddouble cdf = Exp(Eta * (-Expm1(u)));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Log1p(-Log1p(-p) / Eta) * Theta;

                return x;
            }
            else {
                ddouble x = Log1p(-Log(p) / Eta) * Theta;

                return x;
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => Theta * Exp(Eta) * -Ei(-Eta);

        public override ddouble Median => Quantile(0.5d);

        public override ddouble Mode => (Eta >= 1d) ? 0d : -Log(Eta) * Theta;

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
            IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 2048);

        public static GompertzDistribution operator *(GompertzDistribution dist, ddouble k) {
            return new(dist.Eta, dist.Theta * k);
        }

        public static GompertzDistribution operator /(GompertzDistribution dist, ddouble k) {
            return new(dist.Eta, dist.Theta / k);
        }

        public override string ToString() {
            return $"{typeof(GompertzDistribution).Name}[eta={Eta},theta={Theta}]";
        }

        public override string Formula => "p(x; eta, theta) := exp(u + (1 - exp(u)) * eta) * eta / theta, u = x / theta";
    }
}
