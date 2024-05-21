using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class CauchyDistribution : StableDistribution<CauchyDistribution>,
        IAdditionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution>,
        ISubtractionOperators<CauchyDistribution, CauchyDistribution, CauchyDistribution>,
        IAdditionOperators<CauchyDistribution, ddouble, CauchyDistribution>,
        ISubtractionOperators<CauchyDistribution, ddouble, CauchyDistribution>,
        IMultiplyOperators<CauchyDistribution, ddouble, CauchyDistribution>,
        IDivisionOperators<CauchyDistribution, ddouble, CauchyDistribution>,
        IFittableContinuousDistribution<CauchyDistribution> {

        public override ddouble Mu { get; }
        public ddouble Gamma { get; }

        private readonly ddouble pdf_norm, gamma_inv;

        public CauchyDistribution() : this(mu: 0d, gamma: 1d) { }

        public CauchyDistribution(ddouble gamma) : this(mu: 0d, gamma: gamma) { }

        public CauchyDistribution(ddouble mu, ddouble gamma) {
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(gamma), ParamAssert.IsFinitePositive(gamma));

            Mu = mu;
            Gamma = gamma;

            pdf_norm = RcpPI / gamma;
            gamma_inv = 1d / gamma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * gamma_inv;

            ddouble pdf = pdf_norm / (Square(u) + 1d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * gamma_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsInfinity(u)) {
                    return Sign(u) < 0 ? 0d : 1d;
                }

                ddouble cdf = (u > -2d)
                    ? Atan(u) * RcpPI + 0.5d
                    : -Atan(1d / u) * RcpPI;

                return cdf;
            }
            else {
                if (IsInfinity(u)) {
                    return Sign(u) < 0 ? 1d : 0d;
                }

                ddouble cdf = (u < 2d)
                    ? Atan(-u) * RcpPI + 0.5d
                    : Atan(1d / u) * RcpPI;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (p == 0d) {
                    return NegativeInfinity;
                }
                if (p == 1d) {
                    return PositiveInfinity;
                }

                ddouble x = Mu - Gamma / TanPI(p);

                return x;
            }
            else {
                if (p == 1d) {
                    return NegativeInfinity;
                }
                if (p == 0d) {
                    return PositiveInfinity;
                }

                ddouble x = Mu + Gamma / TanPI(p);

                return x;
            }
        }

        public override double Sample(Random random) {
            double u = random.NextUniformOpenInterval01() - 0.5d;
            double r = double.TanPi(u);

            double v = r * (double)Gamma + (double)Mu;

            return v;
        }

        public override bool Symmetric => true;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => Log(4d * PI * Gamma);

        public override ddouble Alpha => 1d;

        public override ddouble Beta => 0d;

        public override ddouble C => Gamma;

        public static CauchyDistribution operator +(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public static CauchyDistribution operator -(CauchyDistribution dist1, CauchyDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, dist1.Gamma + dist2.Gamma);
        }

        public static CauchyDistribution operator +(CauchyDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.Gamma);
        }

        public static CauchyDistribution operator -(CauchyDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.Gamma);
        }

        public static CauchyDistribution operator *(CauchyDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Gamma * k);
        }

        public static CauchyDistribution operator /(CauchyDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Gamma / k);
        }

        public static (CauchyDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (CauchyDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileLinearFitter<CauchyDistribution>.Fit(new CauchyDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(CauchyDistribution).Name}[mu={Mu},gamma={Gamma}]";
        }

        public override string Formula => "p(x; mu, gamma) := 1 / (1 + u^2) / (pi * gamma), u = (x - mu) / gamma";
    }
}