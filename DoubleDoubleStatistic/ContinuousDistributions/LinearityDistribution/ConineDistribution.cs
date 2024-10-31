using DoubleDouble;
using DoubleDoubleStatistic.Misc;
using DoubleDoubleStatistic.RandomGeneration;
using DoubleDoubleStatistic.Utils;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class CosineDistribution : LinearityDistribution<CosineDistribution>,
        IAdditionOperators<CosineDistribution, ddouble, CosineDistribution>,
        ISubtractionOperators<CosineDistribution, ddouble, CosineDistribution>,
        IMultiplyOperators<CosineDistribution, ddouble, CosineDistribution>,
        IDivisionOperators<CosineDistribution, ddouble, CosineDistribution>,
        IFittableContinuousDistribution<CosineDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble sigma_inv;

        public CosineDistribution() : this(mu: 0d, sigma: 1d) { }

        public CosineDistribution(ddouble sigma) : this(mu: 0d, sigma: sigma) { }

        public CosineDistribution(ddouble mu, ddouble sigma) {
            ParamAssert.ValidateLocation(nameof(mu), IsFinite(mu));
            ParamAssert.ValidateScale(nameof(sigma), ParamAssert.IsFinitePositive(sigma));

            Mu = mu;
            Sigma = sigma;

            sigma_inv = 1d / sigma;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv;

            if (Abs(u) >= 1d) {
                return 0d;
            }

            ddouble pdf = (1d + CosPi(u)) * sigma_inv * 0.5d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (interval == Interval.Lower) {
                if (u <= -1d) {
                    return 0d;
                }
                if (u >= 1d) {
                    return 1d;
                }

                ddouble cdf = (1d + u + SinPi(u) * RcpPi) * 0.5d;

                return cdf;
            }
            else {
                if (u <= -1d) {
                    return 1d;
                }
                if (u >= 1d) {
                    return 0d;
                }

                ddouble cdf = (1d - u - SinPi(u) * RcpPi) * 0.5d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (p < 1e-28) {
                return (interval == Interval.Lower) ? Mu - Sigma : Mu + Sigma;
            }
            if ((1d - p) < 1e-28) {
                return (interval == Interval.Lower) ? Mu + Sigma : Mu - Sigma;
            }
            if (Abs(p - 0.5d) < 1e-28) {
                return Mu;
            }

            ddouble u = KeplerE(Pi * p * 2d, 1d) * RcpPi - 1d;

            ddouble x = (interval == Interval.Lower) ? (Mu + Sigma * u) : (Mu - Sigma * u);

            return x;
        }

        public override double Sample(Random random) {
            double r, thr;

            do {
                r = random.NextUniform() * 2d - 1d;
                thr = (double.CosPi(r) + 1d) * 0.5d;
            } while (random.NextUniform() > thr);

            double v = r * (double)Sigma + (double)Mu;

            return v;
        }

        public override (ddouble min, ddouble max) Support => (Mu - Sigma, Mu + Sigma);

        public override bool Symmetric => true;

        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance =>
            Square(Sigma) * (1d - 6d / (Pi * Pi)) / 3d;
        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis =>
            6d * (90d - Square(Square(Pi))) / (5d * Square(Pi * Pi - 6d));

        public override ddouble Entropy => 2d * Ln2 + Log(Sigma) - 1d;

        public static CosineDistribution operator +(CosineDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.Sigma);
        }

        public static CosineDistribution operator -(CosineDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.Sigma);
        }

        public static CosineDistribution operator *(CosineDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Sigma * k);
        }

        public static CosineDistribution operator /(CosineDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Sigma / k);
        }

        public static (CosineDistribution? dist, ddouble error) Fit(IEnumerable<double> samples, (double min, double max) fitting_quantile_range, int quantile_partitions = 100)
            => Fit(samples.Select(v => (ddouble)v), fitting_quantile_range, quantile_partitions);

        public static (CosineDistribution? dist, ddouble error) Fit(IEnumerable<ddouble> samples, (ddouble min, ddouble max) fitting_quantile_range, int quantile_partitions = 100) {
            return QuantileLinearFitter<CosineDistribution>.Fit(new CosineDistribution(), samples, fitting_quantile_range, quantile_partitions);
        }

        public override string ToString() {
            return $"{typeof(CosineDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := (1 + cos(u * pi)) / (2 * sigma), u = (x - mu) / sigma";
    }
}
