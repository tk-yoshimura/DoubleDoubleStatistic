using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class NormalDistribution : StableDistribution<NormalDistribution>,
        IAdditionOperators<NormalDistribution, NormalDistribution, NormalDistribution>,
        ISubtractionOperators<NormalDistribution, NormalDistribution, NormalDistribution>,
        IAdditionOperators<NormalDistribution, ddouble, NormalDistribution>,
        ISubtractionOperators<NormalDistribution, ddouble, NormalDistribution>,
        IMultiplyOperators<NormalDistribution, ddouble, NormalDistribution>,
        IDivisionOperators<NormalDistribution, ddouble, NormalDistribution> {

        public override ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_sq, exp_scale, erf_scale;

        public NormalDistribution() : this(mu: 0d, sigma: 1d) { }

        public NormalDistribution(ddouble sigma) : this(mu: 0d, sigma: sigma) { }

        public NormalDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;

            sigma_sq = sigma * sigma;
            pdf_norm = 1d / (sigma * Sqrt(2d * PI));
            exp_scale = -1d / (2 * sigma_sq);
            erf_scale = -1d / (Sqrt2 * sigma);
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = pdf_norm * Exp(Square(x - Mu) * exp_scale);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Lower) {
                ddouble cdf = Erfc((x - Mu) * erf_scale) * 0.5d;

                return cdf;
            }
            else {
                ddouble cdf = Erfc((Mu - x) * erf_scale) * 0.5d;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Mu - Sigma * Sqrt2 * InverseErfc(2d * p);

                return x;
            }
            else {
                ddouble x = Mu + Sigma * Sqrt2 * InverseErfc(2d * p);

                return x;
            }
        }

        public override bool Symmetric => true;

        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance => sigma_sq;

        public override ddouble Skewness => 0d;

        public override ddouble Kurtosis => 0d;

        public override ddouble Entropy => Log(Sigma * Sqrt(2d * PI * E));

        public override ddouble Alpha => 2d;

        public override ddouble Beta => 0d;

        public override ddouble C => Sigma / Sqrt2;

        public static NormalDistribution operator +(NormalDistribution dist1, NormalDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static NormalDistribution operator -(NormalDistribution dist1, NormalDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, Hypot(dist1.Sigma, dist2.Sigma));
        }

        public static NormalDistribution operator +(NormalDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.Sigma);
        }

        public static NormalDistribution operator -(NormalDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.Sigma);
        }

        public static NormalDistribution operator *(NormalDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Sigma * k);
        }

        public static NormalDistribution operator /(NormalDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(NormalDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-u^2) / (sqrt(2 * pi) * sigma), u = (x - mu) / sigma";
    }
}
