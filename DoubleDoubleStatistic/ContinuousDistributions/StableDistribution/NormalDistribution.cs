using DoubleDouble;
using DoubleDoubleStatistic.RandomGeneration;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class NormalDistribution : StableDistribution<NormalDistribution>,
        IAdditionOperators<NormalDistribution, NormalDistribution, NormalDistribution>,
        ISubtractionOperators<NormalDistribution, NormalDistribution, NormalDistribution>,
        IAdditionOperators<NormalDistribution, ddouble, NormalDistribution>,
        ISubtractionOperators<NormalDistribution, ddouble, NormalDistribution>,
        IMultiplyOperators<NormalDistribution, ddouble, NormalDistribution>,
        IDivisionOperators<NormalDistribution, ddouble, NormalDistribution> {

        public override ddouble Mu { get; }
        public ddouble Sigma { get; }

        private static readonly ddouble sqrt2_inv = 1d / Sqrt2;

        private readonly ddouble pdf_norm, sigma_inv;

        public NormalDistribution() : this(mu: 0d, sigma: 1d) { }

        public NormalDistribution(ddouble sigma) : this(mu: 0d, sigma: sigma) { }

        public NormalDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;

            sigma_inv = 1d / sigma;
            pdf_norm = 1d / (sigma * Sqrt(2d * PI));
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * sigma_inv;

            ddouble pdf = pdf_norm * Exp(u * u * -0.5d);

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * sigma_inv;

            if (interval == Interval.Lower) {
                ddouble cdf = Erfc(-u * sqrt2_inv) * 0.5d;

                return cdf;
            }
            else {
                ddouble cdf = Erfc(u * sqrt2_inv) * 0.5d;

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

        public override double Sample(Random random) {
            double u = random.NextGaussian();
            double v = u * (double)Sigma + (double)Mu;

            return v;
        }

        public override bool Symmetric => true;

        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance => Sigma * Sigma;

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

        public static ddouble KLDivergence(NormalDistribution dist_p, NormalDistribution dist_q) {
            ddouble d = (Log(dist_q.Sigma / dist_p.Sigma) + (dist_p.Sigma + Square(dist_p.Mu - dist_q.Mu)) / dist_q.Sigma - 1d) * 0.5d;

            return d;
        }

        public static ddouble JSDivergence(NormalDistribution dist_p, NormalDistribution dist_q) {
            ddouble d = ((dist_p.Sigma + dist_q.Sigma) * Square(dist_p.Mu - dist_q.Mu) + Square(dist_p.Sigma - dist_q.Sigma))
                / (4d * dist_p.Sigma * dist_q.Sigma);

            return d;
        }

        public override string ToString() {
            return $"{typeof(NormalDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := exp(-u^2) / (sqrt(2 * pi) * sigma), u = (x - mu) / sigma";
    }
}