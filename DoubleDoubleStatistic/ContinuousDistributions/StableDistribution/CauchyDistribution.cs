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
                    : CDFPade.Value(-u);

                return cdf;
            }
            else {
                if (IsInfinity(u)) {
                    return Sign(u) < 0 ? 1d : 0d;
                }

                ddouble cdf = (u < 2d)
                    ? Atan(-u) * RcpPI + 0.5d
                    : CDFPade.Value(u);

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

                ddouble x = Mu + Gamma * TanPI(p - 0.5d);

                return x;
            }
            else {
                if (p == 1d) {
                    return NegativeInfinity;
                }
                if (p == 0d) {
                    return PositiveInfinity;
                }

                ddouble x = Mu - Gamma * TanPI(p - 0.5d);

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

        internal static class CDFPade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_limit = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 30, 0xD4136876A9E77C4CuL, 0x7E49F52100000000uL),(+1, 30, 0xD4136876A9E77C4CuL, 0x7E49F52100000000uL)),
                ((+1, 33, 0xC13E731F2FCCFC1AuL, 0x55B29A7740000000uL),(+1, 33, 0xCA14977976E14BF2uL, 0xDAF5AF58A0000000uL)),
                ((+1, 35, 0x9E4F11AFBB7A4F38uL, 0x78F35E7058000000uL),(+1, 35, 0xADD2CC0E72306FF9uL, 0x08465051A8000000uL)),
                ((+1, 36, 0x9A017B023AD8BCE2uL, 0x532DB7EB10000000uL),(+1, 36, 0xB265D16D245A2213uL, 0xCBDC5FE804000000uL)),
                ((+1, 36, 0xC5F2175F5831C968uL, 0x2D38CF7B64000000uL),(+1, 36, 0xF345064F03068B8FuL, 0x5BCF6B8234000000uL)),
                ((+1, 36, 0xB0EECAE780E36066uL, 0xBCF1D16998000000uL),(+1, 36, 0xE840EED68B1991CCuL, 0x7D5408EF44000000uL)),
                ((+1, 35, 0xE182F7E973E9B265uL, 0x0277E60BD8000000uL),(+1, 36, 0x9F63C70BB4CB45F0uL, 0xBA66D8F484000000uL)),
                ((+1, 34, 0xCEB1C72AC3E57700uL, 0x2810704380000000uL),(+1, 35, 0x9EECD0EA160A1EEDuL, 0x0A76F7AB68000000uL)),
                ((+1, 33, 0x87B138D333371C4DuL, 0x637DDDC760000000uL),(+1, 33, 0xE5EF235D92451C67uL, 0x4B0E2A6A60000000uL)),
                ((+1, 30, 0xFB2839B2A867E347uL, 0xB7B2E36600000000uL),(+1, 31, 0xEE73411EAAA67C45uL, 0x3164060F80000000uL)),
                ((+1, 28, 0x9EFBD3FEA5A6C524uL, 0xBF0A2E6400000000uL),(+1, 29, 0xAD03D984FC6A82A6uL, 0xE4BAE1DE00000000uL)),
                ((+1, 25, 0x8306ED31E5AD6714uL, 0x32BCD88000000000uL),(+1, 26, 0xA8CB8F87FC80E35EuL, 0x23D58B3000000000uL)),
                ((+1, 21, 0x81B3C9D8DECBF538uL, 0x806B960000000000uL),(+1, 22, 0xCFBF7593ACED669BuL, 0x3FCBBF0000000000uL)),
                ((+1, 16, 0x8612B6BF46605DA6uL, 0x422E800000000000uL),(+1, 18, 0x911F06064B06B32FuL, 0x3E20B00000000000uL)),
                ((+1, 9, 0xDA1674E373DAF8FDuL, 0x9C20000000000000uL),(+1, 12, 0xBD8BC92CCA7BB091uL, 0x4BEC000000000000uL)),
                ((+1, 0, 0xEA9C000000000000uL, 0x0000000000000000uL),(+1, 5, 0x930AB1C8C2507026uL, 0x3E00000000000000uL)),
            }));

            public static ddouble Value(ddouble u) {
                Debug.Assert(u >= 2d);

                ddouble v = 1d / (u * u);

                ddouble y = ApproxUtil.Pade(v, pade_plus_limit) / (u * PI);

                return y;
            }
        }
    }
}