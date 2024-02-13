using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class DeltaDistribution : StableDistribution,
        IAdditionOperators<DeltaDistribution, DeltaDistribution, DeltaDistribution>,
        ISubtractionOperators<DeltaDistribution, DeltaDistribution, DeltaDistribution>,
        IAdditionOperators<DeltaDistribution, ddouble, DeltaDistribution>,
        IMultiplyOperators<DeltaDistribution, ddouble, DeltaDistribution> {

        public override ddouble Mu { get; }

        public DeltaDistribution() : this(mu: 0d) { }

        public DeltaDistribution(ddouble mu) {
            ValidateLocation(mu);

            Mu = mu;
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = x == Mu ? PositiveInfinity : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (interval == Interval.Upper) {
                return CDF(-x, Interval.Lower);
            }

            if (IsZero(Mu)) {
                ddouble cdf = IsNegative(x) ? 0 : 1;

                return cdf;
            }
            else {
                ddouble cdf = x < Mu ? 0 : 1;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) => InRangeUnit(p) ? Mu : NaN;

        public override bool Symmetric => true;

        public override ddouble Mean => Mu;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Variance => 0d;

        public override ddouble Entropy => NegativeInfinity;

        public override ddouble Alpha => 0d;

        public override ddouble Beta => 0d;

        public override ddouble C => 0d;

        public static DeltaDistribution operator +(DeltaDistribution dist1, DeltaDistribution dist2) {
            return new(dist1.Mu + dist2.Mu);
        }

        public static DeltaDistribution operator -(DeltaDistribution dist1, DeltaDistribution dist2) {
            return new(dist1.Mu - dist2.Mu);
        }
        public static DeltaDistribution operator +(DeltaDistribution dist, ddouble s) {
            return new(s + dist.Mu);
        }

        public static DeltaDistribution operator *(DeltaDistribution dist, ddouble k) {
            return new(k * dist.Mu);
        }

        public override string ToString() {
            return $"{typeof(DeltaDistribution).Name}[mu={Mu}]";
        }
    }
}
