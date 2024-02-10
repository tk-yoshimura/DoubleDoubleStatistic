using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class DeltaDistribution : ContinuousDistribution,
        IAdditionOperators<DeltaDistribution, DeltaDistribution, DeltaDistribution>,
        ISubtractionOperators<DeltaDistribution, DeltaDistribution, DeltaDistribution> {

        public ddouble Mu { get; }

        public DeltaDistribution() : this(mu: 0) { }

        public DeltaDistribution(ddouble mu) {
            ValidateLocation(mu);

            Mu = mu;
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = x == Mu ? PositiveInfinity : Zero;

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

        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;
        public override bool Scalable => true;
        public override bool Symmetric => true;

        public override ddouble Mean => Mu;
        public override ddouble Median => Mu;
        public override ddouble Mode => Mu;
        public override ddouble Variance => 0;

        public override ddouble Entropy => NegativeInfinity;

        public static DeltaDistribution operator +(DeltaDistribution dist1, DeltaDistribution dist2) {
            return new(dist1.Mu + dist2.Mu);
        }

        public static DeltaDistribution operator -(DeltaDistribution dist1, DeltaDistribution dist2) {
            return new(dist1.Mu - dist2.Mu);
        }

        public static DeltaDistribution operator *(ddouble k, DeltaDistribution dist) {
            return new(k * dist.Mu);
        }

        public override string ToString() {
            return $"{typeof(DeltaDistribution).Name}[mu={Mu}]";
        }
    }
}
