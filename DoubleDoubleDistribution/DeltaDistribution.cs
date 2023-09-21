using DoubleDouble;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class DeltaDistribution : Distribution {

        public ddouble Mu { get; }

        public DeltaDistribution() : this(mu: 0) { }

        public DeltaDistribution(ddouble mu) {
            ValidateLocation(mu);

            this.Mu = mu;
        }

        public override ddouble PDF(ddouble x) {
            ddouble pdf = x == Mu ? PositiveInfinity : Zero;

            return pdf;
        }

        public override ddouble CDF(ddouble x) {
            if (IsZero(Mu)) {
                ddouble cdf = IsNegative(x) ? 0 : 1;

                return cdf;
            }
            else {
                ddouble cdf = x < Mu ? 0 : 1;

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p) => Mu;

        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;

        public override ddouble Mean => Mu;
        public override ddouble Median => Mu;
        public override ddouble Mode => Mu;
        public override ddouble Variance => 0;

        public override ddouble Entropy => NegativeInfinity;

        public override string ToString() {
            return $"{typeof(DeltaDistribution).Name}[mu={Mu}]";
        }
    }
}
