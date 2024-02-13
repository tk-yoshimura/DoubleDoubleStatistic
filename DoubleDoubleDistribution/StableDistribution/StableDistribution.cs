using DoubleDouble;

namespace DoubleDoubleDistribution {
    public abstract class StableDistribution : ContinuousDistribution {
        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;
        public override bool Scalable => true;
        public override bool Shiftable => true;

        public virtual ddouble Alpha => throw new NotImplementedException();
        public virtual ddouble Beta => throw new NotImplementedException();
        public virtual ddouble C => throw new NotImplementedException();
        public virtual ddouble Mu => throw new NotImplementedException();
    }
}