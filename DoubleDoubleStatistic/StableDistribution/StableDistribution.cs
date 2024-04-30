using DoubleDouble;
using System.Numerics;

namespace DoubleDoubleStatistic {
    public abstract class StableDistribution<TSelf> : LinearityDistribution<TSelf>
        where TSelf :
            IAdditionOperators<TSelf, TSelf, TSelf>,
            ISubtractionOperators<TSelf, TSelf, TSelf>,
            IMultiplyOperators<TSelf, ddouble, TSelf>,
            IDivisionOperators<TSelf, ddouble, TSelf>,
            IAdditionOperators<TSelf, ddouble, TSelf>,
            ISubtractionOperators<TSelf, ddouble, TSelf> {

        public sealed override bool AdditiveClosed => true;
        public sealed override bool SubtractiveClosed => true;

        public virtual ddouble Alpha => throw new NotImplementedException();
        public virtual ddouble Beta => throw new NotImplementedException();
        public virtual ddouble C => throw new NotImplementedException();
        public virtual ddouble Mu => throw new NotImplementedException();
    }
}