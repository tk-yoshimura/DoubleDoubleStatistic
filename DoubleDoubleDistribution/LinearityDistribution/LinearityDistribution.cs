using DoubleDouble;
using System.Numerics;

namespace DoubleDoubleDistribution {
    public abstract class LinearityDistribution<TSelf> : ScalableDistribution<TSelf>
        where TSelf :
            IMultiplyOperators<TSelf, ddouble, TSelf>,
            IAdditionOperators<TSelf, ddouble, TSelf>,
            ISubtractionOperators<TSelf, ddouble, TSelf> {

        public sealed override bool Shiftable => true;
    }
}