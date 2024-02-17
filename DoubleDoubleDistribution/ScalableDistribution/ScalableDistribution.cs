using DoubleDouble;
using System.Numerics;

namespace DoubleDoubleDistribution {
    public abstract class ScalableDistribution<TSelf> : ContinuousDistribution
        where TSelf :
            IMultiplyOperators<TSelf, ddouble, TSelf> {

        public sealed override bool Scalable => true;
    }
}