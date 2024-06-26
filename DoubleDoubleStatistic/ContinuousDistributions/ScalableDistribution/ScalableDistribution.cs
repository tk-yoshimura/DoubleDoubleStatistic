﻿using DoubleDouble;
using System.Numerics;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    public abstract class ScalableDistribution<TSelf> : ContinuousDistribution
        where TSelf :
            IMultiplyOperators<TSelf, ddouble, TSelf>,
            IDivisionOperators<TSelf, ddouble, TSelf> {

        public sealed override bool Scalable => true;
    }
}