using DoubleDouble;

namespace DoubleDoubleStatistic.Optimizer {
    internal interface IConvexSearch<T> {
        T BestParam { get; }

        bool Finished { get; }

        void Iteration();

        static abstract T Search(Func<T, ddouble> cost_func, (T min, T max) range, int iter);
    }
}
