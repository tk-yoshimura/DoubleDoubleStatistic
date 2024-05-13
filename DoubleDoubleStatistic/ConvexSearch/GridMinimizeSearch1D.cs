using DoubleDouble;

namespace DoubleDoubleStatistic.Optimizer {
    internal class GridMinimizeSearch1D : IConvexSearch<ddouble> {

        private readonly Func<ddouble, ddouble> cost_func;

        private int iter;

        private (Int128 position, ddouble param) range_min, range_max;
        private (Int128 position, ddouble param) next_min, next_max;
        private readonly Dictionary<Int128, (ddouble param, ddouble cost)> searched = [];

        public GridMinimizeSearch1D(Func<ddouble, ddouble> cost_func, (ddouble min, ddouble max) range, int depth, int max_iter) {
            ArgumentOutOfRangeException.ThrowIfLessThan(max_iter, 1, nameof(max_iter));
            ArgumentOutOfRangeException.ThrowIfLessThan(depth, 1, nameof(depth));
            ArgumentOutOfRangeException.ThrowIfGreaterThan(depth, 120, nameof(depth));

            if (range.min > range.max) {
                (range.min, range.max) = (range.max, range.min);
            }

            this.cost_func = cost_func;
            this.iter = max_iter;

            next_min = range_min = (Int128.NegativeOne << depth, range.min);
            next_max = range_max = (Int128.One << depth, range.max);
        }

        public ddouble BestParam => searched.Count >= 1
            ? searched.MinBy(item => item.Value.cost).Value.param
            : ddouble.NaN;

        public bool Finished => iter <= 0 || Int128.Abs(next_max.position - next_min.position) <= 1;

        private ddouble CostCache(Int128 position, ddouble param) {
            ddouble cost;

            if (!searched.TryGetValue(position, out (ddouble param, ddouble cost) value)) {
                cost = cost_func(param);
                cost = ddouble.IsNaN(cost) ? ddouble.PositiveInfinity : cost;
                searched.Add(position, (param, cost));

                return cost;
            }

            return value.cost;
        }

        public void Iteration() {
            if (Finished) {
                return;
            }

            iter--;

            (Int128 position_min, ddouble param_min) = next_min;
            (Int128 position_max, ddouble param_max) = next_max;
            (Int128 position_center, ddouble param_center) = ((position_min + position_max) / 2, (param_min + param_max) * 0.5d);

            ddouble cost_min = CostCache(position_min, param_min);
            ddouble cost_max = CostCache(position_max, param_max);
            ddouble cost_center = CostCache(position_center, param_center);

            if ((ddouble.Abs(param_min - param_max) <= ddouble.Abs(param_center) * 1e-30) ||
                (cost_min == cost_center && cost_center == cost_max && param_center == BestParam)) {

                next_min = next_max = (position_center, param_center);

                searched[position_center] = (param_center, ddouble.BitDecrement(cost_center));

                return;
            }

            if (cost_min >= cost_center && cost_center <= cost_max) {
                next_min = ((position_min + position_center) / 2, (param_min + param_center) * 0.5d);
                next_max = ((position_max + position_center) / 2, (param_max + param_center) * 0.5d);
            }
            else if (cost_min < cost_max) {
                next_min = (
                    Int128.Max(range_min.position, 2 * position_min - position_center),
                    ddouble.Max(range_min.param, 2d * param_min - param_center)
                );
                next_max = (position_center, param_center);
            }
            else {
                next_min = (position_center, param_center);
                next_max = (
                    Int128.Min(range_max.position, 2 * position_max - position_center),
                    ddouble.Min(range_max.param, 2d * param_max - param_center)
                );
            }
        }

        public static ddouble Search(Func<ddouble, ddouble> cost_func, (ddouble min, ddouble max) range, int iter) {
            GridMinimizeSearch1D searcher = new(cost_func, range, depth: 100, max_iter: iter);

            while (!searcher.Finished) {
                searcher.Iteration();
            }

            return searcher.BestParam;
        }
    }
}
