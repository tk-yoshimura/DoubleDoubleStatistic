using DoubleDouble;

namespace DoubleDoubleStatistic.Optimizer {
    internal class GridMinimizeSearch2D : IConvexSearch<(ddouble x, ddouble y)> {

        private readonly Func<(ddouble x, ddouble y), ddouble> cost_func;

        private int iter;

        private (Int128 position, ddouble param) range_xmin, range_xmax, range_ymin, range_ymax;
        private (Int128 position, ddouble param) next_xmin, next_xmax, next_ymin, next_ymax;
        private readonly Dictionary<(Int128 x, Int128 y), ((ddouble x, ddouble y) param, ddouble cost)> searched = [];

        public GridMinimizeSearch2D(Func<(ddouble x, ddouble y), ddouble> cost_func, (ddouble min, ddouble max) range_x, (ddouble min, ddouble max) range_y, int depth, int max_iter) {
            ArgumentOutOfRangeException.ThrowIfLessThan(max_iter, 1, nameof(max_iter));
            ArgumentOutOfRangeException.ThrowIfLessThan(depth, 1, nameof(depth));
            ArgumentOutOfRangeException.ThrowIfGreaterThan(depth, 120, nameof(depth));

            if (range_x.min > range_x.max) {
                (range_x.min, range_x.max) = (range_x.max, range_x.min);
            }
            if (range_y.min > range_y.max) {
                (range_y.min, range_y.max) = (range_y.max, range_y.min);
            }

            this.cost_func = cost_func;
            this.iter = max_iter;

            next_xmin = range_xmin = (Int128.NegativeOne << depth, range_x.min);
            next_xmax = range_xmax = (Int128.One << depth, range_x.max);

            next_ymin = range_ymin = (Int128.NegativeOne << depth, range_y.min);
            next_ymax = range_ymax = (Int128.One << depth, range_y.max);
        }

        public (ddouble x, ddouble y) BestParam => searched.Count >= 1
            ? searched.MinBy(item => item.Value.cost).Value.param
            : (ddouble.NaN, ddouble.NaN);

        public bool Finished => iter <= 0 ||
            (Int128.Abs(next_xmax.position - next_xmin.position) <= 1 && Int128.Abs(next_ymax.position - next_ymin.position) <= 1);

        private ddouble CostCache((Int128 x, Int128 y) position, (ddouble x, ddouble y) param) {
            ddouble cost;

            if (!searched.TryGetValue(position, out ((ddouble x, ddouble y) param, ddouble cost) value)) {
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

            (Int128 position_xmin, ddouble param_xmin) = next_xmin;
            (Int128 position_xmax, ddouble param_xmax) = next_xmax;
            (Int128 position_xcenter, ddouble param_xcenter) = ((position_xmin + position_xmax) / 2, (param_xmin + param_xmax) * 0.5d);

            (Int128 position_ymin, ddouble param_ymin) = next_ymin;
            (Int128 position_ymax, ddouble param_ymax) = next_ymax;
            (Int128 position_ycenter, ddouble param_ycenter) = ((position_ymin + position_ymax) / 2, (param_ymin + param_ymax) * 0.5d);

            ddouble cost_xmin = CostCache((position_xmin, position_ycenter), (param_xmin, param_ycenter));
            ddouble cost_xmax = CostCache((position_xmax, position_ycenter), (param_xmax, param_ycenter));
            ddouble cost_ymin = CostCache((position_xcenter, position_ymin), (param_xcenter, param_ymin));
            ddouble cost_ymax = CostCache((position_xcenter, position_ymax), (param_xcenter, param_ymax));
            ddouble cost_center = CostCache((position_xcenter, position_ycenter), (param_xcenter, param_ycenter));

            if ((ddouble.Abs(param_xmin - param_xmax) <= ddouble.Abs(param_xcenter) * 1e-30 &&
                 ddouble.Abs(param_ymin - param_ymax) <= ddouble.Abs(param_ycenter) * 1e-30) ||
                (cost_xmin == cost_center && cost_center == cost_xmax &&
                 cost_ymin == cost_center && cost_center == cost_ymax && (param_xcenter, param_ycenter) == BestParam)) {

                next_xmin = next_xmax = (position_xcenter, param_xcenter);
                next_ymin = next_ymax = (position_ycenter, param_ycenter);

                searched[(position_xcenter, position_ycenter)] = ((param_xcenter, param_ycenter), ddouble.BitDecrement(cost_center));

                return;
            }

            if (cost_xmin >= cost_center && cost_center <= cost_xmax) {
                next_xmin = ((position_xmin + position_xcenter) / 2, (param_xmin + param_xcenter) * 0.5d);
                next_xmax = ((position_xmax + position_xcenter) / 2, (param_xmax + param_xcenter) * 0.5d);
            }
            else if (cost_xmin < cost_xmax) {
                next_xmin = (
                    Int128.Max(range_xmin.position, 2 * position_xmin - position_xcenter),
                    ddouble.Max(range_xmin.param, 2d * param_xmin - param_xcenter)
                );
                next_xmax = (position_xcenter, param_xcenter);
            }
            else {
                next_xmin = (position_xcenter, param_xcenter);
                next_xmax = (
                    Int128.Min(range_xmax.position, 2 * position_xmax - position_xcenter),
                    ddouble.Min(range_xmax.param, 2d * param_xmax - param_xcenter)
                );
            }

            if (cost_ymin >= cost_center && cost_center <= cost_ymax) {
                next_ymin = ((position_ymin + position_ycenter) / 2, (param_ymin + param_ycenter) * 0.5d);
                next_ymax = ((position_ymax + position_ycenter) / 2, (param_ymax + param_ycenter) * 0.5d);
            }
            else if (cost_ymin < cost_ymax) {
                next_ymin = (
                    Int128.Max(range_ymin.position, 2 * position_ymin - position_ycenter),
                    ddouble.Max(range_ymin.param, 2d * param_ymin - param_ycenter)
                );
                next_ymax = (position_ycenter, param_ycenter);
            }
            else {
                next_ymin = (position_ycenter, param_ycenter);
                next_ymax = (
                    Int128.Min(range_ymax.position, 2 * position_ymax - position_ycenter),
                    ddouble.Min(range_ymax.param, 2d * param_ymax - param_ycenter)
                );
            }
        }

        public static (ddouble x, ddouble y) Search(Func<(ddouble x, ddouble y), ddouble> cost_func, ((ddouble x, ddouble y) min, (ddouble x, ddouble y) max) range, int iter) {
            GridMinimizeSearch2D searcher = new(cost_func, (range.min.x, range.max.x), (range.min.y, range.max.y), depth: 100, max_iter: iter);

            while (!searcher.Finished) {
                searcher.Iteration();
            }

            return searcher.BestParam;
        }
    }
}
