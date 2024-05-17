using DoubleDouble;

namespace DoubleDoubleStatistic.Optimizer {
    internal class BisectionMinimizeSearch2D : IConvexSearch<(ddouble x, ddouble y)> {

        private readonly Func<(ddouble x, ddouble y), ddouble> cost_func;

        private int iter;

        private (Int128 position, ddouble param) range_xmin, range_xmax, range_ymin, range_ymax;
        private (Int128 position, ddouble param) next_xmin, next_xmax, next_ymin, next_ymax;
        private readonly Dictionary<(Int128 x, Int128 y), ((ddouble x, ddouble y) param, ddouble cost)> searched = [];

        public BisectionMinimizeSearch2D(Func<(ddouble x, ddouble y), ddouble> cost_func, (ddouble min, ddouble max) range_x, (ddouble min, ddouble max) range_y, int depth, int max_iter) {
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
            (Int128 position_xcen, ddouble param_xcen) = ((position_xmin + position_xmax) / 2, (param_xmin + param_xmax) * 0.5d);

            (Int128 position_ymin, ddouble param_ymin) = next_ymin;
            (Int128 position_ymax, ddouble param_ymax) = next_ymax;
            (Int128 position_ycen, ddouble param_ycen) = ((position_ymin + position_ymax) / 2, (param_ymin + param_ymax) * 0.5d);

            ddouble cost_xmin_ymin = CostCache((position_xmin, position_ymin), (param_xmin, param_ymin));
            ddouble cost_xcen_ymin = CostCache((position_xcen, position_ymin), (param_xcen, param_ymin));
            ddouble cost_xmax_ymin = CostCache((position_xmax, position_ymin), (param_xmax, param_ymin));
            ddouble cost_xmin_ycen = CostCache((position_xmin, position_ycen), (param_xmin, param_ycen));
            ddouble cost_xcen_ycen = CostCache((position_xcen, position_ycen), (param_xcen, param_ycen));
            ddouble cost_xmax_ycen = CostCache((position_xmax, position_ycen), (param_xmax, param_ycen));
            ddouble cost_xmin_ymax = CostCache((position_xmin, position_ymax), (param_xmin, param_ymax));
            ddouble cost_xcen_ymax = CostCache((position_xcen, position_ymax), (param_xcen, param_ymax));
            ddouble cost_xmax_ymax = CostCache((position_xmax, position_ymax), (param_xmax, param_ymax));

            if ((ddouble.Abs(param_xmin - param_xmax) <= ddouble.Abs(param_xcen) * 1e-30 &&
                 ddouble.Abs(param_ymin - param_ymax) <= ddouble.Abs(param_ycen) * 1e-30) ||
                (cost_xmin_ycen == cost_xcen_ycen && cost_xcen_ycen == cost_xmax_ycen &&
                 cost_xcen_ymin == cost_xcen_ycen && cost_xcen_ycen == cost_xcen_ymax && (param_xcen, param_ycen) == BestParam)) {

                next_xmin = next_xmax = (position_xcen, param_xcen);
                next_ymin = next_ymax = (position_ycen, param_ycen);

                searched[(position_xcen, position_ycen)] = ((param_xcen, param_ycen), ddouble.BitDecrement(cost_xcen_ycen));

                return;
            }

            (ddouble cost, (int x, int y) index, (Int128 x, Int128 y) position)[] pts = [
                (cost_xmin_ymin, (-1, -1), (position_xmin, position_ymin)),
                (cost_xcen_ymin, ( 0, -1), (position_xcen, position_ymin)),
                (cost_xmax_ymin, (+1, -1), (position_xmax, position_ymin)),
                (cost_xmin_ycen, (-1,  0), (position_xmin, position_ycen)),
                (cost_xcen_ycen, ( 0,  0), (position_xcen, position_ycen)),
                (cost_xmax_ycen, (+1,  0), (position_xmax, position_ycen)),
                (cost_xmin_ymax, (-1, +1), (position_xmin, position_ycen)),
                (cost_xcen_ymax, ( 0, +1), (position_xcen, position_ycen)),
                (cost_xmax_ymax, (+1, +1), (position_xmax, position_ycen)),
            ];

            (ddouble cost, (int x, int y) index, (Int128 x, Int128 y) position) pts_mincost = pts.MinBy(item => item.cost);

            if (pts_mincost.index == (0, 0)) {
                next_xmin = ((position_xmin + position_xcen) / 2, (param_xmin + param_xcen) * 0.5d);
                next_xmax = ((position_xmax + position_xcen) / 2, (param_xmax + param_xcen) * 0.5d);
                next_ymin = ((position_ymin + position_ycen) / 2, (param_ymin + param_ycen) * 0.5d);
                next_ymax = ((position_ymax + position_ycen) / 2, (param_ymax + param_ycen) * 0.5d);
            }
            else {
                Int128 position_xsft = pts_mincost.index.x * (position_xmax - position_xmin) / 2;
                ddouble param_xsft = pts_mincost.index.x * (param_xmax - param_xmin) * 0.5d;

                Int128 position_ysft = pts_mincost.index.y * (position_ymax - position_ymin) / 2;
                ddouble param_ysft = pts_mincost.index.y * (param_ymax - param_ymin) * 0.5d;

                next_xmin = (
                    Int128.Clamp(next_xmin.position + position_xsft, range_xmin.position, range_xmax.position),
                    ddouble.Clamp(param_xmin + param_xsft, range_xmin.param, range_xmax.param)
                );
                next_xmax = (
                    Int128.Clamp(next_xmax.position + position_xsft, range_xmin.position, range_xmax.position),
                    ddouble.Clamp(param_xmax + param_xsft, range_xmin.param, range_xmax.param)
                );

                next_ymin = (
                    Int128.Clamp(next_ymin.position + position_ysft, range_ymin.position, range_ymax.position),
                    ddouble.Clamp(param_ymin + param_ysft, range_ymin.param, range_ymax.param)
                );
                next_ymax = (
                    Int128.Clamp(next_ymax.position + position_ysft, range_ymin.position, range_ymax.position),
                    ddouble.Clamp(param_ymax + param_ysft, range_ymin.param, range_ymax.param)
                );
            }
        }

        public static (ddouble x, ddouble y) Search(Func<(ddouble x, ddouble y), ddouble> cost_func, ((ddouble x, ddouble y) min, (ddouble x, ddouble y) max) range, int iter) {
            BisectionMinimizeSearch2D searcher = new(cost_func, (range.min.x, range.max.x), (range.min.y, range.max.y), depth: 100, max_iter: iter);

            while (!searcher.Finished) {
                searcher.Iteration();
            }

            return searcher.BestParam;
        }
    }
}
