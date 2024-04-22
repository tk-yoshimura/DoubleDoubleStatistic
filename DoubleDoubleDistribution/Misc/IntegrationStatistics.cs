using DoubleDouble;
using DoubleDoubleIntegrate;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution.Misc {
    public static class IntegrationStatistics {
        public static ddouble Mean(ContinuousDistribution dist, ddouble eps, long discontinue_eval_points) {
            ddouble f(ddouble x) {
                ddouble pdf = dist.PDF(x);

                ddouble y = x * pdf;

                return y;
            }

            (ddouble mean, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(
                f, dist.Support.min, dist.Support.max, eps, discontinue_eval_points: discontinue_eval_points
            );

            Debug.WriteLine($"Mean integrate err: {err}");
            Debug.WriteLine($"Mean integrate eval_points: {eval_points}");

            return mean;
        }

        public static ddouble Variance(ContinuousDistribution dist, ddouble eps, long discontinue_eval_points) {
            ddouble mu = dist.Mean;

            ddouble f(ddouble x) {
                ddouble pdf = dist.PDF(x);

                ddouble y = Square(x - mu) * pdf;

                return y;
            }

            (ddouble variance, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(
                f, dist.Support.min, dist.Support.max, eps, discontinue_eval_points: discontinue_eval_points
            );

            Debug.WriteLine($"Variance integrate err: {err}");
            Debug.WriteLine($"Variance integrate eval_points: {eval_points}");

            return variance;
        }

        public static ddouble Skewness(ContinuousDistribution dist, ddouble eps, long discontinue_eval_points) {
            ddouble mu = dist.Mean, var = dist.Variance;

            ddouble f(ddouble x) {
                ddouble pdf = dist.PDF(x);

                ddouble y = Cube(x - mu) * pdf;

                return y;
            }

            (ddouble value, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(
                f, dist.Support.min, dist.Support.max, eps, discontinue_eval_points: discontinue_eval_points
            );

            Debug.WriteLine($"Skewness integrate err: {err}");
            Debug.WriteLine($"Skewness integrate eval_points: {eval_points}");

            ddouble skewness = value / ExMath.Pow3d2(var);

            return skewness;
        }

        public static ddouble Kurtosis(ContinuousDistribution dist, ddouble eps, long discontinue_eval_points) {
            ddouble mu = dist.Mean, var = dist.Variance;

            ddouble f(ddouble x) {
                ddouble pdf = dist.PDF(x);

                ddouble y = Square(Square(x - mu)) * pdf;

                return y;
            }

            (ddouble value, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(
                f, dist.Support.min, dist.Support.max, eps, discontinue_eval_points: discontinue_eval_points
            );

            Debug.WriteLine($"Kurtosis integrate err: {err}");
            Debug.WriteLine($"Kurtosis integrate eval_points: {eval_points}");

            ddouble kurtosis = value / Square(var) - 3d;

            return kurtosis;
        }

        public static ddouble Entropy(ContinuousDistribution dist, ddouble eps, long discontinue_eval_points) {
            ddouble f(ddouble x) {
                ddouble pdf = dist.PDF(x);

                if (pdf == 0d) {
                    return 0d;
                }

                ddouble y = -pdf * Log(pdf);

                return y;
            }

            (ddouble entropy, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(
                f, 0, PositiveInfinity, eps, discontinue_eval_points: discontinue_eval_points
            );

            Debug.WriteLine($"Entropy integrate err: {err}");
            Debug.WriteLine($"Entropy integrate eval_points: {eval_points}");

            return entropy;
        }
    }
}
