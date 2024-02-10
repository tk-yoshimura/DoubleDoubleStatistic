﻿using DoubleDouble;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class LevyDistribution : StableDistribution,
        IAdditionOperators<LevyDistribution, LevyDistribution, LevyDistribution>,
        ISubtractionOperators<LevyDistribution, LevyDistribution, LevyDistribution>,
        IMultiplyOperators<LevyDistribution, ddouble, LevyDistribution> {

        public override ddouble Mu { get; }
        public override ddouble C { get; }

        private readonly ddouble pdf_norm;

        public LevyDistribution() : this(mu: 0, c: 1) { }

        public LevyDistribution(ddouble mu, ddouble c) {
            ValidateLocation(mu);
            ValidateScale(c);

            Mu = mu;
            C = c;

            pdf_norm = Sqrt(C / (2 * PI));
        }

        public override ddouble PDF(ddouble x) {
            ddouble v = x - Mu;
            if (v <= 0d) {
                return 0d;
            }

            ddouble pdf = pdf_norm * Exp(-C / (2 * v)) / Cube(Sqrt(v));

            if (IsNaN(pdf)) {
                return 0d;
            }

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble v = x - Mu;

            if (interval == Interval.Lower) {
                if (v <= 0d) {
                    return 0d;
                }

                ddouble cdf = Erfc(Sqrt(C / (2 * v)));

                return cdf;
            }
            else {
                if (v <= 0d) {
                    return 1d;
                }

                ddouble cdf = Erf(Sqrt(C / (2 * v)));

                return cdf;
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble quantile = Mu + C / (2 * Square(InverseErfc(p)));

                return quantile;
            }
            else {
                ddouble quantile = Mu + C / (2 * Square(InverseErf(p)));

                return quantile;
            }
        }

        public override bool AdditiveClosed => true;
        public override bool SubtractiveClosed => true;
        public override bool Scalable => true;

        public override ddouble Mean => PositiveInfinity;
        public override ddouble Median => Mu + C / (2 * Square(Erfc(0.5d)));
        public override ddouble Mode => Mu + C / 3d;
        public override ddouble Variance => PositiveInfinity;

        public override ddouble Entropy => (1d + 3d * EulerGamma + Log(16 * PI * C * C)) / 2;

        public override ddouble Alpha => 0.5d;
        public override ddouble Beta => 1d;

        public static LevyDistribution operator +(LevyDistribution dist1, LevyDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, Hypot(dist1.C, dist2.C));
        }

        public static LevyDistribution operator -(LevyDistribution dist1, LevyDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, Hypot(dist1.C, dist2.C));
        }

        public static LevyDistribution operator *(LevyDistribution dist, ddouble k) {
            return new(k * dist.Mu, k * dist.C);
        }

        public override string ToString() {
            return $"{typeof(LevyDistribution).Name}[mu={Mu},c={C}]";
        }
    }
}
