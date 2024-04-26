using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class IrwinHallDistribution : ContinuousDistribution,
        IAdditionOperators<IrwinHallDistribution, IrwinHallDistribution, IrwinHallDistribution> {

        public int N { get; }

        private readonly ReadOnlyCollection<ddouble> cdf_intway;
        private readonly List<ReadOnlyCollection<ddouble>> pdf_table = [], cdf_table = [];
        private readonly ddouble pdf_norm, cdf_norm, quantile_norm;

        public IrwinHallDistribution(int n) {
            ValidateShape(n, n => n > 0 && n <= 64);

            N = n;

            pdf_norm = TaylorSequence[n - 1];
            cdf_norm = TaylorSequence[n];
            quantile_norm = Factorial[n];

            for (int k = 0; k < (n + 1) / 2; k++) {
                Polynomial pdf = PolynomialCoef.PDF(n, k);
                Polynomial cdf = PolynomialCoef.CDF(n, k);

                pdf_table.Add(new ReadOnlyCollection<ddouble>(pdf.Coef.Select(c => (ddouble)c).Reverse().ToArray()));
                cdf_table.Add(new ReadOnlyCollection<ddouble>(cdf.Coef.Select(c => (ddouble)c).Reverse().ToArray()));
            }

            List<ddouble> intway = [0];

            for (int k = 1; k < cdf_table.Count; k++) {
                intway.Add(cdf_table[k][^1]);
            }

            if (n % 2 != 0) {
                intway.Add(quantile_norm - intway[^1]);
            }
            else {
                intway.Add(quantile_norm / 2);
            }

            cdf_intway = Array.AsReadOnly(intway.ToArray());
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x) || x > N) {
                return 0d;
            }

            if (x * 2 > N) {
                x = N - x;
            }

            if (x <= 1d) {
                ddouble pdf = Pow(x, N - 1) * pdf_norm;
                return pdf;
            }

            int index = int.Clamp((int)Floor(x), 0, pdf_table.Count - 1);
            ddouble c = x - index;

            ReadOnlyCollection<ddouble> coef = pdf_table[index];

            ddouble y = coef[0];

            for (int i = 1; i < coef.Count; i++) {
                y = coef[i] + y * c;
            }

            y *= pdf_norm;

            return y;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                x = N - x;
            }

            if (x * 2 > N) {
                return 1d - CDF(x: N - x, interval: Interval.Lower);
            }

            if (IsNegative(x)) {
                return 0d;
            }

            if (x <= 1d) {
                ddouble cdf = Pow(x, N) * cdf_norm;
                return cdf;
            }

            int index = int.Clamp((int)Floor(x), 0, cdf_table.Count - 1);
            ddouble c = x - index;

            ReadOnlyCollection<ddouble> coef = cdf_table[index];

            ddouble y = coef[0];

            for (int i = 1; i < coef.Count; i++) {
                y = coef[i] + y * c;
            }

            y *= cdf_norm;

            return y;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Upper) {
                p = 1d - p;
            }

            if (p > 0.5d) {
                return N - Quantile(1d - p, interval: Interval.Lower);
            }

            ddouble v = p * quantile_norm;

            if (v <= 1d) {
                ddouble x = RootN(v, N);
                return x;
            }
            else {
                int index;
                for (index = cdf_intway.Count - 2; index > 0; index--) {
                    if (v >= cdf_intway[index]) {
                        break;
                    }
                }

                Debug.Assert(index >= 0 && index < cdf_intway.Count, nameof(index));
                Debug.Assert(v >= cdf_intway[index], nameof(index));

                ddouble c = (v - cdf_intway[index]) / (cdf_intway[index + 1] - cdf_intway[index]);

                ReadOnlyCollection<ddouble> pdf_coef = pdf_table[index];
                ReadOnlyCollection<ddouble> cdf_coef = cdf_table[index];

                for (int iter = 0; iter < 16; iter++) {
                    ddouble pdf = pdf_coef[0], cdf = cdf_coef[0];

                    for (int i = 1; i < pdf_coef.Count; i++) {
                        pdf = pdf_coef[i] + pdf * c;
                    }
                    for (int i = 1; i < cdf_coef.Count; i++) {
                        cdf = cdf_coef[i] + cdf * c;
                    }

                    pdf *= N;

                    ddouble dc = (cdf - v) / pdf;

                    c -= dc;
                    c = Clamp(c, 0d, 1d);

                    if (Abs(dc) < 1e-30) {
                        break;
                    }
                }

                ddouble x = index + c;

                return x;
            }
        }

        public override bool AdditiveClosed => true;

        public override (ddouble min, ddouble max) Support => (0d, N);

        public override ddouble Mean => N * 0.5d;

        public override ddouble Median => N * 0.5d;

        public override ddouble Mode => (N > 1) ? (N * 0.5d) : NaN;

        public override ddouble Variance => (ddouble)N / 12d;

        public override ddouble Skewness => 0;

        public override ddouble Kurtosis => -(ddouble)6 / (5 * N);

        public override ddouble Entropy => entropy_table[N - 1];

        public static IrwinHallDistribution operator +(IrwinHallDistribution dist1, IrwinHallDistribution dist2) {
            return new(dist1.N + dist2.N);
        }

        public override string ToString() {
            return $"{typeof(IrwinHallDistribution).Name}[n={N}]";
        }

        public override string Formula => "p(x; n) := sum((-1)^k * binom(n, k) * (x - k)^(n - 1), k, 0, floor(x)) / (n - 1)!";

        private static class PolynomialCoef {
            static readonly Dictionary<(int n, int k, int j), BigInteger> table = [];

            private static BigInteger ACoef(int n, int k, int j) {
                if (n <= 0 || k < 0 || j < 0 || j >= n || k >= n) {
                    return BigInteger.Zero;
                }

                if (k == 0) {
                    return (j < n - 1) ? BigInteger.Zero : BigInteger.One;
                }

                if (table.TryGetValue((n, k, j), out BigInteger value)) {
                    return value;
                }

                BigInteger c = Binom.Value(n, k) * Binom.Value(n - 1, j) * BigInteger.Pow(k, n - j - 1);

                int s = n + k - j - 1;

                BigInteger new_value = ACoef(n, k - 1, j) + (((s & 1) == 0) ? c : -c);

                table.Add((n, k, j), new_value);

                return new_value;
            }

            public static Polynomial PDF(int n, int k) {
                if (n <= 0 || k < 0 || k >= n) {
                    return Polynomial.Zero;
                }

                List<BigInteger> coef = [];

                for (int j = 0; j < n; j++) {
                    coef.Add(ACoef(n, k, j));
                }

                Polynomial p = new([.. coef]);
                Polynomial p_sft = Polynomial.XShift(p, k);

                return p_sft;
            }

            public static Polynomial CDF(int n, int k) {
                if (n <= 0 || k < 0 || k >= n) {
                    return Polynomial.Zero;
                }

                List<BigInteger> coef = [.. Polynomial.XMul(n * PDF(n, k), 1).Coef];

                for (int i = 1; i < coef.Count; i++) {
                    if ((coef[i] % i) != 0) {
                        throw new ArithmeticException();
                    }

                    coef[i] /= i;
                }

                BigInteger c0 = BigInteger.Zero;

                for (int i = 0; i <= k; i++) {
                    c0 += (((i & 1) == 0) ? +1 : -1) * Binom.Value(n, i) * BigInteger.Pow(k - i, n);
                }

                coef[0] = c0;

                Polynomial p = new([.. coef]);

                return p;
            }
        }

        private static readonly ReadOnlyCollection<ddouble> entropy_table = new([
            0,
            0.5,
            "0.71929485277519245478725740580423449843973",
            "0.86673511712544196694891114667595801960822",
            "0.97955793901712421984811590313404104529006",
            "1.0712972171474640187858180659188173938623",
            "1.1486893990383667482223436097512365523907",
            "1.2156482856684973838857166733356307947329",
            "1.2746665990034460625452336438279804400602",
            "1.3274346769608224188284261532418828582042",
            "1.3751531514119085978864996710110577405451",
            "1.4187061074463438524411014027392511242367",
            "1.4587636599363741437056996762596012084364",
            "1.4958459872410572415849866428619781141359",
            "1.5303650305066239407542294249196999393678",
            "1.5626526155878752900412291932424853707171",
            "1.5929799866799252076196889672449469197520",
            "1.6215717219910005063933763887311778846628",
            "1.6486158670259684368739325751786184985394",
            "1.6742714567294664369875389837176607895850",
            "1.6986741949778706333898739436942732947029",
            "1.7219408081656692249866983412248674350146",
            "1.7441724279967231014524362653719308262294",
            "1.7654572523051930575537036100069649844674",
            "1.7858726613326554362939282002129069575579",
            "1.8054869179942024500748791212944163340533",
            "1.8243605465953447367672900069169316641070",
            "1.8425474603400621889258729539538594248151",
            "1.8600958906434253812340247107989542217212",
            "1.8770491586491007512168591592802234512813",
            "1.8934463200569933111807761531835803442125",
            "1.9093227074382482886594860625103280648043",
            "1.9247103889964243719030693940922058215781",
            "1.9396385587643197328584495124719567156560",
            "1.9541338701790117074652393864121228594609",
            "1.9682207226187634102116031627481747205866",
            "1.9819215086445376985448920643152518940796",
            "1.9952568282413676382008091917756377148833",
            "2.0082456752086177486558388486755138454090",
            "2.0209055999344627430330765889818394846396",
            "2.0332528520569541548950610247308517866159",
            "2.0453025059225782762261416279560266051290",
            "2.0570685712732233824905472743069253615108",
            "2.0685640912008720802666064477519034562271",
            "2.0798012290882162080359766649803452951801",
            "2.0907913459887972221641325305220315058806",
            "2.1015450696812508755307592031041209821646",
            "2.1120723564501367756096760791006313253860",
            "2.1223825464938012358073820016577427976644",
            "2.1324844137322833488275341131780524057760",
            "2.1423862106810419337335319834177625172351",
            "2.1520957089657203316972384649094654148840",
            "2.1616202359764163066560651975999954492065",
            "2.1709667080946601404976639529837606033675",
            "2.1801416608706270829868321818753285253451",
            "2.1891512764804646641267766497717271596182",
            "2.1980014087527206761016471050280668345326",
            "2.2066976060176564134010262092011645059322",
            "2.2152451320028454654169541501409935365572",
            "2.2236489849721618647725779482963504365809",
            "2.2319139152824441933203372248569536022086",
            "2.2400444415122748402885195155964058011487",
            "2.2480448653000080744928293155776943266176",
            "2.2559192850130564506229632421139753581484",
        ]);
    }
}
