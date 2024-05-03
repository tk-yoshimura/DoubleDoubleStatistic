using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.Utils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class FoledNormalDistribution : ScalableDistribution<FoledNormalDistribution>,
        IMultiplyOperators<FoledNormalDistribution, ddouble, FoledNormalDistribution>,
        IDivisionOperators<FoledNormalDistribution, ddouble, FoledNormalDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }

        private readonly ddouble pdf_norm, sigma_sq, exp_scale, erf_scale;

        private readonly bool needs_fold;

        private QuantileBuilder quantile_lower_builder = null, quantile_upper_builder = null;

        public FoledNormalDistribution() : this(mu: 0d, sigma: 1d) { }

        public FoledNormalDistribution(ddouble sigma) : this(mu: 0d, sigma: sigma) { }

        public FoledNormalDistribution(ddouble mu, ddouble sigma) {
            ValidateLocation(mu);
            ValidateScale(sigma);

            mu = Abs(mu);

            Mu = mu;
            Sigma = sigma;

            sigma_sq = sigma * sigma;
            pdf_norm = 1d / (sigma * Sqrt(2d * PI));
            exp_scale = -1d / (2 * sigma_sq);
            erf_scale = 1d / (Sqrt2 * sigma);

            needs_fold = Exp(Square(mu) * exp_scale) > 0d;

            Debug.WriteLine($"{nameof(needs_fold)} : {needs_fold}");
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (IsNegative(x)) {
                return 0d;
            }

            if (needs_fold) {
                ddouble pdf = pdf_norm * (Exp(Square(x - Mu) * exp_scale) + Exp(Square(x + Mu) * exp_scale));

                return pdf;
            }
            else {
                ddouble pdf = pdf_norm * Exp(Square(x - Mu) * exp_scale);

                return pdf;
            }
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                if (IsNegative(x)) {
                    return 0d;
                }

                if (needs_fold) {
                    ddouble cdf = (Erf((x - Mu) * erf_scale) + Erf((x + Mu) * erf_scale)) * 0.5d;

                    return cdf;
                }
                else {
                    ddouble cdf = Erfc((Mu - x) * erf_scale) * 0.5d;

                    return cdf;
                }
            }
            else {
                if (IsNegative(x)) {
                    return 1d;
                }

                if (needs_fold) {
                    ddouble cdf = (Erfc((x - Mu) * erf_scale) + Erfc((Mu + x) * erf_scale)) * 0.5d;

                    return cdf;
                }
                else {
                    ddouble cdf = Erfc((x - Mu) * erf_scale) * 0.5d;

                    return cdf;
                }
            }
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (!needs_fold) {
                if (interval == Interval.Lower) {
                    ddouble x = Mu - Sigma * Sqrt2 * InverseErfc(2d * p);
                    x = Max(0d, x);

                    return x;
                }
                else {
                    ddouble x = Mu + Sigma * Sqrt2 * InverseErfc(2d * p);
                    x = Max(0d, x);

                    return x;
                }
            }
            else {
                ddouble bias = Mu / Sigma, c = 1d / Sqrt(2d * PI);

                ddouble df(ddouble u) {
                    return c * (Exp(-Square(u - bias) * 0.5d) + Exp(-Square(u + bias) * 0.5d));
                }

                if (interval == Interval.Lower) {
                    if (p == 0d) {
                        return 0d;
                    }
                    if (p == 1d) {
                        return PositiveInfinity;
                    }

                    ddouble f(ddouble u) {
                        return (Erf((u - bias) / Sqrt2) + Erf((u + bias) / Sqrt2)) * 0.5d;
                    }

                    this.quantile_lower_builder ??= new QuantileBuilder(0d, 64d, f, samples: 256);

                    (ddouble x, ddouble x0, ddouble x1) = quantile_lower_builder.Estimate(p);

                    for (int i = 0; i < 8; i++) {
                        ddouble y = f(x), dx = (y - p) / df(x);

                        if (!IsFinite(dx)) {
                            break;
                        }

                        x = Clamp(x - dx, x0, x1);

                        if (Abs(dx / x) < 1e-29 || Abs(dx) < Epsilon) {
                            break;
                        }
                    }

                    x *= Sigma;

                    return x;
                }
                else {
                    if (p == 0d) {
                        return PositiveInfinity;
                    }
                    if (p == 1d) {
                        return 0d;
                    }

                    ddouble f(ddouble u) {
                        return (Erfc((u - bias) / Sqrt2) + Erfc((bias + u) / Sqrt2)) * 0.5d;
                    }

                    this.quantile_upper_builder ??= new QuantileBuilder(64d, 0d, f, samples: 256);

                    (ddouble x, ddouble x0, ddouble x1) = quantile_upper_builder.Estimate(p);

                    for (int i = 0; i < 8; i++) {
                        ddouble y = f(x), dx = (y - p) / df(x);

                        if (!IsFinite(dx)) {
                            break;
                        }

                        x = Clamp(x + dx, x1, x0);

                        if (Abs(dx / x) < 1e-29 || Abs(dx) < Epsilon) {
                            break;
                        }
                    }

                    x *= Sigma;

                    return x;
                }
            }
        }

        public override (ddouble min, ddouble max) Support => (0d, PositiveInfinity);

        public override ddouble Mean => needs_fold
            ? Sigma * Sqrt(2d * RcpPI) * Exp(-Square(Mu / Sigma) * 0.5d) + Mu * (1d - Erfc(Mu / (Sigma * Sqrt2)))
            : Mu;

        public override ddouble Median => needs_fold
            ? Quantile(0.5d)
            : Mu;

        public override ddouble Mode {
            get {
                if (!needs_fold) {
                    return Mu;
                }

                ddouble x = ModePade.Value(Mu / Sigma) * Sigma;

                return x;
            }
        }

        public override ddouble Variance => needs_fold
            ? Square(Mu) + Square(Sigma) - Square(Mean)
            : sigma_sq;

        private ddouble? skewness = null;
        public override ddouble Skewness => needs_fold
            ? skewness ??= IntegrationStatistics.Skewness(this, eps: 1e-28, discontinue_eval_points: 8192)
            : 0d;

        private ddouble? kurtosis = null;
        public override ddouble Kurtosis => needs_fold
            ? kurtosis ??= IntegrationStatistics.Kurtosis(this, eps: 1e-28, discontinue_eval_points: 8192)
            : 0d;

        private ddouble? entropy = null;
        public override ddouble Entropy => needs_fold
            ? entropy ??= IntegrationStatistics.Entropy(this, eps: 1e-28, discontinue_eval_points: 16384)
            : Log(Sigma * Sqrt(2d * PI * E));

        public static FoledNormalDistribution operator *(FoledNormalDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.Sigma * k);
        }

        public static FoledNormalDistribution operator /(FoledNormalDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.Sigma / k);
        }

        public override string ToString() {
            return $"{typeof(FoledNormalDistribution).Name}[mu={Mu},sigma={Sigma}]";
        }

        public override string Formula => "p(x; mu, sigma) := (exp(-(x - mu)^2) + exp(-(x + mu)^2)) / (sqrt(2 * pi) * sigma)";

        private static class ModePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                (Zero, (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0x9CC470A0490973E8uL, 0x190AEA051EC10A4EuL), (-1, -1, 0x9A3E86EBA335EDF5uL, 0xA57431F81BE77138uL)),
                ((-1, 0, 0xBCE8EBE5A3F586F6uL, 0xAF26082AA16A1E69uL), (+1, -2, 0xF0798539AFD14C3EuL, 0x2B2E882A5F6322EEuL)),
                ((-1, -3, 0xC978D88EA4BB4A43uL, 0xA7A6AB37C4BAFBD2uL), (-1, -2, 0x872D3C2452EB4AB0uL, 0x4CE202EB48D9CA89uL)),
                ((+1, -3, 0xA8F9B7ED8F1B3BE0uL, 0x81DE7E33F1441533uL), (-1, -2, 0xB0C077E7E44CAD9DuL, 0xDC59DB30F9E68368uL)),
                ((+1, -3, 0xE40FAD6B58343687uL, 0xCCB6C78E34DC4EF0uL), (+1, -2, 0x8B65057D7BEDCF4AuL, 0xFFB19B6769427C8BuL)),
                ((-1, -9, 0x9A65BA1CAE3402CDuL, 0x72213CEC2C3DEAC0uL), (+1, -2, 0x8F83687BA8D34287uL, 0x31B0B4FA1DCC55CCuL)),
                ((+1, -2, 0xC2F20FAE37F33B41uL, 0xE2577338C5475B30uL), (-1, -3, 0x9422DEBE28771CAAuL, 0x6E28885EDF62D155uL)),
                ((-1, -3, 0xE5C2BDAFAFDD7416uL, 0x67EC37E0B9276ECCuL), (+1, -4, 0x94341778237D7CDBuL, 0x70BEE0F128D20053uL)),
                ((+1, -5, 0x9332CA2DE47FAD50uL, 0xB50A5AB4D87589B8uL), (-1, -5, 0xE0CC2803A2EFCE01uL, 0xDE2D4039E0A733E8uL)),
                ((-1, -9, 0x9CFCB55A34A7E043uL, 0x3A5299B4A61F30F7uL), (-1, -6, 0xBB177C0EDC365BC8uL, 0x72C287D397507D2EuL)),
                ((+1, -7, 0xF0717084FAFAD988uL, 0x851B24EF359515B9uL), (+1, -6, 0xBFFFF2D1514AC2ACuL, 0xAA4AEB3A767D854DuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x97C7D8127BB4003DuL, 0x2FA29C9908C74E40uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xA5DD01C17C5F9C1BuL, 0xA31859314DDCF52BuL), (-1, 0, 0xC92BC2928000EE2BuL, 0x111EB4778CE144C8uL)),
                ((-1, 1, 0xDE281A11FADB828DuL, 0x6181F44658F061C8uL), (+1, 0, 0xAD79E4AF5AD242EDuL, 0xAE40F8A13027EB87uL)),
                ((+1, 1, 0xC40A54494A14216AuL, 0x9F6F892217CF3930uL), (-1, 0, 0x977A2C1C4CD84705uL, 0x1D60EDF391FE4EB3uL)),
                ((-1, -1, 0xE2EDA68214505D6CuL, 0xE4FACB9994636F76uL), (+1, 0, 0x87759FF44A19894AuL, 0xFD0DE9D0BFB8647CuL)),
                ((+1, -1, 0x8ECA6DBC7EB3369DuL, 0xC44E4A015007CE79uL), (-1, -3, 0x9246B22E7C41F145uL, 0x259ECE59CD7B3A03uL)),
                ((-1, -2, 0xDF5E3220824F9646uL, 0x43C65671FABE8C42uL), (-1, -3, 0xCB06E182F6D8D587uL, 0x9515340810D46B49uL)),
                ((+1, -3, 0xE4FE2B6C6B864723uL, 0x11767B54FF76BF84uL), (-1, -5, 0xE994C6C7B8A42FD1uL, 0x694CDE362B12A810uL)),
                ((-1, -3, 0x8A77ACBB75433602uL, 0x4B4BE41D4ABE8F19uL), (+1, -4, 0xBFBD85609FB72C34uL, 0xE1CB70F931ABE6D1uL)),
                ((+1, -4, 0xF3131E35B1033CACuL, 0x407BE9AE3A09F44AuL), (-1, -6, 0x921F95FF83F1E138uL, 0x5F16AB5C17EBD19FuL)),
                ((-1, -4, 0x803C515172AF6E42uL, 0xFD94B8E32B18908BuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xDBD33A19C542A1DAuL, 0x8D65B27C0B9D255DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xA391B78853FA295DuL, 0xF9063E8A937E0E08uL), (-1, 0, 0xCE60FBCF69FDE883uL, 0x64A4D73F112E858CuL)),
                ((-1, 1, 0x993F1C526AC2C32AuL, 0xC951C640643AADCBuL), (+1, 0, 0xFAE33F98B1BC798EuL, 0x3C918EE5AAE489D6uL)),
                ((+1, 1, 0xF05D4ADE2507EB96uL, 0x5F7076AD53E73D75uL), (-1, 1, 0x8395D72FE12DC44AuL, 0xE34DB167ADBA007FuL)),
                ((-1, 1, 0x9EF89B319B67E083uL, 0xDA82E63D6EA79ADBuL), (+1, 1, 0xAD7E420FCD193DBFuL, 0xC68C40CCF66049F3uL)),
                ((+1, 1, 0xE79EFDC3DE8F9642uL, 0x08AE44E41F4911A2uL), (-1, 0, 0xCC7C789BB37EB92AuL, 0xF5779D87FA19ECA7uL)),
                ((-1, 1, 0x9DBE650026B63F55uL, 0xAD56F23A816075F8uL), (+1, -1, 0xDBB6E25F9563DE52uL, 0x8B49CC88DB6BA187uL)),
                ((+1, 0, 0xB2737FBEB082DC1FuL, 0xE1F2500A40DED22CuL), (-1, -1, 0x8156A6BBD234B1A6uL, 0xA216076FE43D87B1uL)),
                ((-1, -6, 0xCB2B6CEFCEAC0D13uL, 0x9E64D71451E21681uL), (+1, -3, 0xFB5E6529758C0D08uL, 0xE6575CF735EAF077uL)),
                ((+1, -3, 0xEFBACCE48A68B1E3uL, 0xB63BFBF2AE42BEDBuL), (+1, -3, 0xD3B0E4183D4DFE39uL, 0x2528887536DED25CuL)),
                (Zero, (-1, -3, 0x82429A5D9D0461C8uL, 0x8A167A68FBDE16BFuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p5_0p625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0x8CC1D3B4333579C0uL, 0x76AF5D87EEDE2AC1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -5, 0xBB2E97B299B45089uL, 0x7A1C18970BA1BE34uL), (-1, 0, 0xDB92AD271BB7F073uL, 0xD0C6FE4786028AEAuL)),
                ((-1, 1, 0x84E577B89E7F08ABuL, 0xDB2F8DF9F2543836uL), (+1, 0, 0xC4FA1D99F07789ECuL, 0xDACD9D101BA8403CuL)),
                ((+1, 2, 0xB98B446DC372BA8AuL, 0x00E1B5D536A902D5uL), (+1, 0, 0x80B9E5E012E0D1CFuL, 0xB52029BA94734E3EuL)),
                ((-1, -1, 0xC9870C310A1638DDuL, 0xD07570886A407B9DuL), (-1, 0, 0x8E4B8CCD3510EAECuL, 0x82086EAC8CB491CBuL)),
                ((-1, 1, 0x901F7DE270A65443uL, 0x0A76F4C8ED9C7446uL), (+1, 0, 0xB1030DB5F2F39147uL, 0x02A525657ED9DD6DuL)),
                ((+1, 2, 0xB62CA65E3A972E7AuL, 0x66577F0946194C81uL), (-1, -1, 0xD3DC73C1AEB3E7EDuL, 0xD594A48523F8E98EuL)),
                ((-1, -1, 0xE80206BA77462E6EuL, 0x6E23F386C5684E0CuL), (+1, 0, 0xF08055FFCE916758uL, 0x3FB4548B4E4CEAA0uL)),
                ((+1, 0, 0xA98AE25287E01077uL, 0xA628FB1DB0CF60AFuL), (-1, -1, 0x882E6260F10E0A84uL, 0xEDAAAC4E20C0D888uL)),
                ((+1, -1, 0x9C2E821888AEE134uL, 0x5041D8369C9BE55FuL), (-1, -4, 0xE5BD5A64C084758BuL, 0x6ECA0252C8725268uL)),
                ((+1, -2, 0xBF144168F1ECBD85uL, 0xCE7607ED0E712D21uL), (+1, -2, 0xC6FDD6E23F0BB7CBuL, 0x702E61A5D5C6C6EEuL)),
                ((+1, -2, 0xA58B22A854F10CDCuL, 0xD629BCCDDC4C3B1FuL), (-1, -4, 0xA19FFC016D3008A7uL, 0x81A6DF2E59CCE700uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p625_0p75 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xA936BFEDCA455287uL, 0x838CA773D477FD83uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0x8B3FEF4F04635544uL, 0xA89492BDEA737579uL), (-1, -2, 0xF97CDE1A5B404264uL, 0x77153CC994013474uL)),
                ((-1, 2, 0x800733B221F7272AuL, 0x65D4FBA4BC0196CCuL), (-1, 1, 0x8D59C468FDB93157uL, 0x0992BF16FDA5D2AEuL)),
                ((+1, 3, 0xD61D67DC1E18D2E0uL, 0x9C7025D4EF69F974uL), (+1, 3, 0xC112DC707335A994uL, 0x039CA4E79FC61094uL)),
                ((-1, -2, 0xAD76B1FFDBCCB9E2uL, 0x43740632EC9E80C2uL), (-1, 4, 0x82AA88619C868B9EuL, 0x4238F48AB3B0E733uL)),
                ((-1, 3, 0xAD779AE9935711D7uL, 0xA477B5AE5F97CC35uL), (+1, 4, 0x90D08B27B30C8620uL, 0xC1F21CE5255702E7uL)),
                ((+1, 5, 0x88EB2F1B05063332uL, 0xE48D833245A20EC1uL), (-1, 3, 0xAB44D34405A7AA63uL, 0x76D2EE08E10AA144uL)),
                ((-1, 4, 0xA269071CB7F1316BuL, 0x7F95AC56CDFEA028uL), (+1, 3, 0xAF72AEDAD6AF6BB7uL, 0x43A431DBDE7834C0uL)),
                ((+1, 3, 0xFEDC4A7A504152BAuL, 0x380B76B2103A2B3BuL), (-1, 2, 0xE6F7A3067F30F79CuL, 0x809CB5AD0936F4F6uL)),
                ((-1, -1, 0xDAA84B192AB445CCuL, 0xACF1A7EB9AA87F7AuL), (+1, 2, 0x841D54515A7677FBuL, 0x7BBDB6CFE7C83C11uL)),
                ((+1, 1, 0xB5D9CAFC880522A3uL, 0xAF8CFE28B840A584uL), (+1, -2, 0xAF1050B3219E28FEuL, 0x0C3F558EBA47AB70uL)),
                (Zero, (-1, 0, 0x805C41D6919468B7uL, 0xE1D029514A2CFAE1uL)),
                (Zero, (+1, -2, 0x8640932EA5825BCAuL, 0x0356D333E1216216uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p75_0p875 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xC4BEF02EC2EA7006uL, 0x91CAC96D130F4E35uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0x8039C9BE6610EC42uL, 0x7E9577B75C2D05CFuL), (-1, 0, 0xE3530ED9B95F1F48uL, 0xD819FF8BB9647B0CuL)),
                ((+1, 2, 0xC6ED6BDF62189F1CuL, 0xE3CCEC79E3308CDBuL), (+1, 2, 0xBCFD699CE56E7079uL, 0xC2D16FB418CB5454uL)),
                ((+1, 2, 0xB43BE1C8E1CCEEC0uL, 0x5B01682889398DA9uL), (-1, 1, 0xE03BA9DC8E8F85AEuL, 0x0028A6CAB8F74D9EuL)),
                ((+1, 2, 0x9F9A318FF498E35FuL, 0xEFB59FB62053A08CuL), (+1, 2, 0xFCDFE0CD17C2BA70uL, 0xCD6F4DB5B46B21F3uL)),
                ((+1, 4, 0xAB7FF1F49F4ED501uL, 0xD6BFDF6D12CF31F5uL), (+1, 0, 0xB3D2AF9AFE297364uL, 0xF6635E38D0134D28uL)),
                ((+1, 3, 0x9E46CFD3D8697EF8uL, 0x36AC66B5186B1686uL), (+1, 2, 0xC595CE919D5E9494uL, 0xB673AE7FD813E66EuL)),
                ((+1, 4, 0x9A03E3CB8AFDBBFAuL, 0x68DEB123BE156E3EuL), (+1, 1, 0x90CD11E1C9EAA100uL, 0x62A223097E2368B3uL)),
                ((+1, 3, 0xE330FDAB4DC77E10uL, 0x5F3F21DC330ECE6AuL), (+1, 1, 0xE3A9F47E35A4F295uL, 0xEFA50DE87D27E437uL)),
                ((+1, 3, 0x9A82FDAFD71C4148uL, 0xE40CD7033F01F875uL), (+1, 0, 0xA5D13FBFEC55683FuL, 0x391621FC43126D28uL)),
                ((+1, 2, 0x9F938A60FE94787CuL, 0x2994AA339A77E641uL), Zero),
                ((+1, 0, 0xCF78981370CF2BC1uL, 0xEDAC63E2EB48F030uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p875_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xE117F0965CBE0650uL, 0xB4535458CFCAD9A1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x83EC6AEF55D448F7uL, 0x090CD9893ED25755uL), (-1, -1, 0xBFB49596CDF8D458uL, 0x034AA97111ECCFA9uL)),
                ((+1, 3, 0x99DAB2A12747F48AuL, 0x6FCDBC518A05E6E5uL), (+1, 2, 0xBCDE2B062902AD88uL, 0x8F1AADB123EF7A27uL)),
                ((+1, 3, 0xD904B8A4D87C7EC9uL, 0x561CFE2FB1E73E37uL), (+1, 0, 0xA70242A2B9296EF4uL, 0xCF1C1DFF81A2090FuL)),
                ((+1, 4, 0xADBF6DC91152A507uL, 0xFB69C21A50F9C10EuL), (+1, 3, 0x9CAC08D14A6B4404uL, 0xED275ABB2A67A21EuL)),
                ((+1, 5, 0xB5E7B3250EB7CE19uL, 0xFA1E7CB21440BD31uL), (+1, 3, 0xC2A7A7E3AECDA2EDuL, 0x827E05859DFF2FC0uL)),
                ((+1, 5, 0xB5EBF788EAA1E9F4uL, 0xA3DEE37D4DBFEF9DuL), (+1, 3, 0xAE0D187890618D38uL, 0xAD80A9D7DFE844D4uL)),
                ((+1, 5, 0xDD07AE44D7C58517uL, 0xA0E4BA0EA08D364BuL), (+1, 3, 0xD3B3167808F17B5FuL, 0x472F1934558F7C4BuL)),
                ((+1, 5, 0xDAE89514395ACA5BuL, 0xF4E096EF39B421F7uL), (+1, 3, 0x997DA6200BAA7D79uL, 0x4AD41EDA08B40875uL)),
                ((+1, 5, 0x89BF550C1F0851F5uL, 0x314AA410DA454EF2uL), (+1, 2, 0x9DC5BAA2B8921128uL, 0x58B8D9AD398682F0uL)),
                ((+1, 4, 0xA41D81BB1D67AB58uL, 0x648984406E9154B8uL), (+1, -1, 0x93FA550EA206CFEBuL, 0xCC9C8F9848023CD3uL)),
                ((+1, 2, 0xDA001C6CCE1B5C27uL, 0xB7B2DF1D54F21E26uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_1_1p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xFFD3CEABA73039DAuL, 0x708244915B2729ADuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0xC66B40B77C915223uL, 0x0FC5B69B91E3261BuL), (+1, 1, 0x85DBDBF58C3915C4uL, 0x7A2B92466F28C237uL)),
                ((+1, 4, 0x919FF6F07D03FB64uL, 0xCD6FDD0D0A32FA9EuL), (+1, 2, 0xD1E55032B8BF62CBuL, 0x148A1B967E06312CuL)),
                ((+1, 5, 0xE301ED5015ACFBBFuL, 0xEFFCEC1FF64DD4FAuL), (+1, 4, 0xA4FD74B317785DF3uL, 0x8FB70560D4B0023DuL)),
                ((+1, 6, 0xA425D04FE5A59A17uL, 0x7F1EFE5B423E5E18uL), (+1, 4, 0x8A8589B8C75BC48DuL, 0x6DB45F770F530C16uL)),
                ((+1, 7, 0xC6A71B004934164DuL, 0xDB52B04E592F1FB1uL), (+1, 6, 0x902ED55D1CFE2E09uL, 0xAF96A00B86D04C0AuL)),
                ((+1, 7, 0xF33A10E3167078C0uL, 0xE198FF3E9743F7A4uL), (+1, 5, 0x9900F7ED6C4398C8uL, 0xB6E863460D58F2F8uL)),
                ((+1, 8, 0xB1CA96E601679FEDuL, 0x8DA99ECDF3796554uL), (+1, 6, 0xDD9FD798DBC3DF12uL, 0x03912ACD1251AA9CuL)),
                ((+1, 8, 0xD653376AD2812778uL, 0x6BE574382C0FDFB1uL), (+1, 6, 0x8D8BD871F48D6B98uL, 0x9D8A02DA15BEEF4CuL)),
                ((+1, 8, 0xBBF0DBC4376310B0uL, 0x52018DE752EC6BBEuL), (+1, 6, 0x9F5B1D09C19D1C80uL, 0xB628955F4B844F49uL)),
                ((+1, 8, 0xB9754DF02FC15FA9uL, 0x6B576B008EEEB2DAuL), (+1, 5, 0xD99F80A4AE157675uL, 0x91EA22BBB3654112uL)),
                ((+1, 7, 0xE3D86299851FF9E1uL, 0x94FF85AD41B0475FuL), (+1, 4, 0xDBDF12D5CC92E025uL, 0x6FA98E3221E17CABuL)),
                ((+1, 6, 0xF584EFBBEDE69CC5uL, 0x4E314495ADF71960uL), (+1, 2, 0xEF7BD252D9BD6E31uL, 0x55462E850E61FFA2uL)),
                ((+1, 5, 0xC0CC8EC8E0B4E405uL, 0x820C1DCBE9D175B0uL), (-1, 0, 0x9721E6A482C46F4FuL, 0x71CED1DE162D215FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_1p25_1p75 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 1, 0xA3FFD574012002EBuL, 0xEF5F62BA425B2F41uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 3, 0x8BC962687C4C58B7uL, 0xC58FCA1D105644D0uL), (+1, 1, 0x9BC247F39BD452C6uL, 0xA9D95AB1CF7B1C18uL)),
                ((+1, 5, 0x8D3E6B01B5E43FF2uL, 0x7819F137E1B23778uL), (+1, 3, 0xB0424C1DC7A7578BuL, 0x077C79EA493D84E2uL)),
                ((+1, 6, 0xB42C4D30A17EEDEBuL, 0xCFF21C471B27D647uL), (+1, 4, 0xBB9EB7DB839678CFuL, 0x19B9CC3CC6D83A69uL)),
                ((+1, 7, 0xD6BC7CC5D3034B9AuL, 0x25E957948E8A065DuL), (+1, 5, 0xE2967E75259942D4uL, 0xDE748BD8CFABCBABuL)),
                ((+1, 8, 0xD6BE45F8652229CCuL, 0x858F7C0F35B385B5uL), (+1, 6, 0xCE32484F04B968E8uL, 0x7E9F8F03BAA9002EuL)),
                ((+1, 9, 0xBDD9E571FA1EA927uL, 0x3B0CB77D9D2B1BBDuL), (+1, 7, 0xADE1B46EA41BEE3EuL, 0xF9E5AC516FDF9DE2uL)),
                ((+1, 10, 0x984CB8940E0DB72DuL, 0xDFC2F01C68463B5BuL), (+1, 8, 0x84997EE6AC2A71C4uL, 0x7610359DD02D2F56uL)),
                ((+1, 10, 0xD67B615E87E88A6AuL, 0x2E39B6C681F8B8B9uL), (+1, 8, 0xABC3A24114006FF3uL, 0x8305916A1A620225uL)),
                ((+1, 11, 0x89C2FDF5929D9433uL, 0x9E005133652A4AA5uL), (+1, 8, 0xD2948BF6DB8248B1uL, 0x6C9E6A9617575B52uL)),
                ((+1, 11, 0x9C79E853E032563CuL, 0x4CD4BB4FBD23DBFBuL), (+1, 8, 0xD7F73E5EF65F5FAAuL, 0xEC24C324090FEC8EuL)),
                ((+1, 11, 0x9DD289741387BB46uL, 0x28D226DFE31D86C3uL), (+1, 8, 0xC843807E7B741A1FuL, 0xD103CE22C324188DuL)),
                ((+1, 11, 0x8A1DF669B044CF7FuL, 0x6DC0EBBEBDAFA712uL), (+1, 8, 0x96D9BFEAA1DB2B72uL, 0xE4AE21D6D5C79EA9uL)),
                ((+1, 10, 0xCC6019361E2CFEF0uL, 0x9E522AE0427EB4CAuL), (+1, 7, 0xBCCF93C4E1DB9617uL, 0x2C238DD75105C110uL)),
                ((+1, 9, 0xF2DE322CC162A8C5uL, 0xA947CF413CD92BBEuL), (+1, 6, 0x984C380FA39B4952uL, 0x06E4150F5A49A9CBuL)),
                ((+1, 8, 0xBD15957788573B71uL, 0x2A79F84D81A5FB33uL), (+1, -3, 0x84E8746D0067097EuL, 0xB4B2784876EEF50FuL)),
                ((+1, 6, 0x99D7CB70A4ABB449uL, 0xD1C78B8D56A7390AuL), (-1, -7, 0xACA07E309D647D2FuL, 0x2EDD9AEB641D5439uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_1p75_2p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 2, 0x81FFFFFFFFFEADB3uL, 0x2330B71727BE7AD8uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 6, 0xBB31E00256D9B4A5uL, 0x1A007483E882D562uL), (+1, 4, 0xB16C2F4541CFAA91uL, 0x2A472B8030651CA9uL)),
                ((+1, 10, 0x8C9F98C0CAED6E0EuL, 0xFEC42878B33B51DAuL), (+1, 8, 0x80C88BFB1502D537uL, 0x2BE0ED8D3FC2EAE2uL)),
                ((+1, 13, 0x91595CC3AC63C2DCuL, 0x3B5818EE941B278AuL), (+1, 11, 0x80E71AF0D0323735uL, 0xDB6F49C11C9D80B1uL)),
                ((+1, 15, 0xE687F36A47B98CB2uL, 0x85E06CFE1913F680uL), (+1, 13, 0xC63AE8CA7A22554DuL, 0xB7D88C41BDD3B0F7uL)),
                ((+1, 18, 0x947A9056EDAC0BA2uL, 0x3A5FB29F5F6B3212uL), (+1, 15, 0xF7B5CA2D18F3895CuL, 0x10C86B167EB3EAADuL)),
                ((+1, 20, 0xA0A41ED35C89CD12uL, 0x9B940A50E46D77F5uL), (+1, 18, 0x81F7E747B9F768C0uL, 0xE3E6F1B76D178BF5uL)),
                ((+1, 22, 0x9509A0811086E66CuL, 0x68702D34A83B69D5uL), (+1, 19, 0xE9B12FEA61BF9390uL, 0x4D2354A42DCC74E1uL)),
                ((+1, 23, 0xF0232150149EB54FuL, 0xF58C620E884D0D0EuL), (+1, 21, 0xB61C2EDF151D556BuL, 0x1086B8E6FC0947A1uL)),
                ((+1, 25, 0xA90DA870F77782A9uL, 0x2928790C7767C04FuL), (+1, 22, 0xF744615B944DFA47uL, 0x7207F0886DD86644uL)),
                ((+1, 26, 0xD030DC2283D57957uL, 0xCCC31CFAEEEA4B58uL), (+1, 24, 0x9220835239D649C1uL, 0xCC5A222209F9968CuL)),
                ((+1, 27, 0xDEA8C45B381E2D88uL, 0x2CAE7B4D3B39269BuL), (+1, 25, 0x94ADB1180223F76BuL, 0x7980A757878C6978uL)),
                ((+1, 28, 0xCD72B3F101742DCBuL, 0x712F5D36AD5415B0uL), (+1, 26, 0x813FAF834FA9329EuL, 0xDC77A23F68585FCAuL)),
                ((+1, 29, 0x99176B6BFD184279uL, 0x64BD2C12DF8741DBuL), (+1, 26, 0xABD21EA7765DD8A8uL, 0x234784D70076D09DuL)),
                ((+1, 29, 0xC6175E21AF383242uL, 0xC66AB3DA9144B1BEuL), (+1, 26, 0xD23DFD999C825B2FuL, 0xAFAF547D97BECD17uL)),
                ((+1, 28, 0xE2EAC5904689ADEDuL, 0x5EF1F2603D9E6F72uL), (-1, -27, 0xDAD03579BF02DE90uL, 0x875676C3E2A79EBFuL)),
                ((+1, 26, 0xD23DFD999C820F2EuL, 0x1FE34C874D6CD370uL), Zero),
            }));

            public static ddouble Value(ddouble x) {
                if (x < 1d) {
                    return 0d;
                }

                ddouble y;

                if (x < 6d) {
                    ddouble u = Sqrt(x - 1d);

                    if (u <= 0.25d) {
                        y = ApproxUtil.Pade(u, pade_0_0p25);
                    }
                    else if (u <= 0.375d) {
                        y = ApproxUtil.Pade(u - 0.25d, pade_0p25_0p375);
                    }
                    else if (u <= 0.5d) {
                        y = ApproxUtil.Pade(u - 0.375d, pade_0p375_0p5);
                    }
                    else if (u <= 0.625d) {
                        y = ApproxUtil.Pade(u - 0.5d, pade_0p5_0p625);
                    }
                    else if (u <= 0.75d) {
                        y = ApproxUtil.Pade(u - 0.625d, pade_0p625_0p75);
                    }
                    else if (u <= 0.875d) {
                        y = ApproxUtil.Pade(u - 0.75d, pade_0p75_0p875);
                    }
                    else if (u <= 1d) {
                        y = ApproxUtil.Pade(u - 0.875d, pade_0p875_1);
                    }
                    else if (u <= 1.25d) {
                        y = ApproxUtil.Pade(u - 1d, pade_1_1p25);
                    }
                    else if (u <= 1.75d) {
                        y = ApproxUtil.Pade(u - 1.25d, pade_1p25_1p75);
                    }
                    else {
                        y = ApproxUtil.Pade(u - 1.75d, pade_1p75_2p25);
                    }
                }
                else {
                    y = x;
                }

                return y;
            }
        }
    }
}
