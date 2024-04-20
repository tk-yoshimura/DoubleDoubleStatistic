using DoubleDouble;
using DoubleDoubleIntegrate;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleDistribution {
    public class SkewNormalDistribution : LinearityDistribution<SkewNormalDistribution>,
        IAdditionOperators<SkewNormalDistribution, ddouble, SkewNormalDistribution>,
        ISubtractionOperators<SkewNormalDistribution, ddouble, SkewNormalDistribution>,
        IMultiplyOperators<SkewNormalDistribution, ddouble, SkewNormalDistribution> {

        public ddouble Mu { get; }
        public ddouble Sigma { get; }
        public ddouble Alpha { get; }

        private readonly ddouble pdf_norm, cdf_norm, erfc_scale, sigma_inv, s;

        public SkewNormalDistribution() : this(alpha: 0, mu: 0, sigma: 1) { }

        public SkewNormalDistribution(ddouble alpha) : this(alpha, mu: 0, sigma: 1) { }

        public SkewNormalDistribution(ddouble alpha, ddouble mu, ddouble sigma) {
            ValidateShape(alpha, alpha => Abs(alpha) <= 16d);
            ValidateLocation(mu);
            ValidateScale(sigma);

            Mu = mu;
            Sigma = sigma;
            Alpha = alpha;

            pdf_norm = 1d / (sigma * Sqrt(2d * PI));
            cdf_norm = 1d / Sqrt(2d * PI);
            erfc_scale = Alpha / Sqrt2;
            sigma_inv = 1d / sigma;

            s = Alpha / Hypot(1, alpha);
        }

        public override ddouble PDF(ddouble x) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv;

            ddouble pdf = pdf_norm * Exp(-u * u * 0.5d) * Erfc(-u * erfc_scale);
            pdf = IsFinite(pdf) ? pdf : 0d;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            if (IsNaN(x)) {
                return NaN;
            }

            ddouble u = (x - Mu) * sigma_inv;

            ddouble f(ddouble u) {
                ddouble c = Exp(-u * u * 0.5d) * Erfc(-u * erfc_scale);
                return c;
            }

            if (interval == Interval.Lower) {
                ddouble cdf = Erfc(-u / Sqrt2) / 2 - 2 * OwenT(u, Alpha);

                if (cdf < 1e-5) {
                    ddouble eps = Ldexp(Erfc(-u / Sqrt2), -94);

                    cdf = cdf_norm * GaussKronrodIntegral.AdaptiveIntegrate(f, NegativeInfinity, u, eps, discontinue_eval_points: 2048).value;
                }

                cdf = IsFinite(cdf) ? Clamp(cdf, 0d, 1d) : (x < Mu) ? 0d : 1d;

                return cdf;
            }
            else {
                ddouble cdf = Erfc(u / Sqrt2) / 2 + 2 * OwenT(u, Alpha);

                if (cdf < 1e-5) {
                    ddouble eps = Ldexp(Erfc(u / Sqrt2), -94);

                    cdf = cdf_norm * GaussKronrodIntegral.AdaptiveIntegrate(f, u, PositiveInfinity, eps, discontinue_eval_points: 2048).value;
                }

                cdf = IsFinite(cdf) ? Clamp(cdf, 0d, 1d) : (x < Mu) ? 0d : 1d;

                return cdf;
            }
        }

        public override ddouble Mean =>
            Mu + Sigma * s * Sqrt(2d / PI);
        public override ddouble Mode => ModePade.Value(Alpha) * Sigma + Mu;

        public override ddouble Variance =>
            Sigma * Sigma * (1d - 2d * s * s / PI);

        public override ddouble Skewness =>
            (4d - PI) / 2 * Cube(s * Sqrt(2 / PI)) / Cube(Sqrt(1 - 2 * s * s / PI));

        public override ddouble Kurtosis =>
            2 * (PI - 3d) * Square(Square(s * Sqrt(2 / PI))) / Square(1 - 2 * s * s / PI);

        public override ddouble Entropy {
            get {
                ddouble f(ddouble x) {
                    ddouble pdf = pdf_norm * Exp(-x * x * 0.5d) * Erfc(-x * erfc_scale);

                    if (pdf == 0d) {
                        return 0d;
                    }

                    ddouble y = -pdf * Log(pdf);

                    return y;
                }

                (ddouble value, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(
                    f, NegativeInfinity, PositiveInfinity, 1e-28, discontinue_eval_points: 1024
                );

                value += Log(Sigma);

                Debug.WriteLine($"Entropy integrate err: {err}");

                return value;
            }
        }

        public static SkewNormalDistribution operator +(SkewNormalDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu + s, dist.Sigma);
        }

        public static SkewNormalDistribution operator -(SkewNormalDistribution dist, ddouble s) {
            return new(dist.Alpha, dist.Mu - s, dist.Sigma);
        }

        public static SkewNormalDistribution operator *(SkewNormalDistribution dist, ddouble k) {
            return new(dist.Alpha, dist.Mu * k, dist.Sigma * k);
        }

        public override string ToString() {
            return $"{typeof(SkewNormalDistribution).Name}[alpha={Alpha},mu={Mu},sigma={Sigma}]";
        }

        private static class ModePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                (Zero, (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xCC42299EA1B28468uL, 0x7E59E2805EBEB97CuL), (+1, -1, 0xD08DE3E8DC210445uL, 0x04E3DC0B3295F7FEuL)),
                ((+1, -1, 0xA666F81A4640AFBFuL, 0x7432113B758FD719uL), (+1, 1, 0xFAFA61A725E6362BuL, 0xD4D9F1C60FE01B1DuL)),
                ((+1, 1, 0xA7BE26484F0B0BF3uL, 0x75022EB9A6A0CE86uL), (+1, 1, 0xAFF02C8863796476uL, 0x254F458923473FCDuL)),
                ((+1, 0, 0xE3CA0BBA22CD69E3uL, 0x0A1EEB9A1F3A92ACuL), (+1, 2, 0xC384CC6D591069C9uL, 0x7D0B50B2DED64CC9uL)),
                ((+1, 1, 0xD1A7DF03118BB1A9uL, 0xD3EA71DE32A0E40DuL), (+1, 1, 0xE5EE23D0FD2545FEuL, 0xD43098858C46373CuL)),
                ((+1, 0, 0xE52307B92B2E7773uL, 0x1ABF82207DB6F030uL), (+1, 2, 0x995812C09819C6DBuL, 0xCF8F58573A818B14uL)),
                ((+1, 0, 0xF73A9DCC2074D81BuL, 0xAEDBDB64F35B7978uL), (+1, 1, 0x919A7F77294A0CF7uL, 0xB0860AA9EFD3C58DuL)),
                ((+1, -1, 0xCD63BFCFCEA71ADCuL, 0x00DC19CD0775B2A3uL), (+1, 0, 0xFC5662CB04B19DD0uL, 0xAC4D360F0CB1B350uL)),
                ((+1, -1, 0x896868B2E972EDF5uL, 0x0328122BECC08AC2uL), (+1, -1, 0xB5C61FB31E49E4A0uL, 0xFD1D49F2EE5A31D0uL)),
                ((+1, -3, 0x9BF4B2BB0606D6C0uL, 0x5F51BB1126BE81C2uL), (+1, -2, 0xCCD809E3493BCC3CuL, 0x7D98DECF9D1D6ECEuL)),
                ((+1, -5, 0xF5BB8F77E2137EF3uL, 0x61CDBEDAAD0CC686uL), (+1, -4, 0xC705BE05DC5557E5uL, 0xCD34945AFECCF7C3uL)),
                ((+1, -7, 0x92E6A6EEF3CAFDDAuL, 0x10BA8E153369E0E5uL), (+1, -5, 0x86852D34C872CF79uL, 0x3BA4CBD524B36EB7uL)),
                ((+1, -10, 0xD41A8059E7871661uL, 0xD3E454012847A502uL), (+1, -8, 0x82FC6BE0007448F0uL, 0xA62B9809A8B1D1C6uL)),
                ((+1, -16, 0xBAC98DCA1A38256FuL, 0x66ED521574716CECuL), (+1, -11, 0xA09AF4EFE6EF43A5uL, 0x143ADA4747A912B1uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_0p5_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xB0F4306E16B4CA1EuL, 0x47C1CB7DDBED83DEuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xECF298F7F4C11159uL, 0xFBD9224B4D1F298BuL), (+1, 1, 0xF80D10C4F18CADF9uL, 0xF067C5D891BFB9A7uL)),
                ((+1, 2, 0x9101AB465720BDC6uL, 0x704120F516E05024uL), (+1, 3, 0x89BCFF653DC6D412uL, 0x42727DBE63E898C2uL)),
                ((+1, 2, 0xE5399A78B27CB7FFuL, 0x709FBF5607950ED5uL), (+1, 3, 0xCBAA4BA9259E6620uL, 0x8D3B5F584499D006uL)),
                ((+1, 2, 0xFCFEA68DB670E20BuL, 0x2EB4E7028A028854uL), (+1, 3, 0xDD36C4057595D938uL, 0x3F50528042DAFB2DuL)),
                ((+1, 2, 0xCCBCAF9DC4801C8BuL, 0xE1B86D680C912311uL), (+1, 3, 0xB47F13A5DECAA238uL, 0xF4AE7D55EF87BE27uL)),
                ((+1, 1, 0xF358BEBB57EFB4F4uL, 0xBAAFB54827B3A5B0uL), (+1, 2, 0xE0E3BA50B7657D72uL, 0xEF2EA56F79169617uL)),
                ((+1, 0, 0xD23F81B4898D128DuL, 0xE71A99916432EAD4uL), (+1, 1, 0xD2EDE0F9B99D614DuL, 0x454CB30BD57F62C3uL)),
                ((+1, -2, 0xFA7B24E554F7575DuL, 0x515FF570B431D252uL), (+1, 0, 0x918878D935FD4494uL, 0xB1A09BCF16D6D35EuL)),
                ((+1, -4, 0xBB0699D608467CF3uL, 0x823EF19366CC885BuL), (+1, -2, 0x8A52A18495C90A1CuL, 0xCA03BD9446EA2D9FuL)),
                ((+1, -7, 0x82ACFB0A9FE0526AuL, 0xBEFC1F8AD7133748uL), (+1, -5, 0xA1D29DBCE3054F04uL, 0xF50D4207C4A52752uL)),
                ((+1, -15, 0xD8706F7A88778B04uL, 0x9C091D9D010877F6uL), (+1, -9, 0xA1DBFD2EF2934C3FuL, 0x8C84B3057AC5406DuL)),
                ((-1, -21, 0xF7F2F8EDAABB9440uL, 0xDDE488B1548AB0D6uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x818CC9224D9D0B57uL, 0x9BB66FAAC33F193FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0xA60BB7E8E067AD7DuL, 0xA700D48449B03FDCuL), (+1, 2, 0x99BC7AFB03ADB01BuL, 0xC8808E69C4EA9519uL)),
                ((+1, 2, 0xC5C8481C6831F6E4uL, 0x68AE07B00E1F77BFuL), (+1, 3, 0xB2435DFAF052DAC6uL, 0xE69002A695EC1F31uL)),
                ((+1, 3, 0x9063DD2A56233FF1uL, 0x19973353EF6C3974uL), (+1, 4, 0x82063AECCD7091AAuL, 0x77E1BA882BE98628uL)),
                ((+1, 3, 0x8F081D1433DDF273uL, 0xAC88860855DB0202uL), (+1, 4, 0x83B336E116AEAE33uL, 0x49EB9DEA28EA1E47uL)),
                ((+1, 2, 0xC8F878A0D2F505A2uL, 0xEC60C4BD9DF9E59FuL), (+1, 3, 0xC1B61EF0429E8162uL, 0x519B4AEE03526ADDuL)),
                ((+1, 1, 0xCB4770EA6A8D0537uL, 0xA02467CCF97FC099uL), (+1, 2, 0xD2AD240FE901D176uL, 0x51DC17A4F0DEE3CCuL)),
                ((+1, 0, 0x931BDB4F10260AFDuL, 0xED7209EE0F9BF6EFuL), (+1, 1, 0xA97F8D7F11635BF0uL, 0x58EE46F73A3907BCuL)),
                ((+1, -2, 0x943E90B9E0AA6BEFuL, 0x6BD26C3AAD9FB1C4uL), (+1, -1, 0xC6AFAA2F87D38B75uL, 0xF23CA241030F5656uL)),
                ((+1, -5, 0xC4ECC5C7B672C541uL, 0xB686F1E879145E1AuL), (+1, -3, 0xA40948E4B0A033B9uL, 0x83E2A38E17BCE0A8uL)),
                ((+1, -8, 0x9C1815D535099ED0uL, 0xF2269688743C8C75uL), (+1, -6, 0xB3A9FD8B6088D614uL, 0x13BB7BBF3AF26E97uL)),
                ((+1, -13, 0xF5E72DB4B9C1A830uL, 0x8C37A9E0DD4BD185uL), (+1, -10, 0xEB39F044A54242E2uL, 0x72691EC873539BCDuL)),
                ((+1, -18, 0x8297B607FFAADC90uL, 0x6F864124041B6069uL), (+1, -14, 0x982D120202C16099uL, 0x5732DA6CCC7C008CuL)),
                ((+1, -29, 0xAF5B2065006FAA8FuL, 0x04A1F4710E48C4A9uL), (+1, -20, 0x803887900796D69FuL, 0x155F35142AD4AC30uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x87DFC3EB8C2CEB75uL, 0x32366FABEB75DFC7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xDDF77C59996BCA05uL, 0xEF25D350D1FEB347uL), (+1, 1, 0xD687F97A08AB5FABuL, 0xC9DA857B5EEDD422uL)),
                ((+1, 1, 0xA2CCCC0815BBC54CuL, 0x83CE9C83675F4281uL), (+1, 2, 0xA41E2E436955F630uL, 0x6060B26CAB7CD5BDuL)),
                ((+1, 1, 0x8CE2E86800EA19C0uL, 0xEF55225136856087uL), (+1, 2, 0x96C484FFFBBC45D2uL, 0xA5C77EDA0AEE979CuL)),
                ((+1, 0, 0x9EEF2D648F4CFBFFuL, 0x6DC93E7FABB511CDuL), (+1, 1, 0xB8347D9268AD923EuL, 0x7613B915231D9D8CuL)),
                ((+1, -2, 0xF3AF5BAFCA7CC6A5uL, 0xA5BDBA7DB11A0A8AuL), (+1, 0, 0x9CAA9FCFF9C548A4uL, 0xDF1E5B0760771E22uL)),
                ((+1, -3, 0x806D6248990A7D4CuL, 0x1FC4BA48A6D8C9D8uL), (+1, -2, 0xBCCDDCF8D4469174uL, 0x1885DBBFDD9FADE6uL)),
                ((+1, -6, 0xB8352C8D8DE9246DuL, 0x1AFD8E5265CF72F4uL), (+1, -4, 0xA0FE6A3BF408760AuL, 0x1738714FEB07322AuL)),
                ((+1, -9, 0xAE35F4A2E3DB4CF5uL, 0x384B212A79B3CF62uL), (+1, -7, 0xBED092A581B55A64uL, 0x1ECDFACF3AC906DBuL)),
                ((+1, -13, 0xCD3455A09770DEA5uL, 0x00BD8A2EE23F65A4uL), (+1, -10, 0x9790D7237670B8B3uL, 0xB7572373439D8A2BuL)),
                ((+1, -17, 0x88CAEB865ED401E6uL, 0xF26E1DBB98F580CFuL), (+1, -14, 0x97CE8DEEC9CE62B9uL, 0xDD19A0ED1C4F064EuL)),
                ((+1, -23, 0xAE48D94B8A58488FuL, 0xF0E9BD7F0ED51B42uL), (+1, -19, 0xAD635264599C8AFDuL, 0x2E6CF41F8290BB7CuL)),
                ((+1, -30, 0x92B1479DA8E9F0E2uL, 0x7BC535A7482EE711uL), (+1, -25, 0xBCDC664BB7E96DB7uL, 0x579C864C39109A2BuL)),
                ((+1, -42, 0x8B389F4BB0465CE5uL, 0x0550D857F1EFE559uL), (+1, -32, 0x838EE98BE1779F48uL, 0xBD210F7AC12D0E35uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xD57D7C0174BF927BuL, 0x0A23451D6D6CF5ADuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xAA0C371225032B05uL, 0x2DAEE9B1699C46F1uL), (+1, 0, 0xDBBA710355252BBAuL, 0x36622EE7CB2A0A81uL)),
                ((+1, -2, 0xECF76BD618FEFA83uL, 0xF451C2C2F6F3C3E8uL), (+1, 0, 0xA7A105E74219E2ECuL, 0xB7AB09EB6E69D9E0uL)),
                ((+1, -3, 0xBD4FA5D95064D1E4uL, 0x691472F3A7993746uL), (+1, -1, 0x956FD2EDFA783A47uL, 0xA53675DFFB9F2A7FuL)),
                ((+1, -5, 0xBF010EEB2428F288uL, 0x2F3A177AB4F4752FuL), (+1, -3, 0xAC2BD56090F37648uL, 0xF1FB6D3BE4392AFAuL)),
                ((+1, -8, 0xFCCEEB8F24B5D228uL, 0x693DA2542DD07928uL), (+1, -5, 0x85E2CE8247733522uL, 0x5E74E94ED7A97FD0uL)),
                ((+1, -11, 0xDCFA4AB376876A11uL, 0xBD2F90B9D9126A46uL), (+1, -8, 0x8E95E717DC60B486uL, 0xF05EB51ABCF60A87uL)),
                ((+1, -15, 0xFAE4B61F181B8864uL, 0xB40B21BAF7DA47C8uL), (+1, -12, 0xCEC08CCDA7C6F83EuL, 0x8521E359D608B1FDuL)),
                ((+1, -19, 0xB186E24C33C53574uL, 0xF88A918074920773uL), (+1, -16, 0xC720EC016D0291DAuL, 0xC93C230CBD2B9D89uL)),
                ((+1, -24, 0x91230870F0A503D9uL, 0xB10D82149E8941CEuL), (+1, -21, 0xF2FF0FA4609EA619uL, 0xAB5B33039C6E0E33uL)),
                ((+1, -31, 0xEE498EF152395F4BuL, 0xFD3A90829E08271DuL), (+1, -26, 0xAD0BD634C32B7084uL, 0x6EE0F690C5727EF9uL)),
                ((+1, -38, 0x8F3AD0C3167B9BE7uL, 0xEE22FFB7C41BFB57uL), (+1, -33, 0xF788C6E5B470BC0FuL, 0x5CF6307542E503D1uL)),
                ((+1, -50, 0x823CD272BA6FCDCAuL, 0xC04DF90E07A3A6C9uL), (+1, -41, 0xFBEEE97223EAD7B0uL, 0xEDAAA1920753237CuL)),
                ((-1, -60, 0xA3866DD932267101uL, 0x05AB0DCB48BBE77AuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8DE12EACFE1CBC9CuL, 0xBD88458214875998uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xC9AFEDF332D2D3F9uL, 0x9EABC88C3B04BF0CuL), (+1, -1, 0xCB46505B7804C9BEuL, 0x7A83F4A9B1ACA1E6uL)),
                ((+1, -5, 0xF7D6C658E36AF37EuL, 0xE588E76BBDC48322uL), (+1, -2, 0x8DEB7CE7E3340284uL, 0x0DDD0DCBF6DE68FFuL)),
                ((+1, -7, 0xAC5CDEC250E5F3D3uL, 0x012D882271FC7771uL), (+1, -5, 0xE4E600FAFC6FCE65uL, 0xF7E9793778D520FDuL)),
                ((+1, -10, 0x95480BED6F207C9EuL, 0x9AF3C648DC8E71C5uL), (+1, -8, 0xEB95C03831C04796uL, 0x20EE87CE5931AC33uL)),
                ((+1, -14, 0xA701377009F20DA3uL, 0xAB79F4BB8B9C8655uL), (+1, -11, 0xA1673AED82E69734uL, 0x8C6F95AF72A8368BuL)),
                ((+1, -19, 0xF278FFF0B74A25EFuL, 0x6905E76CA53070EAuL), (+1, -15, 0x9519D06A09CD8E6EuL, 0x177976FF9720324CuL)),
                ((+1, -24, 0xDFE6FFFCB7542598uL, 0x1E60F8AA1604C860uL), (+1, -20, 0xB82D4CA880A9399DuL, 0x19FF623AE7B55C2DuL)),
                ((+1, -30, 0xFAE1FE63F76276DAuL, 0x2FF884CE87EE5053uL), (+1, -25, 0x93D65453DC8CCDABuL, 0xBFA395A91E2A090DuL)),
                ((+1, -36, 0x9C64F39117329655uL, 0x46CCC641203BFAB2uL), (+1, -31, 0x922AE9AA788E7510uL, 0x833C476475F7787AuL)),
                ((+1, -44, 0xB87C675A48CBE08AuL, 0x66A610E7A03D27BBuL), (+1, -38, 0xA226EDA9EF26BB45uL, 0x006D9351220224FCuL)),
                ((+1, -53, 0x8EA2C6211E066DBDuL, 0x0017DF1312D03B32uL), (+1, -46, 0xA9F1B6EC50A9AA2AuL, 0xE9D258BA3091EE59uL)),
                ((+1, -68, 0xD1BD5A5F7DF7DDF0uL, 0x61B6F3926C953E63uL), (+1, -56, 0xE2E9DE00AE2EDD8BuL, 0xD19C905E9736C7E9uL)),
            }));

            public static ddouble Value(ddouble x) {
                if (IsNegative(x)) {
                    return -Value(-x);
                }

                ddouble y;
                if (x <= 0.5d) {
                    y = ApproxUtil.Pade(x, pade_0_0p5);
                }
                else if (x <= 1d) {
                    y = ApproxUtil.Pade(x - 0.5d, pade_0p5_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_1_2);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 2d, pade_2_4);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 4d, pade_4_8);
                }
                else {
                    y = ApproxUtil.Pade(x - 8d, pade_8_16);
                }

                return y;
            }
        }
    }
}
