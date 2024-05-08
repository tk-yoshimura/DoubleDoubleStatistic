using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using DoubleDoubleStatistic.RandomGeneration;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic.ContinuousDistributions {
    [DebuggerDisplay("{ToString(),nq}")]
    public class LandauDistribution : StableDistribution<LandauDistribution>,
        IAdditionOperators<LandauDistribution, LandauDistribution, LandauDistribution>,
        ISubtractionOperators<LandauDistribution, LandauDistribution, LandauDistribution>,
        IAdditionOperators<LandauDistribution, ddouble, LandauDistribution>,
        ISubtractionOperators<LandauDistribution, ddouble, LandauDistribution>,
        IMultiplyOperators<LandauDistribution, ddouble, LandauDistribution>,
        IDivisionOperators<LandauDistribution, ddouble, LandauDistribution> {

        public override ddouble Mu { get; }

        public override ddouble C { get; }

        private readonly ddouble c_inv;

        private static readonly ddouble mode_base = "-0.14182805081482592285930203871252083765";
        private static readonly ddouble median_base = "0.86311662299158754170254190531890461851";
        private static readonly ddouble entropy_base = "2.3726364400044818244844049010588577710";

        public LandauDistribution() : this(mu: 0d, c: 1d) { }

        public LandauDistribution(ddouble c) : this(mu: 0d, c: c) { }

        public LandauDistribution(ddouble mu, ddouble c) {
            ValidateLocation(mu);
            ValidateScale(c);

            Mu = mu - Log(2d * RcpPI * c);
            C = c;

            c_inv = 1d / c;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsInfinity(u)) {
                return 0d;
            }

            ddouble pdf = PDFPade.Value(u) * c_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            ddouble cdf = CDFPade.Value(u, interval != Interval.Lower);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            ddouble x = Mu + C * QuantilePade.Value(p, interval != Interval.Lower);

            return x;
        }

        public override double Sample(Random random) {
            double z = random.NextUniformOpenInterval01(), u = z - 0.5d;
            double w = random.NextUniformOpenInterval01();

            double r = 2d / double.Pi * (z * double.TanPi(u) * double.Pi - double.Log(double.Log(w) * double.CosPi(u) / (-2d * z * (double)C)));
            double v = r * (double)C + (double)Mu;

            return v;
        }

        public override ddouble Median => Mu + median_base * C;

        public override ddouble Mode => Mu + mode_base * C;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => entropy_base + Log(C);

        public override ddouble Alpha => 1d;

        public override ddouble Beta => 1d;

        public static LandauDistribution operator +(LandauDistribution dist1, LandauDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, dist1.C + dist2.C);
        }

        public static LandauDistribution operator -(LandauDistribution dist1, LandauDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, dist1.C + dist2.C);
        }

        public static LandauDistribution operator +(LandauDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.C);
        }

        public static LandauDistribution operator -(LandauDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.C);
        }

        public static LandauDistribution operator *(LandauDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public static LandauDistribution operator /(LandauDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.C / k);
        }

        public override string ToString() {
            return $"{typeof(LandauDistribution).Name}[mu={Mu},c={C}]";
        }

        public override string Formula => "p(x; mu, c) := stable_distribution(x; alpha = 1, beta = 1, mu, c)";

        private static class PDFPade {
            private static readonly ddouble pi_half = Ldexp(PI, -1);

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8FD7D22D9A0E40DEuL, 0x1F5E03328F7B0D03uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xE2039E22F7805780uL, 0x7A25A10DBFCA1548uL), (+1, 1, 0xD1D11CF3DD01D9C6uL, 0x9ABB1F2D93855B6AuL)),
                ((+1, 0, 0x9D5ECA5E835A4BB5uL, 0xF9595E99F30547B6uL), (+1, 2, 0xA78C49E2477E53D5uL, 0x2A453A28FE23B7CDuL)),
                ((+1, 0, 0x807ACA05940D54F4uL, 0xB8076A575A7C2BE1uL), (+1, 2, 0xAAF87690E1EE8A35uL, 0x861B2CD7F2A6439CuL)),
                ((+1, -1, 0x8949FCC7AFC9D57EuL, 0xB230E3583909A02FuL), (+1, 1, 0xF73C1D2EE6B2B26BuL, 0xE7E3F83C501EFFA2uL)),
                ((+1, -3, 0xCACADFC8013E692EuL, 0xEA198529DE9DCCA7uL), (+1, 1, 0x851E13293DDB5593uL, 0x200A72300D38A7C0uL)),
                ((+1, -5, 0xD3EA3B4AE1C3DDC9uL, 0x123CFA2C5F7E88E1uL), (+1, -1, 0xDAA4135C63D4FFC4uL, 0x00678F3B85F1DA6FuL)),
                ((+1, -7, 0x9BB03832BB6EF7B7uL, 0xFACCDC3A55D9ADD4uL), (+1, -2, 0x89EBF04EEE3A5EF4uL, 0x83DDE6237D6B517BuL)),
                ((+1, -10, 0x9A1C9C9E6AC26DB1uL, 0xEACC482E865D1390uL), (+1, -4, 0x84AF97178C0EAFCCuL, 0xAA799EF01C369EEAuL)),
                ((+1, -14, 0xB959B664279F82B4uL, 0xD06B05B2BD49B5B6uL), (+1, -7, 0xBE026AD2579A16ECuL, 0xE2AF2557336DF77EuL)),
                ((+1, -19, 0xC71B33A576EE5E8AuL, 0x256C48663246121DuL), (+1, -10, 0xC05FC7E3C276BF2CuL, 0x6854A9D7CB143961uL)),
                ((-1, -30, 0xEB7F7D1147AC2E7FuL, 0xC897571D47512A6DuL), (+1, -14, 0xF759D394E2BF7809uL, 0xA491FCB848E296DEuL)),
                ((+1, -36, 0xE13E732FD0872D04uL, 0x194AFDE298F32B71uL), (+1, -18, 0x9733A98E92B570B7uL, 0xDC74295BEC6D9868uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xC330E2A239EDEC0EuL, 0xF0984A589B77614DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xA781511D0F41909DuL, 0x4CF64FABA242FF19uL), (+1, 1, 0x8F364E15FA3C01D7uL, 0xF1DF82DC4673E363uL)),
                ((+1, -3, 0xFC061BB4C7FEC4F4uL, 0x8955583CE3C245F7uL), (+1, 1, 0x980E714F5E63FBB3uL, 0x990819D5720C5427uL)),
                ((+1, -4, 0xDC044CE30E1701BBuL, 0x8F71352774FAA287uL), (+1, 0, 0xC81886AE1F2404E3uL, 0x431AF715ACED7A60uL)),
                ((+1, -6, 0xF6E392604B14AAFBuL, 0xDF3B947C43011F70uL), (+1, -1, 0xB419B22E994CFAC3uL, 0xBBA376885CEAC3DBuL)),
                ((+1, -8, 0xB7B1A529FFD33864uL, 0xAA3656B6C92A8ED9uL), (+1, -3, 0xE75919642228C2AAuL, 0xF38D4D824CFCECABuL)),
                ((+1, -11, 0xB1A24705E7347082uL, 0xA7669633725EACB0uL), (+1, -5, 0xD6CAA4978CD6060DuL, 0xC274A84809234C9FuL)),
                ((+1, -15, 0xCD67ADD12FAC827BuL, 0x08A89CD3C1EFCAFAuL), (+1, -7, 0x8E8CDAA2F5952FD5uL, 0x48C11091CBD90A5AuL)),
                ((+1, -20, 0xD732AADEF79EF327uL, 0x346966923AAE35CCuL), (+1, -10, 0x81F4F983EFDB8629uL, 0x70AFA9DFF16FD9FEuL)),
                ((-1, -30, 0xA98A97F404E39DE9uL, 0x03E234AA1CA45765uL), (+1, -14, 0x940B0C7C40D313B7uL, 0x5DDB4A876223B7BFuL)),
                ((+1, -35, 0x8BCA1034A70DEC18uL, 0xF765CC1A7798664AuL), (+1, -19, 0xA309F1D1F5F05700uL, 0xC4DA5106EFA850BDuL)),
                ((-1, -41, 0x94AD41040331C786uL, 0x3684F5D46EC725CCuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xE39A2D43E77886E5uL, 0xCEC269E2B632FCEAuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xB405BE3737025177uL, 0xFFACCC92D31EC392uL), (+1, 1, 0x876186A4796FACFEuL, 0x901DF5FCA32E814DuL)),
                ((+1, -4, 0xFC16688087BED437uL, 0xD2E3A3E1F5BFEB4CuL), (+1, 1, 0x8499D0E79799FAC2uL, 0xB49A7A4F4B82F6CBuL)),
                ((+1, -5, 0xCD84894BA4810106uL, 0x6789F72B506D7EDDuL), (+1, 0, 0x9E5AD5D7364F03C8uL, 0x5B853DDAB55DD103uL)),
                ((+1, -7, 0xD716FAEE2F81260FuL, 0xD89466CE180CDF78uL), (+1, -2, 0xFF8DD33DC1D4A382uL, 0x88F697ADCEBC6CC7uL)),
                ((+1, -9, 0x9572B7DB6DEB12F8uL, 0x6D2E2AC8AB2E5608uL), (+1, -3, 0x91E3F84BF71EF923uL, 0x64A59FC1304B42EDuL)),
                ((+1, -12, 0x88F63322CBBB903EuL, 0x6D0905A0A00F356BuL), (+1, -6, 0xEFAECF9C87C3FF89uL, 0x5B90A0A0CACA9AFFuL)),
                ((+1, -16, 0x9DAA2CF12DF29907uL, 0xE1F74F15C68B0064uL), (+1, -8, 0x8CFCFE4DFB737439uL, 0x4588D08127AAFDDAuL)),
                ((+1, -21, 0xC9FCFB8B957B9419uL, 0x4D59C744DA343B4EuL), (+1, -12, 0xE72B60AE89B334F3uL, 0xB123486068EF346EuL)),
                ((+1, -27, 0xD14149F822E56C87uL, 0x59D439C0FB432F76uL), (+1, -16, 0xF8AF68325B148A28uL, 0x6D53EA1D170A2B1AuL)),
                ((-1, -40, 0xF89A690AE1E45F10uL, 0xD496FC9AE2C9B1F1uL), (+1, -20, 0x9ABB52A642004DC5uL, 0xC7637F6967927AE6uL)),
                ((+1, -46, 0xB215479171E9CE2EuL, 0xFE7280E6D29466AFuL), (+1, -26, 0xA1BBA60DED09C237uL, 0x489DCEDF59CA709BuL)),
                ((-1, -53, 0xADE16CF44CF958A3uL, 0xF76B0E5FA9DAD36FuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xB06C90E0B97BFF54uL, 0x6228C103930E83EEuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0xBC7AA0657EFAD7B9uL, 0x5108C49EAF9DCA55uL), (+1, 0, 0xBD5C4962FB728BBFuL, 0xF9618F2C96BA0A8DuL)),
                ((+1, -6, 0xB2B4816CBBDF90AEuL, 0x9F2735D7A2CFFC96uL), (+1, 0, 0x800A486E4134AEEEuL, 0x9E2A203D1BD715A2uL)),
                ((+1, -8, 0xC5218D783147C1A4uL, 0xABAEC027E86A1E0BuL), (+1, -2, 0xD0858E8B62C3C8ACuL, 0x6F9DB99C4D44B610uL)),
                ((+1, -10, 0x8B21D82A69ED18F9uL, 0x40D096AA869F3E06uL), (+1, -4, 0xE2BE2388ED31D0B4uL, 0x5FFA2ED119F010CDuL)),
                ((+1, -13, 0x82055FB091D7435FuL, 0xF165813BC21925C9uL), (+1, -6, 0xAC7F3C5BEE5F321CuL, 0x191E2048E2253D20uL)),
                ((+1, -17, 0xA0B539A6201D80A5uL, 0x9B8615B2177082EAuL), (+1, -9, 0xBAFB20223636AC82uL, 0x57B39F72518D94BDuL)),
                ((+1, -22, 0xFE6EC2DD011067D3uL, 0x345F7E037F63A550uL), (+1, -12, 0x9025890E0EFCFF5CuL, 0xCEACA5D4BA5F443FuL)),
                ((+1, -27, 0xF03F01D75293E54DuL, 0x7A51BD579B5162ADuL), (+1, -16, 0x9AE26B67F1D07A65uL, 0xFB43A44A34B176CBuL)),
                ((+1, -33, 0xEB4E66221E5C5902uL, 0x02A48044B8B863B9uL), (+1, -21, 0xDE28DAEB4A494438uL, 0x2F1431EAF4A06887uL)),
                ((+1, -40, 0xAC44CC1AB8B1BF03uL, 0xFE915D198D1F609AuL), (+1, -26, 0xC4B212A3E2170481uL, 0x7BC390467BDAF2DAuL)),
                ((-1, -56, 0x92A97A51FF762FB1uL, 0x6F9D8932F3C4DD80uL), (+1, -32, 0xB9F3C4E61D539B8CuL, 0x80CE24F250D5FA49uL)),
                ((+1, -65, 0xB05D7C014F0BC69FuL, 0x0969F555D90279E4uL), (+1, -39, 0x869ABD9721504AA5uL, 0x3E0785912E9971C9uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -7, 0xC5EFA14FB084C1EDuL, 0x8730E036DAD13627uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -8, 0xF147183B26CD7A11uL, 0x195F73951948BE2EuL), (+1, -1, 0xDB3767047C173569uL, 0x3E63A737F1292C21uL)),
                ((+1, -9, 0x80938800B9364447uL, 0xA2C75EAF47643844uL), (+1, -2, 0xA8F7183AEA9CFA01uL, 0x7D8DFF23E5E74090uL)),
                ((+1, -12, 0x9CB9DE056D629875uL, 0x02BBF7BF4E7C8537uL), (+1, -4, 0x9A5F5C181D353F74uL, 0x93DF9150172CAC14uL)),
                ((+1, -16, 0xEFE31BF064E77E13uL, 0xBE27385C62E38A01uL), (+1, -7, 0xB946010298058346uL, 0x3404DD68238F3692uL)),
                ((+1, -20, 0xEE3F51C19F16E0E7uL, 0x2590F5CF6DE8A631uL), (+1, -10, 0x98E9757A7D80E200uL, 0x42E8754BF03140B4uL)),
                ((+1, -24, 0x994D53BB5A55803BuL, 0xB1F6AC629B19733AuL), (+1, -14, 0xB097E162C1137EE3uL, 0x338AF8568924747AuL)),
                ((+1, -30, 0xF7D2048866C466B2uL, 0x4208490BE761C38EuL), (+1, -18, 0x8E53D866BB369A6BuL, 0x57E210BF07B03688uL)),
                ((+1, -36, 0xEB0BC35E62A10216uL, 0x0C60AB9B7BA98F7DuL), (+1, -23, 0x9CD38AA9E118D9A7uL, 0xC301F9143AC740AFuL)),
                ((+1, -43, 0xE4A845B5EAB1D71EuL, 0x990EF155D0C41C89uL), (+1, -29, 0xE266E8D468F62868uL, 0x09C981A4F6ECFBD5uL)),
                ((+1, -51, 0xA5788F25B1AC0CC5uL, 0xE21EA01E9B327D32uL), (+1, -35, 0xC68017747E2C6297uL, 0xD3455261DA71FF83uL)),
                ((-1, -69, 0x898F2650A5793A73uL, 0x8FD304086184C474uL), (+1, -42, 0xB7BC9FC39DEFA17FuL, 0xF9232E0EF85C7049uL)),
                ((+1, -79, 0xA40A759B6D2A0AF8uL, 0x2C40575144E173B5uL), (+1, -50, 0x81A14CF5AB198706uL, 0xB31D4EFC4A3FCE71uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -9, 0xC103A81288484B91uL, 0x73C6853F5E2A054DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -11, 0xFE63B5AB323A55F4uL, 0x4C393C27DC4540D0uL), (+1, -2, 0xEAFDCEA8CE7D17ECuL, 0xFDEF5C52F56199E6uL)),
                ((+1, -13, 0x915560E79216B072uL, 0x59C3687E726A5F35uL), (+1, -4, 0xC0CCF7C420FCC7DAuL, 0x48C46C182C0198D9uL)),
                ((+1, -17, 0xBC653703A35ED989uL, 0x97A9E763696F39CBuL), (+1, -7, 0xBA3BE2F974460A20uL, 0x977C824C5BDEC544uL)),
                ((+1, -21, 0x98437518AE261C2CuL, 0x45653C8007CE16F5uL), (+1, -11, 0xEAC8BFF3F02C1EE0uL, 0x837D3A2226B8FCB3uL)),
                ((+1, -26, 0x9EE4B57D04244721uL, 0x9F24BE0E3BCFA421uL), (+1, -15, 0xCA58B889317DDFE5uL, 0x3BBEE46A300606FCuL)),
                ((+1, -32, 0xD670E904C72125B0uL, 0xCFEBC55B21664F1BuL), (+1, -20, 0xF2D0563D3F980AC7uL, 0x3283D5106CC347C1uL)),
                ((+1, -38, 0xB64C22634E3E3DF3uL, 0x0AE391992223BD48uL), (+1, -25, 0xCAA9A4FADA4976A4uL, 0xEE8AB3D5C5EA63B2uL)),
                ((+1, -45, 0xB7AD416042126932uL, 0x1BBF4EA42AC1E835uL), (+1, -31, 0xE715BFDF00DD3CCDuL, 0x98AA00A76E395328uL)),
                ((+1, -53, 0xC1B7E1534F046BEDuL, 0x039E7DF048C7FF7AuL), (+1, -37, 0xAD36A6161F5D94EFuL, 0x2542F7035A83501FuL)),
                ((+1, -62, 0x9D8FA57B3DC290D9uL, 0xD8AB30844E7DFE0AuL), (+1, -44, 0x9F4C165DBBF0B4D8uL, 0xB95A7A0A0F1746FFuL)),
                ((-1, -82, 0xD4180CD96D07F242uL, 0xDBE7E6FCFD55D205uL), (+1, -52, 0x9DD13A3FBB53BE8CuL, 0x279EB9BAA3509818uL)),
                ((+1, -92, 0xE3F51542C2140CA7uL, 0x0A68DFCC81DB4C6EuL), (+1, -62, 0xF71E539E894930B8uL, 0x2EE58044330A5914uL)),
                ((-1, -102, 0xACAB9E813F02BA5AuL, 0x3F36A89BB9F63355uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -11, 0xB72219FFE913ED20uL, 0x4353D7383EC6EEC3uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -14, 0xD65DD33B765D774BuL, 0xD53DBEEB718D7A4BuL), (+1, -3, 0xD8272D629AA37FA8uL, 0x73ECFDBE1FCE0E5FuL)),
                ((+1, -18, 0xD543ED0EB3667119uL, 0xE75618AD3AA733FDuL), (+1, -6, 0xA0F671D912F4A353uL, 0x1DC87A194CA53664uL)),
                ((+1, -23, 0xEB0B4C4DF305826DuL, 0x88032AD19007EC72uL), (+1, -10, 0x8AEFED7E14CB2BD9uL, 0x63C94043381434C3uL)),
                ((+1, -28, 0x9CD07501A8F242A9uL, 0x5247FE16DA4AF066uL), (+1, -15, 0x99A5A4456AC0E90CuL, 0xF3B24725B833E01EuL)),
                ((+1, -34, 0x820516909B20DD62uL, 0x6D71EC2CBBA98A69uL), (+1, -21, 0xE31F5A91BEC03C63uL, 0xEBA96557102DDDA3uL)),
                ((+1, -41, 0x8445EDA329D9AF79uL, 0x169D4E8BA4ECC6E3uL), (+1, -27, 0xE32C1408FA9A77BCuL, 0x3FF1725B5ECE5DE2uL)),
                ((+1, -49, 0x9CA1B76A20643896uL, 0x87BAD4A3CCCD20D5uL), (+1, -33, 0x983EBA25C82E5985uL, 0xB845872309AD2F38uL)),
                ((+1, -58, 0xBFB9FB884A4E5E90uL, 0xF7263662B5091680uL), (+1, -40, 0x8454DC5BCDFB9868uL, 0x0336C620CD85888EuL)),
                ((+1, -68, 0xB370AFF5FBEB8938uL, 0x65307970B3320972uL), (+1, -48, 0x8BCCC49D14B4E61DuL, 0x118045CBFCDA6A2AuL)),
                ((-1, -89, 0x9D97869552FC614EuL, 0x67B86A023AEB49FDuL), (+1, -57, 0x9E1DEA2A0A6D2BBAuL, 0xFE4E4AC4402BBA50uL)),
                ((+1, -100, 0xBFD0E30E53455574uL, 0x41433968504275D2uL), (+1, -67, 0x8CCEFEF85F28237BuL, 0x45F0A33F70537E75uL)),
                ((-1, -111, 0xA40030BC82F844B6uL, 0xDEBF0FB686AE13A1uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xAF54925E7E1322A2uL, 0x20A11BFCB8F91D68uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xB8887971BD5715CEuL, 0xCD262FF181410707uL), (+1, 0, 0x8B5C8778B1BD05B2uL, 0xCDE0D3DA4BFEE3DFuL)),
                ((+1, -2, 0xB448B7693093C293uL, 0x1AD3B856FEBF3FF1uL), (+1, -1, 0x8BB105B5F779CD8FuL, 0xE094196B1DE0CE22uL)),
                ((+1, -4, 0xD54CAB4A72752284uL, 0x56DDC3FD27C9EB65uL), (+1, -3, 0xA7BAC0E762134EDDuL, 0x8F6C27E9C235602DuL)),
                ((+1, -6, 0xA944B0371713613EuL, 0x7182549AFC218B42uL), (+1, -5, 0x858E9731B21938BFuL, 0x4D826873D6FC2936uL)),
                ((+1, -9, 0xBDB9054BC31A437BuL, 0x854C06049B722B11uL), (+1, -8, 0x94F88A8FDC3D5322uL, 0xDAC59537EBEBBAD4uL)),
                ((+1, -12, 0x9AC4BD49757B6885uL, 0x8BA51052FA28609FuL), (+1, -12, 0xF22EA5B76CAB68E2uL, 0x7ECBF891477D55EDuL)),
                ((+1, -16, 0xBAD567104F7CC77AuL, 0xF003C56D383C9389uL), (+1, -15, 0x930B32857187C366uL, 0x4DF995058B0A93D1uL)),
                ((+1, -20, 0xA7333FF5EAF52A97uL, 0x27A94D83E4ABE341uL), (+1, -19, 0x842109FDF7459BEDuL, 0x81B40E28F23C1921uL)),
                ((+1, -25, 0xD859562288A59DD7uL, 0xFF0FF91A20D048C0uL), (+1, -24, 0xA6919971AA1EC9FEuL, 0x65F9B81F1632FEE6uL)),
                ((+1, -30, 0xB6E8846F3A55BDD1uL, 0x1789C60BEEA6AEB9uL), (+1, -29, 0x95502EAF7D36DD7EuL, 0xCA1EB5AA59D276B4uL)),
                ((+1, -36, 0x8EAC32202F8BD2AAuL, 0x3F64CCCECB8EC1FDuL), (+1, -36, 0xCC801C3A1B553D52uL, 0xAB5ECC83AF6D66C9uL)),
                ((-1, -46, 0x947DB7676FAA0988uL, 0xB363F21852B86CEDuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA7259682B06C78CEuL, 0xEA86611597374BA0uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x95444046E901FD61uL, 0x7C0469EC65F50182uL), (+1, -1, 0xE83755B5F690F26FuL, 0xAB556943302BD3B5uL)),
                ((+1, -3, 0xFB8348FC5AE51A54uL, 0x8E62FF8347A4D74AuL), (+1, -2, 0xC5365C165607A597uL, 0x74FEC4E7F343A219uL)),
                ((+1, -4, 0x84EB65A7F63E5891uL, 0x33B569E207353159uL), (+1, -4, 0xD0EF50DDA48EBE4EuL, 0x47ADF2E7989800D3uL)),
                ((+1, -7, 0xC5FFEEB7068F3220uL, 0xD50E24EF45D48DACuL), (+1, -6, 0x9B91407CBA41CB61uL, 0x1D7B00E47BA112F3uL)),
                ((+1, -10, 0xDD5238F7001442EDuL, 0xEAB5811124228125uL), (+1, -9, 0xADCC64B196D8B5C9uL, 0xFB829C119E4E1138uL)),
                ((+1, -13, 0xC0DD006853441C72uL, 0x9E32B858A485279AuL), (+1, -12, 0x97754C675F0EE4B0uL, 0x3EBE22056A9A0138uL)),
                ((+1, -16, 0x86456AC96317F786uL, 0xD007A11F30214F72uL), (+1, -16, 0xD2EC151FB68DFC04uL, 0x7294291A6D730ADBuL)),
                ((+1, -20, 0x97B8E294BD22A4C1uL, 0x3C2070D8327F98B7uL), (+1, -20, 0xEE5849C6D3D7F473uL, 0x203DC8A30DB58497uL)),
                ((+1, -24, 0x8C5EF8E50E17D7DDuL, 0xAA299D1DFAA7709DuL), (+1, -24, 0xDC7CB6C3B6C3701CuL, 0x403C800383B7DBB5uL)),
                ((+1, -29, 0xD5387003631DF476uL, 0x0F0BCEFC3A123B9CuL), (+1, -28, 0xA76C200B4D563FD6uL, 0xAA03B3D6E30DBF60uL)),
                ((+1, -33, 0x8459E8FB9C84681CuL, 0x212BF425E23CBB9FuL), (+1, -33, 0xD006CF3929464B87uL, 0xD0D59115BBF1A0CCuL)),
                ((+1, -38, 0x846034B73A35EC7EuL, 0x355DE6A01A8373D5uL), (+1, -38, 0xCFBC1E60A1E7077BuL, 0x7FCC45E6CA6730BFuL)),
                ((+1, -44, 0xCEBF6C00539E0B28uL, 0xE9A3BC00B434A749uL), (+1, -43, 0xA291814C2DDF8F66uL, 0xDC6E32265F8B5190uL)),
                ((+1, -50, 0xE9F3AA94596465FAuL, 0x5644C2E901CB8E02uL), (+1, -49, 0xB7856527592E287DuL, 0x143A164EEEF0C03DuL)),
                ((+1, -56, 0xA5EEA66ECB8B0EABuL, 0xE0B3182BBC723303uL), (+1, -55, 0x827AA5417C1977A2uL, 0xEAAF97CAA74CAB68uL)),
                ((+1, -73, 0x80D6C404215CA8AEuL, 0x759DD087B93518DFuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA3021F470665549AuL, 0x889E3F319ADA964BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8666B479B5AAF870uL, 0xECB7FABA547D20E7uL), (+1, -2, 0xD323C1962DCF2EE9uL, 0x1B0A7C4BF67DCF24uL)),
                ((+1, -5, 0xD79595BEFDDF3B54uL, 0xCF8AC04357677739uL), (+1, -4, 0xA94FC47972146155uL, 0x71D41E36F9BA97B5uL)),
                ((+1, -8, 0xDFE1A3E6F4CD187AuL, 0x6D96A41A933428D2uL), (+1, -7, 0xAFD68CFB56161388uL, 0x113A98CD00C63622uL)),
                ((+1, -11, 0xA8E71198D9EB7A24uL, 0x3E6200048728EE7BuL), (+1, -10, 0x84A81B09C61CFBF2uL, 0x4EFA41CCE2070C4EuL)),
                ((+1, -15, 0xC4C5990D717E1F5CuL, 0x2915D3D76A225483uL), (+1, -14, 0x9A8AC8C57EDCCC99uL, 0x68FBED5660A39FA5uL)),
                ((+1, -19, 0xB7889212D6E3919AuL, 0xCBDB332BFDD52179uL), (+1, -18, 0x902625877810B568uL, 0xFB1AAF0F68B5D785uL)),
                ((+1, -23, 0x8C27D445356AFBE4uL, 0x207444C7011A31DBuL), (+1, -23, 0xDC272CDA6B0568D5uL, 0x1D44122EAF6A9F77uL)),
                ((+1, -28, 0xB19384D24CA27005uL, 0xA47DCCC641C8B2B7uL), (+1, -27, 0x8B782C15CC0E02F9uL, 0x0CC7052082D456C8uL)),
                ((+1, -33, 0xBBBF8E96BF4C0B7FuL, 0xC4B78294AC6A2A9BuL), (+1, -32, 0x9374C711FABB0148uL, 0x6017F336217FD56FuL)),
                ((+1, -38, 0xA57F4EE5C1B078DBuL, 0xD9D8CCAE80F9BEEAuL), (+1, -37, 0x81FB57C75EB179E8uL, 0x16B3C6BBD6E45B8CuL)),
                ((+1, -44, 0xF1119D20BF4FFAEDuL, 0xFE0D513CD75D48E6uL), (+1, -43, 0xBD55A5F0C8192EA3uL, 0x6CD558E8583A8EC6uL)),
                ((+1, -49, 0x8E2E59D2C0C76B19uL, 0x42388EF7AA9AD90DuL), (+1, -49, 0xDF5675CD904AB664uL, 0xE496EF3AFD68F3F2uL)),
                ((+1, -55, 0x8256FDC0269DE47CuL, 0x775742E47DB48EB4uL), (+1, -55, 0xCCBCB8E4684E05E9uL, 0xB5A02CD4A4C61944uL)),
                ((+1, -62, 0xAADB745F8FF1D239uL, 0x7E6A48334A9726C6uL), (+1, -61, 0x8630E496D5F25BC8uL, 0x26A24210BDC871B7uL)),
                ((+1, -70, 0xFD285FA06AFF7BF4uL, 0xF2A9D4D2A97F2ABDuL), (+1, -69, 0xC6D460229632CCC7uL, 0xD11922BF7312C482uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9837FE71D4E88uL, 0x706A4BD196CFA50FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xF4D7CDBBF2F9E5BFuL, 0x6675271C8812A8D0uL), (+1, -2, 0xC04C9A583FE97058uL, 0x35EC202FF34DCA91uL)),
                ((+1, -5, 0xB4ADE306B60838AFuL, 0x7F15B3F668D1E5EFuL), (+1, -4, 0x8DE7B7B67121A57CuL, 0x906F9634ED62B96FuL)),
                ((+1, -8, 0xAE1C9BE3F306D8DAuL, 0xB18F2E75ED61E1FFuL), (+1, -7, 0x88BF3ECBBD04C620uL, 0x0630A14A919A3CEEuL)),
                ((+1, -12, 0xF5B1592B74FE7C66uL, 0x2C864BD7EFF3E9BCuL), (+1, -11, 0xC0F776360BE1E82BuL, 0x5DAACCF0DE6DAAF7uL)),
                ((+1, -15, 0x86E5AF6303FB18E8uL, 0x094D8E8F6ED67079uL), (+1, -15, 0xD3E562C28065EE34uL, 0x5F028AF917AA6465uL)),
                ((+1, -20, 0xEF26E740146A71F3uL, 0xECD6951F64C481CDuL), (+1, -19, 0xBBD45B14E2E809AEuL, 0xF74CFA90BFC23C66uL)),
                ((+1, -24, 0xAF2669E6A080CBBEuL, 0x7A793FCF06EB8E04uL), (+1, -23, 0x899001FD1B841DF8uL, 0x7FEB0DE507F6F02CuL)),
                ((+1, -29, 0xD716F94A83CC57A5uL, 0x1A1751A9CE4F19CCuL), (+1, -28, 0xA8EE5BC503EDDF14uL, 0x1FCBF08AE1EABD5CuL)),
                ((+1, -34, 0xDF46DA5703A3DC5DuL, 0xF55C66021D8BEB53uL), (+1, -33, 0xAF5C75456259FF05uL, 0x828F4917BAF665B7uL)),
                ((+1, -39, 0xC472BFD9A4AE1CF6uL, 0xA51AB52FCBE3CAF8uL), (+1, -38, 0x9A4A432100509D63uL, 0x9E78787B455A6F1BuL)),
                ((+1, -44, 0x920E0E4EB6B9A016uL, 0x5348BD49162F70A4uL), (+1, -44, 0xE56C29966A99F5A2uL, 0x83F5703A3199809FuL)),
                ((+1, -50, 0xB5AA7833FF611519uL, 0x258F586BDF9228BDuL), (+1, -49, 0x8EAE1884077EEFB7uL, 0xE152E009352390BEuL)),
                ((+1, -56, 0xB9218964EC7F0A84uL, 0x170EECB510ECB131uL), (+1, -55, 0x9166CBF31FBFBE67uL, 0x5D1EEA614B3F3C96uL)),
                ((+1, -62, 0x948D365B460C05FAuL, 0xA664C3BCC0CF9401uL), (+1, -62, 0xE95825A4586DB093uL, 0xC7ED3D8AD2300D4DuL)),
                ((+1, -69, 0xACC3CE4D7DA11E1DuL, 0xC68EBCF622936B21uL), (+1, -68, 0x87B070001731DAFCuL, 0x5B63EA9BD3ED0A42uL)),
                ((+1, -77, 0xEAA465EDD2B54DD0uL, 0x18C787226AC5856CuL), (+1, -76, 0xB8499C0E4C7C5527uL, 0xA397ACCD383E35EDuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9836E4E44154DuL, 0x8FA0B5750D899B7AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8EAF20FBAA9C0782uL, 0xF2DB6C3A02950E86uL), (+1, -2, 0xE020AE080E774C8CuL, 0xE37FB9B37C7C9C19uL)),
                ((+1, -5, 0xF5FF1F040A9B51E6uL, 0x082900723C572C9BuL), (+1, -4, 0xC1348B6409863134uL, 0x6E0340143426C120uL)),
                ((+1, -7, 0x8AF910EFB9647CCBuL, 0xC8B131161D6D0447uL), (+1, -7, 0xDA4C5328A58F3081uL, 0x9FADA3C51E30B31CuL)),
                ((+1, -11, 0xE71AAF0E96B8314DuL, 0x0B19F798DF65CB79uL), (+1, -10, 0xB5824361CC4FA9CBuL, 0xBDC2FF611EA1731EuL)),
                ((+1, -14, 0x96807C9292C69721uL, 0x1B4BB1C935663ED2uL), (+1, -14, 0xEC6867BE39B8A4F1uL, 0x12706DBF5636299DuL)),
                ((+1, -18, 0x9F8112F9530C09CDuL, 0x3507CAB17A71B3F0uL), (+1, -18, 0xFA8C715D9AEF198FuL, 0x5B707FEABB215269uL)),
                ((+1, -22, 0x8D08399019779010uL, 0x2B24C50EA8D28779uL), (+1, -22, 0xDD88623FAC228D7CuL, 0xE9DBB0FF2C223069uL)),
                ((+1, -27, 0xD384D323709AF594uL, 0x833EF57B57C59F41uL), (+1, -26, 0xA620633443F1BBDAuL, 0xBDB2B4A61D9E6F7BuL)),
                ((+1, -31, 0x881E362D7F7545E9uL, 0x9D63671F05E0909BuL), (+1, -31, 0xD5D04D1DA9845739uL, 0xC51FFDF7BAA6A5F7uL)),
                ((+1, -36, 0x96632F1519A55311uL, 0x40B02E42294A2CF2uL), (+1, -36, 0xEC3A606E0FB6CB18uL, 0xDC5B2509297E41E3uL)),
                ((+1, -41, 0x92268A8F645B47A6uL, 0xFE3DF18B22F856B0uL), (+1, -41, 0xE5929FBC2B4E3A65uL, 0x51BDE9FBD069B42BuL)),
                ((+1, -47, 0xE3CF4074453A9D10uL, 0xA4CFB3CC17E4A3C4uL), (+1, -46, 0xB2EBD54ED4EFD09CuL, 0x2E4DC547C7EE29D0uL)),
                ((+1, -52, 0xB94EC04876D70FD2uL, 0x8F2CDF82CA9822B1uL), (+1, -51, 0x918A4ED89550306FuL, 0xEC3DF804B98BE051uL)),
                ((+1, -58, 0x874C543B4E29B723uL, 0x3B0B76F0C8AD592CuL), (+1, -58, 0xD4869E505F10AFA9uL, 0xF28C493CC2648114uL)),
                ((+1, -63, 0x87F173332818B37EuL, 0xA7204C7E74B90E58uL), (+1, -63, 0xD589FD653E39EBEDuL, 0x83D2EEB1B207CC23uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_0p5_0 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x856EA0D7D42F1804uL, 0x7DA6626665D2AEFAuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xDB59AF6FE7B1AA91uL, 0x718020CAC20F0A43uL), (+1, 1, 0xB0B130EF6EE6FF88uL, 0x8B326EE70B5EF7CFuL)),
                ((+1, 0, 0x99E3B3CC4EE0EDB1uL, 0xEDCD629D0B1707CEuL), (+1, 1, 0xFF7D4726F60976BDuL, 0xAC8D523E43E46E39uL)),
                ((+1, -1, 0xF0D088387FEE232EuL, 0x7D33209A482932F4uL), (+1, 1, 0xF13B3BA2399D5813uL, 0xD9212D5A7E91DD6AuL)),
                ((+1, -2, 0xE7C83DB051CB3A1DuL, 0xD49E1C11C7930C28uL), (+1, 1, 0xA3CEB6D30ACB7118uL, 0x3906414FA726169BuL)),
                ((+1, -3, 0x900987CDCB978CC6uL, 0x6FE01EB6EA693C06uL), (+1, 0, 0xA5C0829E56AECA8FuL, 0x2F89B16D0FFD5DCDuL)),
                ((+1, -6, 0xEBE7614E31744062uL, 0x057EF06FB615EF26uL), (+1, -2, 0xFE10E099BF6F86B8uL, 0x856CDFC9BF6ABFDDuL)),
                ((+1, -9, 0xEFEB9BE17145BE86uL, 0xBD0F9B1E675C4142uL), (+1, -3, 0x91CD3E9CC24C57A3uL, 0x28F511E5A9B36D50uL)),
                ((+1, -13, 0xD974839F8C5400C9uL, 0x97562EAAC93F0DF5uL), (+1, -6, 0xF3C3B75AE69449B3uL, 0x286E9D97E5815821uL)),
                ((+1, -21, 0xCC3A1D33B929DDBDuL, 0x11948FDF3B781882uL), (+1, -8, 0x878146E5E2FABF60uL, 0x3D4D4570F0C77186uL)),
                ((-1, -24, 0xDA389A993E25C4DCuL, 0x5BC7258DBA1F30EDuL), (+1, -12, 0xA66AFFF9C61AE2E2uL, 0xBE809886D3DF82C4uL)),
                ((+1, -28, 0xAF5E8CC7B70B5F26uL, 0xFFFD9FF6AE36BC58uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0x944B88F0B94725AAuL, 0x34A814E2F301D43EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xBEDF4845081BA12FuL, 0xE3B933FF44086C1DuL), (+1, 1, 0xC8327ED082C11EFDuL, 0xB23A97E2351EB365uL)),
                ((+1, 0, 0xD981DF71B91115A2uL, 0x0038C7CE130E67F4uL), (+1, 2, 0xB14EF061E2628B52uL, 0xAAE5E46857BF6E4AuL)),
                ((+1, 1, 0x904C9FD701AA2180uL, 0xDE82BF00671506EDuL), (+1, 2, 0xD42254275892592DuL, 0xF5879B83A792DC50uL)),
                ((+1, 0, 0xF71091B3D61FEDB7uL, 0x00B3F0EAB0BA13D9uL), (+1, 2, 0xBEB2D36E8C2A8285uL, 0xC80975C3A2B40DEAuL)),
                ((+1, 0, 0x8F5C827519D66235uL, 0xA17093718912C1B1uL), (+1, 2, 0x850DE56CC194C693uL, 0x0E9010079152800BuL)),
                ((+1, -2, 0xE7AA62B94FCA9AF8uL, 0xC84652990117BB68uL), (+1, 1, 0x947C216C9332F785uL, 0x9E5F61B1DC07CADCuL)),
                ((+1, -3, 0x842F33C38AB3DFBCuL, 0xC165A5946DB272E5uL), (+1, 0, 0x84997C0BD049E7E2uL, 0x449D298B8C362728uL)),
                ((+1, -6, 0xD35BED1CE692A8A7uL, 0x31AA875EF36AC976uL), (+1, -2, 0xBE913D3ACB8767F5uL, 0x815B80008E10DA2AuL)),
                ((+1, -9, 0xDB8A29F2E1DC1B1AuL, 0x66B12C45CF91C9CCuL), (+1, -4, 0xD673F1C99BBB5B52uL, 0x2BF51B0317AC7E91uL)),
                ((+1, -13, 0xE11FF3B8A1BD56D9uL, 0xB1414F9559110F4CuL), (+1, -6, 0xB9DD3DCF5DAACDC1uL, 0x1C54D4B7D141B56BuL)),
                ((-1, -23, 0xECE08FB97F21C5D9uL, 0xA2B62C5A09834E3EuL), (+1, -9, 0xDDDA0150459C4EC3uL, 0xD0660154C7C6638CuL)),
                (Zero, (+1, -12, 0xA722385F8D3A6AE7uL, 0xAA1CB5249F168406uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA38FA9BBA5400807uL, 0x59F23E253384D94FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -4, 0xC413E9EDEF9E59EDuL, 0x2312AE4ACD95E300uL), (-1, -4, 0xFDBA571DECF8DFFFuL, 0xAA192E3D3A6D5A77uL)),
                ((+1, -3, 0xC7D2C68217FD39EEuL, 0x3166E6055BF2F1F2uL), (+1, -2, 0x92B529D1C8EA5BB3uL, 0x2D0896C95D0F4629uL)),
                ((-1, -7, 0xC3FD19E7D6E6C1F3uL, 0x9223512C5EB8E1DEuL), (-1, -8, 0x8BE1DF9B7326F0ACuL, 0xB565238555A4DFF1uL)),
                ((+1, -6, 0xB81AC4B3370C1981uL, 0x6D3D2244D1637717uL), (+1, -6, 0xEF6C9CC81FCEB7CCuL, 0x750C9A57F5F3B083uL)),
                ((+1, -14, 0x9E921CBC290B9201uL, 0xE4EB84BE2D8AE888uL), (+1, -9, 0x9F5CC06044CF3957uL, 0x6DA88B1CDE039F33uL)),
                ((+1, -10, 0xAAB12C68A6D675CCuL, 0x3594FDA851FA9E5EuL), (+1, -10, 0xBD482BC1225EF49BuL, 0xE5DE0E783758A519uL)),
                ((+1, -15, 0xEE5AEA332597E238uL, 0x52B45D09A1BF9FF2uL), (+1, -13, 0xEEBC92C2E9B6DD2CuL, 0x25A832C5B5E5E7ADuL)),
                ((+1, -15, 0xA6E44EA66E52044EuL, 0xE53CEAB570917A6BuL), (+1, -15, 0xA9D85A5CF5EED235uL, 0x6D34D54A774B8789uL)),
                ((+1, -19, 0x8FF4182DCC66C9F1uL, 0xC27D50D527FF89DDuL), (+1, -18, 0xD498CB4C6D2AE644uL, 0x7EE4D76D6337200AuL)),
                ((+1, -21, 0x91356BB79D640EC7uL, 0x5456DCD2E09339BEuL), (+1, -21, 0xA11C43816F927E40uL, 0xC2A8CC78D678F333uL)),
                ((+1, -26, 0x98E36DB07D3770F6uL, 0x45C03489EEF7A6F8uL), (+1, -25, 0xBB1B1DBC76B1C94FuL, 0xA1BB764BB736CE8FuL)),
                ((+1, -30, 0xDC3408C445B0115BuL, 0x3DDC00782B7CCCD4uL), (+1, -29, 0x92C65667CD825E4FuL, 0xD200249E901E82BEuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA12BBCB3BE24EDE3uL, 0x68DDC63161A2AD92uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xA00507F3553905EEuL, 0x469A55C8A06DC697uL), (+1, -1, 0xFFF3180EEB5FFFC9uL, 0x38DBD2B1984A889FuL)),
                ((+1, -1, 0x9352B2AA3E975820uL, 0x367800093E020B3BuL), (+1, -1, 0xEA8108857A7C1415uL, 0x4CF98C5ECBFD1EFBuL)),
                ((+1, -2, 0xAF4D0AA9EAF56F95uL, 0xA58D8F2904511A9AuL), (+1, -1, 0x8C23369344CCA1F5uL, 0x2558ECD5CEEC1532uL)),
                ((+1, -3, 0xB4166564859AE3B8uL, 0xD0B4CED034179472uL), (+1, -2, 0x8F713113A181629EuL, 0x626326A3CE27EE00uL)),
                ((+1, -4, 0x918BEC57848EBEB6uL, 0x0360E6C5049A7BA1uL), (+1, -4, 0xE89F3A5B43140FB4uL, 0x840E1688EB2E5319uL)),
                ((+1, -6, 0xC8D6A235F0C03746uL, 0xD8F803A617FFD6B3uL), (+1, -5, 0xA00680808F06AC76uL, 0xE630D963D5C76169uL)),
                ((+1, -8, 0xE4313E18963088A7uL, 0xBCA96E9AC54ED823uL), (+1, -7, 0xB64E1262BCD163D4uL, 0x2832FFA13E255BE6uL)),
                ((+1, -10, 0xDAAC619D02C1ACB0uL, 0xE5B3B915C1EE4404uL), (+1, -9, 0xAE44CB642B5E2B48uL, 0x4F08DD22D85110ADuL)),
                ((+1, -12, 0xAC9B00BCEAE4E3FFuL, 0x1623408FB48A2ED8uL), (+1, -11, 0x89DE7686922C6CC9uL, 0x3FD5A84BBE147669uL)),
                ((+1, -15, 0xE095FC1E8E226728uL, 0x1F274DD8EC6D43D2uL), (+1, -14, 0xB3042A03D9C0D7A3uL, 0xC3059BB91550D0F3uL)),
                ((+1, -18, 0xEAB55EFA6A5E3E05uL, 0x61ABFDFF24D9BA83uL), (+1, -17, 0xBB710A25BDCF6231uL, 0x00EEFD39C8A9D958uL)),
                ((+1, -21, 0xC0BCAE0C1C90E6ECuL, 0x6FCF547D631D5063uL), (+1, -20, 0x99A716246717FC47uL, 0x1CE7466756781ACAuL)),
                ((+1, -25, 0xE66FE8DEC21FA48AuL, 0x9AD7495EDD7ABECAuL), (+1, -24, 0xB800B650F5535921uL, 0xB8D750907D53E3EFuL)),
                ((+1, -29, 0x9F29B54C6FC93243uL, 0xF4CB18717E44A37BuL), (+1, -29, 0xFDCB3767CF7F0AF5uL, 0x4E1E3D30BCBC6132uL)),
                ((-1, -45, 0xA190A257FC821A5AuL, 0x34A58C5500350FACuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA07543591DEF1503uL, 0x52FB9925049D89A1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0x8C99FEC34E5ED0C8uL, 0x7501239A8376AA1CuL), (+1, 0, 0xE05CD2429462F300uL, 0xBD0CF2AC2DEC050EuL)),
                ((+1, -1, 0xF20C6016E8167D7AuL, 0x5EEE462F6897C17BuL), (+1, 0, 0xC1208242A527C81EuL, 0x64AF6303FCF89D96uL)),
                ((+1, -1, 0x84C78D1FD90D83F9uL, 0x5A33A67FAE207714uL), (+1, -1, 0xD3E30331529E5826uL, 0xFC34FE432EA2CF91uL)),
                ((+1, -3, 0xCE18E2A32479A2EBuL, 0x82D4E11D9841F24DuL), (+1, -2, 0xA470FF98EDA37515uL, 0x06F4E0057C668A31uL)),
                ((+1, -5, 0xEF29A2CF3508964BuL, 0x7B6B0E7361902B67uL), (+1, -4, 0xBED2F6D249A22775uL, 0xEA04A49ECE6A3E13uL)),
                ((+1, -7, 0xD6ECC337273512C6uL, 0x6240C3CCE5859FB9uL), (+1, -6, 0xAB7C3548C71DAA9EuL, 0x0A0D19653CA829CAuL)),
                ((+1, -9, 0x99371810160CB01AuL, 0xB2BFD1576A8394F5uL), (+1, -9, 0xF47F2A1EEC22F6ACuL, 0xC536BBF9BF32D796uL)),
                ((+1, -12, 0xB02A3F609D44DB65uL, 0x3DFE3C4CF8A336BEuL), (+1, -11, 0x8C8F060995BBEFD1uL, 0xACD5159E99F208B5uL)),
                ((+1, -15, 0xA4F0CCC8B300903AuL, 0x1101AE7A5F8A539BuL), (+1, -14, 0x839AB6C91A8E8F19uL, 0x99BCADF0BBC58CC4uL)),
                ((+1, -19, 0xFC6C7492FF105AABuL, 0x9B6215C43DF0734DuL), (+1, -18, 0xC9675D08CB4958E4uL, 0x547D3BEFC3A1746BuL)),
                ((+1, -22, 0x9D37BE66B9A6ECEFuL, 0xE1A59318DFDB2409uL), (+1, -22, 0xFAE26C86711B18D3uL, 0x5AF87AFC1D74EA8AuL)),
                ((+1, -26, 0x9CF69EA89186D99DuL, 0x9E94625291EFDAB7uL), (+1, -26, 0xFA79F2DC4C4A275EuL, 0x4B8CC2395229AB44uL)),
                ((+1, -31, 0xF3C6B33265FEF631uL, 0xFC2A269E04B0D868uL), (+1, -30, 0xC28175EAAECD5B5BuL, 0x2A19D5C4426452ECuL)),
                ((+1, -35, 0x87EA2BE870FD0668uL, 0x60F2272E21E1EE71uL), (+1, -35, 0xD8E3642FE1C8A1D2uL, 0x2421405DF424976AuL)),
                ((+1, -41, 0xBBA421AB0D0081C5uL, 0x8E9EBE92247BD8BCuL), (+1, -40, 0x95B74F16F5F9D4F5uL, 0x9D6A64E12924F9A0uL)),
            }));

            private static ddouble PlusValue(ddouble x) {
                Debug.Assert(x >= 0);

                if (x <= 64d) {
                    ddouble y;

                    if (x <= 1d) {
                        Debug.WriteLine("pade minimum segment passed");

                        y = ApproxUtil.Pade(x, pade_plus_0_1);
                    }
                    else if (x <= 2d) {
                        y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                    }
                    else if (x <= 4d) {
                        y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                    }
                    else if (x <= 8d) {
                        y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                    }
                    else if (x <= 16d) {
                        y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                    }
                    else if (x <= 32d) {
                        y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                    }

                    return y;
                }
                else {
                    int exponent = double.ILogB((double)x);

                    ddouble v;
                    if (exponent < 8) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -6)), pade_plus_expp6_8);
                    }
                    else if (exponent < 16) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -8)), pade_plus_expp8_16);
                    }
                    else if (exponent < 32) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -16)), pade_plus_expp16_32);
                    }
                    else if (exponent < 64) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -32)), pade_plus_expp32_64);
                    }
                    else if (exponent < 128) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -64)), pade_plus_expp64_128);
                    }
                    else {
                        v = 2 * RcpPI;
                    }

                    ddouble y = v * Square(1d / x);

                    return y;
                }
            }

            private static ddouble MinusValue(ddouble x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 1d) {
                    ddouble y;
                    if (x <= 0.5d) {
                        y = ApproxUtil.Pade(0.5d - x, pade_minus_0p5_0);
                    }
                    else {
                        y = ApproxUtil.Pade(1d - x, pade_minus_1_0p5);
                    }

                    return y;
                }
                else if (x <= 8d) {
                    ddouble v;
                    if (x <= 2d) {
                        v = ApproxUtil.Pade(x - 1d, pade_minus_1_2);
                    }
                    else if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }

                    ddouble sigma = Exp(x * pi_half - 1);

                    ddouble y = v * Sqrt(sigma) * Exp(-sigma);

                    return y;
                }
                else {
                    return 0d;
                }
            }

            public static ddouble Value(ddouble x) {
                return (x >= 0d) ? PlusValue(x) : (x <= 0d) ? MinusValue(x) : NaN;
            }
        }

        private static class CDFPade {
            private static readonly ddouble pi_half = Ldexp(PI, -1);

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xB6921ED0719FA5CCuL, 0xBB19AD4402C94B1EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xCF699217F9B5AF4EuL, 0xCCBBE91C4E5437D4uL), (+1, 1, 0xAAA0DFB1B164845EuL, 0x97166E026D6D2FBFuL)),
                ((+1, 0, 0xDF2BF973ED2DF013uL, 0xB8C52EE9E06B7621uL), (+1, 1, 0xDDF82885AE738C4FuL, 0x893F55D518C628ECuL)),
                ((+1, 0, 0x9335B70AFAB8FA6DuL, 0x4CB5F76952DDA7D6uL), (+1, 1, 0xB69A45965CB8E7D1uL, 0xDDC39D0ED6AC2D73uL)),
                ((+1, -1, 0x833C25BF4680EA60uL, 0x53D7A64032197A4FuL), (+1, 0, 0xD13503FD6F3ABC71uL, 0x7D82C0E99BEBD406uL)),
                ((+1, -3, 0xA44C1C5E95859881uL, 0xB775D9552B5EBB29uL), (+1, -1, 0xADE13F5275843608uL, 0x7FA546BE8727D5F8uL)),
                ((+1, -5, 0x90DF0A515CF1D92AuL, 0x045BFC85DAA4C49AuL), (+1, -3, 0xD43EDAB141EB3848uL, 0x13037DA580D6C13DuL)),
                ((+1, -8, 0xADC08C50167F35E4uL, 0x3A66E390D0181C88uL), (+1, -5, 0xBC0EAA954729E9D4uL, 0x3DFD2BADEBDE748AuL)),
                ((+1, -11, 0x812AA0D8FC7EE252uL, 0x3DE3CBBEB425E816uL), (+1, -8, 0xE8469C37F3CDC97DuL, 0xE946182C5CEDCDB3uL)),
                ((+1, -16, 0xB3F0B81A19AFA49FuL, 0xEA5EF5B539C54147uL), (+1, -11, 0xB59869AAFEF183C0uL, 0x01A8E20611ABECC4uL)),
                ((-1, -27, 0xD55B0611A8837E1EuL, 0xE360A699E5DD3BDBuL), (+1, -15, 0x8A2AA58125308A26uL, 0x2D82F789C5023BEEuL)),
                ((+1, -32, 0xC9FC423D68B2664EuL, 0x6D60BDDE0149984BuL), Zero),
                ((-1, -38, 0xF01CF7F5A0E261E2uL, 0xB167FA0345085900uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xF22864D7D5E62B95uL, 0x4EC3CA6B75DE36FCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xD6B50A1CF8713F73uL, 0xC07402E2CFCAD97CuL), (+1, 1, 0x8B489FEAD1625695uL, 0xA670ECF0F0A1C0D1uL)),
                ((+1, -1, 0xAEFC3163D4EC0632uL, 0xE00AC4132567A096uL), (+1, 1, 0x8DE790682E6691A5uL, 0x095A753771BCBA94uL)),
                ((+1, -2, 0xAB395A47089EE02DuL, 0x554178C5A030C0D6uL), (+1, 0, 0xB107AB69B5E86F5BuL, 0x361A84BB667CB8F7uL)),
                ((+1, -4, 0xDD02296BEF4993FFuL, 0xF46BFCB6A2A7DBA7uL), (+1, -1, 0x950A28FCDEA9C230uL, 0x470C5C2A02C30D83uL)),
                ((+1, -6, 0xC230B80EDC21F2E8uL, 0x8BD7A3B54C22D099uL), (+1, -3, 0xB031A1A3E1E2B28EuL, 0x7524CFE585F5EB33uL)),
                ((+1, -9, 0xE68EADE1B73CDEBEuL, 0xB17C469BAC12B952uL), (+1, -5, 0x9362AE3F1FE86073uL, 0x182D8D0E82900C2CuL)),
                ((+1, -12, 0xB00A240523850C10uL, 0x72078D6D68E3EC06uL), (+1, -8, 0xAB1832C18FFD2AA4uL, 0x882056133BE02BE9uL)),
                ((+1, -16, 0x98B9C2214325C309uL, 0xDBF15B6B393E6C52uL), (+1, -11, 0x825CCF2EEF286B4EuL, 0x7FD7822BE9F6225CuL)),
                ((+1, -22, 0xD6E9AFC84BFB0581uL, 0xDA2651DD87195250uL), (+1, -16, 0xE5E4F3321056ED99uL, 0xC02C3B89C2EAD5BCuL)),
                ((-1, -35, 0xB9209CEBB1F2AAE1uL, 0x870525D248018278uL), (+1, -21, 0xA7125F35869196DCuL, 0xE78C7EA6C95B8EEEuL)),
                ((+1, -42, 0xDFDD62B3D4254104uL, 0x1B3664D3F3990797uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xA6ABFEE84696F96BuL, 0xCD6DC16AB79D59BBuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x8A4A113FCEE0D6D1uL, 0x485D01BEB6DB3BB6uL), (+1, 1, 0x800D51D0582FD8E1uL, 0x174EB7017C7C6919uL)),
                ((+1, -2, 0xCFE0180D0F23C9CEuL, 0x4E5FFF9D36EF85DFuL), (+1, 0, 0xEB6BA7D48A0B668FuL, 0xD2F22259727B9A3BuL)),
                ((+1, -3, 0xB9D77E07BC1089DDuL, 0xF5326F4ECAA799C2uL), (+1, 0, 0x82CD687559734AC7uL, 0x8FFF80571525FE10uL)),
                ((+1, -5, 0xD9EF9117F21DC5C2uL, 0xE10EA9E054EA0C61uL), (+1, -2, 0xC27C09E8B523B6DBuL, 0xA00C163014DA4FEDuL)),
                ((+1, -7, 0xAE0CBDA71B298DFFuL, 0xA010E52A63B7470BuL), (+1, -4, 0xCA2A92257325011EuL, 0x5EAE295FFA67D3B9uL)),
                ((+1, -10, 0xBDED1C0CCA997580uL, 0x99710CB12CF9395AuL), (+1, -6, 0x94F66D7428DA2070uL, 0x80C2E349BD0AFBB3uL)),
                ((+1, -13, 0x89C69FC095308009uL, 0x4D1950256F029145uL), (+1, -9, 0x9A2E47F5F88DEB57uL, 0xD905F249E842C98DuL)),
                ((+1, -18, 0xF85C9B94ED423278uL, 0x916F4921FBAD0E86uL), (+1, -13, 0xD89EA00DBDE20B3CuL, 0xDC7C27E69E93EAE6uL)),
                ((+1, -23, 0xF1545A5F70C1E6F7uL, 0x6D5A789E529F4A3DuL), (+1, -17, 0xC04839CFDFB92044uL, 0x9F09882AD72CBE16uL)),
                ((+1, -29, 0xB44A1002C8916092uL, 0x4CBC4436885C3159uL), (+1, -22, 0xBAA141C4DBF9D76CuL, 0xF50FDB3145E53B58uL)),
                ((-1, -44, 0xA32CFF511DC50887uL, 0xD99F3DCAD55090DBuL), (+1, -28, 0x8CD9F1057C3A5F66uL, 0x9F418A0BC10752F6uL)),
                ((+1, -52, 0xCA320E4FA053F697uL, 0x739A6E8D08285787uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xBD48ADBA6BE348A5uL, 0xE43A799DAE026D8FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xCE1FA0B259E3D2ACuL, 0x5ECFF843E7E75B47uL), (+1, 0, 0xA936A5986F859CF0uL, 0xB1FF8E996391B529uL)),
                ((+1, -4, 0xC87FFB04FB9F7189uL, 0x334EA86412DBC9F5uL), (+1, -1, 0xCA2EED8230DAF0F8uL, 0x59D6B1F7D2688ABCuL)),
                ((+1, -6, 0xE3F3FEA604B4D545uL, 0x75A067983B1CBED2uL), (+1, -2, 0x8F6D48D58D612A71uL, 0x03536DA206532087uL)),
                ((+1, -8, 0xA66C1B60F7EFDF5BuL, 0x0F92E90C445D7DBDuL), (+1, -4, 0x857B7588A02F4AD1uL, 0xBC88CC2350EED9F7uL)),
                ((+1, -11, 0xA167AE484D494B8DuL, 0x0B7CC2DCE966D7ACuL), (+1, -7, 0xA9E873B8D00E2B75uL, 0x5D387A20B23FD826uL)),
                ((+1, -15, 0xCFBFBB8442D66BB7uL, 0xD4309FADEA49926DuL), (+1, -10, 0x958214008D50695FuL, 0xC4892F3CD75B517BuL)),
                ((+1, -19, 0xABEE5D0E6939213AuL, 0x6C24B3C734FF5BC6uL), (+1, -14, 0xB392AD2CBA87D141uL, 0x0053B4157E5B1269uL)),
                ((+1, -24, 0xAA868BB345883AAEuL, 0x309E0F379603D56FuL), (+1, -18, 0x8DAF2B27DE6CC9DFuL, 0xFE0463C4308E8492uL)),
                ((+1, -30, 0xB05E71C398B44E47uL, 0xBA3C9F16CC81BD21uL), (+1, -23, 0x885192ABBC64DB73uL, 0x6F26ECCAC395F1A3uL)),
                ((+1, -37, 0x891F798181529CE8uL, 0x7A558BE0D1E271AAuL), (+1, -29, 0x8ABDD8D9F47F5DC8uL, 0xFA19161819482A00uL)),
                ((-1, -54, 0x84363F38A5004BF3uL, 0xA9D5C4AB81EF6758uL), (+1, -37, 0xD6CD9B8F7CA9F13BuL, 0x1519F475E9C5070EuL)),
                ((+1, -63, 0xA9A7E9F1EE44A1ADuL, 0x489F7AD20C2E3307uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xBD3CAC5BD395F849uL, 0xC936ED939BFCE479uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x80013825C7D059BFuL, 0xA90B762DA3E886D1uL), (+1, -1, 0xCEA2CF475197ADB8uL, 0x3C94F3004A3253ECuL)),
                ((+1, -6, 0x989BE62C54F492C0uL, 0x5F6BD3AD890E1050uL), (+1, -2, 0x94FF311D7AA06E55uL, 0xEB96BCFADC9E6F9EuL)),
                ((+1, -9, 0xD2478934EB9B2C2AuL, 0x3F5F835028DCB68FuL), (+1, -5, 0xFC8A615C2A6B15D4uL, 0xE8099F36A902957BuL)),
                ((+1, -12, 0xB880D3DEE7F1F2CAuL, 0x4C6D7564405A7EC5uL), (+1, -7, 0x8B32BFFFD4CB6DBAuL, 0x53C7595E8162F9ABuL)),
                ((+1, -16, 0xD631790C47D66DBDuL, 0x546E272FDEB3A472uL), (+1, -11, 0xD09E56645B4F5AA4uL, 0x3285470047D5579FuL)),
                ((+1, -20, 0xA59195F9F83A440FuL, 0xAF99DDAA44AFD715uL), (+1, -15, 0xD7B7D74BD8D7A074uL, 0x72FE2A4442715F89uL)),
                ((+1, -25, 0xA766C1D204B928A7uL, 0x7AC5A84BECECADF6uL), (+1, -19, 0x9905A795A3A59A23uL, 0xA1EC613F58590BB9uL)),
                ((+1, -31, 0xD3E03E341906196BuL, 0x2ECE1570E6B33CE6uL), (+1, -24, 0x913FC27F62C21BE2uL, 0x75CCB7C85D555999uL)),
                ((+1, -37, 0x9AD21D1B6B511DFAuL, 0x3CA56C5260395F08uL), (+1, -30, 0xAFC2D51854AADDFBuL, 0x4DC34906B6693C8AuL)),
                ((+1, -45, 0xE0BE57BC353FDA0FuL, 0x57B67F025292125CuL), (+1, -37, 0xF945F87565828C64uL, 0xEF6B00714D09A64AuL)),
                ((+1, -54, 0xE514DA975436CE61uL, 0x2C73160290E92F78uL), (+1, -44, 0xB1C052BA9174D697uL, 0x6841975DC0AE51A9uL)),
                ((-1, -75, 0xE056A9F319514F64uL, 0x349C37AEDC8D47C5uL), (+1, -53, 0xB3C6E2A9086C530AuL, 0xCD0756C72363430AuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xB5A524D5C35E8565uL, 0x1A077112BC03FABEuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -7, 0xF1DB00CA813483FFuL, 0xBB177A135BEF9538uL), (+1, -2, 0xCC6E6F7780C365D0uL, 0xF66642E1B5075B47uL)),
                ((+1, -9, 0x8B37A7A9ED605423uL, 0x38F10C6B513141BDuL), (+1, -4, 0x8F9AB330FDE8B552uL, 0x6780726F3A117CACuL)),
                ((+1, -13, 0xB5351EF7D3E56987uL, 0x4A420E9ABD389E55uL), (+1, -8, 0xE8FCF5F25749739FuL, 0xFDE1B30C0708D9D9uL)),
                ((+1, -17, 0x92635FBB07A21D16uL, 0x6F050B16A12E38EEuL), (+1, -12, 0xF0E782D60274E8BBuL, 0xF06BD1915C1A1C05uL)),
                ((+1, -22, 0x97C29D2A13BF471CuL, 0xBC2B8DC96D56A1CFuL), (+1, -16, 0xA5501EB6734323A0uL, 0x8236ACBD7A5E6CF3uL)),
                ((+1, -28, 0xC9C3E094DBB96E6FuL, 0x40BE99FA68DAD32CuL), (+1, -21, 0x98129FDB2D217854uL, 0x06A68EEB3FD45EEFuL)),
                ((+1, -34, 0xA6FE1503002A2ABDuL, 0x5BFCB39E8559F1A4uL), (+1, -27, 0xB91B284AEA579FB8uL, 0xFCB18B52EDE8965FuL)),
                ((+1, -41, 0xA116F1166FA3A679uL, 0xDF0239EE454D17F7uL), (+1, -33, 0x8FA8B9D4FBD1E42FuL, 0x2EF3292717ED8943uL)),
                ((+1, -49, 0x9EB67FA98A87675BuL, 0xF9C5E28E0320713CuL), (+1, -40, 0x84721ECF978BF4BEuL, 0x0810F941F64AF292uL)),
                ((+1, -59, 0xE8515C72DF62619BuL, 0xE41F63AC2C4DBE07uL), (+1, -49, 0xFD58C2B249110242uL, 0xBC736C08A75EE5FEuL)),
                ((-1, -80, 0xC55B64A6DFB6B96EuL, 0xC0414ADF87B0EECAuL), (+1, -58, 0xB658345EB2623D7BuL, 0x9FE9985750FCA545uL)),
                ((+1, -91, 0xEE44F70EBFBC1E7AuL, 0xACF58DB8CDC6DA8FuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0xAEB045E7371BF9AFuL, 0xED82F240A0E24C5DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -9, 0xE3C0FAED2FA55367uL, 0x0C09E6B39DDF2E61uL), (+1, -3, 0xC86DF7BBBBF1B047uL, 0x503F313F3866BD51uL)),
                ((+1, -13, 0xFFBE81141F34F647uL, 0xD99452030306C7A7uL), (+1, -6, 0x8988D9C6AF8EFAD4uL, 0x469A17F85D789035uL)),
                ((+1, -17, 0xA1C0F7A1AEC71236uL, 0x49CB3FA17F710430uL), (+1, -11, 0xD93418DB3CDB7AFFuL, 0x5765CB7333F5E013uL)),
                ((+1, -23, 0xFD2C92EACD347B1CuL, 0x2697D51CA22E32AAuL), (+1, -16, 0xD9DFD4A82E29D58FuL, 0x6F63A69B783F5A8EuL)),
                ((+1, -29, 0xFD9702477DA9A418uL, 0xBA437B4478743046uL), (+1, -21, 0x90998E3373E170F6uL, 0xDDEAF59EE080CEE7uL)),
                ((+1, -35, 0xA28E1F0DD52FB43EuL, 0xF5805AE597C8EAE1uL), (+1, -27, 0x80502C0F6CE525D9uL, 0x83E0BB00F562E9B3uL)),
                ((+1, -42, 0x81943D7F7B6BFC4EuL, 0x7B44CE8811E84E5CuL), (+1, -34, 0x9659E9E9D17B3F58uL, 0x8B5F3E1FCA77EE38uL)),
                ((+1, -51, 0xF0AB902747428935uL, 0x81D451C3592FB8BBuL), (+1, -42, 0xE0595CF284116A14uL, 0x5C6275DFB09E14CBuL)),
                ((+1, -60, 0xE456D0CB47F8538DuL, 0x59370A561034D951uL), (+1, -50, 0xC6B9BEF374B65939uL, 0xE21DBAC1CEE432E7uL)),
                ((+1, -70, 0xA117D46584DDB1A1uL, 0x4F59EFCED48B9F0EuL), (+1, -59, 0xB6A8623025C418F1uL, 0x5B935AE4C35AD035uL)),
                ((-1, -94, 0xFF99B508C72C77B7uL, 0xB36AAE061C98F099uL), (+1, -70, 0xFCF73BCD36B1BF2EuL, 0xBFF0254E243BF812uL)),
                ((+1, -105, 0x95722855C2D431D3uL, 0x8290E2E3713DE58BuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA9EFEB701EE7F2B2uL, 0x64A7B0D08954FF9DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xA3F34C2D25D879CBuL, 0x58EB22A74A166970uL), (+1, -1, 0xFC9CDA5F87AA8C8AuL, 0x3DDE7390C2F76832uL)),
                ((+1, -2, 0x92531B8790531C1EuL, 0x7E7FFB5E9AD36CE3uL), (+1, -2, 0xE4CE407885907D81uL, 0xB7C4EA25D18CC378uL)),
                ((+1, -4, 0x9EACD45F988F3258uL, 0x2E7CD00E5E9EADA9uL), (+1, -4, 0xF9B6C9400F94FF30uL, 0xBE162BC99347FB5DuL)),
                ((+1, -7, 0xE93B4BE96799FDC6uL, 0x38BC5E14622F7E22uL), (+1, -6, 0xB77ABDD2B0AA3DD7uL, 0xCAC4F6FC2EA7D6CBuL)),
                ((+1, -10, 0xF6E0A809DA4449F1uL, 0x72480A14D1D38E6BuL), (+1, -9, 0xC1B36E3BA879B982uL, 0x0967A0642B537CD9uL)),
                ((+1, -13, 0xC35C8136930045A3uL, 0x62867C082A77E222uL), (+1, -12, 0x99515E53B82687DCuL, 0x2DD5AEE8824A1802uL)),
                ((+1, -17, 0xEBE483ACBE224681uL, 0x31125F50FB2DE0E1uL), (+1, -16, 0xB98A4DE26332175EuL, 0xDA2530FDCC850ED8uL)),
                ((+1, -21, 0xD9B2C312157436DDuL, 0x7AC2CDBB5581F6C7uL), (+1, -20, 0xAADA28192287DC9BuL, 0xE68DC5397075E2A5uL)),
                ((+1, -25, 0x96CCD9C4052BF880uL, 0x7A557F795F1B4508uL), (+1, -25, 0xECBD3BE4EB6264F9uL, 0xA6654B22E4BBC9F0uL)),
                ((+1, -30, 0x92501E059B538D73uL, 0xC9AE389EE0585F0EuL), (+1, -30, 0xE6309D590E457E3FuL, 0xEBBA1D4DE915536CuL)),
                ((+1, -36, 0xAEA42669D34B9280uL, 0x5DD63064A5F84B06uL), (+1, -35, 0x88F1E395BA665140uL, 0x3E4D080E3D899644uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA541C42790A2B0ABuL, 0x10D583663B8F9005uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x86A97A3D0AA65334uL, 0xDD040F913BD6E794uL), (+1, -1, 0xD2A27B978E1FBF75uL, 0x2C3D98F3714A5597uL)),
                ((+1, -3, 0xD21AEFEC7EAFAD62uL, 0xA30FAD4E545C8046uL), (+1, -2, 0xA4F568237E749E52uL, 0xC9DE1A954E486D20uL)),
                ((+1, -5, 0xD11925BE318FB86DuL, 0x61ADEDF574E9A37CuL), (+1, -4, 0xA446DABD14CE0E1AuL, 0x89D5235D9800851DuL)),
                ((+1, -7, 0x95500351BB07BA3FuL, 0xE1D2B21B89CF32F5uL), (+1, -7, 0xEA8E9432570B2ABEuL, 0x38EC45B946C85703uL)),
                ((+1, -10, 0xA2A31E396383FB87uL, 0x53499FCF4C95BF2BuL), (+1, -10, 0xFF776B0E954F4B0FuL, 0x9DAFE2D3C6F2453CuL)),
                ((+1, -13, 0x8BD11FFCEF9894F6uL, 0xB7E59B4F4E5150BDuL), (+1, -13, 0xDB9DC9AF676ED181uL, 0x137D54B28020101CuL)),
                ((+1, -17, 0xC1983133108A9C6FuL, 0x7437CE029D613469uL), (+1, -16, 0x980A8F16487B40ABuL, 0x7958FFB8380195A8uL)),
                ((+1, -21, 0xDA8608A3DC60BF0CuL, 0x415CD7EB75FD45BCuL), (+1, -20, 0xABA6017ECC8A5C1BuL, 0x46ABE080A28DBBDEuL)),
                ((+1, -25, 0xCA5D562CDE306445uL, 0x91A28770B110F298uL), (+1, -24, 0x9EED7B8DCA3C57F5uL, 0xFBF57CA363EF7CC9uL)),
                ((+1, -29, 0x99E5F49B6B65DD46uL, 0xF4DE2738387CD37FuL), (+1, -29, 0xF1B5C0D997AEE040uL, 0xA0490A10DFBB36FFuL)),
                ((+1, -34, 0xBF1B9BE57E51C0B2uL, 0x0DD1CA9C23954F3DuL), (+1, -33, 0x9627425DA5DCC38CuL, 0x675A68597D889DBFuL)),
                ((+1, -39, 0xBED740E15F63BAD2uL, 0xC1D6FB7BDBE76F26uL), (+1, -38, 0x95CBDD688905A778uL, 0xB93DE10D63AFF0F6uL)),
                ((+1, -44, 0x93FFFE7E6B39E251uL, 0x4469255DBF55C6D1uL), (+1, -44, 0xE8A596E381FD2E9BuL, 0x8D0DF93671E9E3E4uL)),
                ((+1, -50, 0xA62EE9C89580FEA4uL, 0xFEA44F96A86A9492uL), (+1, -49, 0x826B7A6AEFAA3AB7uL, 0x754A6854AEE97E49uL)),
                ((+1, -57, 0xE06D91C5C93F5D91uL, 0xF55B095F0479D267uL), (+1, -56, 0xB067B211F7C7D992uL, 0xCDBA6831693261F0uL)),
                ((+1, -75, 0xE7068AF22FDCE1F8uL, 0xF7C4B36E13944976uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2FE052F99657C0EuL, 0x682DA8F9F9648154uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x869AF8BE4B86CB17uL, 0xCFD8EB63001489EFuL), (+1, -2, 0xD3731C6B59D3901EuL, 0xBBD02684E265131AuL)),
                ((+1, -5, 0xD6FE32F5B12DC00EuL, 0xF5E12E9F9ED4865BuL), (+1, -4, 0xA8D9D35AD4634CDAuL, 0x9E07A422676ED056uL)),
                ((+1, -8, 0xDCE175012281DA89uL, 0x97C9BD52F464DA16uL), (+1, -7, 0xAD7AF75B4363B6DAuL, 0x74574F21FF5DA537uL)),
                ((+1, -11, 0xA3AB94E31BD956D6uL, 0xADE8063EA0009B32uL), (+1, -10, 0x808BF989D805C8ADuL, 0xD0E11057FAB842B9uL)),
                ((+1, -15, 0xB9C5EC24CE2FF505uL, 0xBBEDCD8115D45CE4uL), (+1, -14, 0x91E7A17BBAD73FBEuL, 0x17E111F4C91E4EA1uL)),
                ((+1, -19, 0xA74472ABBDD53B01uL, 0x7013505034AB3E85uL), (+1, -18, 0x835F5901D06234C1uL, 0x69BE9CCDD1BA994BuL)),
                ((+1, -24, 0xF3F00163C5301B40uL, 0xF5409C1E911BD486uL), (+1, -23, 0xBF9649CF5D66A35CuL, 0xB376FEE019963BFFuL)),
                ((+1, -28, 0x919D3B99FC9A468AuL, 0xADE7474A5D1D33C3uL), (+1, -28, 0xE4BB1573BECD1F37uL, 0x41EBA0BD08623502uL)),
                ((+1, -33, 0x8EA4BB1A1AD41FA5uL, 0xF8AC710AD0051E1FuL), (+1, -33, 0xE010600C7EE96B4EuL, 0xEE51FCC8BB55FCD0uL)),
                ((+1, -39, 0xE3CD6148789B3AA6uL, 0xFB4EFF77ADE5A01BuL), (+1, -38, 0xB2EA3892C29C6D20uL, 0x69A7BA531C5A007FuL)),
                ((+1, -44, 0x918469FC4FC19EF9uL, 0x43006A5CB0EF446DuL), (+1, -44, 0xE49442CEDEA7CB50uL, 0xD402FFF92D369079uL)),
                ((+1, -50, 0x8F16D64FDAEBD41FuL, 0x0BEDEA8503F5950AuL), (+1, -50, 0xE0C32EBD004753D8uL, 0x8CB7CAB185BA577AuL)),
                ((+1, -57, 0xC6F17F100E337285uL, 0xD31323312E8012EAuL), (+1, -56, 0x9C404A5E04F94985uL, 0xBA85B40575F0FCA5uL)),
                ((+1, -64, 0x9E08C00FE1618027uL, 0xE9688C9943B5C5D1uL), (+1, -64, 0xF83C761992FF57F1uL, 0x129B5322161C804AuL)),
                ((-1, -88, 0xCB70426D4FD4CB7AuL, 0x23C2715379D9BCF3uL), Zero),
                ((+1, -96, 0x8B3AF45407CCA6B2uL, 0xF73EBD357C17C0A7uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F983774E911066uL, 0xA857A2F016ABE3A0uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xF5488508174795AAuL, 0x35E3004DA31E60E5uL), (+1, -2, 0xC0A521365302BFB3uL, 0x3DE5E5EBF9DBC267uL)),
                ((+1, -5, 0xB547F613F9A68A85uL, 0x6A7B441A54A9CCC7uL), (+1, -4, 0x8E60BA3AC611911FuL, 0xA6AE9485C890CC8CuL)),
                ((+1, -8, 0xAEEC6A75355BD6DDuL, 0x8249E7DE879E7E1BuL), (+1, -7, 0x896274D77088DB97uL, 0x0DD620CE7C66CB84uL)),
                ((+1, -12, 0xF7204D1DB870A7B6uL, 0x6CC28F612A280172uL), (+1, -11, 0xC217AA7BD696C888uL, 0x23A3734B6D100C8DuL)),
                ((+1, -15, 0x87D32E5133BDFEE3uL, 0x14139DE6E4005C60uL), (+1, -15, 0xD55A717CE8D0C5FDuL, 0x117C6E6065155EC9uL)),
                ((+1, -20, 0xF1054CAF12E2D80DuL, 0x7562255186AEA738uL), (+1, -19, 0xBD4C165ACBFF006FuL, 0xA429650D9D15721FuL)),
                ((+1, -24, 0xB0AACCDB81BCF528uL, 0xC02417A66BED0B64uL), (+1, -23, 0x8AC10BBC4ACE3BB3uL, 0x1CE06EAC42E617FFuL)),
                ((+1, -29, 0xD91E79EB7256FF37uL, 0x83C1AB7909851644uL), (+1, -28, 0xAA865FF0035AE340uL, 0xDA0EB70BC4CBD07DuL)),
                ((+1, -34, 0xE18AA4902FC1D6C0uL, 0x04327E3616797393uL), (+1, -33, 0xB123D2F406A87817uL, 0x5F4B31A107AF321CuL)),
                ((+1, -39, 0xC691C3A1F992CD2AuL, 0x530B8D823E189150uL), (+1, -38, 0x9BF4BEBA64B4703BuL, 0x1C5CADF95C50FD57uL)),
                ((+1, -44, 0x93B879E4A4042B28uL, 0xD7DEC992C7FF2A7CuL), (+1, -44, 0xE809FB58B7B0BB36uL, 0x6328FBF0B8C1134BuL)),
                ((+1, -50, 0xB7D73D45D8A17072uL, 0x199AC578D7387141uL), (+1, -49, 0x906361BA1EDD8188uL, 0x8F321C6BAFE66DF2uL)),
                ((+1, -56, 0xBB723B46861F5C35uL, 0xBB491E2F91FFA52DuL), (+1, -55, 0x93384C51D9227698uL, 0xF9B4DA637D2FD839uL)),
                ((+1, -62, 0x967C1B729160604FuL, 0x5D52B4B5E4AA4D08uL), (+1, -62, 0xEC6186AF18967299uL, 0x22E94958DAC6B326uL)),
                ((+1, -69, 0xAF190B4946261C00uL, 0x4BD9EB14596F1791uL), (+1, -68, 0x898581DE6A221A65uL, 0x2EB8C86B78DF5E09uL)),
                ((+1, -77, 0xEDED8E1D7840EB07uL, 0x6A20021E25166C52uL), (+1, -76, 0xBADE40B3715E1924uL, 0xCF50C2491CEAA622uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expp64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9836E4E44153BuL, 0xF9C46567C6ED7B1EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8EBC9265AD941DCFuL, 0x700AB625BA462FCFuL), (+1, -2, 0xE035CBCAC7BCD20EuL, 0x88D82F65CC0EB38AuL)),
                ((+1, -5, 0xF62C2F3E4CE5D6A1uL, 0x1C6A70BC0DA39770uL), (+1, -4, 0xC157EFEC375858A8uL, 0xC7224E0A235DD131uL)),
                ((+1, -7, 0x8B1E3B4FDBA672B1uL, 0x35F2C12DEB459709uL), (+1, -7, 0xDA86B44E0DDDD7B6uL, 0x9E295350116036E9uL)),
                ((+1, -11, 0xE76AF84E84C875E2uL, 0xB1EE0464B4387424uL), (+1, -10, 0xB5C151DDDA32A815uL, 0x0368775ABE5358CCuL)),
                ((+1, -14, 0x96C03A0977CD84CFuL, 0xFF2C13C3F6FCE8F8uL), (+1, -14, 0xECCC872800B74889uL, 0xBB0010B83B443B9EuL)),
                ((+1, -18, 0x9FD0365110A08B3EuL, 0xE2CDA082A17C14B0uL), (+1, -18, 0xFB08C0AAB3595609uL, 0xA96B54E533147359uL)),
                ((+1, -22, 0x8D58008AD368AB19uL, 0x0F246E186F89D70FuL), (+1, -22, 0xDE05B296F970A0CBuL, 0x3C60B39901EB483AuL)),
                ((+1, -27, 0xD40A932A49A66220uL, 0x3ED3AF5E67E6D228uL), (+1, -26, 0xA6896F421EFD3965uL, 0x486EE7463D8E9F8CuL)),
                ((+1, -31, 0x887D119B58655F81uL, 0xAA1B864D1929D29AuL), (+1, -31, 0xD6654D6FE4EA5A41uL, 0x8777C2AB5A98AA4EuL)),
                ((+1, -36, 0x96D528AF77C03679uL, 0xC6B915F1CBB58C7EuL), (+1, -36, 0xECED6880650C9AA0uL, 0x12A840D55BAB45E7uL)),
                ((+1, -41, 0x929EBFD6D4AD7BA1uL, 0x29D3A93DC855EDE5uL), (+1, -41, 0xE64F7249EB0EF37AuL, 0x810C664E017198F2uL)),
                ((+1, -47, 0xE49529F8E1AE8AA6uL, 0xBE96E13A316DDE99uL), (+1, -46, 0xB38745E9A8BE882EuL, 0xC7E103504AF4B71FuL)),
                ((+1, -52, 0xBA002D4366BF28D4uL, 0x4EA6EB306E66F4A8uL), (+1, -51, 0x9215A866858D496AuL, 0x70DF54B266172A1EuL)),
                ((+1, -58, 0x87CB321A35872F56uL, 0x205512A8BCC1DA53uL), (+1, -58, 0xD54DE66F1912208EuL, 0x30E950238F6F9A89uL)),
                ((+1, -63, 0x88863DB2E8E6F0B2uL, 0x772800E2C26C0671uL), (+1, -63, 0xD673B5CFBB4F143EuL, 0x89958A476D6CEC70uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_0 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xB021B0B7012C7211uL, 0x60F3DE9F6A4A6B44uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x96E5D9D12DAB1849uL, 0xDB90AE7A5E43C34FuL), (+1, 1, 0xDF1B5A69B848EB50uL, 0x0B53FE452E7A8045uL)),
                ((+1, -1, 0xEE3976D5F81A337CuL, 0x703A2243B127984CuL), (+1, 2, 0xCFF05FBC77820002uL, 0xF966E8F6CC9A6B44uL)),
                ((+1, 0, 0xE5CE2E97B45A6E15uL, 0x8F8FA5ABAE64F695uL), (+1, 3, 0x826589FA01563D4DuL, 0x0FDDC93B584F3362uL)),
                ((+1, 1, 0x9799D284FF4FACEAuL, 0x3D658DC01E8FD807uL), (+1, 2, 0xF2FAFE6E2F2971ABuL, 0x28EAC843E219B032uL)),
                ((+1, 1, 0x912CC0246331012CuL, 0xCA700C9E992F53E7uL), (+1, 2, 0xB009B14A9D1E8F87uL, 0x3CE4991737FC5993uL)),
                ((+1, 0, 0xD130E1CD1BF04C21uL, 0xC59C481C1B19C2FCuL), (+1, 1, 0xCBACF7351FF23F16uL, 0x29258ABB837625C0uL)),
                ((+1, -1, 0xE80956E477CE4A67uL, 0x12453748F4A5EE55uL), (+1, 0, 0xBE631E2A266C1F91uL, 0x6AFCA99753AC8C60uL)),
                ((+1, -2, 0xC8C22A98A2549A47uL, 0x6733600A2AE0D1F9uL), (+1, -1, 0x904A6D26CF01BC5FuL, 0xC5F6A17698EA27B8uL)),
                ((+1, -3, 0x87EB25B346E49B6CuL, 0x14FC1CE45F28B16CuL), (+1, -3, 0xB028CEC47D6B3430uL, 0x18CFD864AE7455D8uL)),
                ((+1, -5, 0x8EA307789A6F98E7uL, 0xE606DFF4964AC094uL), (+1, -5, 0xAA83A75A9D11D6A1uL, 0x00E6F5929DE3166AuL)),
                ((+1, -8, 0xE1D9D82274B61496uL, 0xB9AC9EE5EEB6E67DuL), (+1, -8, 0xFD4B7DE8DFE67ECBuL, 0x4B23E4347A9564FBuL)),
                ((+1, -11, 0xFF4E94453329C690uL, 0xE001A3D9451E46F3uL), (+1, -10, 0x8825D36D25AA295CuL, 0x40F9CCE1D871DDE9uL)),
                ((+1, -14, 0xB6568EBD779C8489uL, 0x2781CBC08DFC9073uL), (+1, -14, 0xBAFD559E225F56D5uL, 0xD606700C577CF43BuL)),
                ((+1, -19, 0xE6A050A648FB3DDDuL, 0xF1763A3DD0F28A17uL), (+1, -19, 0xE69B667C796D45C0uL, 0xAAE29B2C798E8002uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xABE43D404DEF0C20uL, 0x5C0CE74E6382627FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -8, 0xDBD6EFDDB775CC69uL, 0xB3A30B52EC1B0B84uL), (-1, -3, 0xDE930CB7F8EAF496uL, 0xB7E65A8D80FF068AuL)),
                ((+1, -4, 0x999294C6FDD0FEFEuL, 0xEF4463B2406CB42BuL), (+1, -2, 0xAA42482A22E6A4B8uL, 0x8E34F4830C3A1070uL)),
                ((+1, -6, 0xBD703FFBFD06F01AuL, 0x7B08FBD9C60BFF77uL), (-1, -8, 0xC5FC557209F06856uL, 0x9346CC45FECBDDB5uL)),
                ((+1, -8, 0x999EC6C60F394F92uL, 0x87002B6BBE9A7DBDuL), (+1, -6, 0xE73CE27A39C017F5uL, 0xA4037D224CA3AF1CuL)),
                ((+1, -8, 0x89E7A974321B972AuL, 0x648900E124F4DD1AuL), (+1, -7, 0xA561B8D9DF050AB4uL, 0xBF371EBEEB479492uL)),
                ((+1, -11, 0x8890779886C83935uL, 0xB41076A770DC10F5uL), (-1, -12, 0xC439527C98B4584CuL, 0x0D187C3130F8C0CBuL)),
                ((+1, -13, 0xEC10582BBB64DB62uL, 0x17614A9CB6D087ABuL), (+1, -10, 0xB9DE14E6CDD46ABCuL, 0x3322B205D5D436AEuL)),
                ((+1, -15, 0xD4266B0B8BD9D831uL, 0xA9B0131CCC9BEC3BuL), (-1, -13, 0x80E639358967A255uL, 0xF98268BF308E76B9uL)),
                ((+1, -18, 0xF9955AECB2B8CA2CuL, 0x7B6B223768B7AA76uL), (+1, -14, 0x8F6636BB02C95A19uL, 0x4C4D400635DF8550uL)),
                ((+1, -20, 0xB826AC778D23703CuL, 0xD288F935CFF5C6FEuL), (-1, -19, 0xE40182658FD444B7uL, 0x996F40741A399038uL)),
                ((+1, -23, 0xCDBF000844E37A4FuL, 0x0CF79BF5F791050AuL), (+1, -20, 0x92531E9E7D7E3604uL, 0x96F8820DD26756B1uL)),
                ((+1, -27, 0xCAE9F71B497BC688uL, 0x43A388EC4300F640uL), (+1, -30, 0x964BFEFC5BB2B45DuL, 0x8479C5FACCD6D127uL)),
                ((-1, -39, 0xC10A96339EAF4004uL, 0x9B1E875C0BAEC40BuL), Zero),
                ((-1, -36, 0xA2A72772D2AC5410uL, 0x2E1843C8267335ADuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xC2BF6F4B055A8AB2uL, 0x01F950C3658B0A1FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xF6C46735EAED9285uL, 0x587643DDF14D6D68uL), (+1, -1, 0x91045E8137F42CBBuL, 0xD4D0D6DE2D60DA24uL)),
                ((+1, -2, 0x90CD18B7FC1838CBuL, 0x11A4114B5915E063uL), (+1, -1, 0xBEF2F1EE8606F377uL, 0x1AB78B53B71BBBD5uL)),
                ((+1, -3, 0x914037EB689CA625uL, 0x2884D5A195D72BFFuL), (+1, -2, 0xAB14B61A744B36DEuL, 0x91341972CC5B8250uL)),
                ((+1, -4, 0xACD959E81F51B0EAuL, 0x4AA392DE61886711uL), (+1, -3, 0xE43C96ABC2CFF78CuL, 0xE519F743F75F71BEuL)),
                ((+1, -5, 0x87D4B0CFAC7492BEuL, 0x3689626FE9127CB7uL), (+1, -4, 0xA0658AB7437484CFuL, 0xE9B76F99F679BF61uL)),
                ((+1, -7, 0xD1C40FD6EEDEF57DuL, 0x3C0094029A7E26E4uL), (+1, -5, 0x8A8DFD1BED76CCA1uL, 0xA09A1E7A15C27A89uL)),
                ((+1, -9, 0xFC4A0432A2858DFDuL, 0x10CADED05B067A7EuL), (+1, -7, 0x954BD410D0E91D6DuL, 0x8AC2F6C18703A184uL)),
                ((+1, -10, 0x869B6294172AD29EuL, 0x7361205294F69528uL), (+1, -9, 0xB1EBB081F90D2A0BuL, 0x1ED871FA45850578uL)),
                ((+1, -13, 0xEC4DBCDB618F0BA2uL, 0xD2C72B706FAE1BC4uL), (+1, -11, 0x8BEA59CDCC385F91uL, 0x30D49DC138E239BBuL)),
                ((+1, -15, 0xAAB41EA12A8F8236uL, 0xCC6869890B6F545BuL), (+1, -14, 0xE21522F970F7228CuL, 0x202C20C8CF6E462EuL)),
                ((+1, -18, 0xC5A145BE5C9233CEuL, 0x489EE530E6B5BC27uL), (+1, -17, 0xE908C7E6514E89C9uL, 0x2AF03549186B22A7uL)),
                ((+1, -21, 0xAAA65901D70F099FuL, 0x29CFE2259ECB8CF3uL), (+1, -20, 0xE3AABB261C58F506uL, 0x8DB2857559AD3BFCuL)),
                ((+1, -25, 0xC5A4BAA94449D836uL, 0xE161F5EFB5AD83F8uL), (+1, -24, 0xE46D7882E5925049uL, 0x3BBB4346C68E79B2uL)),
                ((+1, -30, 0xF21F48A13E3CE628uL, 0x1C4DDB28759B4C77uL), (+1, -28, 0xAA8D7F3DBF6150EAuL, 0x5929CA6D281DFDFCuL)),
                ((+1, -37, 0x9297D555B5380DC0uL, 0x1DBBF5008D69A0B6uL), Zero),
                ((-1, -43, 0xA8FF66070765BCF0uL, 0xDF76887EFA8219FBuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_minus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xCBC9646C7FEB635BuL, 0x0BEDC57221640E92uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xDB025F483498DC8DuL, 0x7C9CB25AD49567F8uL), (+1, 1, 0x89549416C6476EA3uL, 0xF8FAEF5933770F3FuL)),
                ((+1, -1, 0xE574BAC1DD665F47uL, 0x892D4D6EFC7FB9D2uL), (+1, 1, 0x8FCE604A8368EAB0uL, 0xC12AEF445E3D8B9CuL)),
                ((+1, -1, 0x999C10D1AA7053C4uL, 0x914B563F85D56B41uL), (+1, 0, 0xC0839E5700BDBBEEuL, 0xAA0A4B8EF3E51EC4uL)),
                ((+1, -2, 0x921CCA45BD3B3557uL, 0xBB292F850B191BACuL), (+1, -1, 0xB71EC526BEBFC957uL, 0x06EEC94DC99CEB28uL)),
                ((+1, -4, 0xD0DAC7976BAAD524uL, 0x4336AFFDA2C95D11uL), (+1, -2, 0x82E1932946F62FCDuL, 0xD5FE6DD75BF2C33EuL)),
                ((+1, -6, 0xE839F0F7D8E7EC22uL, 0x56EE0CA9E5258B2DuL), (+1, -4, 0x91870F84CF4762CEuL, 0xF9381443040A8AEAuL)),
                ((+1, -8, 0xCDB7DDD5EF7B5541uL, 0xAE9FEAB38E9FDC0EuL), (+1, -6, 0x80EA092C4574C753uL, 0x20FD77BD5FD75559uL)),
                ((+1, -10, 0x93B360ABC80FAC01uL, 0x3DEBB38A7D45CF52uL), (+1, -9, 0xB91D2EA30BD62850uL, 0x6B7D1FBBF4FE07A6uL)),
                ((+1, -13, 0xADF3D7CEC2EB26BDuL, 0x6E66E0CB60FCC994uL), (+1, -12, 0xDA053AB8D0FB7D5AuL, 0x464A4027475B005FuL)),
                ((+1, -16, 0xA92E6F38B69C3C81uL, 0xB330C196CE799AD6uL), (+1, -15, 0xD408BCF8C375B2AEuL, 0x42614C1E98A5C724uL)),
                ((+1, -19, 0x881CC74303E6F04BuL, 0xEDD11D430005CC19uL), (+1, -18, 0xAA97F126E3AEAFD4uL, 0x60757BFC26236C4AuL)),
                ((+1, -23, 0xB4453BFEF6F2BA4AuL, 0xE70E6080D2494385uL), (+1, -22, 0xE1EF286365C84DE0uL, 0xC870BA25FA43DA1CuL)),
                ((+1, -27, 0xC18A7B8CF628607DuL, 0x629FD3E70EAD2FD7uL), (+1, -26, 0xF2916653204E74E3uL, 0xAC577F363333AB9AuL)),
                ((+1, -31, 0xA3749A324B7EF567uL, 0xB57115074587B431uL), (+1, -30, 0xCCDC777C12E31100uL, 0x47D844B15A01483EuL)),
                ((+1, -36, 0xC93A0341FEFCFECFuL, 0x1B03F1F26B5737C1uL), (+1, -35, 0xFC332A7A8267DB4FuL, 0xB04E389E58C5A789uL)),
                ((+1, -41, 0x9D4F9E4A2794E83FuL, 0x32D31A11C61BC877uL), (+1, -40, 0xC52906A3AEF41539uL, 0x196F84BCD0C75538uL)),
            }));

            private static ddouble PlusValue(ddouble x) {
                Debug.Assert(x >= 0);

                if (x <= 64d) {
                    ddouble y;

                    if (x <= 1d) {
                        Debug.WriteLine("pade minimum segment passed");

                        y = ApproxUtil.Pade(x, pade_plus_0_1);
                    }
                    else if (x <= 2d) {
                        y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                    }
                    else if (x <= 4d) {
                        y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                    }
                    else if (x <= 8d) {
                        y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                    }
                    else if (x <= 16d) {
                        y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                    }
                    else if (x <= 32d) {
                        y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                    }

                    return y;
                }
                else {
                    int exponent = double.ILogB((double)x);

                    ddouble v;
                    if (exponent < 8) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -6)), pade_plus_expp6_8);
                    }
                    else if (exponent < 16) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -8)), pade_plus_expp8_16);
                    }
                    else if (exponent < 32) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -16)), pade_plus_expp16_32);
                    }
                    else if (exponent < 64) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -32)), pade_plus_expp32_64);
                    }
                    else if (exponent < 128) {
                        v = ApproxUtil.Pade(Log2(Ldexp(x, -64)), pade_plus_expp64_128);
                    }
                    else {
                        v = 2 * RcpPI;
                    }

                    ddouble y = v / x;

                    return y;
                }
            }

            private static ddouble MinusValue(ddouble x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 1d) {
                    ddouble y = ApproxUtil.Pade(1d - x, pade_minus_1_0);

                    return y;
                }
                else if (x <= 8d) {
                    ddouble v;
                    if (x <= 2d) {
                        v = ApproxUtil.Pade(x - 1d, pade_minus_1_2);
                    }
                    else if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }

                    ddouble sigma = Exp(x * pi_half - 1);

                    ddouble y = v / Sqrt(sigma) * Exp(-sigma);

                    return y;
                }
                else {
                    return 0d;
                }
            }

            public static ddouble Value(ddouble x, bool complementary) {
                if (x >= 0d) {
                    return complementary ? PlusValue(x) : 1d - PlusValue(x);
                }
                else if (x <= 0d) {
                    return complementary ? 1d - MinusValue(x) : MinusValue(x);
                }
                else {
                    return NaN;
                }
            }
        }

        private static class QuantilePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 2, 0xBF02D3DA2465C3F1uL, 0x777A04EF9EF5BE52uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 7, 0xB3432F8B37090313uL, 0xB1D7634232D7FA20uL), (+1, 5, 0x97AD03F603B08CB8uL, 0x4F8637C6F287FF7BuL)),
                ((+1, 11, 0x856DE6D10707BAF5uL, 0xA351A62FCAA77EB6uL), (+1, 9, 0x95A2106DCFB6EDC8uL, 0x9902FE027EA9FF3EuL)),
                ((+1, 13, 0xC191987E4935A863uL, 0x00C0F8E186E11E4DuL), (+1, 12, 0x9DA772E89EC9E41CuL, 0x41913E5DC4881C81uL)),
                ((+1, 14, 0xF7D9BA79634782B7uL, 0x305D6559D8683AF3uL), (+1, 14, 0xBB080F92AA21EB70uL, 0xC92F011A93090660uL)),
                ((-1, 11, 0xCDEE8326869036E4uL, 0xF47053B69CE8D755uL), (+1, 15, 0xE78334802B4841C8uL, 0x18BECD671818FC86uL)),
                ((-1, 17, 0xAEDD57FAE428D477uL, 0x9F851ADF91051090uL), (+1, 15, 0xA7C51A2BDB85BFE6uL, 0x7DB24436E59EA244uL)),
                ((-1, 17, 0xD913D61A0AC28409uL, 0x5C851711AE284153uL), (-1, 16, 0xDBEED9EDC85F4FB8uL, 0x4A62126838721F6FuL)),
                ((+1, 17, 0xF010221BF6A7D217uL, 0x357DE23C3A7C8514uL), (-1, 17, 0xCA3B19B8B72EEEB7uL, 0x6CAA72148740723AuL)),
                ((+1, 18, 0xD55B42225740A710uL, 0xB2F511E312E605C3uL), (+1, 14, 0xE4290C4306C230EAuL, 0xAE190AF9B1E74BE6uL)),
                ((-1, 17, 0x91F08A205649181DuL, 0x1C6F72920FCEB8A0uL), (+1, 17, 0xB498E3AFCBEBCB76uL, 0x7714D9A5C6DBF5B8uL)),
                ((-1, 17, 0xDEFF25B2843386C5uL, 0x8E47077398D5D7CEuL), (+1, 13, 0x85711EFA71DE3311uL, 0x3DE81F11A76E4085uL)),
                ((+1, 15, 0xB52D3842CB0514FFuL, 0xFA61D51F2160CBB9uL), (-1, 15, 0xA14153D8F7B5D920uL, 0x0A596D2FBABD4C09uL)),
                ((+1, 13, 0xEFAE712565091A75uL, 0x982519362433AE68uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 1, 0xB5A6BE1C589F63C8uL, 0x19656AC98B6723B5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 4, 0xB8F640732FF224B6uL, 0xCE6A448A20A99B24uL), (+1, 3, 0xD0289F5B91A79EE5uL, 0x3FDF3F1BB1926220uL)),
                ((+1, 5, 0x80AB4A35C987344AuL, 0xE237114946D65560uL), (+1, 5, 0xE952BBB04B5704A1uL, 0xDFB5D2E0CEC6B301uL)),
                ((-1, 7, 0x9F751FC0EBB7210EuL, 0x02C95140CA1E82F2uL), (+1, 6, 0x9CAB7EDB0CEF768BuL, 0x96473D2DA2BAA516uL)),
                ((-1, 8, 0xC266817848CEB241uL, 0x7C3FD123A723E50EuL), (-1, 7, 0x875791395537ECE2uL, 0x7E1A4E8AE013DE1FuL)),
                ((+1, 8, 0x8A02AEF983A6DE92uL, 0xA71F83C96D7F1A4FuL), (-1, 8, 0xC70A272B38DD8DB0uL, 0x524B16D52A7701C5uL)),
                ((+1, 9, 0xFC25D31119CB9B1FuL, 0x2CFD2F29E02EF376uL), (+1, -1, 0xF9EAFA364522E58CuL, 0x68EBF03DF4560A48uL)),
                ((-1, 7, 0xB4B4E60E10B3680DuL, 0xA90C40E18D7EADDEuL), (+1, 9, 0x831E16024AD76680uL, 0x94880510BA94251BuL)),
                ((-1, 9, 0xE15E177A8E8C55B9uL, 0xD7704138CBC3FBB3uL), (+1, 6, 0xD7AD00BAE8297A2EuL, 0x2715639AEB4CFCE3uL)),
                ((+1, 6, 0xBDD0C7A158E529BBuL, 0x90460C62016F8008uL), (-1, 7, 0xE5A24152C73923BFuL, 0x05D342BFF593B038uL)),
                ((+1, 7, 0xCF7677E3A3DF7E57uL, 0xD9DACB605042A452uL), (-1, 4, 0xD478018A79E4C07CuL, 0x422B233959E945CBuL)),
                ((-1, 4, 0x857BE42C3E6950E7uL, 0x8B37E86F5D69A58DuL), (+1, 4, 0x80DC416A9A8D4322uL, 0x1F8DAC4C99F4B483uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_0p375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xCCECC552E5E5B8D0uL, 0x6F753D7B798FACB7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x83E2AC3F43AAC803uL, 0xF76564CA3114ACE3uL), (+1, 2, 0x9B4555E2DEF2B2EFuL, 0x246D124E06B69494uL)),
                ((-1, 4, 0xB6CF461E01516894uL, 0x0C41BC3A3A2B61A8uL), (-1, 0, 0x9F4E59A961889F2DuL, 0x8E9B88FBE8679B16uL)),
                ((-1, 3, 0x9960A99FA6BBA5D9uL, 0x93702B40852903C6uL), (-1, 4, 0xFF22512E5351FC1EuL, 0x3C30C71C717E7FEDuL)),
                ((+1, 6, 0xCE69FEAAF10944A2uL, 0x32A29FAA7540332DuL), (-1, 3, 0xF2A609F5EF949171uL, 0x0077750083A3543BuL)),
                ((+1, 4, 0xFA4A36A2805E17DDuL, 0xC34910D9C69E67D5uL), (+1, 6, 0x99EE36091B6EA866uL, 0xEB6CFD50A1CC0E38uL)),
                ((-1, 7, 0xC22E0F980A1ACFEAuL, 0xF1E9C5E1317D5743uL), (+1, 5, 0xA6B880E277834524uL, 0xD8A8ACADD4AD1F11uL)),
                ((-1, 4, 0xCE8469757719AC2EuL, 0xDCF9251C85385108uL), (-1, 6, 0xA11893ABB1A53CF2uL, 0xF9C1357675B37849uL)),
                ((+1, 7, 0x90A8E06F533E7326uL, 0x694CAC562B16440DuL), (-1, 4, 0xF9645497B6126211uL, 0x559569DF0CC1ED3DuL)),
                ((-1, 0, 0x9B1272A9828396B7uL, 0xADA21A1839B3D5EEuL), (+1, 4, 0xF91B0F56C7316DA8uL, 0x679C8D5B315B60C5uL)),
                ((-1, 4, 0xE694F72D1F78339AuL, 0xAC3201893C23CAC2uL), (+1, 2, 0xA4DB18ABE1790314uL, 0x372D7935ACF44E45uL)),
                ((+1, 0, 0xBA1A8B3024505DF3uL, 0x23236538177B3F77uL), (-1, 0, 0xF748FD8A49FC7FC3uL, 0x4AD069ABCF49B863uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xBF02D3DA2465C3F1uL, 0x777A04EF9EA91152uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0x89312D97675D46F4uL, 0x46DB00DE9223ABF6uL), (+1, 0, 0xB91EF28C54DD6F8DuL, 0xD5C68FC6D5AD0606uL)),
                ((+1, -1, 0xB17EEB48A015C824uL, 0xFE8DAB3EDFE5FE4FuL), (+1, -1, 0xFA0EDED2C4C569F8uL, 0xA25300F9AB5BB42EuL)),
                ((+1, -2, 0x940E5DC905B502A6uL, 0xCF50999BB60C31C0uL), (+1, -2, 0xDAC16FFF5EFD0E78uL, 0x014033ED676146E4uL)),
                ((+1, -4, 0xB6B23FCBF547CB18uL, 0xFFA1B77D263B8F17uL), (+1, -3, 0x8AEA147B6501F736uL, 0x44DA9B54803088E6uL)),
                ((+1, -6, 0xAA02D4F929B45170uL, 0x88B8A9C5C2337638uL), (+1, -5, 0x84D6894633B091C1uL, 0x1525838AB962572BuL)),
                ((+1, -9, 0xF5C37B98A7CDC6E5uL, 0x66130F7CD67C6C3CuL), (+1, -8, 0xC31C852F26AE46D6uL, 0x65AB2EFFDCEEC4B2uL)),
                ((+1, -11, 0x8BE26FD2461E5F8FuL, 0x69B793C2C460943FuL), (+1, -11, 0xDC6516DD89542638uL, 0x3ABE74B96E5E74B7uL)),
                ((+1, -15, 0xEDF5F8672C136A46uL, 0xD4ABF9D15941DB5CuL), (+1, -14, 0xBCCE0949E74E924CuL, 0xA908114CF027F2C3uL)),
                ((+1, -18, 0x9FA643406C260279uL, 0xDB6D85AC13C36E06uL), (+1, -18, 0xF0D11F8F991FB4D5uL, 0xF20781B948660DCDuL)),
                ((+1, -22, 0x819A39D3C6340746uL, 0xC711B16298E65589uL), (+1, -22, 0xD8C5CA65172FD78BuL, 0x8556AC1EFD0C9E04uL)),
                ((+1, -27, 0xA22F07EDF5AB22BFuL, 0x5E6DBAB440C47222uL), (+1, -27, 0xF2A37CE127B9D35DuL, 0xF6AA800F73E42CE0uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm4_6 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xB95F9EDB18B12DE9uL, 0xDD968503D7F212D5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xBA579859D668C1EBuL, 0x63EF70659BE92E5FuL), (+1, 0, 0x85D6BFA7E38CB8A8uL, 0xA7DEBBF9567AC042uL)),
                ((+1, -2, 0xBFFFB6CB8B317ABDuL, 0xAB3CC54481DECFABuL), (+1, -1, 0x900C5B7210D3F231uL, 0x505E99B6DB09ED07uL)),
                ((+1, -3, 0x8586835D3FAEBA90uL, 0x70CAEF76D0ECBBF6uL), (+1, -3, 0xCD441DCE350238E5uL, 0xFE4E9BD75168D49AuL)),
                ((+1, -5, 0x85A5C75E44172AE7uL, 0xCF7AA697B6972520uL), (+1, -5, 0xD18865AC7ADE2F80uL, 0x030AC551F200D54EuL)),
                ((+1, -8, 0xCA6AFF843CECB0CFuL, 0x158FD0DBFFFFE237uL), (+1, -7, 0x9FDC817279D97722uL, 0xFE08DF7B3417CC3DuL)),
                ((+1, -11, 0xEBAC8711CD96E91EuL, 0xC3E4CB1E4D752A4CuL), (+1, -10, 0xB9845C3CF11AD0C0uL, 0x7B26F7E7B7CF0443uL)),
                ((+1, -14, 0xD1F0F95588DEDB28uL, 0x9473EEB4C96C6E35uL), (+1, -13, 0xA513C5C362FBC885uL, 0x1EBB16BD7CE2FDC6uL)),
                ((+1, -17, 0x92BF26930C3FF28EuL, 0x0A74521EDF64737BuL), (+1, -17, 0xE33EE22F646C379EuL, 0xB34012BBB672AC76uL)),
                ((+1, -21, 0x95DFD20D05D663A5uL, 0x726CB537A800509FuL), (+1, -21, 0xF22482D03096E7E9uL, 0x5F0A2CC2633274A5uL)),
                ((+1, -25, 0x825647AB1A4635E9uL, 0x97288BD4F48262DBuL), (+1, -25, 0xC5C57561425CA8CCuL, 0x3A0206819A87DE5CuL)),
                ((+1, -30, 0x8DC15F2D2A413186uL, 0x1FB578650368DE3CuL), (+1, -30, 0xE74AE798EFE93FACuL, 0xFA55F4959990D0E4uL)),
                ((+1, -36, 0xFE88F3D082C1D60EuL, 0x239EB7E2FB1B1F29uL), (+1, -35, 0xC17CBD3E1CA83D90uL, 0xDC2E8F00D1EBFF73uL)),
                ((-1, -47, 0xBFA2747451C818D4uL, 0xF52F88F89B6C8911uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xAC62CBF54A8361E5uL, 0x30A7DB042EF67DC3uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x87B02B6D6E88CC0EuL, 0x4815ABAEBA93F98CuL), (+1, -1, 0xD0772E4097FD4016uL, 0x358459E6BCEEED9EuL)),
                ((+1, -3, 0xD9D6265144629373uL, 0x8D96FB180B6EF40EuL), (+1, -2, 0xAA37E6526159ED57uL, 0x13B4133CF47C26C0uL)),
                ((+1, -5, 0xE3802E08F55E0DC6uL, 0x380174D662B62189uL), (+1, -4, 0xB2D1C24D4FC2E7ABuL, 0xAD521D64B1BA4E5DuL)),
                ((+1, -7, 0xA94F458AD64D2C32uL, 0x25D86855D5E8D6AFuL), (+1, -6, 0x8540D3497A554AF9uL, 0xE1093AB25F647AF8uL)),
                ((+1, -10, 0xBBFECC6B87057CE2uL, 0xE6856FF80AA071C3uL), (+1, -9, 0x9399D082A6713519uL, 0xAB4872418D674DCAuL)),
                ((+1, -13, 0x9E4FBF7111CFFD8AuL, 0xDE85740BDA5A1027uL), (+1, -13, 0xF8978DD654B8C70AuL, 0xB75A68734734F53BuL)),
                ((+1, -17, 0xCD2443B8360F00CCuL, 0xC7062C7488CB1802uL), (+1, -16, 0xA0D94F1D5659F831uL, 0x3C04BEA41F975693uL)),
                ((+1, -21, 0xCA2E7CE0E4835AFFuL, 0xA58C8F10D3456E90uL), (+1, -20, 0x9F8BB92068D84024uL, 0x94CA277055437787uL)),
                ((+1, -25, 0x9872C8E9F4105336uL, 0x7E5332B086E66C9BuL), (+1, -25, 0xEDDE166BF5A39417uL, 0xA69221B3D825960CuL)),
                ((+1, -30, 0x9D98562581AF782EuL, 0xA654D6391D32C052uL), (+1, -30, 0xF935CD428458C948uL, 0x7357543C2FDC6ADAuL)),
                ((+1, -36, 0xCED7491207BC2B59uL, 0x68440DF9EC5FDA76uL), (+1, -35, 0xA1BBD042E63C5B1EuL, 0xCDFC33E889E64A62uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm8_12 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA63AA91E97EC0EE3uL, 0x10E64EFC9DB2650CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x8642F31AFDB49424uL, 0xC4CA7DC16A283331uL), (+1, -1, 0xD18E5A72D55A775EuL, 0x575D65FCD8C90004uL)),
                ((+1, -3, 0xCCF9B3A889C230AAuL, 0x6015BDDA7F0E9748uL), (+1, -2, 0xA0F291CD2E1046CCuL, 0x2DD8D1D2B800E0BBuL)),
                ((+1, -5, 0xC3BF64EEDFA2D6A2uL, 0x42D1AC42271929B5uL), (+1, -4, 0x99E035550D5D5DEFuL, 0x8258D698CF7888D8uL)),
                ((+1, -7, 0x82A2EF5B814BA325uL, 0x22C47CFB9ED79F7AuL), (+1, -7, 0xCD3A383F4B119ADDuL, 0xC129A1EECC54502DuL)),
                ((+1, -10, 0x8147A63E97808F38uL, 0x66239E7B338924DBuL), (+1, -10, 0xCB005A1F28981795uL, 0x97D2D1CA3B1F79C8uL)),
                ((+1, -14, 0xC4C57A027025706FuL, 0xD66F021C5536E4DBuL), (+1, -13, 0x9A8C162361B11C92uL, 0x378BCC163B0127A8uL)),
                ((+1, -18, 0xEBCA0CD6E1A0FA00uL, 0xE4ABADE0A8368D5DuL), (+1, -17, 0xB93C41F99F18DD00uL, 0x0907D1E36E2FDD43uL)),
                ((+1, -22, 0xE12EC4B92C8E1548uL, 0xB86722C56B972224uL), (+1, -21, 0xB0D1FF8EC275848EuL, 0x1FFD3CC25F4B9FC0uL)),
                ((+1, -26, 0xAB68B0ACDF71E42AuL, 0x96348DEC8C4DFB12uL), (+1, -25, 0x869D5B7886828325uL, 0xEA83AC535E360A40uL)),
                ((+1, -31, 0xCDB8A9238348C461uL, 0x2EA0D35DB2ABF8C2uL), (+1, -30, 0xA1A354BAACEF9D2CuL, 0x987570F21504C2BDuL)),
                ((+1, -36, 0xBD136DBFB92E546FuL, 0x01E4AFAA161F780AuL), (+1, -35, 0x946904CA1C6605ACuL, 0x0D958F8FB753C873uL)),
                ((+1, -42, 0xF5C77AE8014FEA3FuL, 0xECAB9D31C82E1927uL), (+1, -41, 0xC126F4D3DF157A15uL, 0xEAA2AC40B4881AADuL)),
                ((+1, -48, 0xBC6D1F9D886FAF99uL, 0xD721D092CF9ADA86uL), (+1, -47, 0x93ECE1E378FC8C8DuL, 0xC85E0622C6C6F835uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm12_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA349ECA8DDF65726uL, 0x0AC9CC99D8481200uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x9D01575E97FDBDE9uL, 0x9EBC7027014D6BC5uL), (+1, -2, 0xF6BEB781AECBE22BuL, 0xFADD69A443FC7A79uL)),
                ((+1, -4, 0x8ECCC13A92262E20uL, 0xCDE4DF0A669BB5D4uL), (+1, -4, 0xE0521C1CC9AC5269uL, 0x9AE8B4EB4B180B88uL)),
                ((+1, -7, 0xA36082B626304FDDuL, 0xBFAC88CC8F820FB6uL), (+1, -6, 0x804B40F0ABB5DFE6uL, 0x4A13E0DDE6A0D5A5uL)),
                ((+1, -10, 0x83C9C5D14141238FuL, 0x30D6C5C08736C20EuL), (+1, -10, 0xCF0957EAB39A8C52uL, 0x9FAFDD3CFCC91610uL)),
                ((+1, -14, 0x9EEEF76D58C29D1EuL, 0x2A9B428B29D3DAE7uL), (+1, -14, 0xF9A6181E02F2BCB5uL, 0x579DC0D0E6AA4685uL)),
                ((+1, -18, 0x938925DD6D881AEBuL, 0x248A1B38EC9F9188uL), (+1, -18, 0xE7BB22A5BBD7ACA8uL, 0xF33A7CAF66735B4CuL)),
                ((+1, -23, 0xD52C7B4DEB82FC5AuL, 0xD410D919D492E115uL), (+1, -22, 0xA7729393633D65DDuL, 0x59D427D1E7240B9BuL)),
                ((+1, -28, 0xEE5F372498EEF6C2uL, 0x3A0224A9369F55EEuL), (+1, -27, 0xBB30DECEC2C5678EuL, 0x40952098901E658DuL)),
                ((+1, -33, 0xC836CA47BF3E18E6uL, 0x466B5A2296EE023FuL), (+1, -32, 0x9D444592BB53C555uL, 0xF5D7F86D6897C481uL)),
                ((+1, -39, 0xE9F25492306ACAFDuL, 0xE5ACA052A795D359uL), (+1, -38, 0xB7B9AE07A54CD0DAuL, 0x779BD1FC1E5B83B2uL)),
                ((+1, -45, 0x9B76171586D02CD2uL, 0xAFA617D37F9EFEEFuL), (+1, -45, 0xF435C38F44BAAB90uL, 0xB8D7D5AC7416731AuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA3004DF7CF4DBC90uL, 0x7ECA581E8272CB84uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x9940605507FB0842uL, 0x853C0FE8949182BAuL), (+1, -2, 0xF0BD7163396BEDDFuL, 0x15E0F9E43F7DD111uL)),
                ((+1, -4, 0x892E3A76B045B7EBuL, 0xEA22D0734565681BuL), (+1, -4, 0xD77B27734A86AB4AuL, 0xBB1F40C748BCF928uL)),
                ((+1, -7, 0x9C4BDE1202D013ADuL, 0x16333E69EF743EF4uL), (+1, -7, 0xF581F765D4CFB595uL, 0x1D421F64DB41B408uL)),
                ((+1, -11, 0xFF282CB95C56DEA4uL, 0x045E7961A0C67AEFuL), (+1, -10, 0xC866DBDCEF77D98FuL, 0xEC23BF4D327B7E05uL)),
                ((+1, -14, 0x9EF86FC081072D8AuL, 0x74CF86BBAD76369FuL), (+1, -14, 0xF9B573440AC2CA82uL, 0xEE2AF2252AB9D956uL)),
                ((+1, -18, 0x9D0881E98C9D773DuL, 0x8E90B47E945F9D31uL), (+1, -18, 0xF6AACE081FB420E3uL, 0xE5FDE148690B543BuL)),
                ((+1, -23, 0xFBD2566752DAB539uL, 0x1097C9C608C66075uL), (+1, -22, 0xC5C7FF920B556B8BuL, 0x066D01B5D084D4E6uL)),
                ((+1, -27, 0xA6374B6AE35353C9uL, 0x6E172BB7328713DCuL), (+1, -26, 0x828B6DBE072C2D10uL, 0x0391B441CD2F19D8uL)),
                ((+1, -32, 0xB5D0A46F9C80DA45uL, 0x9771FCFC739DB700uL), (+1, -31, 0x8ECC7103E720BE29uL, 0x10ED08F377E5B497uL)),
                ((+1, -37, 0xA4B94AC22E5A64F8uL, 0x55EBEBE6E301E86BuL), (+1, -36, 0x815F59473B7EFB4FuL, 0x83F2CC244D0518F4uL)),
                ((+1, -43, 0xF517F0E9F588EB93uL, 0x834AFF5209E5579AuL), (+1, -42, 0xC07F6F112E3D45A2uL, 0xD8BD3659D87A85D1uL)),
                ((+1, -48, 0x92BD90729A90FBE2uL, 0x980180F1B92E1EADuL), (+1, -48, 0xE67F626A4F43DCF6uL, 0x2FA69C107F670066uL)),
                ((+1, -54, 0x87E072E1F12D4005uL, 0xDEBDDBA6D3704BE2uL), (+1, -54, 0xD56FA1C160916356uL, 0x941ABBD2F387E549uL)),
                ((+1, -61, 0xB2941225D6EA2B0AuL, 0x1218118026E1349DuL), (+1, -60, 0x8C412320DA266D53uL, 0xCEE6979B68CF5B36uL)),
                ((+1, -68, 0x86318F49ACFA5F91uL, 0x130AC33BDA2FC6D6uL), (+1, -68, 0xD2CAAE88E3C24D05uL, 0x08A46D6CA006A309uL)),
                ((+1, -95, 0xBF355F9121676AEDuL, 0xDE4F1462FF43270DuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9837C2841A91DuL, 0x78ACC9B4917FA491uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xF4E9D6B776318209uL, 0x71E8F84CCDBD63EBuL), (+1, -2, 0xC05AC4804AE7BAE0uL, 0xB06E73F16DAD28A4uL)),
                ((+1, -5, 0xB4C79ED5ED337F57uL, 0x6FB1BD73B8610196uL), (+1, -4, 0x8DFBEDC5F5B4D6A2uL, 0x45FBD9284BFBC090uL)),
                ((+1, -8, 0xAE40E0C2C6DA0FE3uL, 0x7AB9D93FCA5ED6BFuL), (+1, -7, 0x88DBBB1C6CBBC484uL, 0xF123DF60D99C6233uL)),
                ((+1, -12, 0xF5F45F11C1300D76uL, 0x7DF02CD46E56A2D3uL), (+1, -11, 0xC12C19FFB72CB08CuL, 0x8A89D68BE763F3D2uL)),
                ((+1, -15, 0x871328FE350EFC55uL, 0x919B7176434CD32BuL), (+1, -15, 0xD42CD1592F399EB1uL, 0x47B8048B059EA94DuL)),
                ((+1, -20, 0xEF8722C0C55F50B5uL, 0x22E448FD6E2360A8uL), (+1, -19, 0xBC1FEFC3263770A8uL, 0x1178CA97AD811E2CuL)),
                ((+1, -24, 0xAF78B180FBDE0103uL, 0x0B546BF9742C3D76uL), (+1, -23, 0x89D0A14D83807FCEuL, 0x4980C0794284004CuL)),
                ((+1, -29, 0xD78B390BEC4BD7B1uL, 0xA568317AE1EBB31FuL), (+1, -28, 0xA949A90713813C61uL, 0x0CB89BAD3A1E5717uL)),
                ((+1, -34, 0xDFD0602B05D46073uL, 0x69BB5AA34C23B306uL), (+1, -33, 0xAFC877DD0590606AuL, 0xE97961506040549FuL)),
                ((+1, -39, 0xC4FBDC7B1567463DuL, 0xB08A401FFAF19736uL), (+1, -38, 0x9AB5F319967D97C1uL, 0xACAAD9F0A1DC9A98uL)),
                ((+1, -44, 0x9281493BFFA9A55BuL, 0xE5D9B5129CBD2D9FuL), (+1, -44, 0xE6212A64B682BF54uL, 0x1B0FE924FB545F95uL)),
                ((+1, -50, 0xB64C95CE04E3973BuL, 0x61E8AAF8C9876B27uL), (+1, -49, 0x8F2D6BCC28F18F4BuL, 0x08C1EAE8F1C618EBuL)),
                ((+1, -56, 0xB9DD2F54DA5C4306uL, 0xE6CC3D067987C096uL), (+1, -55, 0x91FA2CDB7B437B96uL, 0xC4832600ECD82CCDuL)),
                ((+1, -62, 0x9539DDF74DA5370CuL, 0x536795203A4432B4uL), (+1, -62, 0xEA675A39E9BA6071uL, 0xCE250497522F0777uL)),
                ((+1, -69, 0xADADEA679A7E5266uL, 0x954FF41ABFE90165uL), (+1, -68, 0x88684E902D6553C2uL, 0x9A98C404E422B7F5uL)),
                ((+1, -77, 0xEC2DD51A2A9AC8A5uL, 0xA777135BD4E3DC01uL), (+1, -76, 0xB97E9CB5A302F67CuL, 0x5AD07DF3240FE317uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_upper_expm64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9836E4E441545uL, 0xF509B3866D4C1747uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0x8EB60EA72C04D067uL, 0xA6404783BA63CC6AuL), (+1, -2, 0xE02B901B0F55B6ADuL, 0xE9CCBC583CFD01C9uL)),
                ((+1, -5, 0xF617AF66ADE609BCuL, 0x810F7371D8D86D60uL), (+1, -4, 0xC147D646EC05FF12uL, 0x285D4E6B669E42F7uL)),
                ((+1, -7, 0x8B0E876FD43C7228uL, 0x8C122001291A7F2BuL), (+1, -7, 0xDA6E09E65BBB3712uL, 0x39242FEE16DD4204uL)),
                ((+1, -11, 0xE74BE6B67D7D5432uL, 0x2A73A62DA473AE00uL), (+1, -10, 0xB5A8EB20EBE796D3uL, 0xC4280E8790905942uL)),
                ((+1, -14, 0x96AA1054D9A8C6FAuL, 0x7176F92A23FC042EuL), (+1, -14, 0xECA9B6EBA21CDF4FuL, 0xC82D38B2A8722836uL)),
                ((+1, -18, 0x9FB8236AAF10EF82uL, 0x9CB2BD8ED993E5CDuL), (+1, -18, 0xFAE2F0019877D258uL, 0x39FA3DC736A8CD8FuL)),
                ((+1, -22, 0x8D43B927BE091134uL, 0x133F728B2D8FB618uL), (+1, -22, 0xDDE5D7FA60FF9B1FuL, 0x824CB1646C050FC1uL)),
                ((+1, -27, 0xD3EF6E3622AE283EuL, 0x6DA02BD1C73586B0uL), (+1, -26, 0xA6741D90136D34D8uL, 0x73D4C9BD3E7D28F4uL)),
                ((+1, -31, 0x8871A9DF318AE5B7uL, 0xA104A5D0385DB862uL), (+1, -31, 0xD6536320CCBC0306uL, 0xA4E5B9A3CBEED24DuL)),
                ((+1, -36, 0x96C8794EBE0DAD96uL, 0x286E97E73F31E7A1uL), (+1, -36, 0xECD97B88495A7A76uL, 0xD41B78F5FD2FDCE4uL)),
                ((+1, -41, 0x92B454875CAEAFDFuL, 0x0AA82DDF8F2A1405uL), (+1, -41, 0xE671587369E340B5uL, 0x0AC870CDF6440BFEuL)),
                ((+1, -47, 0xE461CA36F3341EDAuL, 0xD5629BCDF92872B5uL), (+1, -46, 0xB35EEC8BE80EFE34uL, 0x53D91BFBEDE7D681uL)),
                ((+1, -52, 0xBAAE3DB5A521C35FuL, 0x8D3D3B80296AE7CDuL), (+1, -51, 0x929E5E17D778C54DuL, 0xB9F94F18C50EE535uL)),
                ((+1, -58, 0x870FA2AC222D799AuL, 0xF3304E1E5A8B7406uL), (+1, -58, 0xD42747F911EF07EAuL, 0x33E6D4B36EB000CBuL)),
                ((+1, -63, 0x8982E82AC6B37E71uL, 0x5CE140B582C65D6FuL), (+1, -63, 0xD80098CB820420FEuL, 0xC6183101B6BFD5BAuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -1, 0x96F1891554E9D3CCuL, 0x141B8F8D032F164FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 3, 0xFDD2308700FB1AEEuL, 0x157E7BC01EF647F7uL), (+1, 5, 0x872297D4C8704A00uL, 0xE64EEE2CFA35BE4AuL)),
                ((-1, 7, 0x8BF71A17638C9D45uL, 0xAF4656BCBF74808DuL), (+1, 8, 0xE5E6C2D4A8BB3601uL, 0xAB56966C3A4FDE7BuL)),
                ((-1, 7, 0x9B22F5974770BDA8uL, 0x8A08EE105138C0DEuL), (+1, 11, 0xC4EE0B7B61945B99uL, 0xC474FF3B4CDED141uL)),
                ((+1, 12, 0x97C221CBC3982461uL, 0xA6ADE18D458EF799uL), (+1, 13, 0xA5CDB0601E4C2DB4uL, 0x44B0118FBDA01923uL)),
                ((+1, 14, 0xEB169592E87A24AEuL, 0x6AFF0010DCE09170uL), (+1, 13, 0x9A4A77E254F87D22uL, 0x3A83C5F75E494E37uL)),
                ((+1, 15, 0xDE66CE4598DB9DCBuL, 0xA4564EDE3D2C904FuL), (-1, 15, 0x8A9122FE87055E6AuL, 0x81D21E666E2B15ABuL)),
                ((-1, 15, 0x95D1B8F056CEA080uL, 0x243BD779BCBDFCE2uL), (-1, 16, 0x9AF1F9513F815379uL, 0x73B0AE11950D2F56uL)),
                ((-1, 17, 0xD2225419226FFEC3uL, 0xD505317B1FF5ED6CuL), (+1, 14, 0xEEC3D2DE249A24FCuL, 0xD49DD93809E07362uL)),
                ((-1, 15, 0xCDB239B83AF8EE8CuL, 0xCA39C5468542A1C9uL), (+1, 17, 0x84CA7F7B8887DAF4uL, 0xD9200DDECA08E20DuL)),
                ((+1, 17, 0xE115C98FB830731FuL, 0xF2317A89DD151211uL), (-1, 14, 0x83D8F9208399020CuL, 0xBE50361A2F14EDF7uL)),
                ((+1, 15, 0xAE2D3FEF91E88D93uL, 0xD7C7965F70BA942EuL), (-1, 16, 0x8F12A1A5C45CD50EuL, 0x605532834144214FuL)),
                ((-1, 16, 0x865CDAF616E355B7uL, 0x1C54D6534D3E9064uL), (+1, 13, 0xB75155370EC641F9uL, 0x64CFDB67FE5D632EuL)),
                (Zero, (+1, 12, 0x9302F3D6A36D3656uL, 0x4BCA4A314E4A9E40uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p25_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -3, 0x8567ABB5E361E39BuL, 0x88B00F2DC7268B49uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0x80FB349B0FF89497uL, 0x4D58BFC664AB9750uL), (+1, 3, 0xB94F26A06168D45DuL, 0xCD5767023A38664AuL)),
                ((+1, 5, 0x8D7236660261ACDFuL, 0xA19F7CE52AE692ACuL), (+1, 5, 0xA9D6C8F982AB2FC2uL, 0x78E1822E3CB72816uL)),
                ((+1, 7, 0x9B67E3C829837BFBuL, 0xB5E88CBFF6FD98DCuL), (+1, 4, 0xADB733AB8B782186uL, 0x8964FE94DEF81D37uL)),
                ((+1, 7, 0xB4B2E6C4AC78A6FEuL, 0x85898A4D7A067E6FuL), (-1, 7, 0xAC599343D848AE89uL, 0xDC67B9D49EE2C676uL)),
                ((-1, 8, 0xA3E25013D22EC3AAuL, 0x9C5A207B8FD69FC0uL), (-1, 7, 0xE101C166DAAF1BF5uL, 0x1F30B087250560EDuL)),
                ((-1, 9, 0xB39EBA7EBF4CC126uL, 0xE253C11B0E0F6C37uL), (+1, 7, 0xF02E6B6E585B02C6uL, 0xCF681987CB4A4AEFuL)),
                ((+1, 7, 0x97B542FB781024EAuL, 0xB1338A16EF6CA4C7uL), (+1, 8, 0xB444E1F76CD08491uL, 0x4B7409EFAEFDE6BBuL)),
                ((+1, 9, 0xAC9EE810ECC0F5E6uL, 0x20F1ECEB8FEA4690uL), (-1, 7, 0xA1E4BBFB1B8202D8uL, 0x6A9D248DFD40CAD1uL)),
                ((-1, 4, 0x9886382E16F95296uL, 0xD4D94358E3034D35uL), (-1, 7, 0xA37367955E8A0E35uL, 0x7615400623158291uL)),
                ((-1, 7, 0xA621BB0A903B60BDuL, 0xECD6D70B7596A2AEuL), (+1, 5, 0xADD420149C71B3F3uL, 0xAC063AD60A8DA86BuL)),
                ((+1, 2, 0xADE36A611DCBC603uL, 0xA94D9F11995AAD83uL), (+1, 3, 0x8C7E06C079EB082EuL, 0x69266E8EE96CC423uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_0p375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xA65ED3A6076F1105uL, 0x94936A869F2EB30EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0x971251929C8904CAuL, 0x87F593EC06F9D1CBuL), (+1, 1, 0xA93F33FAB1487993uL, 0xA233E18A5BF98383uL)),
                ((+1, 3, 0x9E073E871E671BA0uL, 0xD7FFAAAC69D6E217uL), (-1, 3, 0x8F72A5D7CEDF9A1CuL, 0x869EE69DFE56849EuL)),
                ((-1, 4, 0xDE41080F65A35B59uL, 0x9860500C20419122uL), (-1, 4, 0xA8539057851BA337uL, 0x8A945CCEF4F42092uL)),
                ((-1, 6, 0x8C9CDB6251887CADuL, 0x22FAA27C34C36B74uL), (+1, 5, 0x88A48F442CDDA517uL, 0xAEFED27127DF931BuL)),
                ((+1, 6, 0x83B9168100AB1965uL, 0x8FA1DBDCC4427155uL), (+1, 5, 0xE55B4CDCCE7F3119uL, 0xA6979C2472C19D55uL)),
                ((+1, 7, 0x9ADF32CA5BC26544uL, 0x8F366042C7A9C00FuL), (-1, 6, 0x8446094D694CA875uL, 0x9D3538D07A152994uL)),
                ((-1, 6, 0xA6CDF7CB1A1F5042uL, 0xB289BD0641849A72uL), (-1, 5, 0xE9F6CC0BBF60754DuL, 0xCF8124C80CB9861AuL)),
                ((-1, 6, 0xF2EE3B5F3362F2F7uL, 0x85D0F98B4DBFEACEuL), (+1, 5, 0xE19BC58F4D1DBB8CuL, 0xB263A6513C024788uL)),
                ((+1, 5, 0xBCE38373BF0431A3uL, 0x94E6EE7D3575B7D8uL), (+1, 3, 0xFF59585DE6CF40AFuL, 0x53B8142B8D9768A5uL)),
                ((+1, 4, 0xB90EDCB87AC32131uL, 0x3F915BFAE72E4A79uL), (-1, 3, 0xD157567B5CBC6C24uL, 0xC9F02810EFB43FC0uL)),
                ((-1, 2, 0x971A2545AF7B6F54uL, 0xA0A5547FD8C6E1E7uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -1, 0x96F1891554E9D3CCuL, 0x141B8F8D065DC997uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0x9B79C814C7DBFC9CuL, 0xA69578CCE9FC384CuL), (+1, 0, 0xBB660D1C0CEC7C5AuL, 0xF1871DC76DA48863uL)),
                ((-1, -1, 0xF68A101B83AAC0C0uL, 0xB373944A6127CCD6uL), (+1, -1, 0xE478C15021E762A1uL, 0xDAD71537DA797993uL)),
                ((-1, -2, 0xC9CB9FE25789E50DuL, 0x5D7829A1CA01E587uL), (+1, -2, 0x963A5005B8F05E20uL, 0x05C6E627689EFA98uL)),
                ((-1, -4, 0xBBFFD3BBF06630E8uL, 0x35E9FC2015AEA471uL), (+1, -5, 0xE6ADCCDA45437BABuL, 0x630392208A89FE74uL)),
                ((-1, -7, 0xCC0497417AC1AD78uL, 0xB71D773BBED7DD62uL), (+1, -8, 0xD11A43F2B962819CuL, 0xB278C1025667EC6DuL)),
                ((-1, -11, 0xFAE315231BA2773EuL, 0x4D7953AC884B680FuL), (+1, -12, 0xD77C7227F566F286uL, 0x56E1899F802B7A3FuL)),
                ((-1, -15, 0xA07691170A3828E4uL, 0xFFFDD3E402C9C118uL), (+1, -17, 0xE4EFBFD4F4B05B23uL, 0xEA4410875D52ECAEuL)),
                ((-1, -21, 0xB098C873E91E1AF0uL, 0xD3F6F5A613D14063uL), (+1, -23, 0xCAE231CA56344001uL, 0x47961F18642543ECuL)),
                ((-1, -29, 0xCDA2C5D5B304F6EAuL, 0x5D63628A08B43EECuL), (+1, -31, 0xABE899E95454AAC8uL, 0x1536824C6B1A4EC5uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, -1, 0xE139AFDD8A9C575FuL, 0xBA48EAFF3F5CD7FCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0xB9AE111E5B73DA35uL, 0xDB617AD67BAC802BuL), (+1, 0, 0xAFE4E881F55B763BuL, 0xCCCDBBE0F0FDB6DDuL)),
                ((-1, 0, 0x81EA9931DB762FB2uL, 0xD37F35A473AB8CCDuL), (+1, -1, 0xD1F2BB1A1018D3F1uL, 0x31841826D187E4FFuL)),
                ((-1, -2, 0xCBE29A1164B2C32CuL, 0x15B8281B8A113EDEuL), (+1, -2, 0x8ECBC6E36ADF53F3uL, 0xD475A016056E98F9uL)),
                ((-1, -4, 0xC6E46623A6DFF7C2uL, 0x2B0C2CD0984F6C6CuL), (+1, -5, 0xF43100FB319C3237uL, 0xA09B90E7B230C591uL)),
                ((-1, -7, 0xFBFED3DC6AF321D9uL, 0xBEE7F9C66AD19AAFuL), (+1, -7, 0x8898427FA6EC92EFuL, 0x1550AAD0C5376D55uL)),
                ((-1, -10, 0xD1FAF9D5E3258577uL, 0x6DA49A21A7F45348uL), (+1, -11, 0xC9D93238968536D9uL, 0x94F4F4677360F1DEuL)),
                ((-1, -14, 0xE3D928C4B421C234uL, 0x5E7DA4B204523224uL), (+1, -15, 0xC26AB85916F870A9uL, 0x85C975FEED90C7C5uL)),
                ((-1, -18, 0x9BD2E29B0B5121B2uL, 0x5AA75966C8039EAAuL), (+1, -20, 0xEB765540C77AF5ABuL, 0xB2175B3F95BA5092uL)),
                ((-1, -24, 0xFC92AA1CC0322539uL, 0x558AC1C23F8210C5uL), (+1, -25, 0xA7BD602AAAC17041uL, 0x72C9058568B2845DuL)),
                ((-1, -30, 0xD96516BE5D07D575uL, 0x21C575965220832DuL), (+1, -32, 0xF9EA3E8589791D5FuL, 0x156D3A1CA11F153FuL)),
                ((-1, -37, 0xA1393E252A01C48EuL, 0x6DBED8B13639D09CuL), (+1, -39, 0x9B21466EEB1B2911uL, 0xEF8F667D0A8DC31EuL)),
                ((-1, -47, 0xF8514D926551269AuL, 0x9279BBE28E97BDD3uL), (+1, -49, 0xB5A5FC0D8D61F0D1uL, 0xE2ADE3DD2F859014uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 0, 0xBF7CBA705E988511uL, 0x1A03391F0FAB076CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0x90CA5BB920170ED9uL, 0x1C9A991EA05EA195uL), (+1, -1, 0xAFED2AE845139D06uL, 0x39A712C5DBEAFBFEuL)),
                ((-1, -2, 0xBD72A79783744BD6uL, 0x29FDDE8EFFDA4B3BuL), (+1, -3, 0xD1FE72439F2F1993uL, 0x85B92037ED392D36uL)),
                ((-1, -4, 0x8CDF44BB2A53C1DEuL, 0xF140479C787FAE4DuL), (+1, -5, 0x8ECED3A8CEC2AD58uL, 0x9F11C75DC92A8945uL)),
                ((-1, -7, 0x837E2BB6A6C814ABuL, 0x68BFE0E6EC0F0FF2uL), (+1, -9, 0xF421FC51747635B2uL, 0x1C885DCB5074A68FuL)),
                ((-1, -11, 0xA08E88BD4E65C694uL, 0xEE623D69FBE72F15uL), (+1, -12, 0x887D3EAFF3881840uL, 0x31509A2DB26E6BB7uL)),
                ((-1, -15, 0x81A08669004CA6F8uL, 0xAB0C06EE09020864uL), (+1, -17, 0xC98A081EDE8C9E74uL, 0x23C62BEA18316AC6uL)),
                ((-1, -20, 0x88D552F3700BB272uL, 0x5B6066F3E4F70D5FuL), (+1, -22, 0xC1EC8685A87F7E90uL, 0x136F9103A3E27A79uL)),
                ((-1, -26, 0xB69BFCE558751D6DuL, 0x8019DD44186AEDA4uL), (+1, -28, 0xEA916CB99B098ABEuL, 0x43F8E4F941FDDDA2uL)),
                ((-1, -32, 0x90B0317D27B0D1D8uL, 0x504610C06140C73FuL), (+1, -34, 0xA6D887800A097CFEuL, 0xB81AF9AFD9CA1FACuL)),
                ((-1, -40, 0xF3C144A532C7958EuL, 0x13F5C51D6B6642B6uL), (+1, -42, 0xF821910F7F209765uL, 0x93642D98BB1824E2uL)),
                ((-1, -48, 0xB0D5B9CE2FDC8810uL, 0xAA573EFC2567CD7EuL), (+1, -50, 0x99B2DC45839102ADuL, 0x7E9C87758E55E841uL)),
                ((-1, -58, 0x84A8AF5EF39ABF4CuL, 0x09C50883B0A95D52uL), (+1, -61, 0xB38ADA78AB34123CuL, 0x16539866ED9B95B2uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0x82627C073A8B94D1uL, 0xF6728CD3AF74C6DCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -1, 0xBF265B93D54472CFuL, 0x126216A1BF8DB827uL), (+1, -2, 0xAFF650330A25CD96uL, 0x753CCD556499E978uL)),
                ((-1, -4, 0xF365A42EEA57A341uL, 0x93D59DFDF39F3E57uL), (+1, -5, 0xD2107AD67C2FD0A1uL, 0x790332C3C3DCC114uL)),
                ((-1, -7, 0xB0AF79A81A07154FuL, 0x072AB3487EC63EF7uL), (+1, -8, 0x8EDCE475DCC69F99uL, 0xF6F137A6F83A9152uL)),
                ((-1, -11, 0xA1692BF07A1F6AF6uL, 0x8C89084171A9DEF6uL), (+1, -13, 0xF4372A703A5B77EBuL, 0xA02A4413A9AF2AA0uL)),
                ((-1, -16, 0xC14A24C1D385CF1EuL, 0x4E0FF23E6CAF0E7AuL), (+1, -17, 0x8883C2F21562FF73uL, 0x5A994174DB3842DCuL)),
                ((-1, -21, 0x994C7900F017CD69uL, 0x2AA77B08D115C4FDuL), (+1, -23, 0xC985C0F50D3FB791uL, 0x08CBCFD8FE1D7A37uL)),
                ((-1, -27, 0x9F2A25FBC1ECBA52uL, 0x98E60451EC63B190uL), (+1, -29, 0xC1D4C4EEF4B6648FuL, 0xFA09F2AEC794FDE4uL)),
                ((-1, -34, 0xD11C440A21FF069EuL, 0x5257D746FF7E8B7BuL), (+1, -36, 0xEA54E9D09D3C6345uL, 0xE3725515989B8AD0uL)),
                ((-1, -41, 0xA32FEC1A396E6470uL, 0x36B4BD1D6DF8C5E0uL), (+1, -43, 0xA690F7B43FBD33E8uL, 0x8544872862F88356uL)),
                ((-1, -49, 0x875EFA7A2D5E3288uL, 0x874F534F8E2AD24DuL), (+1, -52, 0xF783CAE6BC360CDCuL, 0x62393805E00E7DD9uL)),
                ((-1, -59, 0xC1329ACF17B3CCDDuL, 0xB733D8B1C5330203uL), (+1, -61, 0x992BD8E2B6DB806FuL, 0x40FAC52101FC43F8uL)),
                ((-1, -70, 0x8DEF6DE189D61050uL, 0x5296D74EA4704080uL), (+1, -73, 0xB2BB6B2B462C50EEuL, 0x81CBD3759FFB0578uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xA248E9007D6242EBuL, 0x1B024AC5EFC0241DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -2, 0xE9EF4FBE4AC86536uL, 0x9EABE847B4C02162uL), (+1, -3, 0xAFB407423E7431EBuL, 0xB86E01904BDD8138uL)),
                ((-1, -5, 0x929AE51D2D71F6BAuL, 0x511F3C03E222F3FAuL), (+1, -7, 0xD1709267300E03A6uL, 0x2A93739B5D398266uL)),
                ((-1, -10, 0xD1B4D9F3D64E572BuL, 0xF289D7E672FB58C2uL), (+1, -11, 0x8E37FE09B09DC47FuL, 0x535E7E24839F57C1uL)),
                ((-1, -15, 0xBCE671A8A8D24D4DuL, 0x1784DA5062742C66uL), (+1, -17, 0xF2BB1B20C43BCF30uL, 0xDFE45F271B4AD1D9uL)),
                ((-1, -21, 0xDF31037FFB5DEDF7uL, 0x827609866D007920uL), (+1, -22, 0x87772BF1FC7F9385uL, 0x0948F59B33E11FE6uL)),
                ((-1, -27, 0xAEBDD5997F691162uL, 0xFD8D120E631E3CBAuL), (+1, -29, 0xC7A4879A1AF64801uL, 0x4CCBED0F1B3AE055uL)),
                ((-1, -34, 0xB3296614ADB2113BuL, 0xB0E11ED7AC040B02uL), (+1, -36, 0xBFB2B0CCEDC8E192uL, 0xD9AE6548DE5985D8uL)),
                ((-1, -42, 0xE87B14C87659AD67uL, 0x28AEF9C7F8948817uL), (+1, -44, 0xE75A37D319B06219uL, 0x4B0A23B3249BF043uL)),
                ((-1, -50, 0xB32C1EC76AC7869BuL, 0x1A4C2FEB5F862CC6uL), (+1, -52, 0xA428D11432E03D46uL, 0x2E5164AB59606503uL)),
                ((-1, -59, 0x92B5F3AAE32D8466uL, 0x6891362B3ED1E36DuL), (+1, -62, 0xF380F3C200C285EEuL, 0x6129034616312EDAuL)),
                ((-1, -70, 0xCE6321E22B4D55D4uL, 0xD1C7F2237893F3C0uL), (+1, -72, 0x966B22727CB89A12uL, 0x55ED4C167298F282uL)),
                ((-1, -82, 0x94D2AA5299ACA4ACuL, 0x20320E73ABAD4C5FuL), (+1, -85, 0xAF344EB6B2FBBD1FuL, 0xA053F32185A83C59uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xC09BDF118016657BuL, 0x22249BC58804BF94uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -2, 0x89492A39FD36F468uL, 0x869605EC72D3C934uL), (+1, -4, 0xAF54438376E8051EuL, 0x3492F5DED0253360uL)),
                ((-1, -7, 0xAA38D0E711260949uL, 0x1D6FEEA426DA3115uL), (+1, -9, 0xD08BE287EF973A7AuL, 0x9E378CD737C8127FuL)),
                ((-1, -13, 0xF0F0237C0EA69910uL, 0xAB8E2B7A4E09EA35uL), (+1, -14, 0x8D4EB9265767EA35uL, 0x2C60D1593B841732uL)),
                ((-1, -19, 0xD6CEE2DFD7C0D885uL, 0x41E4DE311F479E15uL), (+1, -21, 0xF0A7A9BD15D52CECuL, 0xAA7DD00F7524B75BuL)),
                ((-1, -26, 0xFB3E07F87C490AAAuL, 0xBBFEA772C63C0005uL), (+1, -27, 0x860428BB3D2E02F5uL, 0x80EB2629A5EEF7D2uL)),
                ((-1, -33, 0xC2BCD4BB053B2725uL, 0xFBD8518633783DB7uL), (+1, -35, 0xC5143EB7EE62FBCFuL, 0x1C7D5015ECF71E58uL)),
                ((-1, -41, 0xC5AB83A90798E2DFuL, 0x675BD30FD77C768CuL), (+1, -43, 0xBCD3CAE6361B6D64uL, 0x480864EFB30AADB0uL)),
                ((-1, -50, 0xFDE82301FFA8AFCFuL, 0x8AE77A5EDDA02CCCuL), (+1, -52, 0xE365CA00C2529811uL, 0x853D5948DCE3FF92uL)),
                ((-1, -59, 0xC1A64F4B76EC0EEFuL, 0x48C0676664A2D31AuL), (+1, -61, 0xA102555B8A6654C6uL, 0x2EFBC4DBB2333475uL)),
                ((-1, -69, 0x9CD034CD80512426uL, 0xEF8A73ED86DB5B7FuL), (+1, -72, 0xEE53FE08AD90CE16uL, 0x794B371B2194E02DuL)),
                ((-1, -81, 0xD9D9FEE79BC4391FuL, 0xBD9514D39ED3D832uL), (+1, -83, 0x92EAD8781103D2F2uL, 0xAB34E877B3938866uL)),
                ((-1, -94, 0x9A837027B1EC183AuL, 0x12923A8D32BA02C8uL), (+1, -97, 0xAAC82533A1CCE488uL, 0x6A50E63AA637B978uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm128_256 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xDE071A2678E33533uL, 0xA2052E48DECBFEBBuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -3, 0x9CFFC7C8F2A6FE5EuL, 0xB8F9BD051A411E55uL), (+1, -5, 0xAEF7D6FB446AB5CFuL, 0x2A05AD75F44E4BE4uL)),
                ((-1, -9, 0xC12292EF1DFDE8A7uL, 0x14EBB0B294EB82E7uL), (+1, -11, 0xCFB046E87A737009uL, 0x3D43B5E0E3DF6E95uL)),
                ((-1, -15, 0x879E305EE64771E4uL, 0xD3EA51BFA0B7499EuL), (+1, -17, 0x8C6FDEF03C811709uL, 0x9BC33D804B51197AuL)),
                ((-1, -23, 0xEFEFF94A9A124478uL, 0xF872263EEA01F405uL), (+1, -25, 0xEEAEA6CC23519816uL, 0x08212E5D5C5B0CB7uL)),
                ((-1, -30, 0x8B38EB82B3019EBFuL, 0x98044D094811A0B6uL), (+1, -32, 0x84A58665FF0FE40FuL, 0xEE2DD0AB18A574ADuL)),
                ((-1, -39, 0xD62011B3EFAB5FD4uL, 0xE02C146A8C7D9527uL), (+1, -41, 0xC2AB7A0EBEC9B161uL, 0x4B5FCC91D48FF231uL)),
                ((-1, -48, 0xD79CF4612D7913C7uL, 0x3C3360531BBD58C6uL), (+1, -50, 0xBA251574424FED15uL, 0x34CD65B2AA70EBCBuL)),
                ((-1, -57, 0x89570A1D5CAE77F2uL, 0x32F1C788F46AA02EuL), (+1, -60, 0xDFB939D6217073C3uL, 0x36A92FBE3586F413uL)),
                ((-1, -68, 0xCFB1225FFB0A2A22uL, 0xF83FA1BE5AEE6CCAuL), (+1, -70, 0x9E197AA0F02328F9uL, 0x101995AC936E1F6FuL)),
                ((-1, -79, 0xA69FAAAD752AFC24uL, 0x6EE58DE771557D70uL), (+1, -82, 0xE99377902D326418uL, 0x529C7B86DD9F079AuL)),
                ((-1, -92, 0xE5041E671ED84EFFuL, 0x72C12962D49D2976uL), (+1, -94, 0x8FB8BA4BD1C469B2uL, 0x015FCDB4E57D7ADDuL)),
                ((-1, -106, 0xA01AB849E0DCEDB2uL, 0x2BBAC455AD3ACA98uL), (+1, -109, 0xA6C4B7686668E2D4uL, 0xA7A670CDFE2EB074uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm256_512 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 1, 0xFAEDF6FF636F44C3uL, 0xBBEF07D954A33C3AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -4, 0xB063D7007211BAD4uL, 0x59F82C62F14FF69AuL), (+1, -6, 0xAEAB3E73C0DFED32uL, 0x40F584FE6EB07860uL)),
                ((-1, -11, 0xD7B4BBA633B6BCA4uL, 0x8A16296ABF661B25uL), (+1, -13, 0xCEFAEFBF4FC60642uL, 0xF5A5A97AC23D7607uL)),
                ((-1, -18, 0x9690ED706027A84DuL, 0x43C89A6F9255ACE0uL), (+1, -20, 0x8BB886420B1A5C6EuL, 0xB5D91BCDDFA05642uL)),
                ((-1, -26, 0x8464658161187CD4uL, 0x26483B885E4FA23CuL), (+1, -29, 0xED10B33D63D8FDAEuL, 0x64EDD05695D35C64uL)),
                ((-1, -35, 0x98B4353B269D0B73uL, 0xBB0F0043E911603CuL), (+1, -37, 0x83872F2B72A4D04CuL, 0xCA9B607CB5EF0F8AuL)),
                ((-1, -45, 0xE96690F40F626059uL, 0x6D2B175F92BD4093uL), (+1, -47, 0xC0B5AF493CE91466uL, 0x1951917BB9992DDFuL)),
                ((-1, -55, 0xE985181864C4EE8CuL, 0x61A0E41948FB8A98uL), (+1, -57, 0xB7F87B80DFBB367CuL, 0x55517FD868051A27uL)),
                ((-1, -65, 0x93C187EF1226C494uL, 0x45C1B271D06571AAuL), (+1, -68, 0xDCC1B7A8FD5D0DA3uL, 0x04825DE2DD4DC596uL)),
                ((-1, -77, 0xDDDD86873C6EA1DAuL, 0x847CB35E9498ACF2uL), (+1, -79, 0x9BC23DA089E3724CuL, 0xD0A9CE8B7E7AF3B3uL)),
                ((-1, -89, 0xB09DFF94BC6A2DC6uL, 0xC447BB66F112891EuL), (+1, -92, 0xE5C45C7149E78FCDuL, 0x81EA6526595BBB4EuL)),
                ((-1, -103, 0xF08EF5A5CF1C9F09uL, 0xC03B28D1B918BA08uL), (+1, -105, 0x8D2B5471F6B219A0uL, 0x05B89CA964B3C464uL)),
                ((-1, -118, 0xA6141F7B8728E575uL, 0x0C21BD205F3EF38EuL), (+1, -121, 0xA393091C0D92F4F4uL, 0x23E2F72318489285uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_lower_expm512_1024 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((-1, 2, 0x8BC4E6F856417BE0uL, 0xB84D5AADC0D2AD41uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -5, 0xC39EB630558EABE4uL, 0xD3AA13178A6691F1uL), (+1, -7, 0xAE707A062A5B694AuL, 0xDFE856CD4DDF5FA6uL)),
                ((-1, -13, 0xEE2331560014F1D0uL, 0xE66CF1E6BD199D16uL), (+1, -15, 0xCE702DEEB41FF829uL, 0xE9AF8C6FE7AB05C6uL)),
                ((-1, -21, 0xA575661CCCA67D8BuL, 0x75B03A687D4C03B1uL), (+1, -23, 0x8B2C9C3987986798uL, 0x3D1A9030A3E22C7EuL)),
                ((-1, -30, 0x90CE05EBD82A0D40uL, 0x64158491A406EC73uL), (+1, -33, 0xEBD5A872E16BF0BFuL, 0x53E0CE00A612C1C5uL)),
                ((-1, -40, 0xA6380CD75D6E61C4uL, 0xE78D80771F24EC42uL), (+1, -42, 0x82ADD8A184EE1ADDuL, 0x874CF91082ECE649uL)),
                ((-1, -51, 0xFCCCFBC792E8DA41uL, 0x7F32A5E78EBD9A91uL), (+1, -53, 0xBF39D4D79C172B3BuL, 0xFDA18022ACFA2349uL)),
                ((-1, -62, 0xFBA104CBB170DDEDuL, 0xE1C607CAAAE5E648uL), (+1, -64, 0xB6544103032FC059uL, 0xD678FA7CE9D183F2uL)),
                ((-1, -73, 0x9E5A18FDA08B1004uL, 0x44A0927924B2A7F6uL), (+1, -76, 0xDA85CBAAFD0E7AE6uL, 0x00CCE82F9EC5EBB7uL)),
                ((-1, -86, 0xEC648376EE249237uL, 0xC2AFC0046F906767uL), (+1, -88, 0x9A002A6AA751E954uL, 0x9234E9C935AC1FB4uL)),
                ((-1, -99, 0xBAF715B3503C395EuL, 0x644F1EF01EC2DA1BuL), (+1, -102, 0xE2E9D2DD44B70694uL, 0x18023AA550329DF6uL)),
                ((-1, -114, 0xFCB261971FFB5DCDuL, 0xAEB340E645664657uL), (+1, -116, 0x8B4302C6DA128760uL, 0xA25523BEE5FBE64AuL)),
                ((-1, -130, 0xAC90B2860DF715DDuL, 0x96DEC8CC891C694EuL), (+1, -133, 0xA13162694E254DADuL, 0x2DC75CB081A97F78uL)),
            }));

            private static ddouble UpperValue(ddouble x) {
                Debug.Assert(x <= 0.5d);

                if (x >= 0.125d) {
                    ddouble y;
                    if (x <= 0.25d) {
                        y = ApproxUtil.Pade(x - 0.125d, pade_upper_0p125_0p25);
                    }
                    else if (x <= 0.375d) {
                        y = ApproxUtil.Pade(x - 0.25d, pade_upper_0p25_0p375);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 0.375d, pade_upper_0p375_0p5);
                    }

                    return y;
                }
                else {
                    ddouble v;
                    int exponent = double.ILogB((double)x);

                    if (exponent >= -4) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 3)), pade_upper_expm3_4);
                    }
                    else if (exponent >= -6) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 4)), pade_upper_expm4_6);
                    }
                    else if (exponent >= -8) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 6)), pade_upper_expm6_8);
                    }
                    else if (exponent >= -12) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 8)), pade_upper_expm8_12);
                    }
                    else if (exponent >= -16) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 12)), pade_upper_expm12_16);
                    }
                    else if (exponent >= -32) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 16)), pade_upper_expm16_32);
                    }
                    else if (exponent >= -64) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 32)), pade_upper_expm32_64);
                    }
                    else if (exponent >= -128) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 64)), pade_upper_expm64_128);
                    }
                    else {
                        v = Ldexp(RcpPI, 1);
                    }

                    ddouble y = v / x;

                    return y;
                }
            }

            private static ddouble LowerValue(ddouble x) {
                Debug.Assert(x <= 0.5d);

                if (x >= 0.125d) {
                    ddouble y;
                    if (x <= 0.25d) {
                        y = ApproxUtil.Pade(x - 0.125d, pade_lower_0p125_0p25);
                    }
                    else if (x <= 0.375d) {
                        y = ApproxUtil.Pade(x - 0.25d, pade_lower_0p25_0p375);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 0.375d, pade_lower_0p375_0p5);
                    }

                    return y;
                }
                else {
                    ddouble y;
                    int exponent = double.ILogB((double)x);

                    if (exponent >= -4) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 3)), pade_lower_expm3_4);
                    }
                    else if (exponent >= -8) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 4)), pade_lower_expm4_8);
                    }
                    else if (exponent >= -16) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 8)), pade_lower_expm8_16);
                    }
                    else if (exponent >= -32) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 16)), pade_lower_expm16_32);
                    }
                    else if (exponent >= -64) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 32)), pade_lower_expm32_64);
                    }
                    else if (exponent >= -128) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 64)), pade_lower_expm64_128);
                    }
                    else if (exponent >= -256) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 128)), pade_lower_expm128_256);
                    }
                    else if (exponent >= -512) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 256)), pade_lower_expm256_512);
                    }
                    else if (exponent >= -1024) {
                        y = ApproxUtil.Pade(-Log2(Ldexp(x, 512)), pade_lower_expm512_1024);
                    }
                    else {
                        return NegativeInfinity;
                    }

                    return y;
                }
            }

            public static ddouble Value(ddouble x, bool complementary) {
                if (x > 0.5) {
                    return Value(1d - x, !complementary);
                }

                return complementary ? UpperValue(x) : LowerValue(x);
            }
        }
    }
}