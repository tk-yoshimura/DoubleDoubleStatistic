using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class SaSPoint5Distribution : StableDistribution<SaSPoint5Distribution>,
        IAdditionOperators<SaSPoint5Distribution, SaSPoint5Distribution, SaSPoint5Distribution>,
        ISubtractionOperators<SaSPoint5Distribution, SaSPoint5Distribution, SaSPoint5Distribution>,
        IAdditionOperators<SaSPoint5Distribution, ddouble, SaSPoint5Distribution>,
        ISubtractionOperators<SaSPoint5Distribution, ddouble, SaSPoint5Distribution>,
        IMultiplyOperators<SaSPoint5Distribution, ddouble, SaSPoint5Distribution> {

        public override ddouble Mu { get; }

        public override ddouble C { get; }

        private readonly ddouble c_inv;

        private static readonly ddouble entropy_base = "3.6399244456803064957308496039071853510";

        public SaSPoint5Distribution() : this(mu: 0, c: 1) { }

        public SaSPoint5Distribution(ddouble mu, ddouble c) {
            ValidateLocation(mu);
            ValidateScale(c);

            Mu = mu;
            C = c;

            c_inv = 1d / c;
        }

        public override ddouble PDF(ddouble x) {
            ddouble u = (x - Mu) * c_inv;

            ddouble pdf = PDFPade.Value(u) * c_inv;

            return pdf;
        }

        public override ddouble CDF(ddouble x, Interval interval = Interval.Lower) {
            ddouble u = (x - Mu) * c_inv;

            ddouble cdf = (interval == Interval.Lower) ? CDFPade.Value(-u) : CDFPade.Value(u);

            return cdf;
        }

        public override ddouble Quantile(ddouble p, Interval interval = Interval.Lower) {
            if (!InRangeUnit(p)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                ddouble x = Mu - C * QuantilePade.Value(p);

                return x;
            }
            else {
                ddouble x = Mu + C * QuantilePade.Value(p);

                return x;
            }
        }

        public override bool Symmetric => true;

        public override ddouble Median => Mu;

        public override ddouble Mode => Mu;

        public override ddouble Mean => NaN;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => entropy_base + Log(C);

        public override ddouble Alpha => 0.5d;

        public override ddouble Beta => 0d;

        public static SaSPoint5Distribution operator +(SaSPoint5Distribution dist1, SaSPoint5Distribution dist2) {
            return new(dist1.Mu + dist2.Mu, Square(Sqrt(dist1.C) + Sqrt(dist2.C)));
        }

        public static SaSPoint5Distribution operator -(SaSPoint5Distribution dist1, SaSPoint5Distribution dist2) {
            return new(dist1.Mu - dist2.Mu, Square(Sqrt(dist1.C) + Sqrt(dist2.C)));
        }

        public static SaSPoint5Distribution operator +(SaSPoint5Distribution dist, ddouble s) {
            return new(dist.Mu + s, dist.C);
        }

        public static SaSPoint5Distribution operator -(SaSPoint5Distribution dist, ddouble s) {
            return new(dist.Mu - s, dist.C);
        }

        public static SaSPoint5Distribution operator *(SaSPoint5Distribution dist, ddouble k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public override string ToString() {
            return $"{typeof(SaSPoint5Distribution).Name}[mu={Mu},c={C}]";
        }

        private static class PDFPade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_0p015625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 8, 0xAA2293247160F5B9uL, 0xD93D9C44B8B9AA95uL), (+1, 9, 0x859FAEDB6D89FF4EuL, 0x60CA3E2F50EFB314uL)),
                ((+1, 16, 0xC14B26C52C433D9CuL, 0x9E2A35FA67F18BA0uL), (+1, 17, 0x97DEF9EC7F1A9CF1uL, 0xCF680A95D2A5DC3CuL)),
                ((+1, 24, 0x91B0D1296ED4E414uL, 0x0C618E1384E765EFuL), (+1, 24, 0xE51856DF37CABF9FuL, 0x69781EFDF858F748uL)),
                ((+1, 31, 0x9F8325B9C61A67B7uL, 0x42B010C691363213uL), (+1, 31, 0xFB1DD90BDA7057B9uL, 0xD6A31165A0863334uL)),
                ((+1, 38, 0x83C2556780ED7E86uL, 0x1346E9520E1AF414uL), (+1, 38, 0xCFCD4C209F6AA950uL, 0x663B8372EBC86D89uL)),
                ((+1, 44, 0xA79503B7574333B4uL, 0x1BA73B63F4A478F8uL), (+1, 45, 0x84879E327421EC23uL, 0x7E3AFA4176F7AAEAuL)),
                ((+1, 50, 0xA527B4C77A48C923uL, 0x46A8D2C3AC410BCBuL), (+1, 51, 0x83357E1A2D1326D8uL, 0x43BDC94C32DBC005uL)),
                ((+1, 55, 0xFBB522529D5A73B5uL, 0x26F1A7ED0A07033DuL), (+1, 56, 0xC9766B6DD223E451uL, 0xF798A774F25C5FA9uL)),
                ((+1, 61, 0x92A477D5FE8DA462uL, 0x7C287BADAA18DEEEuL), (+1, 61, 0xEDAC4D58C40EBE8FuL, 0xDCFE2A7A8DA71E85uL)),
                ((+1, 65, 0xFFCB28CF12A6C0CAuL, 0x86374C3069E63948uL), (+1, 66, 0xD3CF426EA35AC070uL, 0xA4250E2B0F054938uL)),
                ((+1, 70, 0xA11101E755EB33FEuL, 0xB2C6C230F15DE8D4uL), (+1, 71, 0x8ABB1156118FC9F4uL, 0xE12AF29DFE9F3C6EuL)),
                ((+1, 74, 0x89BCDCBD19640F11uL, 0xA5C39DF964B8121EuL), (+1, 75, 0x8023A7AEE54BCD25uL, 0x08D2846B2B4E04B5uL)),
                ((+1, 77, 0x8E4EA0D3D18495F7uL, 0x233ECD0737BBF0C6uL), (+1, 78, 0x9C57054D7148C874uL, 0xA45DB39E651E2EC6uL)),
                ((+1, 79, 0x86138244BD3B0076uL, 0xAD2C0E3C883F13AFuL), (+1, 80, 0xE1AE61890D7BFAC9uL, 0x5089370D7A40AFE8uL)),
                ((+1, 77, 0xA1B0DB06FCCDE5FCuL, 0xB13DD69528661F61uL), (+1, 82, 0x9929B6B6427D903FuL, 0x0FBE523A7FD0BBA3uL)),
                (Zero, (+1, 81, 0xB394627C1C31162DuL, 0x13A03A6D4EE10473uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p015625_0p03125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xA0B7AD1EC3EC167CuL, 0xA4E1EEDE6283264CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 7, 0xE2E9B2C2A35D291AuL, 0x03AE1E0F0853AD43uL), (+1, 8, 0xB5928B69900306D3uL, 0x5300E67ADA4362E9uL)),
                ((+1, 14, 0xEBD76EA26A172DB1uL, 0x6C04579377F1907EuL), (+1, 15, 0xBE6BCC7CE289D04AuL, 0xB6A36AB61DFB2A09uL)),
                ((+1, 21, 0x830E28936DCBD898uL, 0x9D079B231B7A5EBCuL), (+1, 21, 0xD6C3BBFE357FBE22uL, 0xFC0C6B6553EE720AuL)),
                ((+1, 26, 0xB1A377A4553821A1uL, 0xA56C89D2EB8A1D1EuL), (+1, 27, 0x950BC7902A9DBCC4uL, 0x0CDBEDCC224F4EFDuL)),
                ((+1, 31, 0x9A80537B19A0C11FuL, 0x0A85AEF90EB854BBuL), (+1, 32, 0x86C8B77A4F83E206uL, 0x5FBC4FED94518B8DuL)),
                ((+1, 35, 0xAD8453280A5E996DuL, 0xDCB8EA2FC59D2BBDuL), (+1, 36, 0xA1B8C3721712F1AFuL, 0x8D909E87EAF3F131uL)),
                ((+1, 38, 0xF3081D9DE7BE34B2uL, 0xCF94826C50693525uL), (+1, 39, 0xFEDA885B307052EDuL, 0x1484607C6DF67699uL)),
                ((+1, 41, 0xC0A93D71EB8437A0uL, 0x7B1D3FB971C6D457uL), (+1, 42, 0xFD6E88D157C4DF82uL, 0xA3782E9B58AF3DCDuL)),
                ((+1, 43, 0x86FB0FA0E74F2546uL, 0x853B93B66F0A6054uL), (+1, 45, 0x922E09B16596BDB3uL, 0xCEA767C78398A1A4uL)),
                ((+1, 41, 0xED8A017171D14243uL, 0xE911788BB8BDB609uL), (+1, 46, 0xA3B3848BCBA6B2E4uL, 0xACFDC19A3F45687DuL)),
                ((-1, 38, 0xC44F2E3C7B51A693uL, 0x92E04B0B6B851D4DuL), (+1, 45, 0xD5E77D5A1D02F4CFuL, 0xF467BDDEB626DA60uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p03125_0p0625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x9B0B6FAB0A687BF8uL, 0x120A7812DEC54DA3uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 6, 0x92D67AEB48F0565EuL, 0xC51F0B0D591E684FuL), (+1, 6, 0xF806300FAAE4361FuL, 0x3258E497DFDAC9BAuL)),
                ((+1, 11, 0xF152C2B0C2BE2926uL, 0x010BCC8A918C2BC4uL), (+1, 12, 0xD2B04B7E5272F626uL, 0x03DDEA4D7E3D352FuL)),
                ((+1, 16, 0xE01443E5F6A0A183uL, 0x021EC40545DEEEABuL), (+1, 17, 0xCD7A42698D94608FuL, 0x57E18024F935F0DEuL)),
                ((+1, 21, 0x803F2B2BE1DF682BuL, 0xF7B87008B7B5D2A1uL), (+1, 21, 0xFD356F507C4775CCuL, 0xD2649731669A046BuL)),
                ((+1, 24, 0xB91731C587DB3296uL, 0xFB619F19B78FCA9DuL), (+1, 25, 0xCC9280D20F416941uL, 0x97F6978B2C630E07uL)),
                ((+1, 27, 0xA5333577F49A6F9EuL, 0xCDD677DE1C62515CuL), (+1, 28, 0xD9E279B42877AC71uL, 0xC5E3C6045D7AD459uL)),
                ((+1, 29, 0xAB6908423CB53A5DuL, 0x1E8D19E4ABECD5E0uL), (+1, 31, 0x96073A9F3E332005uL, 0xAE657FE73AB47796uL)),
                ((+1, 30, 0xB6C1BAFF6630818BuL, 0x51686CF54D9C65FCuL), (+1, 32, 0xFEEF5B6045634DCEuL, 0xBA3AAC43CF9DBD1EuL)),
                ((+1, 30, 0x9A2B8ECEE9F1390AuL, 0xAE9188F5AF09C8D0uL), (+1, 33, 0xF49DC017520F75A2uL, 0xC506CA79ADDCCE66uL)),
                ((+1, 27, 0xD8469ABA074B43D4uL, 0x5D43594055B781F5uL), (+1, 33, 0xDF8318AB3073AFAAuL, 0x5E6333C913B2471AuL)),
                ((-1, 23, 0xC1F218B63899E7A4uL, 0x065880FDCC031201uL), (+1, 32, 0x86CFAB57BDDF3DFFuL, 0xC223CCDABBDFB080uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p0625_0p125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8BE1F74611E6BB05uL, 0x9296FF5C38B30916uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 5, 0x9736C4243FB94247uL, 0x02CAA583651E15B2uL), (+1, 6, 0x918D799A9645C95DuL, 0xD7E98FED0CA9140EuL)),
                ((+1, 10, 0x8BAD170114B246BAuL, 0x9925F52204AEA943uL), (+1, 11, 0x9009F21D425DD728uL, 0xC154223DFF70F8CBuL)),
                ((+1, 14, 0x8F2953B1DCEEBBFAuL, 0xFEA7BEE60DB54F83uL), (+1, 15, 0xA25CC53405D9A829uL, 0x20B19422F24CA10FuL)),
                ((+1, 17, 0xB0C505F107EC51BEuL, 0xB88CD37D2D16F74EuL), (+1, 18, 0xE51B08C94BA6480AuL, 0xE1AD7C0387F980D3uL)),
                ((+1, 20, 0x85610C2E629D424BuL, 0x4C87F44A6707A4D9uL), (+1, 21, 0xD15AE86D1531BA6FuL, 0x4CE5AEDA91B154E1uL)),
                ((+1, 21, 0xEDB24A4BDAFB20C2uL, 0xF4023C1BB0743F55uL), (+1, 23, 0xF79D4380BB0B983BuL, 0x91428523A1CDF852uL)),
                ((+1, 22, 0xE36E1C536D9038FBuL, 0x8F45B44A6F10B440uL), (+1, 25, 0xB7AA262DBBC9639DuL, 0x3B462606C5AF5EE1uL)),
                ((+1, 22, 0xBC849C684D67693AuL, 0x5206F7EF78BAB89BuL), (+1, 26, 0x9F051A59CB93640BuL, 0x3A6220A97B90F231uL)),
                ((+1, 20, 0x9D92B09570A637C2uL, 0x4D55F29D18DCFF49uL), (+1, 26, 0x8B229B59608F4AE6uL, 0x4B34454490FBE548uL)),
                ((-1, 16, 0xE5C78C6A5D47C84BuL, 0x6682487A794709A5uL), (+1, 24, 0xB281ACC9D4E68493uL, 0xA4B42CDD5BAD128AuL)),
                ((+1, 14, 0x9DCC4848E09BA44EuL, 0xCC71210CFF71CC03uL), Zero),
                ((-1, 11, 0xAB536C8AC3DC39B5uL, 0x018FB1065C4C23BDuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xDF0FEDC2FE5FDF90uL, 0x51101334A0A6DC0CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 3, 0xFEEDEF0AEF493EC7uL, 0x317013E6F88E15C7uL), (+1, 5, 0xA0523EBEDE3C961AuL, 0x980765798CA82744uL)),
                ((+1, 7, 0xF74B12EA23AC78B9uL, 0x15DD402524513BC6uL), (+1, 9, 0xAEC2350966ACA57AuL, 0x8D918CA6EAD0C66DuL)),
                ((+1, 11, 0x8455162060F7F9DDuL, 0x95679909C1E1FE74uL), (+1, 12, 0xD9483C9C10898F2BuL, 0x99C61C448C165A26uL)),
                ((+1, 13, 0xAA1BB5CFA04A5797uL, 0xC57A2BD014D980BFuL), (+1, 15, 0xA9AD2274EA7A56B8uL, 0x6C12FD9F9438E3D1uL)),
                ((+1, 15, 0x861F7BBDF750C3E9uL, 0x27C11CA038E79FF4uL), (+1, 17, 0xACD1965979A5BC51uL, 0x83137BFE8B5DCB0FuL)),
                ((+1, 15, 0xFED54A8ABA45A24EuL, 0x3D31A456DAAF9383uL), (+1, 18, 0xE6F2407FB076F541uL, 0x8CB03F397E100A13uL)),
                ((+1, 16, 0x8A3F28C100C1F730uL, 0x04C39EE615D52CC8uL), (+1, 19, 0xC6CE600A18DF386CuL, 0xA58D89885BDB8C54uL)),
                ((+1, 15, 0x9B7CC662B2B20F39uL, 0x7F0891592C2D6974uL), (+1, 19, 0xD32C5280A7C8E427uL, 0x3960C1540BB446A1uL)),
                ((+1, 13, 0x9568EA71604C23ADuL, 0x203DA8ED79E18EB0uL), (+1, 19, 0x80310868D6FA4588uL, 0x8DC8B9AFC1196573uL)),
                ((+1, 9, 0x88166FFE7C81DBE1uL, 0xD73EC17FFC043212uL), (+1, 17, 0x9A943FB12EECB170uL, 0x099EB6ED5DF3CFA7uL)),
                ((-1, 3, 0xF19492A58CEDD9BDuL, 0xE7A1E3235CDE26BFuL), (+1, 14, 0x856C7B32EE0F486CuL, 0xCE73B21020AECDC4uL)),
                ((+1, -1, 0x98410AF4C395C5CAuL, 0xB294E615E5607A03uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p25_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x975ED700D432D21FuL, 0x0D72B717C1E0DF48uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0xB1C63271B48D22E7uL, 0x62B135CE5A574D79uL), (+1, 4, 0xAC13D1453FF2F3BCuL, 0x82E604CD7E9216BAuL)),
                ((+1, 5, 0xB0C0BBFD76DEA94BuL, 0x028878F3C2BA8115uL), (+1, 7, 0xC9877A78BD26FAF6uL, 0xC4F3E4ED71A5E22DuL)),
                ((+1, 7, 0xC18BEB3FFA6506F1uL, 0x6AD0D802E6E8ECADuL), (+1, 10, 0x86CC449D9BF426BEuL, 0xFCF8B083810A4E0EuL)),
                ((+1, 8, 0xFE7F9F93408DE3C2uL, 0x79C77A2F45217369uL), (+1, 11, 0xE2FB7C600F914864uL, 0x0FBD650052BC24B1uL)),
                ((+1, 9, 0xCDCEC8D2468CD377uL, 0xA9B37354666F25C3uL), (+1, 12, 0xFA02655F6BAD54A8uL, 0xECDD0ABA28970937uL)),
                ((+1, 9, 0xCA1AD6AA1DC92BDFuL, 0xEA4A70C3EAD32EECuL), (+1, 13, 0xB57C68EB703E818DuL, 0x318B7F685319C2C4uL)),
                ((+1, 8, 0xE63C3596153BC608uL, 0x1D57DD10F9A93F24uL), (+1, 13, 0xAAF20AC1F9D5F9A3uL, 0x123212CF257F7A0BuL)),
                ((+1, 7, 0x8AA9D410B485D445uL, 0xE3CF1E23FA44BEC8uL), (+1, 12, 0xC8CA3473A902DC1EuL, 0xC8FD7163F23F3454uL)),
                ((+1, 4, 0x915CBFFCE7DB47E6uL, 0xB5989364EF3BFBBBuL), (+1, 11, 0x8894DE4F57D4073FuL, 0xBE2F7B5DBBF7E330uL)),
                ((+1, -1, 0x912C33DED4A2C59DuL, 0x93E21C00CE5269BAuL), (+1, 8, 0xBB2E2B7725356282uL, 0x9125ABABC397074CuL)),
                ((-1, -7, 0x96DD906605286043uL, 0x3F7DB59D2327BE7DuL), (+1, 4, 0xBA0F5E33F85CFFE6uL, 0x0EDE811F23ADCA8BuL)),
                ((+1, -13, 0xDE5840F7473E46E7uL, 0x95A649F3C8E55CC5uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p5_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xAEDC56CB2F486C57uL, 0x90E8586797D59176uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xD029A64B82676648uL, 0x55E357D9DFA3D287uL), (+1, 3, 0xB518B0A593477DBDuL, 0xBE79EB44BC444281uL)),
                ((+1, 2, 0xD1A1F8706AE5456EuL, 0x60AA15CFED21CC38uL), (+1, 5, 0xDF5CDD8BC1C6F93CuL, 0xF5D00BC2517FEFCAuL)),
                ((+1, 3, 0xE88DF0CD89D6A111uL, 0xBFBBA7D456F44A03uL), (+1, 7, 0x9D8354A7628D2838uL, 0x86FCA72051E3CF98uL)),
                ((+1, 4, 0x9B2DF161EE9A93D9uL, 0xCAD53D57E91B0FCBuL), (+1, 8, 0x8C0DDEABE0887A72uL, 0xDF6DD450DBE288F3uL)),
                ((+1, 3, 0xFFD62883530785E1uL, 0xD6D67BC033F3271BuL), (+1, 8, 0xA35734565AC68BA2uL, 0xEE3B92DBCF6AB672uL)),
                ((+1, 3, 0x810178A82ED87110uL, 0x5BE89CE5923058B7uL), (+1, 7, 0xFC10B6CF696EC7A9uL, 0xE30241503DCF12A2uL)),
                ((+1, 1, 0x98639E30CA628CFFuL, 0x0FF676040E25FA45uL), (+1, 6, 0xFDB36B5CC8C7079CuL, 0x1C7D6645D0EC2776uL)),
                ((+1, -2, 0xC04688220F411F62uL, 0xDE83914EA0581D6BuL), (+1, 5, 0xA04119007F58F726uL, 0x4B896FBF67F9A782uL)),
                ((+1, -6, 0xD4F02FBFDC2AB943uL, 0x460127D9189959B9uL), (+1, 2, 0xEC33D3F314CDC687uL, 0x92AA79E0682846A3uL)),
                ((+1, -12, 0xE0446085CB19C28CuL, 0x4E972111E490219AuL), (+1, -1, 0xB0A6B214D8176D43uL, 0x2E0D08F674192F08uL)),
                ((-1, -18, 0x813E4B1D6ACB4F18uL, 0xC8266250F7C2CCA6uL), (+1, -6, 0xC0F64C30ECB6F219uL, 0x52D857FB1C9E9EB8uL)),
                ((+1, -25, 0xD310446A8B41E8A7uL, 0x2E3F1D8D63A76A47uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xB058F19F88323B2FuL, 0x8BBA73552E0593DFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xD417AF0AA438EE74uL, 0xBEE68EBA33FC6969uL), (+1, 2, 0xBC30A21CA457F30EuL, 0xAAD03F0CE0F0BA20uL)),
                ((+1, -1, 0xD7C637FED71CC77CuL, 0xD4C5AB538E243DBDuL), (+1, 3, 0xF14F844C8655BB85uL, 0x57B12B6E42DE21D8uL)),
                ((+1, -1, 0xF243282B23410A67uL, 0x5CD46FA2E9D87EE2uL), (+1, 4, 0xB12594A4109217E9uL, 0xC63C0FDC806D6C68uL)),
                ((+1, -1, 0xA462B7726F0C862EuL, 0xD7DF4F55F8FD14A9uL), (+1, 4, 0xA4644FD02D1A87CBuL, 0xBC2C8C1CE5CD28F9uL)),
                ((+1, -2, 0x8B02E4A875D51E8CuL, 0x477678224DAA880AuL), (+1, 3, 0xC8F9A4683D911775uL, 0xD606A0F48850E108uL)),
                ((+1, -4, 0x91EEC58B3140D53EuL, 0x8BCA507E07C6B3D8uL), (+1, 2, 0xA3B178EDC37D2DA0uL, 0xA09D2593797FC41CuL)),
                ((+1, -7, 0xB7DD982B6B92E275uL, 0x301520986C74AB41uL), (+1, 0, 0xAFD8580FCBF90782uL, 0x00BD773F0CCED948uL)),
                ((+1, -10, 0x8177435F358A5076uL, 0x3456E95A7BBDBD2EuL), (+1, -3, 0xF169981876766CC4uL, 0x56F4AF69A84B836BuL)),
                ((+1, -15, 0xB221682A08FCCF1FuL, 0x81793964ACBF5CC6uL), (+1, -6, 0xC7AA2DF66973B32CuL, 0x7A3D9252625AD0F3uL)),
                ((+1, -21, 0xB4D7A5C00C2708ADuL, 0x4E08BAB9DC28108EuL), (+1, -10, 0xB37003DC2367E1C4uL, 0x0682AB2B82E9565AuL)),
                ((+1, -30, 0xE3B80A8D30F4607AuL, 0x78F814F9C7924946uL), (+1, -15, 0x9056C87E7507EC79uL, 0x635932F0678E95A4uL)),
                (Zero, (+1, -22, 0x83D32F27289FBB51uL, 0x9841D53434EC9AD0uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xA05442F31FA20206uL, 0xCBE654E3921531FCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0xB77B66A77DBDA7C2uL, 0x0FF1C7284272BD88uL), (+1, 1, 0xB8CF8198DA9522DBuL, 0x20437D991D4CBA86uL)),
                ((+1, -4, 0xB14E9242FE38F5F3uL, 0xE781BCF7F411BEF6uL), (+1, 1, 0xE86976020422E1ABuL, 0x98B45366B76E7EEDuL)),
                ((+1, -5, 0xBCB7EF8B200610FBuL, 0xEA55B08AFF8B8C3CuL), (+1, 1, 0xA70BF9696288F513uL, 0x140A973E9C36C892uL)),
                ((+1, -7, 0xF21BD78DAE9231E8uL, 0xDD84B2F107B985EEuL), (+1, 0, 0x976EE7CCFB6D2DA6uL, 0xFC153B27C7EEF766uL)),
                ((+1, -9, 0xC099870A7B0D36E6uL, 0xA31A68AA88684FA1uL), (+1, -2, 0xB44076AC247AFD36uL, 0x716477DC18845875uL)),
                ((+1, -12, 0xBC703E992F1AE283uL, 0xE4185FF06A8DD1BBuL), (+1, -4, 0x8E349C6C5DF4A9DEuL, 0x430DACE6486FB3E8uL)),
                ((+1, -16, 0xD936BEDE5BDC072CuL, 0xA8248250ABAE7F4EuL), (+1, -7, 0x92B42955A378147FuL, 0x193F1F17321FDB07uL)),
                ((+1, -20, 0x86740FC88A548644uL, 0x614F4112D091D6C8uL), (+1, -11, 0xBE839604C53E9C87uL, 0x90C45B7065637729uL)),
                ((+1, -26, 0x929835A862276003uL, 0x6DD987EAE3EBD1D4uL), (+1, -15, 0x90C325C28C1659D9uL, 0x9AFCBD62337CE5BEuL)),
                ((+1, -34, 0x96698E92E08969D2uL, 0xB3C0AEC1A72DCBD4uL), (+1, -21, 0xDFE6C75B48D9E847uL, 0x886B9D973E7423DEuL)),
                ((-1, -43, 0xB450545C7402C344uL, 0xC3701DD3910238F3uL), (+1, -28, 0xFD88811AC109B813uL, 0x4AAD5B53D0229A17uL)),
                ((+1, -51, 0x98B4D388962D865CuL, 0x99B6427C1C2DCC3BuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0x87370AD70505A575uL, 0x1B60C2A74CE6DD77uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -6, 0x93BF73E9D752C7CBuL, 0xA49BEA9BDB70B9DAuL), (+1, 0, 0xB51939E5BC81EA13uL, 0x6300AC9A5F79EAAFuL)),
                ((+1, -7, 0x88855F60EB405EE2uL, 0x33362E1EE2761B5CuL), (+1, -1, 0xDF32E5D4451602ADuL, 0xBC1741FFD26D0AF6uL)),
                ((+1, -9, 0x8B338ACA2B08AC60uL, 0xFE70C41C26C52A89uL), (+1, -2, 0x9D438F7E0ACB1E7DuL, 0x480B27F70C8AE4A9uL)),
                ((+1, -12, 0xAB6EABF75C870148uL, 0xC1305EFCB7FCF87BuL), (+1, -4, 0x8BCECB29DC6FA67FuL, 0xCB1248BA834258E9uL)),
                ((+1, -15, 0x83321CEA3ED85660uL, 0xE68B422A13F4D101uL), (+1, -7, 0xA345A90B126D1F7CuL, 0xAA5168B3E8F0D14EuL)),
                ((+1, -20, 0xF77E9457E278FA58uL, 0x6C40127F2CDE057BuL), (+1, -11, 0xFCE48F110E9AEDB4uL, 0x0A30996FA537471BuL)),
                ((+1, -24, 0x89CB7E40FC197608uL, 0x4E404D2A75027828uL), (+1, -14, 0x80212369748294FBuL, 0x1B9299A1CC16061AuL)),
                ((+1, -30, 0xA516403D22EEE197uL, 0xFFF3FD5F14A72495uL), (+1, -19, 0xA38AD4B897FE4A9BuL, 0xAD30BA5A02486DBCuL)),
                ((+1, -37, 0xAE66B5C96058C515uL, 0xD2166D9A0F35360BuL), (+1, -25, 0xF4710BB8680F2F38uL, 0x74D3FB131F00844BuL)),
                ((+1, -46, 0xACDACECE240D6ABCuL, 0x070B92E376F49594uL), (+1, -31, 0xBA0D51A205F99FF6uL, 0x8A024BDD3C8FD000uL)),
                ((-1, -56, 0xCC56BB2D4AAC0173uL, 0xA7E7838658088A3EuL), (+1, -39, 0xCF7F29C24E3AA963uL, 0x8C9E79386AA89CE1uL)),
                ((+1, -65, 0xAAB7F3FE26FB2F10uL, 0x3A92613FB737C688uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -8, 0xD848925FF992312FuL, 0x61432B68EE996EB5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -9, 0xF47960858081E5B5uL, 0x55DD336C91D604C0uL), (+1, -1, 0xBBF4DF8BD1E49792uL, 0x72A7EFDECB4E239FuL)),
                ((+1, -11, 0xEABBC20E54478A32uL, 0xBCE2BF864FA0A087uL), (+1, -3, 0xF12AB5742E84CFEEuL, 0xCDB7C848AEA1DB34uL)),
                ((+1, -14, 0xFA036335F5C9C722uL, 0x58D8089AD5BF8211uL), (+1, -5, 0xB18C624FED8EB21BuL, 0x937A34A9600E4616uL)),
                ((+1, -17, 0xA1E1553A07D6D8E9uL, 0x51FC6F50F3D4EF47uL), (+1, -8, 0xA5A7CB2CAF5251DCuL, 0xF6B22EC48BB6F886uL)),
                ((+1, -21, 0x8371FB4843BBAAA8uL, 0x7249BFFEFC14CCC3uL), (+1, -12, 0xCC2FC311A879807CuL, 0x394CEB6726E57416uL)),
                ((+1, -26, 0x8555918667C8685CuL, 0x5834487D9FB738B3uL), (+1, -16, 0xA82DFFB941D9847AuL, 0x2F4FAE53725B2E4CuL)),
                ((+1, -32, 0xA35012C31E4B2CB4uL, 0x5A8CEA0473F982C1uL), (+1, -21, 0xB74220C5AE8C38E8uL, 0x6F67C8586C445059uL)),
                ((+1, -39, 0xE0DAA8B97C82E55BuL, 0x9CFD6CDD9C159952uL), (+1, -27, 0xFFF98C8308F81A08uL, 0xDBE9F60C40847820uL)),
                ((+1, -46, 0x980481A81DC53270uL, 0x273580FAECABBE28uL), (+1, -33, 0xD8039867BA4BDB40uL, 0xB375F29FF6EDA546uL)),
                ((+1, -55, 0x985DCF5659C61AABuL, 0x27F754C2C58D579FuL), (+1, -40, 0xC695EFB9702FA764uL, 0x903103C5FEE2C302uL)),
                ((+1, -67, 0xBE2C7B7BD8C1CEFCuL, 0x4D221DEB3DA5C42AuL), (+1, -48, 0xA3C26318B07F8F1BuL, 0x983828E21F0C6844uL)),
                (Zero, (+1, -58, 0x99935B9AD74E8EA3uL, 0xCBE15699E13113D5uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -9, 0xA6AF14B445823CF6uL, 0x15D3D95051B5ECD6uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -11, 0xB9F6677322C35FF8uL, 0x7206E65F02B1BC0FuL), (+1, -2, 0xBB7F6F050AEB93F0uL, 0x217D3A9EC5180345uL)),
                ((+1, -14, 0xB065046C599FDBDEuL, 0xD6353D6707799BEAuL), (+1, -5, 0xF00A184F991122F1uL, 0x5F25942CF52AB5F6uL)),
                ((+1, -18, 0xB9C3FBA843F9A3CFuL, 0x5D0620AAA8FD0724uL), (+1, -8, 0xB05C4CE1C5CD874FuL, 0xA81B25E062C57B75uL)),
                ((+1, -23, 0xEE0D632B6EA1CF8AuL, 0x591CDF4200626FE4uL), (+1, -12, 0xA44048625565EE0EuL, 0x6B3CF4F7F45F015AuL)),
                ((+1, -28, 0xBF6DBB82630CE56EuL, 0x4A6840C1B3920385uL), (+1, -17, 0xCA237B506AD2A1DFuL, 0xC803EE3EBE333F9BuL)),
                ((+1, -34, 0xC0721E725A0F152AuL, 0x50EAA45DB167C1F0uL), (+1, -22, 0xA6459EB87EF9877BuL, 0x8AB1C7FD8C8F5CFAuL)),
                ((+1, -41, 0xE9C66D8D80CA3F74uL, 0xD9A21AD33F8F787DuL), (+1, -28, 0xB4FB98594CD513ACuL, 0x40BAB185C5EA7A8AuL)),
                ((+1, -48, 0x9FB7D44F137529C2uL, 0xC01688EF16000840uL), (+1, -35, 0xFC9514F74A4C2F1FuL, 0x7D3FB65895C6C0ABuL)),
                ((+1, -57, 0xD677F649A97AA0EAuL, 0x67512A52C3E93BFFuL), (+1, -42, 0xD5054522AB0943EBuL, 0x048200F1EE5D64FAuL)),
                ((+1, -67, 0xD59B6B7F99BBF65EuL, 0x85D40CB13AD7E5F6uL), (+1, -50, 0xC3C1F0C7F3517306uL, 0x3C66D876E72F36B2uL)),
                ((+1, -79, 0x848B684A3E20DB28uL, 0x80C82E7ADA901D91uL), (+1, -59, 0xA165ED0EB507DD38uL, 0x207594B8A3AFCE09uL)),
                (Zero, (+1, -70, 0x975CDAC1CA62304BuL, 0x2EF1D70BD6DA5AEAuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -11, 0xFA5EB74CFFA725C5uL, 0x465F2E2CFC67A632uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -13, 0x8B9B2AD171E362F4uL, 0x6243D3A7796BBB8EuL), (+1, -3, 0xBC6D956E1943ED46uL, 0xB27D36A258F082E8uL)),
                ((+1, -17, 0x845DF510CA215649uL, 0xA39B221B9AD5AA17uL), (+1, -7, 0xF269DF3C1716084CuL, 0x8C318906FA394D14uL)),
                ((+1, -22, 0x8B56A488F6AAC1C5uL, 0x2315DD4CD8316916uL), (+1, -11, 0xB2F79EB51BE8E286uL, 0x0A085B0BBD5AE5AEuL)),
                ((+1, -28, 0xB27A849147CFE36CuL, 0xF3C948247755CD27uL), (+1, -16, 0xA7796E8E5A9F8DAAuL, 0x8E573F9AAE712A62uL)),
                ((+1, -34, 0x8F756279450C37D0uL, 0xEC6A71A1A5EC1933uL), (+1, -22, 0xCF12B138B4AFBAB6uL, 0x58894A23DDC055F2uL)),
                ((+1, -41, 0x9027EB2AF46B9D25uL, 0x31B2B9D9D3D8F8B8uL), (+1, -28, 0xAB1D9C6C3CCDD646uL, 0x525D5ACCEACC93E8uL)),
                ((+1, -49, 0xAF09366B0DC6CF13uL, 0x31BB438A98B1F438uL), (+1, -35, 0xBB1848F5DF92B39AuL, 0x0475725D0E5C96DBuL)),
                ((+1, -58, 0xEF113054519FA7B1uL, 0x4AE636F49C0422C7uL), (+1, -42, 0x83215642C8E2FD3BuL, 0x4445D35906D99EA7uL)),
                ((+1, -67, 0xA0705C2D36AC27CAuL, 0x639DAE08A2478F65uL), (+1, -51, 0xDE207D8B429F24E2uL, 0x069BBFAB5081A3CFuL)),
                ((+1, -78, 0x9FB9FDC2C2F45629uL, 0x8291F7ACF69E63BBuL), (+1, -60, 0xCCF7A940F819FF5BuL, 0x2979B3ABB0CEAF71uL)),
                ((+1, -92, 0xC624406ABA0548F9uL, 0xCE96F474E915D577uL), (+1, -70, 0xA9A8FB291B4AAA36uL, 0x12E5BEF65BC2C87AuL)),
                (Zero, (+1, -82, 0x9FB5EF08AF969F3CuL, 0x7234E2BCB9B15A2CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_limit = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xCC42299EA1B28468uL, 0x7E59E2805D5C717FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x8806D74957711C79uL, 0x71FD3FC57AF6FD94uL), (+1, 0, 0x90C0122D9C0C55FFuL, 0x4FE584D94D6C3F09uL)),
                ((-1, -7, 0xA72FD1F9FF29E20BuL, 0xCA5774A3FB5D481AuL), (+1, -1, 0x99E45528733CFE1DuL, 0x05225643DA28BA61uL)),
                ((+1, -10, 0xA9FD92F11175ED5AuL, 0x28BC3336F6853873uL), (+1, -3, 0xD04F2C27D0125561uL, 0xAC911E4B4CADABEAuL)),
                ((+1, -11, 0xB8B4982D37226594uL, 0x574F7C39FA319D7DuL), (+1, -5, 0xBFBA486FB2D45706uL, 0xDEECBF2EDDA292CEuL)),
                ((-1, -13, 0xCF3F2BFE77105A4CuL, 0xE7C634F52EE29DDFuL), (+1, -8, 0xF71EF53E8F4F2528uL, 0x3E936E75EF2BCB01uL)),
                ((+1, -16, 0xCEC560A761793618uL, 0x83C390E43295DE95uL), (+1, -11, 0xCF27C7EA9363B65FuL, 0x1672BFC4CCFCEE07uL)),
                ((-1, -20, 0xACDDF652457D7FE5uL, 0x1952A91326D8DFDDuL), (+1, -15, 0xC3EDA1B9D6840B30uL, 0x72190271702F7FC2uL)),
            }));

            public static ddouble Value(ddouble x) {
                x = Abs(x);

                ddouble y;
                if (x <= 0.015625d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_0p015625);
                }
                else if (x <= 0.03125d) {
                    y = ApproxUtil.Pade(x - 0.015625d, pade_plus_0p015625_0p03125);
                }
                else if (x <= 0.0625d) {
                    y = ApproxUtil.Pade(x - 0.03125d, pade_plus_0p03125_0p0625);
                }
                else if (x <= 0.125d) {
                    y = ApproxUtil.Pade(x - 0.0625d, pade_plus_0p0625_0p125);
                }
                else if (x <= 0.25d) {
                    y = ApproxUtil.Pade(x - 0.125d, pade_plus_0p125_0p25);
                }
                else if (x <= 0.5d) {
                    y = ApproxUtil.Pade(x - 0.25d, pade_plus_0p25_0p5);
                }
                else if (x <= 1d) {
                    y = ApproxUtil.Pade(x - 0.5d, pade_plus_0p5_1);
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
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    ddouble v = Sqrt(x);
                    ddouble u = 1d / v;

                    y = ApproxUtil.Pade(u, pade_plus_limit) * (u * u * u);
                }

                return y;
            }
        }

        private static class CDFPade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_0p015625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8000000000000000uL, 0x0000000000000000uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 7, 0xBD6968E010CE2AAAuL, 0xC5892B75119334CFuL), (+1, 8, 0xBE0C62637F1C6EBFuL, 0xEF8552CC53FC3EFBuL)),
                ((+1, 15, 0x984F52DD3C63CF4DuL, 0xA30DF2EFFBE9EF78uL), (+1, 16, 0x99414D00FE29C792uL, 0xAB050EE97879D46EuL)),
                ((+1, 22, 0x9EB26449C2469A9EuL, 0xB0C76090FB59A5F4uL), (+1, 23, 0xA0388D9A79A5D9FFuL, 0xA08A537405EDC7EEuL)),
                ((+1, 28, 0xE9F494E5935413BDuL, 0x1F658C4FAD8FD1D5uL), (+1, 29, 0xED23FD4038AAC814uL, 0x68914059A8FDD417uL)),
                ((+1, 34, 0xFB34B0DDBFC37549uL, 0x8C5B976719A1C035uL), (+1, 35, 0xFFEA893B070DB2CBuL, 0xC15C5BA4AD0A1B4CuL)),
                ((+1, 40, 0xC7032119C2BB1606uL, 0x51085F33494415A2uL), (+1, 41, 0xCC16897A44549C81uL, 0x3B780A65929F9315uL)),
                ((+1, 45, 0xE766A371FA6EDC3BuL, 0xEE28B98074BEE6AAuL), (+1, 46, 0xEF79CEC2EADD3D47uL, 0x3591E830F817A71BuL)),
                ((+1, 50, 0xC2305968C61A3409uL, 0xC2FE4C3B358A99A6uL), (+1, 51, 0xCB9EC0C4BFEBFDC3uL, 0x1D2A2801B5C6349EuL)),
                ((+1, 54, 0xE36840B59E534974uL, 0x306962A984375341uL), (+1, 55, 0xF34E892C3AD63F81uL, 0x91B1D06DF1278AB2uL)),
                ((+1, 58, 0xB01ECEB2B8AECD34uL, 0x3D334D5C74E27A8BuL), (+1, 59, 0xC2CAF786124677D6uL, 0x8ABE817C3347982DuL)),
                ((+1, 61, 0xA4C1BF0C1905F65EuL, 0x6179A1F16BA3A37BuL), (+1, 62, 0xC18FAA958F61E8D6uL, 0xC14CB86D39AA30D1uL)),
                ((+1, 63, 0x9FC72F2BEB50917FuL, 0x1B92C74F19475BB8uL), (+1, 64, 0xD424BD356045E02BuL, 0x03813356CCFFA297uL)),
                ((+1, 63, 0xE9223B4B38ACF2D4uL, 0x33242711B4E6DB0DuL), (+1, 65, 0xCEB475D4B95231BFuL, 0xBE41E392BB7FBAE3uL)),
                ((+1, 61, 0xCEB6C238F706B7BBuL, 0xF7223203F1161BD4uL), (+1, 64, 0xDEEE3BAE40447074uL, 0x9FD5452F6AD191D6uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p015625_0p03125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xFAEE5A93CEA64CCDuL, 0x57959C1D018AAB98uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 6, 0xAA2FBE53FA629985uL, 0x9F1170D52129C1ACuL), (+1, 7, 0xAEE7C2E5F30ADFC2uL, 0x39400C9B82AFCCDDuL)),
                ((+1, 12, 0xC926D7DA5D132261uL, 0xC902173A5178BFF3uL), (+1, 13, 0xD0B2E3A4F03BF539uL, 0x1FB81EB24DA94052uL)),
                ((+1, 18, 0x86A3B521C33B23ABuL, 0xB31EB3AC011D8E75uL), (+1, 19, 0x8D7C35A69D5D287AuL, 0x428CE95DA4053019uL)),
                ((+1, 22, 0xDE7C5444B9B76AEAuL, 0x178E22737EC020B3uL), (+1, 23, 0xEE09986330A65A42uL, 0xAC0CF0990E4F317BuL)),
                ((+1, 26, 0xE8021E39D9FA3054uL, 0x1C427B17985E902AuL), (+1, 27, 0xFEF09BF732BB3ACAuL, 0xB1C0527B5CF14F8AuL)),
                ((+1, 30, 0x95CA5F1AC59A9419uL, 0xC3F9BC153CEBC8D6uL), (+1, 31, 0xABB4AF82B19085CCuL, 0xE2859378FC8EE595uL)),
                ((+1, 32, 0xE0FDBAC2F0CC6C8CuL, 0xD9845DE27E2E7920uL), (+1, 34, 0x8AA88818697CA59FuL, 0x2B098328C2F3C1ECuL)),
                ((+1, 34, 0xAD3E797FEFA4510CuL, 0x0FC572DCBC89782BuL), (+1, 35, 0xF4AEA19956642E64uL, 0x2D1EA92F2B712C04uL)),
                ((+1, 34, 0xD8CD8450B32C87BDuL, 0x79B69F1770E31E61uL), (+1, 36, 0xC8DCE644BDB850E6uL, 0xD42B2FD3A451750FuL)),
                ((+1, 32, 0xE4804B8009D4321EuL, 0xD46DDD6B25215A4EuL), (+1, 35, 0xD662330D4FFE5776uL, 0x4CDB87321AE0EF99uL)),
                ((-1, 29, 0xB0BC5C7D3BE08FF4uL, 0xCF2F169EF98C26DEuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p03125_0p0625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xF5FDAA1B4C8FCCDFuL, 0x8C4C4466C759E96BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 5, 0xD50B619C94B7ADD7uL, 0x1039219DAFED39A2uL), (+1, 6, 0xE03BF288B7324D60uL, 0x3C297C3CE11BDDAEuL)),
                ((+1, 11, 0x9D7F5FD45750D376uL, 0x15494F70135B583DuL), (+1, 12, 0xA844873A13D7ABE8uL, 0xED41B7CEE921CDFAuL)),
                ((+1, 16, 0x81418D9CBA64284AuL, 0x00AE1AB41DACB114uL), (+1, 17, 0x8CF163A5C9FB103FuL, 0x0719BC959919701AuL)),
                ((+1, 20, 0x803292852DABE34FuL, 0x31339BBA00F67DFDuL), (+1, 21, 0x8FE36C7796A170E6uL, 0x945520A208D7B04DuL)),
                ((+1, 23, 0x9CC44FEAAC7A1B56uL, 0xC3F4F583718E573AuL), (+1, 24, 0xB7A0D5ED217D46A7uL, 0x60F32040C16E0F65uL)),
                ((+1, 25, 0xE758F3D6E5AE55CCuL, 0xEB7FA0450DE45A17uL), (+1, 27, 0x90BD445953EC5EDFuL, 0x25F54852753ED988uL)),
                ((+1, 27, 0xC1BC9FF43365EA29uL, 0xD41E807E4D139DBFuL), (+1, 29, 0x86D8BAFF7CDD7B77uL, 0xC0F20E6D368AC13AuL)),
                ((+1, 28, 0xA492DFAFE072FEC6uL, 0xF34F35A7E33731B9uL), (+1, 30, 0x88CC81F756010435uL, 0x810832981E95FD1AuL)),
                ((+1, 27, 0xE310C8CC00EF0D93uL, 0x2C1CB8FC2C2B9692uL), (+1, 30, 0x81AF361833388C4FuL, 0xB990D45A51150843uL)),
                ((+1, 25, 0x830EF65D1DF3B97BuL, 0x2EF82BC880B48CE6uL), (+1, 28, 0xA0B0B72A7011A24EuL, 0x9FFDB56E0A301644uL)),
                ((-1, 20, 0xEA6A3C107DACFFE4uL, 0xB21E52F6870F4383uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p0625_0p125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xECC384415B68DD0DuL, 0x456E28B34CB53B01uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 4, 0xE75C61024732B637uL, 0xB9F11EF3F6281263uL), (+1, 5, 0xFEE277EBB0AF22B4uL, 0x28609DFAEE52AEFCuL)),
                ((+1, 9, 0xC10B0472343287A6uL, 0xAE770336DBABD88FuL), (+1, 10, 0xD9DFA842B723E672uL, 0xFD5BC23033CC1D1FuL)),
                ((+1, 13, 0xB316CEB7607E3417uL, 0xB01C6C5B62B294B6uL), (+1, 14, 0xD0ADC0D739661604uL, 0xFF02C7437CD9F504uL)),
                ((+1, 16, 0xC9BEDAC9B8D64DA6uL, 0x8A18AA90271555EDuL), (+1, 17, 0xF56CCECC36B5C503uL, 0xD9F138046A39F69AuL)),
                ((+1, 19, 0x8DB738FFB33DBE70uL, 0x5811726B102B878DuL), (+1, 20, 0xB6F13F510D3DF1BDuL, 0xA417ABF5680E5FF2uL)),
                ((+1, 20, 0xF64029511567CBCAuL, 0xE8E0447714E838FFuL), (+1, 22, 0xACBEADCC6EE66FECuL, 0xBA4AF8C6E7D2F8CAuL)),
                ((+1, 21, 0xFEEE34EAF9C068BBuL, 0x65DB1F982F63EB8BuL), (+1, 23, 0xC9510A7AFB66FD75uL, 0x0BABB38ED74FDE1EuL)),
                ((+1, 22, 0x9221DE8F9A2E57F9uL, 0x029FCD8F4DCE1F02uL), (+1, 24, 0x890D758DAF2B5B5CuL, 0xB6653FF230865428uL)),
                ((+1, 21, 0xA20C08352AEE5C89uL, 0x5D3CA67476D601EEuL), (+1, 23, 0xC55650D50243D2A1uL, 0xA8B8F606BA79D779uL)),
                ((+1, 19, 0x8346859DE75147B7uL, 0xF4E838BF885C6DF0uL), (+1, 21, 0xF862B266042E1570uL, 0x67EF22869A558920uL)),
                ((+1, 14, 0x82A4F9FAC7000AFEuL, 0x1D889C3215AF769BuL), (+1, 18, 0xADE96497FB5D4BB0uL, 0x1742ABCE8E0B4D34uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p125_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xDD2069E55811F2BDuL, 0xCC6788804DB5C82EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 3, 0xF2DA606A5B3E3BF5uL, 0xDBD5D7A26113C908uL), (+1, 5, 0x909C86FE5FC0FF8BuL, 0xF9F56E1F6B1D8362uL)),
                ((+1, 7, 0xE3D19FB198CBDE86uL, 0x427D06AF0CD19DE8uL), (+1, 9, 0x8C8C790D46C3D505uL, 0xAA4219D6FE679BD8uL)),
                ((+1, 10, 0xEE7AE1DD134F5A77uL, 0x778094ED791CF52EuL), (+1, 12, 0x99DD69736C7F2937uL, 0x0E054FCB3AD9A7E9uL)),
                ((+1, 13, 0x98CF79568B4AF294uL, 0x06E77A5D818FC80DuL), (+1, 14, 0xD0CBC0F1FBEC3C97uL, 0xFCB60714297AC6B4uL)),
                ((+1, 14, 0xF7D27916F5529DE0uL, 0xEB4B157652BE466BuL), (+1, 16, 0xB647466D6B9A4D96uL, 0xAA6819322D2C146FuL)),
                ((+1, 15, 0xFE649CD04EE69D69uL, 0x9D19145FF244F882uL), (+1, 17, 0xCE01346FF236325DuL, 0x2C7F83A809EB2488uL)),
                ((+1, 16, 0xA0F7AC64A01028E1uL, 0x1ED1BB95BD1D4B72uL), (+1, 18, 0x940A0CE384081F07uL, 0x401C2000A3E484D9uL)),
                ((+1, 15, 0xECE593490FCA1B80uL, 0x8ACE742A44FCEBE5uL), (+1, 18, 0x817A9B5C6569A579uL, 0xF8C38CA85C891F80uL)),
                ((+1, 14, 0xB598718495254E32uL, 0x2BFBE770240B3038uL), (+1, 16, 0xFDEBD7E1378EB2BBuL, 0xEDCF4B33F04BB6BFuL)),
                ((+1, 11, 0xE7895120950D2DA3uL, 0xCFA0173EA9871535uL), (+1, 14, 0xEED88232BE5C8C00uL, 0x83A3AF13E0BCD487uL)),
                ((+1, 6, 0xF86CE328F87953D3uL, 0x94E4696F73523A1BuL), (+1, 11, 0x95EB10EEDF377D0DuL, 0x62111B75D9154DBDuL)),
                ((-1, 0, 0xE47E4150753F7A2BuL, 0x27861C7857075AE8uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p25_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xC6348C79A0E25004uL, 0x7226F6F78A1B8963uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0xE3B8746D1F23DBD0uL, 0xC355C8D33A18FBE2uL), (+1, 4, 0x992BA2EF00DBBA67uL, 0xF5538C628AF103AAuL)),
                ((+1, 5, 0xE01A5234D4A69F0AuL, 0x83C8C4001054849EuL), (+1, 7, 0x9E4EF888FF480B7FuL, 0xE8B16B75801B24D5uL)),
                ((+1, 7, 0xF714133862CEA860uL, 0xABA6467A66CD4DD9uL), (+1, 9, 0xB939EAE1CD2499C4uL, 0x03CCE49D3FEBD33FuL)),
                ((+1, 9, 0xA7A334EAB5FF665DuL, 0x2A54D8F69B29DC25uL), (+1, 11, 0x87244174AE54F4CCuL, 0xA2610E28888ACD86uL)),
                ((+1, 10, 0x90E4767FCAA1A3E2uL, 0x1F995C3C49E794D3uL), (+1, 11, 0xFF870CC5B8E48965uL, 0xD20CA49C420A8810uL)),
                ((+1, 10, 0x9FC178775E300569uL, 0x5CE4747461E3BA96uL), (+1, 12, 0x9D97F04D6ACAE6DBuL, 0x9AA8041E8D053257uL)),
                ((+1, 9, 0xDAF5F06191BDB1C4uL, 0xD9B54C78E7FBDDE9uL), (+1, 11, 0xF93D9201AF9704A1uL, 0x60876C7875CC654EuL)),
                ((+1, 8, 0xAFFF10EC2B74D236uL, 0xC7133E46E6640B9EuL), (+1, 10, 0xF1E2040E0B1515ABuL, 0x132CC95F75280CC8uL)),
                ((+1, 6, 0x949A70752DD5F28DuL, 0xB73DBE27B209898EuL), (+1, 9, 0x84AB28519C400EF8uL, 0xD436B73528256D4AuL)),
                ((+1, 2, 0xD24E7A55D586E3A7uL, 0xCC3D5BB46297879EuL), (+1, 6, 0x8CB5D9BADD3F8FEAuL, 0x0B9794549A79ED33uL)),
                ((+1, -4, 0xFB60B413C6642A27uL, 0x126C02AD45127141uL), (+1, 1, 0xC8A5D635A1EC4EBCuL, 0x8319073F84619A9DuL)),
                ((-1, -10, 0x848698923C826FF1uL, 0x7D4E52308B86BEE9uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p5_0p75 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xA9A167C436229A0EuL, 0x90C25139D92AEAC2uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0xAAE52D0E9A995965uL, 0x069D0A38B4218C1FuL), (+1, 3, 0x8933750AE0D8D554uL, 0xF83BDA8EDFF5C879uL)),
                ((+1, 3, 0x8F4E03CF3FDAC43BuL, 0x07406D41E3DD341AuL), (+1, 4, 0xF7ED6C80A43F748CuL, 0xE9F95B219A2285E2uL)),
                ((+1, 4, 0x81A912EDFE95EF09uL, 0xB317B3F59F118B11uL), (+1, 5, 0xF5A7FD78D480239EuL, 0xD0192CB6DC1640B2uL)),
                ((+1, 4, 0x89405BACB7B728BEuL, 0xA6B0DD3A2478D505uL), (+1, 6, 0x917607A156CC7CA1uL, 0x612A9150550965A3uL)),
                ((+1, 3, 0xAC6949CEA6920981uL, 0x2A8AE74A7043BA4BuL), (+1, 5, 0xD27CDB251809304BuL, 0x3D051E510E1E9B74uL)),
                ((+1, 1, 0xF8F9C7BFDACD7314uL, 0xACC8566AA05D11D8uL), (+1, 4, 0xB69F535582C85A12uL, 0x949E85B534FC9284uL)),
                ((+1, -1, 0xBDC7D38A1CE46297uL, 0x5184B7AC2C93A9CCuL), (+1, 2, 0xB2BE4C166210BD6AuL, 0x0E1220CF33105C74uL)),
                ((+1, -5, 0xFEF00428E1D3BB9AuL, 0xFD39F857063B91A7uL), (+1, -1, 0xADD472A0B30C5357uL, 0xED56968C72B78211uL)),
                ((+1, -10, 0xC1618CA55B98CE41uL, 0xBC67AD6B2C09D95AuL), (+1, -6, 0xFF4F41A30F7DB8AEuL, 0xA0D9DC73C387EB7EuL)),
                (Zero, (+1, -13, 0xEF46002A5B6F59ACuL, 0x2800419A1395C3DCuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p75_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x97AB04EC5FB43AE6uL, 0x1B33857AC0B848D7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xC6142537CE740B9AuL, 0x8D1C6D74641FBFDDuL), (+1, 2, 0xB3B4D9D20565CF35uL, 0xD29622576E63572BuL)),
                ((+1, 1, 0xD5DC174A92D0F640uL, 0x693AE21F7CEB7D02uL), (+1, 3, 0xD37EE1B609588966uL, 0x3B5FF125BE8FE898uL)),
                ((+1, 1, 0xF69ED5B241CF7126uL, 0xFCFA6FFE2B3A8DA1uL), (+1, 4, 0x875C11FD6B2B5C35uL, 0xA6F5E3F947D993ADuL)),
                ((+1, 1, 0xA3CF7C2233CD6373uL, 0x8B79A24DBEA79B49uL), (+1, 3, 0xCC85C7057B5CE5C8uL, 0x0D7812218529C004uL)),
                ((+1, -1, 0xFC06292EFC4E5FDCuL, 0xCA4137CF3E31C875uL), (+1, 2, 0xB9326BAC9CB62FD2uL, 0xC2DA584E5070D0B2uL)),
                ((+1, -3, 0xD6069578A11E5882uL, 0xBC044B70E0358920uL), (+1, 0, 0xC2F4C1FF23EA7E6EuL, 0x68FFCEEA8CED7353uL)),
                ((+1, -6, 0xB232BCFD8DE1C457uL, 0xD5045BB21066AD4CuL), (+1, -3, 0xDB6C57449A356536uL, 0x281A3C762EB19645uL)),
                ((+1, -11, 0xDE634D056DC2E2ACuL, 0x414EB996ADB94A78uL), (+1, -7, 0xDCC85E4B227C0D9AuL, 0x79E5748F92D33C3BuL)),
                ((+1, -18, 0xB090CCAE2CD96E3EuL, 0x6D02D8F69E6352ADuL), (+1, -13, 0xFF13BA1DB3BA6291uL, 0xA82A5933869A1D0FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x8AE540CE1948BF1AuL, 0xC96078B9CE4BB3FCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xA76EDB2786C7DE57uL, 0xD2889C8B8CEF9746uL), (+1, 2, 0xA474BF307CD2FF97uL, 0xA635A7DCEBB5B4B8uL)),
                ((+1, 1, 0xAD647EBCDFAD2713uL, 0x08811D0D445917DBuL), (+1, 3, 0xB72C6E4E76D59F8FuL, 0xD3C51CE156E9698AuL)),
                ((+1, 1, 0xC9DD0167EFF72EE7uL, 0x36962F641703B6FDuL), (+1, 3, 0xE7E5AD86295A1507uL, 0x70A7B48C5DDCBBA3uL)),
                ((+1, 1, 0x913232246CEC57EEuL, 0xAB355E5C70508B44uL), (+1, 3, 0xB7DE316CC91DAED6uL, 0xDF1F24F0A6D44527uL)),
                ((+1, 0, 0x85A706DC20CF92F0uL, 0xA1FDD9B11F26A987uL), (+1, 2, 0xBDCD2D5BEFAAEB6EuL, 0x7B4E3BB980E7A9C0uL)),
                ((+1, -2, 0x9DCDBA8F1D41F262uL, 0x719B4CE70ED3299BuL), (+1, 1, 0x807DE5B48A3CF26AuL, 0xF89DFA91EDD0D355uL)),
                ((+1, -5, 0xE947E2F201200BB3uL, 0xACF73BCDF4171E85uL), (+1, -2, 0xE07F5EF63B00A500uL, 0xFF3C4A9E62D7478AuL)),
                ((+1, -8, 0xCC78C840C1C933E4uL, 0x2C871B0080C8CC49uL), (+1, -5, 0xF2D2FC0B184270E9uL, 0xC3D5B7525CE80013uL)),
                ((+1, -12, 0xC03C62E2DD8219C8uL, 0x585C0D1D2244D246uL), (+1, -8, 0x96A944D90E442753uL, 0x578965BF54706994uL)),
                ((+1, -17, 0xA002A44D079CAA69uL, 0x61D47F4D85C32C35uL), (+1, -13, 0xBACABF54C4A8686DuL, 0xE127403FC345E4B4uL)),
                ((+1, -24, 0x95CA3CE2FB574A34uL, 0xCAD075EF894472C5uL), (+1, -19, 0xAE3F24DE580ECB4BuL, 0x571CB25721C14D98uL)),
                (Zero, (+1, -28, 0xCD16AB668394307BuL, 0xBC2C3B0214E6E08AuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xDB0FFBF2AA24474AuL, 0x336069106E35BA66uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x8B486E31FDDE3691uL, 0xA30E8506DAF55120uL), (+1, 1, 0xAE7A8251E1068C2FuL, 0x4BD6072977B29130uL)),
                ((+1, -1, 0x991BDB7690B3B132uL, 0x798C6D8E412C4508uL), (+1, 1, 0xCF5807F34C8AA32BuL, 0x81DBAF07B5DF5492uL)),
                ((+1, -2, 0xBEBCDF1758A7B90AuL, 0x98ED663A96CDE699uL), (+1, 1, 0x8D098F82B9E7B3CCuL, 0x2755DDCCE2174B03uL)),
                ((+1, -3, 0x945AD6BCAD830650uL, 0x24F61343C0DFF7D8uL), (+1, -1, 0xF28CA79DFA69EC00uL, 0x7B661C2799E2C95BuL)),
                ((+1, -5, 0x95BCA696C89BB5A8uL, 0xBF9858A201621BEFuL), (+1, -2, 0x8966C143F6F21789uL, 0x69953D165FA6400FuL)),
                ((+1, -8, 0xC58C16D8BC81190BuL, 0x225BE82A898690C0uL), (+1, -5, 0xCF7A552243CADE27uL, 0xFA828C05EBEE42A4uL)),
                ((+1, -11, 0xA78C36572952C708uL, 0xFAC94D51AF548B19uL), (+1, -8, 0xCEA6070A93D3A333uL, 0xAE2CD5381D21CE77uL)),
                ((+1, -15, 0xAF47AC24DD966CC9uL, 0xF979B3DE85C93DC3uL), (+1, -11, 0x83898E2D05F74E49uL, 0xDC52AB480A736DEEuL)),
                ((+1, -20, 0xD18B89456267DA45uL, 0x2530338300193636uL), (+1, -16, 0xC9C0D03FD1211D03uL, 0xCB3718A83B5A9D6BuL)),
                ((+1, -26, 0xF948C9C6CDF50266uL, 0xE1DDA9220D487F16uL), (+1, -21, 0xA815201D920D10E0uL, 0x507500FBA0919070uL)),
                ((+1, -33, 0xDE3C949585F78D8DuL, 0xAC4BC00456908335uL), (+1, -28, 0xFA8F2A701B0FF977uL, 0x2EE5B5DEC7FF192FuL)),
                ((+1, -43, 0xF9013BA7F7B42257uL, 0x82AABDB15AC01C92uL), (+1, -36, 0xD3E6C9D9A4A8E640uL, 0x5BB7C2B24364AAB6uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xA7B40EC313E5B148uL, 0xB505352E7C6B5A62uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xDDB5252BAEEA7543uL, 0x9273A319281490B7uL), (+1, 0, 0xB61E88E0B639141AuL, 0x33322A84FDDB43C9uL)),
                ((+1, -4, 0xFCFA58315F144789uL, 0x10FCB01F9198EEC4uL), (+1, -1, 0xE1A3C8EF98BBD94FuL, 0xF4AFC0C378E9115BuL)),
                ((+1, -5, 0xA352B52685E35D3CuL, 0xB171B9FAA8711158uL), (+1, -2, 0x9FD6D9CBBDB55214uL, 0x7B825CD6BAEEC7BFuL)),
                ((+1, -7, 0x837AD0E4E62C8C7FuL, 0xB201D024B7A83195uL), (+1, -4, 0x8EFA6140662AA50FuL, 0xEF06A681EB472810uL)),
                ((+1, -10, 0x892A15067B9C169DuL, 0xC96F29C282803B94uL), (+1, -7, 0xA8552CAA0E9A237FuL, 0x9B154F4E2E3B05EFuL)),
                ((+1, -14, 0xBACB12CEB85DCFC9uL, 0xD75547E85FA1D105uL), (+1, -10, 0x83EC71C145F8AF59uL, 0x7F99DB66D1BC0AE1uL)),
                ((+1, -18, 0xA35185A196442328uL, 0xDAC7BADEBB841ECAuL), (+1, -14, 0x883D1C26D6D3890DuL, 0xC2D7B1FBDB093FD6uL)),
                ((+1, -23, 0xAFE62AAB6CB86DC6uL, 0x4875DF67EA01485EuL), (+1, -19, 0xB39F4654E7C95379uL, 0xE4802C3A97E8749EuL)),
                ((+1, -29, 0xD8344CBFD80453BDuL, 0x5BB57F830BC5E8CDuL), (+1, -24, 0x8E7D50D8F35165D5uL, 0x1FA7A22C2B8C22A3uL)),
                ((+1, -35, 0x840B0E1F57217137uL, 0xFFACFDC52B85A3F3uL), (+1, -31, 0xF54340A4F6E7124EuL, 0x18A1B1DDE8E86994uL)),
                ((+1, -44, 0xF16450D84A75F506uL, 0xCC1F2C3A7E6FC713uL), (+1, -38, 0xBC94FB0E1DE481B4uL, 0x316915D80316C0A8uL)),
                ((+1, -54, 0x8A74DFF612CF961CuL, 0x5F12EC9FAE78E0D5uL), (+1, -47, 0xA44793CF58735A67uL, 0x0CAC271E1CCE0D1BuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xFB1B041498F55719uL, 0xA46D4EE7F6ED7BABuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x9E2706F345DD5218uL, 0xC5B52F2503FB1E78uL), (+1, -1, 0xAF042CBE1E689ECDuL, 0x551BBA812FF8CAFEuL)),
                ((+1, -6, 0xAC3C7AA9B21B88C4uL, 0xE6DAF90DB48942A2uL), (+1, -3, 0xD09EF37DABF2887EuL, 0xD49272ACB32287A6uL)),
                ((+1, -9, 0xD498CF63B5E031E8uL, 0x6984CE212489214BuL), (+1, -5, 0x8E55C6EFE3961435uL, 0xD2FD2415A188DE19uL)),
                ((+1, -12, 0xA3DB1D1B0ED44665uL, 0x6B4CECDEA9EB41B8uL), (+1, -9, 0xF580E405B3B376DCuL, 0xFE07163DC11C0D70uL)),
                ((+1, -16, 0xA3E5966B636AB29AuL, 0xA03FE87AA74A73E2uL), (+1, -12, 0x8B789832BDF3975BuL, 0x2C9845C3AA081CF1uL)),
                ((+1, -21, 0xD64E502A0DBA5EB8uL, 0x36321956FD4B5769uL), (+1, -17, 0xD32D078CF7C410A8uL, 0x00034A6BAC0007E5uL)),
                ((+1, -26, 0xB429B82AC74DF660uL, 0xF4BEE0289D6196FEuL), (+1, -22, 0xD2DEF81EA79590A3uL, 0x5536E9C1C8B98240uL)),
                ((+1, -32, 0xBAD71B47A6443ADDuL, 0x5352438FBD465AAFuL), (+1, -27, 0x868B755839E75F79uL, 0x9932C2C36F00D14EuL)),
                ((+1, -39, 0xDD741569F35A6918uL, 0x56C2C164E67CDB22uL), (+1, -34, 0xCED0C848505C8CD7uL, 0x3031E80F5AE5F310uL)),
                ((+1, -46, 0x829E52518D8CE33FuL, 0x782F04CB476DFFA7uL), (+1, -41, 0xACA313165C2C11A7uL, 0xB7265812B05CA6C7uL)),
                ((+1, -56, 0xE6FAA6EFFD878F1FuL, 0xB61F5D7CEAF82748uL), (+1, -49, 0x80E53C3E026F4C54uL, 0xA37EA04F3201A556uL)),
                ((+1, -67, 0x805DF888828A5FC5uL, 0x62B1CECD4F6570C1uL), (+1, -60, 0xDA56B05B1E1AC44CuL, 0x330198B4CF34B353uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xB8F2265ADCB974E8uL, 0x3AEFE6589340AD11uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -6, 0xED99783083E9D950uL, 0xE3893F8B6834E86EuL), (+1, -2, 0xB2DC82F9BBE324D3uL, 0x34AA423BF8428005uL)),
                ((+1, -8, 0x83F866ED5F3AE548uL, 0x967B00BB99A3C91BuL), (+1, -5, 0xD9EFD40A2DB27173uL, 0xC980236C879A48B6uL)),
                ((+1, -12, 0xA62D5D1261C1EB70uL, 0xF5565E1FAF188630uL), (+1, -8, 0x9805BF14F20A4ABCuL, 0x88B95A14D3ACAB48uL)),
                ((+1, -16, 0x82AB977FDD8EDEE2uL, 0xD9B4C0A2CE13EECAuL), (+1, -12, 0x861141DA13597315uL, 0x56524EDC6AA1FCBDuL)),
                ((+1, -21, 0x855AA5E231265691uL, 0xF3CE1DD9D77134D3uL), (+1, -17, 0x9BCA5FB639468211uL, 0x46DC1D8609081577uL)),
                ((+1, -27, 0xB1E8039CEC033687uL, 0xF3077072BD8B9829uL), (+1, -23, 0xF145472F20544859uL, 0xFC8902F40C8278BEuL)),
                ((+1, -33, 0x98966A0D29187ECAuL, 0x1A211FDFE92E13A5uL), (+1, -29, 0xF66E8D9974096DC0uL, 0x887B3E9EFCE8D26DuL)),
                ((+1, -40, 0xA16C3EE09169F274uL, 0xEC10D65F04DCE380uL), (+1, -35, 0xA0D30E94F5DF701BuL, 0x9A145B8E8EF418E8uL)),
                ((+1, -48, 0xC32187FDADBDF9DEuL, 0x2E83E894E073641CuL), (+1, -43, 0xFCD340C19130A978uL, 0x43FC350B11020023uL)),
                ((+1, -57, 0xEAB09037C1473920uL, 0x26441D1730C43313uL), (+1, -51, 0xD7C938FDE98407B0uL, 0xAE00618EB30213A1uL)),
                ((+1, -67, 0xD37AA478EE8B618FuL, 0x7DDFB5FAF9BA0336uL), (+1, -60, 0xA4AC982EC367B8F4uL, 0xBEB2B9AE5AE71278uL)),
                ((+1, -80, 0xEF6E8BE3C016FB98uL, 0x44371E69296F9AC9uL), (+1, -71, 0x8E7C733EB051B7BEuL, 0xAAE471A9A21FC90CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0x869F35DAA439BFE8uL, 0x75F13DA1CFD52817uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -7, 0xADDC403911C39A3AuL, 0x42497BC03A74083EuL), (+1, -3, 0xB42FCEC9F70E22BEuL, 0x74759F295F69CA8AuL)),
                ((+1, -11, 0xC225AF6791A23A31uL, 0xEEC9279AD45A62C8uL), (+1, -7, 0xDD2D1B9DC7DE689FuL, 0xEB254AE7BAA22D5EuL)),
                ((+1, -16, 0xF5BDDF51B0C5CBF9uL, 0x96663A33DB323F84uL), (+1, -11, 0x9B6B9F92722077A6uL, 0xDEF027702FDA22EBuL)),
                ((+1, -21, 0xC23AD0D5F7F0077DuL, 0x64AE16C723D9EC3BuL), (+1, -16, 0x8A1226D4E6CF0EE7uL, 0x6FE7F4653DBE8590uL)),
                ((+1, -27, 0xC73A7F594C9654FEuL, 0x29FD60EC9629E9A5uL), (+1, -22, 0xA19D291820D3B574uL, 0x088348E954A90107uL)),
                ((+1, -33, 0x858FD224FF7A47C8uL, 0xB4C6C67CA6363898uL), (+1, -29, 0xFC19BF9A421876A2uL, 0x767320D370A57C64uL)),
                ((+1, -41, 0xE63D450EC37EAC44uL, 0xEDDA82FD82B86EB1uL), (+1, -35, 0x81AAFA36E6DFECB8uL, 0x4A87194B47ECBE1BuL)),
                ((+1, -49, 0xF4BFDD3306488E1BuL, 0x579E8F6216E038B1uL), (+1, -43, 0xAA704EA30B9D7F4FuL, 0xFD1027887A8EFBAAuL)),
                ((+1, -57, 0x94A097058BF32D08uL, 0xD8C9DA91CD4C4755uL), (+1, -51, 0x86E570BF54541F19uL, 0x8DDA8677031F87BFuL)),
                ((+1, -67, 0xB393EA5137D64657uL, 0xB4F5FF74A3524CBCuL), (+1, -61, 0xE7D290322CF3A05FuL, 0xC15EDCF3A3669296uL)),
                ((+1, -78, 0xA288A69DAE967D53uL, 0x6528D6CCE4A1F58BuL), (+1, -71, 0xB21264911471B927uL, 0xEE98411D01EED810uL)),
                ((+1, -92, 0xB8CC5634B6B6433AuL, 0x0322C24211D18A16uL), (+1, -83, 0x9B0D692A60ADC84CuL, 0x1CB88A3DE0E046E8uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_limit = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xCC42299EA1B28468uL, 0x7E59E2805D5C717FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x8E663707F1B20D26uL, 0x83492925C9488519uL), (+1, -1, 0xBF5D5F6E6A53F53AuL, 0xDBE222A3FDB4B033uL)),
                ((+1, -5, 0x8CE62E4FB9D4C2EFuL, 0x024FD9AA25BA2BD0uL), (+1, -2, 0x9A2B03B1D75CDAA6uL, 0x345E6DE135E7A2C9uL)),
                ((+1, -8, 0xE043CA87E9256985uL, 0x4ECDA2A36D8E24E8uL), (+1, -4, 0x9993098B8F9A59F3uL, 0xA1052CF94CB31C5AuL)),
                ((+1, -11, 0x9E58224230A2A5AFuL, 0x8E4C793C7D237AFEuL), (+1, -7, 0xCE3B0972F10FC25AuL, 0x7245B292AA7A7A00uL)),
                ((+1, -15, 0x9759A6381207D7C0uL, 0xD8CD411A774CFA05uL), (+1, -10, 0xACD9F8B6DAB12F24uL, 0x4A28C65558B7D6A7uL)),
                ((+1, -20, 0x87B0C79664CFC320uL, 0x31D37B86349AF563uL), (+1, -14, 0xA142A97390314077uL, 0xC15180CEB85BFEB8uL)),
                ((-1, -24, 0x9E022A5D2C60944BuL, 0x7F1DB310E49C4C7FuL), Zero),
            }));

            public static ddouble Value(ddouble x) {
                if (IsNegative(x)) {
                    return 1d - Value(-x);
                }

                ddouble y;
                if (x <= 0.015625d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_0p015625);
                }
                else if (x <= 0.03125d) {
                    y = ApproxUtil.Pade(x - 0.015625d, pade_plus_0p015625_0p03125);
                }
                else if (x <= 0.0625d) {
                    y = ApproxUtil.Pade(x - 0.03125d, pade_plus_0p03125_0p0625);
                }
                else if (x <= 0.125d) {
                    y = ApproxUtil.Pade(x - 0.0625d, pade_plus_0p0625_0p125);
                }
                else if (x <= 0.25d) {
                    y = ApproxUtil.Pade(x - 0.125d, pade_plus_0p125_0p25);
                }
                else if (x <= 0.5d) {
                    y = ApproxUtil.Pade(x - 0.25d, pade_plus_0p25_0p5);
                }
                else if (x <= 0.75d) {
                    y = ApproxUtil.Pade(x - 0.5d, pade_plus_0p5_0p75);
                }
                else if (x <= 1d) {
                    y = ApproxUtil.Pade(x - 0.75d, pade_plus_0p75_1);
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
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    ddouble v = Sqrt(x);
                    ddouble u = 1d / v;

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u;
                }

                return y;
            }
        }

        private static class QuantilePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p125_0p1875 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 1, 0xB5B977A899AD27BAuL, 0x1F50DAC4194AAB4AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 4, 0xFE2D87A0FCF22F02uL, 0xEC4D9F94CC8E2D40uL), (-1, 1, 0xA7F4F3098BD742B5uL, 0xF4A8C3034772CA48uL)),
                ((+1, 7, 0x8BC176E6F0C0AF88uL, 0x41F7428D5D53CD34uL), (-1, 5, 0x86432CB4FC7E81AFuL, 0x7BEEA837E167B816uL)),
                ((+1, 8, 0x97697750F691FA15uL, 0x2AEE7FAA5393E0ABuL), (+1, 2, 0xBD4F27D5679CACC3uL, 0x89A1D9BB56B2D0DDuL)),
                ((+1, 8, 0xA18F8A5DADFC8CB6uL, 0xE10F802A67346930uL), (+1, 8, 0xB54B73E307458168uL, 0xE77718890EAFF061uL)),
                ((+1, 6, 0xFD759BB14A470CC9uL, 0xCECEE65080CC9821uL), (+1, 9, 0xB3E18AE03CE8B334uL, 0x705F44333A77CEEFuL)),
                ((-1, 4, 0xF3C9330627C9A7A2uL, 0xFDF8886A641D4A58uL), (+1, 8, 0xD5A75BEB4D1DCA54uL, 0x3C365DB4AA5FA433uL)),
                ((-1, 4, 0xCA27187F47F4F81DuL, 0xFCFA8875DB6F4107uL), (-1, 5, 0x82FE5626F987008CuL, 0x907BE76BAFBFCE12uL)),
                ((+1, -9, 0xB18BDB99E8CEE6B3uL, 0x9A714AE0FD6669D2uL), (-1, 5, 0xE2E00F79CD4A2737uL, 0xC3CFE78E8FA7FB5EuL)),
                ((+1, -2, 0xD641CCE5855297E0uL, 0x716CF9398F74CC6EuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p1875_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xA454A1E37C4C717AuL, 0x6826AFA033EA4020uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 4, 0xAE626F886CA29B21uL, 0x47945E439ED28143uL), (+1, 2, 0xA2AB71EAFC073E8AuL, 0xFAF3F501C85E94FCuL)),
                ((+1, 7, 0x96272F1D35D0B02BuL, 0xCC90CA8BE09F4D9BuL), (-1, 4, 0xB693D7F98C96504BuL, 0x7901BCDBE9DCDB36uL)),
                ((+1, 9, 0x8509FDA6A85E71DFuL, 0x1519ED9548C43CDDuL), (-1, 7, 0xA194F5D545D9F362uL, 0x88D27B7D81E51E80uL)),
                ((+1, 9, 0xFB3AF132DF220126uL, 0xC0BB1B6128E287AFuL), (-1, 4, 0xE75D9E84E563FA64uL, 0xEB545D2D868FF204uL)),
                ((+1, 9, 0xE417C909521F47BDuL, 0xFF05572730C2419DuL), (+1, 10, 0x993DAFA85DEC385CuL, 0x956BFE712D9A559CuL)),
                ((+1, 7, 0xE30614B8CD7C18A6uL, 0xB1E5E145D4583CA6uL), (+1, 11, 0x8BEC78EEE3FDCCFEuL, 0xD7CF9E4C802F13D3uL)),
                ((-1, 6, 0xD907B6D23D4ABAA0uL, 0x2EDC66FD21089F7FuL), (+1, 9, 0xD87C1BC91AB21ECDuL, 0xF88824E98F0A148CuL)),
                ((-1, 4, 0xF714F5054E25D273uL, 0xDAA29E4E1FBF4CD3uL), (-1, 7, 0xD22B812D1E1E9511uL, 0x16BB90FFAEEA88D7uL)),
                ((+1, 0, 0xDE336CBB9ADF0B04uL, 0x658881C1FA5D07D2uL), (-1, 5, 0xD17382FC2D5D51D7uL, 0x07325C069D77CD5AuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p25_0p3125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x9F3B3FF3ECE73CC0uL, 0xB935437275215948uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 3, 0xD20C39155DD7C0D6uL, 0x9435D20044A6C404uL), (+1, 3, 0x993DACC6B0BDC51FuL, 0xAAD8C2D94E41A6DBuL)),
                ((+1, 6, 0xDC4A215533CAE74AuL, 0x34F0B18C501E2533uL), (+1, 1, 0xC45F4353C69F294FuL, 0xF8B4DFC5C2277332uL)),
                ((+1, 8, 0xE67B3B5AD9D824F2uL, 0xDDE6145A3372227FuL), (-1, 7, 0xC1A7E6E026C7796CuL, 0xBE5DD8A5ADB5271CuL)),
                ((+1, 9, 0xF2F2BBBF2A615C6AuL, 0xFC8E7F4038FE8FA4uL), (-1, 8, 0xBE2B72FA0E52059AuL, 0x9AE9AE46ACA20644uL)),
                ((+1, 9, 0xD4E083531591EABEuL, 0xEB17A717A8B99A90uL), (+1, 9, 0xEE0C8AE1BBE21B42uL, 0x8DEEE06E19CEFF6BuL)),
                ((+1, 4, 0xA549B667B21C2FB8uL, 0x3E2EBA328BC402A7uL), (+1, 11, 0x9E927E5D649E0341uL, 0xBBCD787BB40AECDEuL)),
                ((-1, 6, 0xF33F0B80E527A00CuL, 0xB9A2718E4370CBAAuL), (+1, 8, 0x8777CC945261CC5EuL, 0xA07DB19E236DF328uL)),
                ((+1, 6, 0xB78940F98C3D4319uL, 0xF6FAEA6F281BC33CuL), (-1, 8, 0xEEBB908BB490FFD2uL, 0xF3F1E1415C338BA5uL)),
                ((+1, 2, 0xDC6DF99A7505F6ADuL, 0x1F0AE5B160464CC1uL), (+1, 8, 0x9361EAC9329D96C4uL, 0xCFCB50EE76FFB379uL)),
                ((-1, 2, 0x8C21A04766CCACF1uL, 0xD4945513A0A6AE54uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p3125_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x963737D238DC1D3EuL, 0xCE0A32D89ACCABA0uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 3, 0x962B525FCAC03AC1uL, 0x164F2DDD63BE257BuL), (+1, 4, 0x98A5B7309939FB08uL, 0xFCC27C35711CBFD0uL)),
                ((+1, 6, 0xEA35C551A7F8CDB0uL, 0x77331ADA6454D875uL), (+1, 6, 0xB771F066E5E60701uL, 0x15F694A475B71D4FuL)),
                ((+1, 9, 0xA9DA71D3E579A54DuL, 0x0E7BC50F6ED2F54EuL), (-1, 8, 0x899882B135A7F72AuL, 0xF573AF91F48B2391uL)),
                ((+1, 10, 0xAC0B8EAB767B622EuL, 0xE3A03018104726C5uL), (-1, 11, 0xAD2B68FF735F9B34uL, 0x337CBCD3DE03A890uL)),
                ((-1, 11, 0xE2444815D79DEFD4uL, 0xC19EE12214A69327uL), (-1, 7, 0xB0DDE0AD32FD4C6FuL, 0x74D2EF57D50C7143uL)),
                ((-1, 14, 0xAB190DFC9F6EF4E3uL, 0x8FAEF06B25FA2080uL), (+1, 14, 0xD0EA50C8DB481E68uL, 0x15F2D82F3F0032C0uL)),
                ((-1, 14, 0xE4832B0DE55548B6uL, 0x025A3DE1AD4B45ACuL), (+1, 13, 0xC6AD8DEC918B1D08uL, 0x65B983C45E27C82BuL)),
                ((+1, 11, 0xEF93F54A6FD9DA8CuL, 0x39C10C70818688F0uL), (-1, 16, 0xBAD2FFA291866FEEuL, 0x5722F6C33898CBEDuL)),
                ((+1, 13, 0xFB70FA49151F9190uL, 0x26DC75B82C4BFFFEuL), (-1, 14, 0xA618B2BBB0DA3F39uL, 0x4053F4E3C60720D7uL)),
                ((+1, 8, 0xFC3D63D23A9AF2B4uL, 0x3476D8C0204C2A99uL), (+1, 15, 0xAED98B099E4103BCuL, 0x84AB53A4CA8150C7uL)),
                ((-1, 9, 0x8C882B3C89FF5DC4uL, 0x06851787D53B637EuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p375_0p4375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xE6339135C9318CA6uL, 0x7C00C2FDE40B4575uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 3, 0x90A23745DC58751AuL, 0xD3C54A3DAB68A725uL), (+1, 5, 0xF3924EE917E1B73DuL, 0xD56AE7190FD8BB9BuL)),
                ((+1, 8, 0x9B4519B299E42D1FuL, 0x3FD72B5674448F45uL), (+1, 10, 0xBB1AEA3682954845uL, 0x1BB63A21BD9DDC05uL)),
                ((+1, 12, 0xB93570675B95E847uL, 0x76A8CEBA5523AD02uL), (+1, 14, 0x908DA5658C2D9D8FuL, 0x893A29D3A2EDCE11uL)),
                ((+1, 16, 0x8550C2E64DA23E34uL, 0x3CE667934B8C0489uL), (+1, 16, 0xD903872032267FECuL, 0x13F0A9DC8DAA7DB8uL)),
                ((+1, 18, 0xEA1883253D06C117uL, 0xB5EBB3BC894917EAuL), (+1, 17, 0x9F2118B42B562F72uL, 0x9FCD3C56208BDF22uL)),
                ((+1, 20, 0xECD8D118730DBAD9uL, 0x1F4280D5E0D00C6EuL), (-1, 20, 0xA803942AFB6F89F2uL, 0xC33DD92F13E9CA48uL)),
                ((+1, 21, 0xDC1D412C3DF219F3uL, 0x3B2ACC96A75A8B33uL), (-1, 22, 0x96AB67CDAD4FDAB6uL, 0xF7E311647742EC2EuL)),
                ((-1, 19, 0xE991B69AE490E623uL, 0x056C33A3847315D0uL), (+1, 22, 0xA21085A670C86C1DuL, 0x5C2BCE60A5D25CF4uL)),
                ((-1, 23, 0xBB1806B93C82BC6AuL, 0xA6CA12521DE1751FuL), (+1, 24, 0xCECC2FE8AEB8D262uL, 0xE8A28AD24B84BA08uL)),
                ((-1, 23, 0x8EEB147485724AACuL, 0xF61D9A5F10DACFF5uL), (-1, 23, 0xDE832D1DDDDB801FuL, 0x5E44E3811CF5EFE3uL)),
                ((+1, 18, 0x9DCE7C35D28B560CuL, 0x7B9062532D51F6F8uL), (-1, 24, 0xD1A21BA4A13A6DE8uL, 0x3D6A2B22A6BE3795uL)),
                ((+1, 18, 0xAF160DE54D84DD86uL, 0x5FDDDF37A4BC93D8uL), (+1, 19, 0xD2A48C1819968C65uL, 0xB2B67BB7C893998CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p4375_0p453125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xA3FF543A9F9FBA36uL, 0xD0F54FDE10475DDBuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0xF073BB16FC17C231uL, 0x66B5D0B1451BB52FuL), (+1, 6, 0x8AF09D40CB5C55C2uL, 0x088C9104446E97A1uL)),
                ((+1, 8, 0x937FE0DF79DD84AEuL, 0xE7574ECD8138057CuL), (+1, 10, 0xEDFC790CB74DF7EBuL, 0xB9A37C29E49F5F6DuL)),
                ((+1, 12, 0xC3E84F3B3ABA60DAuL, 0x97D65DE383E6CA20uL), (+1, 14, 0xC67075F54B62C497uL, 0x4A261190EC42FB34uL)),
                ((+1, 16, 0x97D932B51FCAED6DuL, 0x7E481B5C3E3A23C1uL), (+1, 17, 0x9851B8F4F6AE86C2uL, 0x66D69721728C411BuL)),
                ((+1, 19, 0x896E4DC1F583D1C1uL, 0x36EF138CB71F4520uL), (+1, 17, 0xB76C4F9E038A35DCuL, 0xD75AA82480F1B726uL)),
                ((+1, 21, 0x87F9907473746CA9uL, 0xC7A9CBF5F5A3EFB7uL), (-1, 20, 0xE59055696EB1FE08uL, 0xF0BFA19EA482329CuL)),
                ((+1, 21, 0xF3DA0A6CE2FB3FA2uL, 0x95071497CD097BA7uL), (-1, 21, 0xF0D3B803F0B78C57uL, 0x92D05163562A228EuL)),
                ((+1, 20, 0xC4BAEE4303F1947AuL, 0xDA6AF9170F84B24FuL), (+1, 23, 0x84171B05DF838C06uL, 0x4E35AFD23EB942B5uL)),
                ((-1, 17, 0x9AEBE95E0C1CD711uL, 0x71B4557DF5847BA5uL), (+1, 22, 0x8186BB24A351DF30uL, 0x1BDD7762A2BE9075uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p453125_0p46875 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0xD190A665C7709EC2uL, 0xE833A664589DFDC5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0xE6716A5C9EEE6888uL, 0x8E8F426C435C54E0uL), (+1, 6, 0xD4C3DD20B721D276uL, 0xE2B23D1B1C099DEEuL)),
                ((+1, 8, 0xD8B4EEC2A37D8AC1uL, 0x9A9FA484D3B049EBuL), (+1, 12, 0x93385B4EAFB96384uL, 0x0070D29BB618A5D7uL)),
                ((+1, 13, 0xE3452E09505552C3uL, 0xF33ECFFB637C6F05uL), (+1, 16, 0xD7F8AF36E395F9E5uL, 0xDB5CD1995BF856B4uL)),
                ((+1, 18, 0x90F23C2AE1F42AE8uL, 0xE819F42319A64218uL), (+1, 20, 0xAF2A4288863BB844uL, 0xC0A4EB97DC4DA932uL)),
                ((+1, 21, 0xE57377D67B79B485uL, 0xD0F8033B3107710CuL), (+1, 23, 0x8EF052E879F04C87uL, 0xC90DB11A6AA5D0FDuL)),
                ((+1, 24, 0xDB8FF0ACD212CC81uL, 0x77ACCA2B933F7F0CuL), (+1, 23, 0xF92B80C5AD41B04CuL, 0xA0BAED1C7250DC93uL)),
                ((+1, 26, 0xE99B98CEE51A873AuL, 0xAD100FF1F6D3024EuL), (-1, 26, 0xC30E958663F8B5FFuL, 0x4B2262FC437719A7uL)),
                ((+1, 27, 0xE08DEC3ECE920494uL, 0x675E2BE06E3FCB76uL), (-1, 28, 0x85715842C7346D0BuL, 0x1F13974256E5BCD6uL)),
                ((+1, 26, 0xAFA2669A2BEC61B5uL, 0x091CF8E60672A66DuL), (+1, 28, 0xFA1FB702E3B51A77uL, 0x434FAE5F8D550864uL)),
                ((-1, 23, 0xA4F561F5417B7AF5uL, 0x80726EE1EE8B735FuL), (+1, 27, 0xE0DA0D0A03DABC75uL, 0xB67CFE4B116C6F96uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p46875_0p484375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0xCB607F5A2CE0C246uL, 0x946918D1E253277FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0xE413902770DA6165uL, 0x6B327811FB15C781uL), (+1, 7, 0xDDB072B30E516C2CuL, 0xE81EFF0360F3C704uL)),
                ((+1, 9, 0xE0C0D71839932306uL, 0xBBE3138D74A6871AuL), (+1, 14, 0xA8A1E486672B703DuL, 0xF8A44D7A71C4588AuL)),
                ((+1, 16, 0x8033A598767C6577uL, 0x2DE13F2A818D4252uL), (+1, 20, 0x92AD39D01AAC8596uL, 0xDBD2E5E7C3722821uL)),
                ((+1, 21, 0xBB2E221E2349C445uL, 0x3F1FAB4DED9B824AuL), (+1, 25, 0x9F10F90830B39C95uL, 0xCF0B0FC4042CE721uL)),
                ((+1, 26, 0xB61A2146C6868B71uL, 0x528821E2E1055401uL), (+1, 29, 0xDB58F0EDDBA98054uL, 0x157F68A1E13EA915uL)),
                ((+1, 30, 0xEDD7DFD201DA83F1uL, 0x5B1480994B7468BBuL), (+1, 33, 0xBA4ED1100383797CuL, 0xEE4DAE8A3C50217AuL)),
                ((+1, 34, 0xCC630ABCC804E19BuL, 0x1F9B8F3A2290D61EuL), (+1, 36, 0xAD73541CF2012D59uL, 0x8921B099E5C57B6AuL)),
                ((+1, 37, 0xDB0BEA1B1360EC12uL, 0x099A8C85BD0C2F2DuL), (+1, 37, 0xD92239891F7E2101uL, 0xAA3D09BFF00B4FB1uL)),
                ((+1, 40, 0x82C4A8D8E5F361B7uL, 0xE2F407743A0490F6uL), (-1, 39, 0xE104D89351BFAC54uL, 0xFA1C151FDBE1E715uL)),
                ((+1, 41, 0x84020F0856CF9A7CuL, 0xD1F4E0012AC69691uL), (-1, 41, 0xD46E51C38708740EuL, 0xE846C556B0AB003EuL)),
                ((+1, 37, 0xD0858E1B04AF48EBuL, 0x5181A047C90D3BC2uL), (+1, 42, 0xD6825850550D9989uL, 0x530D759E35157FD8uL)),
                ((-1, 37, 0xE162EE9AFEA959DBuL, 0x864544C341D195FCuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p484375_0p4921875 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -7, 0xC9A8E38039CA9B9EuL, 0xE47D86FD6AD987DAuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 2, 0xA126D799CFD7D016uL, 0xFAF674ABA6EEB2D6uL), (+1, 8, 0x8C33A48FCC66AA12uL, 0x8466679D1E91E498uL)),
                ((+1, 9, 0xDA7BBCE5C5621BF5uL, 0x0FECF6BF64FA6408uL), (+1, 15, 0x87C7E7BF28F728EFuL, 0x1560383ABD7E618DuL)),
                ((+1, 16, 0xA78CDFEDC5051D57uL, 0x5939D145BE6B5399uL), (+1, 21, 0x957AFB4E0EB0C7E3uL, 0xC4D63D91904AF8C9uL)),
                ((+1, 22, 0xA09B16403C75C2B5uL, 0x5E2F0000B636F743uL), (+1, 26, 0xC9CBAE94A4199866uL, 0x7561EC2B752F9E85uL)),
                ((+1, 27, 0xC71F84685507A3A3uL, 0xF68B715CE4880A89uL), (+1, 31, 0xA790FCE8D6039701uL, 0xECB10126AB7C90A6uL)),
                ((+1, 32, 0x9E9BD9E8E24BCCCEuL, 0x37D07EEE5EE9F051uL), (+1, 35, 0xA078254111EE13F8uL, 0x146FEAD064844F1BuL)),
                ((+1, 36, 0x9AB7EAF7AC90CFB7uL, 0xA7818D9A138E22C7uL), (+1, 38, 0x8E376294E738441FuL, 0x7A751B6D9029668CuL)),
                ((+1, 39, 0xA47F0A9255C1EA98uL, 0xE90001F55775E401uL), (+1, 36, 0xA4EB04A91C01E510uL, 0x5879F6091F14DABAuL)),
                ((+1, 41, 0x8FA42D996D560DAFuL, 0x293BCCC8A02FCE68uL), (-1, 42, 0x9EC4B771A6E994F6uL, 0x00580742FEE32987uL)),
                ((+1, 39, 0xB07E38620D0F78F4uL, 0x979379104121EFA2uL), (+1, 42, 0xCCCFDD8A885AAD5BuL, 0xD6403EF5998A763FuL)),
                ((+1, 39, 0xA7BD8159C6C36B62uL, 0x30837116EDF9231BuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p4921875_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                (Zero, (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0xC90FDAA22168C234uL, 0xC4C6628B81143EE8uL), (+1, 9, 0xC4583A1B6C8C93FDuL, 0xEE5EA19B0ED7041DuL)),
                ((+1, 10, 0x9A356E7A6AC59481uL, 0xEC07AA6EBCA2FDB8uL), (+1, 18, 0xA4EC3E014EA40DF4uL, 0x03B1922511A49109uL)),
                ((+1, 19, 0x818C8BA08E958E0CuL, 0xA5FEDE09E080EFC6uL), (+1, 26, 0xB6FB33FF391E4B5AuL, 0x0078D4CD81595F09uL)),
                ((+1, 27, 0x8FC56DC3D18C173DuL, 0x5795CE64A7925CD2uL), (+1, 34, 0x92D44FBA1C75264EuL, 0x27D21D5C2631329EuL)),
                ((+1, 34, 0xE6D57C760EB3121DuL, 0xA165A4968FB1743FuL), (+1, 41, 0xB08BC578A7B220B9uL, 0x70DD26362A16FF0FuL)),
                ((+1, 42, 0x8ADFFACD3E199640uL, 0x97A29D407AB3C534uL), (+1, 48, 0xA21675BFE296071FuL, 0x63452D5B97A579A5uL)),
                ((+1, 48, 0xFF4C7FB79BE158FEuL, 0xE04443678FE161AFuL), (+1, 54, 0xE4134244774BDC23uL, 0xA7DD8A79A9230C0DuL)),
                ((+1, 55, 0xB3F5915118BFFEDEuL, 0x62F6E9C6BD7FB196uL), (+1, 60, 0xF4A456AF3248F9E4uL, 0xDEB9512486DA0FE6uL)),
                ((+1, 61, 0xC1A80737E0A46A3BuL, 0x54ADF647CAB088BAuL), (+1, 66, 0xC4B049BC38D1964EuL, 0x9D38ABF357A89E78uL)),
                ((+1, 67, 0x9C97BF472A157466uL, 0x2E3B0CD9C328BD06uL), (+1, 71, 0xE61DE43E59330C90uL, 0xB21ED1999822E52CuL)),
                ((+1, 72, 0xB934702989A9EAE5uL, 0xB457ABEC772EDFC4uL), (+1, 76, 0xB9BCE654B860A667uL, 0xD09958E4DC82C399uL)),
                ((+1, 77, 0x98E57E7A2AB9F7FDuL, 0x01FB75261A746B78uL), (+1, 80, 0xBB5A36B027C7CF5AuL, 0x2E02919C8CB8E340uL)),
                ((+1, 81, 0xA2DDA602359052B3uL, 0xB0E8608FD7783EA5uL), (+1, 83, 0xB904BE526F15E4BEuL, 0x21423F9FE872E73CuL)),
                ((+1, 84, 0xBFEF7A4C0D7D1494uL, 0xD4EF9B1C57FEA595uL), (+1, 82, 0xE30A75D770EA38FAuL, 0xF87819D60ECC5645uL)),
                ((+1, 86, 0xB2870C22DA0515C5uL, 0x544C471519B3D795uL), (-1, 87, 0xDC72A2765EDEE9CCuL, 0x1A0271989D76D65CuL)),
                ((+1, 83, 0xC493DBB5C7621066uL, 0xCAF293ADE52FF7DCuL), (+1, 88, 0x928D92A6445FF120uL, 0x988A9C1FF9191633uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xF4C245CF84B23F18uL, 0x1CDE92FFB45B321CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xE8983249DFCF2F6AuL, 0x7C7EE3A5F642650EuL), (+1, 0, 0xD5F8F6AE185D9C94uL, 0xC231ECB143ADA381uL)),
                ((+1, -3, 0xB31630B0B29C9B56uL, 0x1C4F024B9832D55EuL), (+1, 0, 0x9461BE4939E180D7uL, 0x6267539EB4229BD0uL)),
                ((+1, -4, 0x92175A8A1BE0CF0AuL, 0x459B9A86EA9D499EuL), (+1, -2, 0xE3E1A203E766CD77uL, 0x479CF4C87E3A1622uL)),
                ((+1, -6, 0x8E3103C2977EAE35uL, 0xC32DC2857BD39FA8uL), (+1, -4, 0xDC8462F2B05937ECuL, 0xD2596B23F0F81186uL)),
                ((+1, -9, 0xB75D8A26A1764175uL, 0x8AAD141D65391EF9uL), (+1, -6, 0x90CAC0F31F867243uL, 0x98906BC66CBE76ABuL)),
                ((+1, -12, 0xA9C1159683B31082uL, 0x36E8CCF8EAF7C7C6uL), (+1, -9, 0x85A688F7598EAC4DuL, 0x8BFB38C9636A857CuL)),
                ((+1, -16, 0xDDE1AF1B5B0D68DCuL, 0x3ABB56252758F840uL), (+1, -13, 0xAE9DF37F9AB489A1uL, 0x5D10EBE41A024278uL)),
                ((+1, -20, 0xCD19E604A067C3CEuL, 0xC20601ED033EBC7EuL), (+1, -17, 0x9D7BCF081C1554F8uL, 0xF546BCB24B7E2A82uL)),
                ((+1, -25, 0xD968B112DD1F7ED1uL, 0x47568089DE7B4561uL), (+1, -22, 0xB6EC55E597AEF37EuL, 0x69FA1AB0480CEABEuL)),
                ((+1, -30, 0x8E2565F366F124E3uL, 0xAD16D2413F23264FuL), (+1, -28, 0xBA039C0D9B7F5A86uL, 0x2EEEFE0B39DBDCBAuL)),
                ((-1, -39, 0xF0CF3E4C9BC94D1FuL, 0x6A3F20568A6EEDDDuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm4_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0x8EA2F4187EE1430AuL, 0xFD7DA99039CA1E82uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xEA56B3D3CC5E78B7uL, 0xBA74D039017F4B2FuL), (+1, 0, 0xC5AA2658E4093E3CuL, 0xB9CB161422526E76uL)),
                ((+1, -3, 0xA9815029D6CD9662uL, 0xE6FFF8CBADB0D68DuL), (+1, 0, 0x88F85CC7F6E0C05BuL, 0x30F0C51C6A2F2621uL)),
                ((+1, -4, 0x90965EC4631DB5D7uL, 0x5E2E30590E09AC2DuL), (+1, -2, 0xE46248C4F207A680uL, 0x073906A95DEF4CDDuL)),
                ((+1, -6, 0xA53476546C1015AFuL, 0x895B484104FCE2E9uL), (+1, -3, 0x8187D162944397B6uL, 0x63AB52B509A175D8uL)),
                ((+1, -8, 0x8860200F8CC89F73uL, 0xAAE5E50CBD175B02uL), (+1, -6, 0xD5EDB6309C9C34DCuL, 0x23CFC3EB3EDEE336uL)),
                ((+1, -11, 0xAAC7A027695C7C9CuL, 0x92BBC2DD977A7118uL), (+1, -8, 0x862B55EE7A3FA504uL, 0x466A0CF830CC0020uL)),
                ((+1, -14, 0xA7632AEF0D3ACE3BuL, 0x2D89F74F66539DEFuL), (+1, -11, 0x8389AF427B755787uL, 0x4F01CE3FD19B52D8uL)),
                ((+1, -17, 0x82F201DD32496D05uL, 0x251A0D6D13D84243uL), (+1, -15, 0xCD9B290453B14C08uL, 0xB01DAA741FA8C263uL)),
                ((+1, -21, 0xA52A7A7047F95B56uL, 0xBA20D088B0E808EFuL), (+1, -18, 0x81BB95E2A642226EuL, 0x306D4294C7BE9656uL)),
                ((+1, -25, 0xA952A8833749ECD6uL, 0xFFE287077BAACCA3uL), (+1, -22, 0x84F85E9E6693D377uL, 0xE37FF32B902E8C7FuL)),
                ((+1, -29, 0x8CB4175391D4FC19uL, 0x451BB663EF150798uL), (+1, -27, 0xDD1B8BA23FBFD0E3uL, 0xDF535A7DCEA35D24uL)),
                ((+1, -34, 0xBBD178BCE51F9B05uL, 0x8A0CDE382076EF9FuL), (+1, -31, 0x9369145123AF0BE5uL, 0x39B458487F7699D5uL)),
                ((+1, -39, 0xC38C4F950D466FECuL, 0x315B62D4C5645BE5uL), (+1, -36, 0x99B2F864C8D70311uL, 0xC20F46D3B04E4205uL)),
                ((+1, -44, 0x93EC492587F83A98uL, 0xC80A8445DFD05C16uL), (+1, -42, 0xE837814B5D22804CuL, 0xDF4B4065B976972DuL)),
                ((+1, -50, 0x935567DF485DFB2FuL, 0x04BC42AEBFD439A4uL), (+1, -48, 0xE780E034A61128E8uL, 0xA4A88017EF92E3A6uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xA1B39827F5CF4774uL, 0x8C234FB9D4EB8227uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0xF9CFDEDE83C796A8uL, 0x6326FDF430110BB0uL), (+1, -1, 0xC45984F93ED46726uL, 0x1DEED75D0827245BuL)),
                ((+1, -5, 0xB0C6357E21DCC555uL, 0xDD90C0A8236AFAA2uL), (+1, -2, 0x8AC1AA0189416FBFuL, 0xC23AB1242A960CFBuL)),
                ((+1, -7, 0x99EC367C2215DF4FuL, 0x74E2C1AA6ED790ADuL), (+1, -5, 0xF1C7692EDCBD0301uL, 0x93ABB92DA7BF1D55uL)),
                ((+1, -10, 0xBA6A512E413F09CBuL, 0x4EE58BF13C06987FuL), (+1, -7, 0x926EADF6E71DF5D9uL, 0xBBFDC2A11B369BB6uL)),
                ((+1, -13, 0xA7980576195942E8uL, 0xAA89B86BA28055BEuL), (+1, -10, 0x839E53FC1448B87EuL, 0xED240956592B9F4BuL)),
                ((+1, -17, 0xE89AB642CDDE85B6uL, 0x4E7CC6A829DEB22EuL), (+1, -14, 0xB6AEE8C3AEA90219uL, 0xF1E5D393F3CE9A7DuL)),
                ((+1, -21, 0xFF85F343D09D8FFAuL, 0x40475147B443A7E0uL), (+1, -18, 0xC8B2974F0758578EuL, 0x60EBD490D27E629AuL)),
                ((+1, -25, 0xE1736FBD35B09F1AuL, 0x37CB7F7922BCC1DDuL), (+1, -22, 0xB110BE85152F90F9uL, 0x7F3DD21D99F97804uL)),
                ((+1, -29, 0xA0A8D25E9F7FD978uL, 0xABED1BE28D8D27ADuL), (+1, -27, 0xFC5A011DE71A1415uL, 0x5D5D2F66590BC859uL)),
                ((+1, -34, 0xB849C4803AB3CEA0uL, 0xC8F30D5A6B40701FuL), (+1, -31, 0x90C1A62B1025B773uL, 0x5A05D4F86AABFC98uL)),
                ((+1, -39, 0xA7A2E56C7DF689C3uL, 0x8A0770557926E5D9uL), (+1, -36, 0x83A37A5678D7B05AuL, 0x2E910A2DD19F9B39uL)),
                ((+1, -45, 0xE95C00BB34AD3CD3uL, 0xBB82129F0C866B59uL), (+1, -42, 0xB7515A9D6F225188uL, 0x35DB60AF9831F7C0uL)),
                ((+1, -51, 0xE5A03364FD6F5E88uL, 0xE56E61C8772CB8CBuL), (+1, -48, 0xB44EB31D1F2933C3uL, 0x7BF9065D616DCA55uL)),
                ((+1, -57, 0x837821685B5BC395uL, 0x05327D88B389C503uL), (+1, -55, 0xCE90063FE5A696BAuL, 0x9EB45398203A4818uL)),
                ((+1, -77, 0x9C031AC5FEAC276BuL, 0x9F943AC2A9064061uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xA2F83D7B4F18AF40uL, 0x775BC0EC1C3E39B5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0xABCFF1848185D68BuL, 0x516A48267195437DuL), (+1, -1, 0x86F0938915A7B76BuL, 0x04E7BFFA6610548DuL)),
                ((+1, -6, 0xA515EAF0026EA90EuL, 0xB3D8011609084E5BuL), (+1, -3, 0x81A871C551898AF2uL, 0x2CF6A7348CA0D343uL)),
                ((+1, -9, 0xC4F57704B62601E0uL, 0xC7BC2B7A454E009AuL), (+1, -6, 0x9AB111B917E74074uL, 0xD2E4114FEBC07D8FuL)),
                ((+1, -12, 0xA5AFF29621C61A15uL, 0x53EC6FD2E92BF3B4uL), (+1, -9, 0x82213DBDF68D6C80uL, 0xBD715F7F835DB4E9uL)),
                ((+1, -16, 0xD26275AC787F70E0uL, 0xF233F7AA9F13B8B2uL), (+1, -13, 0xA53C8C83B1C2DE90uL, 0x0EFF58F09EF0640BuL)),
                ((+1, -20, 0xD1FB24EF8B098C3AuL, 0x2492F4A8E39B3343uL), (+1, -17, 0xA4EB008CF3444D3FuL, 0x245FA45ED17A2769uL)),
                ((+1, -24, 0xA8EBEE9EEAC078CFuL, 0x47C22EC43FA9470AuL), (+1, -21, 0x84ABD4709AC13370uL, 0x0665066FE56AA4FFuL)),
                ((+1, -29, 0xDE57C3560409E2B0uL, 0x008812C21C1F49EBuL), (+1, -26, 0xAEA08CEDD970EB49uL, 0xEA8D9FF694C6D054uL)),
                ((+1, -34, 0xF10E6BC7194E99C2uL, 0xA542ED4C89DF4708uL), (+1, -31, 0xBD535A077A12759DuL, 0x4E27FB0091849D46uL)),
                ((+1, -39, 0xD71D9B4D20E9205EuL, 0x36CAB917B987BCF2uL), (+1, -36, 0xA8F3818A7C6A166AuL, 0x7CE9E6EF137B61AFuL)),
                ((+1, -44, 0x9C871C0B7814CC80uL, 0x27EC471C70F5CDA0uL), (+1, -42, 0xF5DF9C476A0C27EBuL, 0xBC0A112503B667B5uL)),
                ((+1, -50, 0xB5BDF62756FC7994uL, 0x9DDCCC81FDA05F51uL), (+1, -47, 0x8EBD6245B0266B59uL, 0xE1FF1BF74827693FuL)),
                ((+1, -56, 0xA13791AD29C2E846uL, 0x8984DB084540F6B1uL), (+1, -54, 0xFD3D400287ABC591uL, 0x5B623E8A36F121BEuL)),
                ((+1, -63, 0xC774769F9243512CuL, 0xB5811056387CAB76uL), (+1, -60, 0x9CA6C99785A0C3DAuL, 0x16B3271FED49A3ACuL)),
                ((+1, -70, 0x86B566CCB3DB8F1CuL, 0x5BC665C54E154D58uL), (+1, -68, 0xD3998BD29722E6EBuL, 0x98F8A8E8D74ACA89uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xA2F9836D08510E4DuL, 0x675054BA68E390BFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x833CEB7DFE7A0D16uL, 0x02341C5C04E76CBBuL), (+1, -2, 0xCE25EB184550AF5DuL, 0xA55A593F03B4EFDAuL)),
                ((+1, -7, 0xCEA549715AE7790EuL, 0x25BA9ACFDA3E3EE7uL), (+1, -4, 0xA24C92D5938F18E1uL, 0xF1EE58389447A110uL)),
                ((+1, -10, 0xD39ACB4AE4FCAD64uL, 0x2DF4D4230F3F0805uL), (+1, -7, 0xA631A46716F6D870uL, 0x39E7512B13893DBBuL)),
                ((+1, -13, 0x9E1CCC6E8A8691FAuL, 0x3DCC085DE96BBC33uL), (+1, -11, 0xF85CCE7751DC1C75uL, 0xB05F539E2117B0F7uL)),
                ((+1, -17, 0xB756756172E732B6uL, 0x11E83A97BF9F6AA2uL), (+1, -14, 0x8FFE3CCDF400A2FDuL, 0xF2E359F99995D478uL)),
                ((+1, -21, 0xAB2DD8BA266E7BF0uL, 0xB7CB207B0455244EuL), (+1, -18, 0x8671990B694F6B05uL, 0x4D0B5109E139EF32uL)),
                ((+1, -25, 0x83C45F8A52D540B7uL, 0xA1E62CCAD773223AuL), (+1, -23, 0xCEFAB01BA92A6E54uL, 0x5A1C0E15B140959DuL)),
                ((+1, -30, 0xA9C218A2747970E8uL, 0x90A90EB78C4BCA5DuL), (+1, -27, 0x8553E8B1D44B810EuL, 0xEE6C1BA972896A69uL)),
                ((+1, -35, 0xB894A02195455AFDuL, 0x518033B21EFD3849uL), (+1, -32, 0x90F8201330F450CEuL, 0x57DA6918ADF4023FuL)),
                ((+1, -40, 0xA9E2E87EBB31DFADuL, 0x0CDE14145D9E6F64uL), (+1, -37, 0x856DADEDF4F409AFuL, 0xA367F412FD91B5E8uL)),
                ((+1, -45, 0x83FE92CBF89D90FEuL, 0x20723322D40FCFBBuL), (+1, -43, 0xCF561BCE40C6EC89uL, 0x8CB2B0A626F15098uL)),
                ((+1, -51, 0xAB73F6CF02583E0FuL, 0x58EFC1D2053E663DuL), (+1, -48, 0x86A8AB0155917FF2uL, 0x7443580D8D47A1C4uL)),
                ((+1, -57, 0xB66ADD61591E1B16uL, 0x7ADA38CB25140406uL), (+1, -54, 0x8F4533DEEF521E47uL, 0x4937B22B10133899uL)),
                ((+1, -63, 0x98E024F672F9991CuL, 0xD4681BD465893C12uL), (+1, -61, 0xF022E96EB829BABCuL, 0x8C19E2CBC0E15560uL)),
                ((+1, -70, 0xB9F1F7AD75B6A200uL, 0x71C90584D3A5E4B7uL), (+1, -67, 0x920A7F7284D88298uL, 0x370C0ABBFCCB982DuL)),
                ((+1, -77, 0x84D4E040BA2C82F4uL, 0xA13782C6FB452C65uL), (+1, -75, 0xD0A6BBFAB9F983F8uL, 0xFC047E8EEE3D0DDAuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xA2F9836E4E441528uL, 0xB63450F559051A43uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x908E7693974FC849uL, 0x74B094450EE4BCBFuL), (+1, -2, 0xE3119DCF64E37E01uL, 0x1D09F9F5F9B79C53uL)),
                ((+1, -7, 0xFBE7447CDB01D18DuL, 0x1F0A789A44BB69B4uL), (+1, -4, 0xC5D82E7584D0D106uL, 0xCA36E784BEEA2EE9uL)),
                ((+1, -9, 0x8F789B80E9BE92DEuL, 0x40E0D009B985CE59uL), (+1, -7, 0xE15D2961CB1787F9uL, 0xDE5F12DE37E94EFEuL)),
                ((+1, -13, 0xEFD84B0A718EB9DAuL, 0x5D4F92B4B45E720FuL), (+1, -10, 0xBC5FAD6DAED92017uL, 0x090D2F099101356BuL)),
                ((+1, -16, 0x9C8235C9C3CE9983uL, 0x5909E7851E4D3699uL), (+1, -14, 0xF5D7DB0D2460BD3DuL, 0x739CF12EB4E60C7EuL)),
                ((+1, -20, 0xA593CEFAB1476A4FuL, 0x28E439D56E74A129uL), (+1, -17, 0x820B4E94B147BA98uL, 0xF082B235E0546143uL)),
                ((+1, -24, 0x91762B80913592DAuL, 0x9DCF4103DBBEC91EuL), (+1, -22, 0xE47D949E7EB4AED7uL, 0x13BD9116213625A3uL)),
                ((+1, -29, 0xD7B92BE89438356CuL, 0xED051D80420410B5uL), (+1, -26, 0xA96DBF8F7DB9CE2CuL, 0xEDB6B4160952FC3DuL)),
                ((+1, -33, 0x87B46151B1974BB7uL, 0x2D7D397608AFE16BuL), (+1, -31, 0xD52A0FC0C9E84406uL, 0x6C80D960B04B0FB9uL)),
                ((+1, -38, 0x9402255924947E84uL, 0xE323CC157AC9CD51uL), (+1, -36, 0xE87DB3B57CD83F22uL, 0x4CAD0F8365143A30uL)),
                ((+1, -43, 0x81D9A7A6007AD196uL, 0xFF8D60AEE9E643D3uL), (+1, -41, 0xCBF7DE8F7D183581uL, 0xE27D72AB605013D1uL)),
                ((+1, -49, 0xE58FE60BF83672B4uL, 0x1B900489D85E7305uL), (+1, -46, 0xB44C331BCD818D74uL, 0xF0352E0A0FDE23D6uL)),
                ((+1, -55, 0xC0A066C44F76A764uL, 0xFC5BB04B8802BE20uL), (+1, -52, 0x9749DE98BE0C85D1uL, 0xBD34827043B2F103uL)),
                ((+1, -60, 0xC95A9ADD1657C503uL, 0xD269BBF75E61A08BuL), (+1, -57, 0x9E249BDD56FC9149uL, 0x00F91D173D3B595CuL)),
                ((+1, -153, 0x8A0D0E75F950470EuL, 0x2EA7A6EFF6D6A908uL), Zero),
            }));

            public static ddouble Value(ddouble x) {
                if (x > 0.5) {
                    return -Value(1d - x);
                }

                if (x >= 0.125d) {
                    ddouble y;
                    if (x <= 0.1875) {
                        y = ApproxUtil.Pade(0.1875 - x, pade_plus_0p125_0p1875);
                    }
                    else if (x <= 0.25) {
                        y = ApproxUtil.Pade(0.25 - x, pade_plus_0p1875_0p25);
                    }
                    else if (x <= 0.3125) {
                        y = ApproxUtil.Pade(0.3125 - x, pade_plus_0p25_0p3125);
                    }
                    else if (x <= 0.375) {
                        y = ApproxUtil.Pade(0.375 - x, pade_plus_0p3125_0p375);
                    }
                    else if (x <= 0.4375) {
                        y = ApproxUtil.Pade(0.4375 - x, pade_plus_0p375_0p4375);
                    }
                    else if (x <= 0.453125) {
                        y = ApproxUtil.Pade(0.453125 - x, pade_plus_0p4375_0p453125);
                    }
                    else if (x <= 0.46875) {
                        y = ApproxUtil.Pade(0.46875 - x, pade_plus_0p453125_0p46875);
                    }
                    else if (x <= 0.484375) {
                        y = ApproxUtil.Pade(0.484375 - x, pade_plus_0p46875_0p484375);
                    }
                    else if (x <= 0.4921875) {
                        y = ApproxUtil.Pade(0.4921875 - x, pade_plus_0p484375_0p4921875);
                    }
                    else {
                        y = ApproxUtil.Pade(0.5 - x, pade_plus_0p4921875_0p5);
                    }

                    return y;
                }
                else {
                    ddouble v;
                    int exponent = double.ILogB((double)x);

                    if (exponent >= -4) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 3)), pade_plus_expm3_4);
                    }
                    else if (exponent >= -8) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 4)), pade_plus_expm4_8);
                    }
                    else if (exponent >= -16) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 8)), pade_plus_expm8_16);
                    }
                    else if (exponent >= -32) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 16)), pade_plus_expm16_32);
                    }
                    else if (exponent >= -64) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 32)), pade_plus_expm32_64);
                    }
                    else if (exponent >= -128) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 64)), pade_plus_expm64_128);
                    }
                    else {
                        v = Ldexp(RcpPI, -1);
                    }

                    ddouble y = v * Square(1d / x);

                    return y;
                }
            }
        }
    }
}
