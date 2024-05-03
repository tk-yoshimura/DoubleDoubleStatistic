using DoubleDouble;
using DoubleDoubleStatistic.InternalUtils;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Numerics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    public class HoltsmarkDistribution : StableDistribution<HoltsmarkDistribution>,
        IAdditionOperators<HoltsmarkDistribution, HoltsmarkDistribution, HoltsmarkDistribution>,
        ISubtractionOperators<HoltsmarkDistribution, HoltsmarkDistribution, HoltsmarkDistribution>,
        IAdditionOperators<HoltsmarkDistribution, ddouble, HoltsmarkDistribution>,
        ISubtractionOperators<HoltsmarkDistribution, ddouble, HoltsmarkDistribution>,
        IMultiplyOperators<HoltsmarkDistribution, ddouble, HoltsmarkDistribution>,
        IDivisionOperators<HoltsmarkDistribution, ddouble, HoltsmarkDistribution> {

        public override ddouble Mu { get; }

        public override ddouble C { get; }

        private readonly ddouble c_inv;

        private static readonly ddouble entropy_base = "2.0694485051346244003155800384542166381";

        public HoltsmarkDistribution() : this(mu: 0d, c: 1d) { }

        public HoltsmarkDistribution(ddouble c) : this(mu: 0d, c: c) { }

        public HoltsmarkDistribution(ddouble mu, ddouble c) {
            ValidateLocation(mu);
            ValidateScale(c);

            Mu = mu;
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

        public override ddouble Mean => Mu;

        public override ddouble Variance => NaN;

        public override ddouble Skewness => NaN;

        public override ddouble Kurtosis => NaN;

        public override ddouble Entropy => entropy_base + Log(C);

        public override ddouble Alpha => 1.5d;

        public override ddouble Beta => 0d;

        public static HoltsmarkDistribution operator +(HoltsmarkDistribution dist1, HoltsmarkDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static HoltsmarkDistribution operator -(HoltsmarkDistribution dist1, HoltsmarkDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static HoltsmarkDistribution operator +(HoltsmarkDistribution dist, ddouble s) {
            return new(dist.Mu + s, dist.C);
        }

        public static HoltsmarkDistribution operator -(HoltsmarkDistribution dist, ddouble s) {
            return new(dist.Mu - s, dist.C);
        }

        public static HoltsmarkDistribution operator *(HoltsmarkDistribution dist, ddouble k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public static HoltsmarkDistribution operator /(HoltsmarkDistribution dist, ddouble k) {
            return new(dist.Mu / k, dist.C / k);
        }

        public override string ToString() {
            return $"{typeof(HoltsmarkDistribution).Name}[mu={Mu},c={C}]";
        }

        public override string Formula => "p(x; mu, c) := stable_distribution(x; alpha = 3/2, beta = 0, mu, c)";

        private static class PDFPade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x931FE65BCE29D188uL, 0xDCFEA56837D33013uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -6, 0xDCE08DEF7E9DAEE8uL, 0xBB1CD0867E6F88DBuL), (-1, -4, 0xC02A5FA404B41DBEuL, 0xEED8A475F82D0B97uL)),
                ((+1, -8, 0xD9B5C6A945958AFCuL, 0x4C768047F14063A3uL), (+1, -2, 0xC8E41E5D2755814EuL, 0x48CF3D45631FE4E9uL)),
                ((+1, -12, 0xD6D32A3F03921E35uL, 0xC586781E29D570A5uL), (-1, -5, 0x881242ABB961B0F3uL, 0x64894BC27915B20AuL)),
                ((+1, -9, 0x9997EEEBE72BDCAEuL, 0x3BBAF91A162F372AuL), (+1, -4, 0x8A57828914BDBFCBuL, 0xA86BC3A46D75B49BuL)),
                ((-1, -13, 0xE2615198AB64DB23uL, 0xDC03DDF5C340A203uL), (-1, -8, 0xA3BB408EBF05A080uL, 0x494E3F067D84B085uL)),
                ((+1, -16, 0xA067917BEA592602uL, 0x32509F1EA808DA9DuL), (+1, -8, 0xD6D9DEDB3171A9DDuL, 0xB2DFB0A3C0F974E4uL)),
                ((+1, -18, 0xBC3B63F5AD467BD3uL, 0x8794B60FA56DFA1BuL), (-1, -12, 0xD279B006C9B43BD4uL, 0xC422111997BC1C2EuL)),
                ((+1, -19, 0xD275CB452BC0AEFBuL, 0xE0DFA137C215072FuL), (+1, -12, 0xC741E768DC0DCE3BuL, 0xE675B0B23ACC7CFEuL)),
                ((-1, -22, 0xBBC00AA5D066FC75uL, 0x80C7FE088ECA0036uL), (-1, -16, 0x918AABB09231E0DCuL, 0xB9B33A91BD732BBDuL)),
                ((-1, -27, 0xE63B9619AC980044uL, 0xE5CCE2958C229D9FuL), (+1, -17, 0xD2303636E5ACAAD3uL, 0x9D799CC05C808EAEuL)),
                ((+1, -27, 0x85DB255A805CD3D6uL, 0xD3AF767D4F723551uL), (-1, -22, 0xAEB26E6CD990B158uL, 0x6A6D90673C9B3A48uL)),
                ((+1, -34, 0xA6DF002B9D741B2EuL, 0xBFC4B114E680C541uL), (+1, -23, 0xC5CF078FFB00A80FuL, 0x3BF7EE53F7F3BEF3uL)),
                ((-1, -34, 0x8C2191B8D5465F76uL, 0x13F844400ADB4C52uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xCEE317603D2266E0uL, 0xE86D323A93526F6EuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -5, 0xAC8434AEDF865338uL, 0xF8797E5281A4CA33uL), (+1, -2, 0xECECE768220A0256uL, 0x0B41422B22B300B5uL)),
                ((+1, -6, 0xC42CEB0CD4894312uL, 0x42FE13597E7CEA8AuL), (+1, -2, 0xF2D957CAE29EF3C4uL, 0xFE3DFC1D5A757ECCuL)),
                ((+1, -9, 0xDD8EEB2E4B9C6008uL, 0xCFDA91D473FC221BuL), (+1, -3, 0xA69FFCA36A4BA714uL, 0xACD4A16DD36A6481uL)),
                ((+1, -11, 0x987EB5158B86BA23uL, 0xD0883F9649DEC63DuL), (+1, -4, 0xB89E2F343694BA7FuL, 0x1BC4307D432670D8uL)),
                ((+1, -13, 0x9562273B5C2F0421uL, 0x0ACAAC4258919962uL), (+1, -6, 0xC081002D051361F9uL, 0x5557BBB1D72BFA0CuL)),
                ((+1, -14, 0xDF349FCDB6EC6241uL, 0xABFD176A373F3041uL), (+1, -7, 0x92F08395A2F8E7FDuL, 0x507390668B686FA4uL)),
                ((-1, -19, 0x9FBA3A5EEA17C5C3uL, 0xB52D0C5A0DE8C9ABuL), (+1, -10, 0xE58CF9D3428FAF03uL, 0x95DA7401A2B4BEA9uL)),
                ((+1, -21, 0xECD024998B487AB1uL, 0x754021C676581E88uL), (+1, -11, 0x824121296E39F293uL, 0x8E1D6B908D7D0AE3uL)),
                ((+1, -22, 0xD5D383F526224BB5uL, 0xD2F0A2222813EF50uL), (+1, -14, 0x8E0952CB65C3702FuL, 0xC8A47CEC1EA31DF1uL)),
                ((-1, -27, 0x978202D7B525BC27uL, 0xD2A4196A5FFC6299uL), (+1, -17, 0xF3BB4F2D674787E7uL, 0x81554140DB94529FuL)),
                ((-1, -28, 0xB893F60968C5B61EuL, 0x6BABBD5A2301F82EuL), (+1, -20, 0x92D37C4ECE97C756uL, 0x90FC023BB6744315uL)),
                ((+1, -30, 0x879A9B75A6AC52C7uL, 0xD5E3EF6D92989B16uL), (+1, -23, 0xB9710A9A75B21436uL, 0x188F07757FB8A72AuL)),
                ((-1, -35, 0xE896659AA0BF7408uL, 0x16FD0C1201BF9C0AuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_3 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xAD231C2457E98CFDuL, 0x919659875288259DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -9, 0xD106043C397661F9uL, 0xBFFA16CF33F994D0uL), (+1, 0, 0x845E0DE5F8A4C035uL, 0xD8EA375A96DDC73AuL)),
                ((+1, -6, 0xB138EA32216BC693uL, 0x8C51E98691751B9CuL), (+1, -1, 0xD70689E341BEDB7DuL, 0x51783F03F6AA08E3uL)),
                ((+1, -9, 0xC16B21459DE4FCC4uL, 0xE5A9EB43787C291BuL), (+1, -2, 0xE8996B4C6591E91BuL, 0x526783C14EF2B3C2uL)),
                ((+1, -10, 0xF2BB4B08FC290E30uL, 0x04DF8C1FA85FEC06uL), (+1, -3, 0xD3D9E85ADDE51E40uL, 0xAA6FA0E43FF4505FuL)),
                ((+1, -12, 0xCB73AD8AFCD1C919uL, 0x19ED9BD21BA83D02uL), (+1, -4, 0x9707A917BEC7E7E1uL, 0x7EDFBBE50846214FuL)),
                ((+1, -14, 0xAB0405D1987B5931uL, 0xFCB41AA394A690B1uL), (+1, -6, 0xB964768194B3CDFCuL, 0xCA920E7301F06FCCuL)),
                ((+1, -17, 0xF291322297A8B976uL, 0xADFD5061236DE3F7uL), (+1, -8, 0xB80736EC3331A5A0uL, 0xFE7FEEC38EB196ADuL)),
                ((+1, -19, 0x8FE17D3A014A3D2CuL, 0xF018F79AAF8FF721uL), (+1, -10, 0x9CB90EA4F132AF3DuL, 0xE6C31F7962298D43uL)),
                ((+1, -24, 0xDB9270436A7CDB71uL, 0xFED97E0BAAC702D0uL), (+1, -13, 0xD20670D6C67EA319uL, 0x58C118F9A9A0B4DAuL)),
                ((+1, -26, 0xA83CC27BD7B8C0E4uL, 0xE769C58C5ACAE5A2uL), (+1, -16, 0xE9259954A7356AC5uL, 0xC4B31A788CEBABEAuL)),
                ((-1, -30, 0x8C7C09E6DCEC284FuL, 0xE7C7111492CF1B3BuL), (+1, -19, 0xAF6DC9B89FCDE2CCuL, 0x3083B6F9183791EDuL)),
                ((+1, -35, 0x8F0E3EE383AF1505uL, 0x7F593B4EC44AC916uL), (+1, -23, 0xC0237D7B42AEEFF8uL, 0xD5F2EBA4FE44D67AuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_3_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -5, 0x8110067F2AB415CFuL, 0xB8FC077FC15371DBuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -6, 0x9D87EE3D0DA9DF2AuL, 0xB478A9B07456AD92uL), (+1, 0, 0xC551D933CE6CFF32uL, 0x0EEB15297D5498EFuL)),
                ((+1, -7, 0xDC941EE4432666B0uL, 0xB630815FEC12B68DuL), (+1, 0, 0xABB35CFB69A33675uL, 0xB5264E2A49D0F67FuL)),
                ((+1, -8, 0xA1822FB451D5C69DuL, 0x286DCC59BF01B676uL), (+1, -1, 0xCCB53FEB5CEE95F8uL, 0x92D6D569AD63C39FuL)),
                ((+1, -10, 0xD5A87753A3F7B6AFuL, 0xAF8C12EE05992900uL), (+1, -2, 0xB7677DE1E93F5BB3uL, 0x8A6E1F55EB2D2DCDuL)),
                ((+1, -12, 0xC9375FC0CC3B4851uL, 0x6D61F5314E263FD9uL), (+1, -4, 0xFFAD08609DD4709DuL, 0x15D93A51DAC5E4FDuL)),
                ((+1, -14, 0x8CB345EBD6BE93C1uL, 0x8DBDBFA49564B237uL), (+1, -5, 0x8CFC2833E3D420A3uL, 0x3B897CCCE9E72A5CuL)),
                ((+1, -17, 0xA1554AF90D70C751uL, 0xF8CB12EB6E835514uL), (+1, -8, 0xF552A446B68645BCuL, 0x8DFABB7CFEB08AA5uL)),
                ((+1, -21, 0xA94AFAA96AEED045uL, 0x0E565E36DD55742CuL), (+1, -10, 0xA55F471630E80AA6uL, 0xBFB2A797B33FEE80uL)),
                ((+1, -25, 0xB63011A67678B52FuL, 0x169FBEB4D980999BuL), (+1, -13, 0xA548B1F7F08BB9FAuL, 0xB537F836AA29DEEAuL)),
                ((-1, -29, 0x9794025D86975189uL, 0xF2C96E41F6BADD81uL), (+1, -17, 0xDEBDA214CE0B25EDuL, 0x7DD1079E029589DFuL)),
                ((+1, -33, 0x81023B0E6CD0EF2BuL, 0xA564B9B184314420uL), (+1, -21, 0x9DE4063A6629D09EuL, 0x10C9B2F790E46B28uL)),
                ((-1, -39, 0xE37FC47FF9E350B3uL, 0x70E59D8D2E9852B2uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_6 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -7, 0xE00479757F8AD0D3uL, 0x671738978DCC39A4uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -6, 0x8426C18F73C11C38uL, 0x01AE269BF312B4A6uL), (+1, 0, 0xF56C0F068BE22F27uL, 0x9388A921E9040A1FuL)),
                ((+1, -7, 0xC4444AA741153FF8uL, 0xB052F821A84807A7uL), (+1, 0, 0xF674559F5EACFAA9uL, 0x8EE992D346B573FCuL)),
                ((+1, -8, 0xC266FDADE5A14FF3uL, 0xC51D22E8057310E4uL), (+1, 0, 0xA61FACEE1CD54078uL, 0xE477726EEB5D9C42uL)),
                ((+1, -9, 0x9273F53533C57D7BuL, 0xA39424BF2738A587uL), (+1, -1, 0xA5AE3BB0D1D45CF0uL, 0xD5DF88B43E7AA933uL)),
                ((+1, -11, 0xA755F8C262D86A17uL, 0x01D1121B48F7D04EuL), (+1, -2, 0x80043167A2438EAAuL, 0xF6B47209CF956028uL)),
                ((+1, -13, 0x95CCD7ED5108907CuL, 0xD575B872FC5814F2uL), (+1, -4, 0x9CDF1D2127E56CDBuL, 0xADE3C156EFEBB045uL)),
                ((+1, -16, 0xCBC94D6A70E8985EuL, 0xD15006E345539336uL), (+1, -6, 0x99BE69C231479C45uL, 0x83A30753AABB8773uL)),
                ((+1, -19, 0xD030CA4181FBD4BEuL, 0x433565B630DC31A9uL), (+1, -9, 0xF05E4FADC7A64BE0uL, 0x5F4E9AAA6DECC021uL)),
                ((+1, -22, 0x917794706CAD3A21uL, 0x9D4BFE3D38F6AC0CuL), (+1, -11, 0x93AD0C53270658C8uL, 0x412F3C04E39031B6uL)),
                ((+1, -27, 0xF6673BE4CA451DD7uL, 0x1EFF5C6A794703ADuL), (+1, -14, 0x8A55A63D25929884uL, 0x36552D67023BED12uL)),
                ((+1, -33, 0xF698782498986B23uL, 0x775F1BF1633BC2CCuL), (+1, -18, 0xBA930A1B0ECDB7B3uL, 0x4EDED6523709E3FCuL)),
                ((-1, -39, 0x9DE52B1A5B260CD7uL, 0xB0E10FF6CB23FE55uL), (+1, -22, 0xA1B0C247EA8693DAuL, 0x2794879055AC42CDuL)),
                ((+1, -45, 0x9EB7F42B2F124773uL, 0x7C0EFB887DFE1ABCuL), (+1, -27, 0x863A5EB850FC7938uL, 0x20F0EFC4EB6A7081uL)),
                ((-1, -52, 0xD21D1C48AAE37B6AuL, 0xAA7FDB867E5961DCuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -8, 0x8A64FA35566AF627uL, 0xD755C16EF31356E1uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -8, 0xC0E3BCEBB38025C8uL, 0x737C89F9B995BD92uL), (+1, 0, 0xEE910C6315461E2DuL, 0xB016454F9CAD32F0uL)),
                ((+1, -9, 0xFB6FAE4B06BE0759uL, 0xD4B6C446CB381DECuL), (+1, 0, 0xD085D70FB7254C79uL, 0x45694A5688ADF88CuL)),
                ((+1, -10, 0xC52E27C67A2FF7EBuL, 0x297E2AEF543DE931uL), (+1, -1, 0xE04E9D82B28EEFBBuL, 0xDD6E9BCBD8F9DB04uL)),
                ((+1, -12, 0xC9AF4CAD93C22D84uL, 0x99820C1B6DA42A4DuL), (+1, -2, 0xA404EAD8E99E38F4uL, 0x5810F37950F5E0C9uL)),
                ((+1, -14, 0x88B025B518BEBE6EuL, 0x68D57E2648E92864uL), (+1, -4, 0xAA2585D0BCD06089uL, 0x5773411103915330uL)),
                ((+1, -18, 0xEE004AB1441D5FC7uL, 0xD7B6C8D6A9E33671uL), (+1, -7, 0xFDA23FF0F157AB48uL, 0x45CACF9CAB347945uL)),
                ((+1, -22, 0xF2478480846F2031uL, 0xE6722DE758700C69uL), (+1, -9, 0x868E2DFB361FC827uL, 0x260DF1A5397BB49CuL)),
                ((+1, -27, 0xE47EA9E63DC1B0F9uL, 0x8BE8EA6CE692B443uL), (+1, -13, 0xC49C8B9BF5062B7DuL, 0x805DA4C220DACB13uL)),
                ((+1, -34, 0xCD79DFE19C6E0C0FuL, 0x57B4E455B48FE05FuL), (+1, -17, 0xB8DC2C3DA9EB8E71uL, 0x1D7A6B73CBF447E2uL)),
                ((-1, -42, 0xF2C8A9E6545632C5uL, 0x946D06BC6562B7CFuL), (+1, -22, 0xC3784EB076A5B817uL, 0x7E20E7E30C0EBB02uL)),
                ((+1, -49, 0xFC825E561B310FC4uL, 0xF89A749DA9EF69B7uL), (+1, -28, 0xA8449B319A12A927uL, 0xC54CBEA80E30ADDCuL)),
                ((-1, -56, 0xC2A689701058EDF6uL, 0x211A3887A9093D3BuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_10 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -10, 0xF9E372F1FE3E3026uL, 0x86781AABBECDFC8AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -10, 0xF545301CC604F9CEuL, 0x6969BF0A61A409B9uL), (+1, 0, 0xA91058B86EEF4640uL, 0x95BFE58349DE6853uL)),
                ((+1, -11, 0xCB2F0753417D1F45uL, 0xBF2D18BA29CD0B6AuL), (+1, -1, 0xC602776C0D01B2D3uL, 0x5043F708331E08B2uL)),
                ((+1, -13, 0xB6A3B2BDD0710176uL, 0xF99E46BD8C394B8DuL), (+1, -2, 0x86CA8EFD1CFCEE48uL, 0x08115F298D4FCE70uL)),
                ((+1, -16, 0xBDDA23CAB36D8CCCuL, 0x2676E5BD87CD6CAFuL), (+1, -5, 0xEAC6EF2D2CB6C0DAuL, 0xA33AF2E7EC300DF9uL)),
                ((+1, -20, 0xE017D661C4762DD2uL, 0x34615C2E7914032AuL), (+1, -7, 0x8787B6AD778CDE91uL, 0xE45FFF36E9723FC5uL)),
                ((+1, -24, 0x8891D93E2A1CFE37uL, 0xDCBCFB03278A46D8uL), (+1, -11, 0xCFB7C13D64EDC143uL, 0x910F71140E5C3CC4uL)),
                ((+1, -30, 0x8C016FEAF57C32C6uL, 0x22F7FECCC0258607uL), (+1, -15, 0xCD63FFA010BEB1D7uL, 0xEE5DDBA27DB50DCDuL)),
                ((+1, -39, 0xFEA9962BAB220226uL, 0xECFD89B6708A06A3uL), (+1, -20, 0xF5347EFF058FA698uL, 0x98260073F436843CuL)),
                ((-1, -47, 0x8AE2B3026C655742uL, 0x9BDDB7DF19A2334FuL), (+1, -25, 0x9A5ECFE1CB04DE95uL, 0x17E84AE03769C318uL)),
                ((+1, -56, 0xCFD86F41EF92DA50uL, 0x0BCE447824250795uL), (+1, -32, 0x93D94DB47FEA0044uL, 0x958C98E6AFFB2EBDuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_10_12 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -10, 0x895587856DFEF67BuL, 0xDAC00F47EB71EE60uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -11, 0x9362880484B89509uL, 0x7FED3D34BD46A6A9uL), (+1, -1, 0xCD4E99B95C445A41uL, 0x851D256CC6931CD6uL)),
                ((+1, -13, 0x80DC592F545F99EAuL, 0xA615876C10203586uL), (+1, -2, 0x8F9BABDEFA252823uL, 0x23DB9E80CC14903DuL)),
                ((+1, -17, 0xE825AA0C8FC7BA96uL, 0x3CEA8886380ACD67uL), (+1, -5, 0xE4B079CA7804802FuL, 0xDF390B01607641D3uL)),
                ((+1, -21, 0xE07DD643D4018FCBuL, 0xC2FA3E21A7844759uL), (+1, -8, 0xE2CE17000E91F2A1uL, 0xA31817569E514288uL)),
                ((+1, -26, 0xDB87F3D8D16D8227uL, 0xF05E91F7FDE5748CuL), (+1, -11, 0x8FE8CCCDBCD37253uL, 0x19939974BA1A3FA0uL)),
                ((+1, -32, 0xB4951B594AC12D4CuL, 0xC410C1798EA4D2F0uL), (+1, -16, 0xE69E97F40095B9A7uL, 0xC4B97E429A86C73DuL)),
                ((+1, -40, 0x8530EBF48E9E926DuL, 0x79218C49EF62E188uL), (+1, -21, 0xDCBC95CFA53FD345uL, 0x86A4CF949C6C3AC1uL)),
                ((-1, -50, 0xEF5B02C1D9475D2FuL, 0x88A98FC5CEADF6A1uL), (+1, -27, 0xDE66B5D8868C58ACuL, 0xFCD1CC0F439D4A15uL)),
                ((+1, -58, 0x9629C355D2420976uL, 0x47BD5199FF4078E0uL), (+1, -34, 0xAB6CFB042086F12EuL, 0x540FED64775D0F6FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_12_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -11, 0xA9E6F72532B98220uL, 0x085906F5598DDA6AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -12, 0x9B64262DE279D00EuL, 0xF5FC1D3BFE2F0ECCuL), (+1, -1, 0xACE545162ED83608uL, 0xA82D496CF34CD5AAuL)),
                ((+1, -15, 0xE954BE4672B79278uL, 0xEE0D8C53F5BF023AuL), (+1, -3, 0xCC6DDE917DF2C518uL, 0x983D13442D55E9F0uL)),
                ((+1, -18, 0xB874EE9D6E795D63uL, 0x61E6156B66EE78F7uL), (+1, -5, 0x8A944C559E1C4FBAuL, 0x1995E55E1F1FE948uL)),
                ((+1, -22, 0xA4088F2F7998D7E8uL, 0x24B5BA543C01C6B1uL), (+1, -9, 0xED03C98D4FBE80D2uL, 0xE6F4F255D7967A3AuL)),
                ((+1, -27, 0xA2AC4F50FE4F39FEuL, 0x1C4B204CF9460A14uL), (+1, -12, 0x8479ED167B2570DAuL, 0xB090B910E631E7A8uL)),
                ((+1, -33, 0xA71A1C0E2C39F195uL, 0x2ABD081532518C83uL), (+1, -17, 0xC1EC283598790AECuL, 0xEE6BFBEA602F7D56uL)),
                ((+1, -40, 0x9360ED8064BDA355uL, 0xF0A6AE5C11493095uL), (+1, -22, 0xB4D38936183F975FuL, 0xC2E344768513F6DFuL)),
                ((+1, -50, 0xED53D4AA71512526uL, 0x8071F57D768D9C6AuL), (+1, -28, 0xC9565DD02B9FE235uL, 0x87966AD4FCB5FD19uL)),
                ((-1, -60, 0xEC3CF594413FF711uL, 0x3235A4E1441853E8uL), (+1, -35, 0xEA5C100294D1ED0DuL, 0x3938AB9DE33A0546uL)),
                ((+1, -69, 0xA5ACCA8C89B011B2uL, 0xF9E750DE6F9E01D1uL), (+1, -43, 0xCE4690BF320724DAuL, 0xBAA8E6C9BAA68A84uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -12, 0xA113ECCF5311D562uL, 0x678D98BA12775F49uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -14, 0xF0CCAC229E1BB742uL, 0x3163EA6BF2725081uL), (+1, -1, 0x88E2D9E80793D2B2uL, 0x579B5680CEB148F4uL)),
                ((+1, -16, 0x984221F4421221A5uL, 0x9C3C3A14E53B69FEuL), (+1, -3, 0x81FDAB06219615E8uL, 0xFC52291A6A684389uL)),
                ((+1, -20, 0xD44C9334A9A8516AuL, 0x4CD6EBAAFEB9C4A7uL), (+1, -6, 0x90530ED085EC2E1CuL, 0x1F5F9AE7C142D713uL)),
                ((+1, -24, 0xB28EAFA8F5F1C8ACuL, 0x2B584F8DA549418DuL), (+1, -10, 0xCF983133FAE16C78uL, 0x6DB43042996E3B2FuL)),
                ((+1, -29, 0xBA94891B797EC702uL, 0x11ACEED0AC3FBD4DuL), (+1, -14, 0xCA88025BA16E5247uL, 0x922D1F82C09DC57EuL)),
                ((+1, -35, 0xF0C2A766A4B27AEDuL, 0xE79DA599D028E419uL), (+1, -18, 0x886868A495211D25uL, 0x55EFF4E59208E659uL)),
                ((+1, -41, 0xB88055A0989C3F86uL, 0xCD71F7AA26750CACuL), (+1, -24, 0xFD4CECD501DFBC0DuL, 0x34BF11ECAFEE46F2uL)),
                ((+1, -48, 0x9A1B8FA11B419957uL, 0x730CA589BF6FB9F3uL), (+1, -29, 0x9F0EFED32DCB7EADuL, 0xA4573C743337F315uL)),
                ((+1, -57, 0xE892550545C88F66uL, 0x5D69171F05BCA425uL), (+1, -35, 0x81DA60BBE0B9BE47uL, 0xBBD2F7252DF07864uL)),
                ((+1, -67, 0xAB028752411F140AuL, 0xBF38C544151B032FuL), (+1, -42, 0x805C48D8EEB1265BuL, 0x52880C9821470BEFuL)),
                ((-1, -78, 0xADBDB846207A04D6uL, 0x7E64F70B0C496C8DuL), (+1, -50, 0x8675A05DA8E74DDFuL, 0x082A2779947C3064uL)),
                ((+1, -88, 0xA3A34B50F24B9E68uL, 0xF1DF05F69B4AA78CuL), (+1, -60, 0xD9C837BB1FABE12DuL, 0xF06FA42DDADDCFF0uL)),
                ((-1, -99, 0xEF3B307E85EC6D08uL, 0x6581707C3B66EB0CuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -15, 0xDC830B27A5755540uL, 0x40DC86ED094FDA9BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -17, 0x84F3D896B0110CAFuL, 0xDDBCEFCDF38EA2DFuL), (+1, -3, 0xEB331831669C181EuL, 0x9BACC92B9D052CE8uL)),
                ((+1, -21, 0x8725187C268F5B64uL, 0xCA28EA3FBF4D691DuL), (+1, -6, 0xBF356FFBE8C44BF0uL, 0xE0AC92FAE64DEA70uL)),
                ((+1, -26, 0x96F6C99394689E65uL, 0x225A39ED366BBA2CuL), (+1, -10, 0xB50287293D2F1E25uL, 0x48ACA8FB8FBCBBDDuL)),
                ((+1, -32, 0xCAAF215B8720494BuL, 0x80CAC5B7AADE2DF0uL), (+1, -15, 0xDCF2EC23B244E51CuL, 0xD4011F17E6D2C0DAuL)),
                ((+1, -38, 0xA8294490B6960FAEuL, 0x49058EB543734D4BuL), (+1, -20, 0xB5E400F2B2841800uL, 0xECFA1EE93E1E65A4uL)),
                ((+1, -45, 0xAAE689596C2ABE56uL, 0x32B0160909CBD5F1uL), (+1, -26, 0xCD4BDDE842B5F866uL, 0xCEB7F66977A9F442uL)),
                ((+1, -53, 0xCBA87E2AB9C16A05uL, 0x5E9AF77A3D4D4DBCuL), (+1, -32, 0x9E4D5CCA2E82FA9AuL, 0xD6122DA37EDE511AuL)),
                ((+1, -61, 0x8194283030532A36uL, 0xB5D5815390781BCAuL), (+1, -39, 0xA3389D6A4CF76528uL, 0xF32FC1D744B263C5uL)),
                ((+1, -71, 0x9011930312F29963uL, 0x6163216C534D6931uL), (+1, -47, 0xD7757C517241FEDBuL, 0x5EA371140AFB1234uL)),
                ((+1, -83, 0x935F53E35253AC48uL, 0x353F98F3355778BDuL), (+1, -55, 0xA891CB6176191E90uL, 0xD5E2A5C39C7697C8uL)),
                ((-1, -96, 0xBA6BDF3491323D19uL, 0xA8F993A6FFD5BCC6uL), (+1, -64, 0x8785B9BA90AE65D1uL, 0x70C0F2E124EE73B0uL)),
                ((+1, -108, 0xA5003A15EDA2AC1BuL, 0x83F80819384182C0uL), (+1, -75, 0xA0A0840F7E17EA03uL, 0x2133605D3CCD8025uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_limit = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0x99319F36F945E34EuL, 0x5EC369E04605551FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 3, 0xC7F44BC1A41A5A09uL, 0x1ADF803B8F7F35CFuL), (-1, 5, 0xB3D62526EA1B0CDCuL, 0xFB7FBFB79C374F01uL)),
                ((+1, 7, 0x8B9D2B9406A2D6D4uL, 0xC7B37BB574585EEDuL), (+1, 9, 0x96E28EA838245BA9uL, 0x3995C87311F18741uL)),
                ((-1, 9, 0xA691874D484AAF81uL, 0x1A598AD24FD7D7EAuL), (-1, 11, 0xF11FA5523F1A33DFuL, 0x426984BEB2F31761uL)),
                ((+1, 10, 0x805DC56D2585DAC7uL, 0x1128712F1F942D46uL), (+1, 13, 0xB9813BB393C374CEuL, 0x8C207DD9C23BBC63uL)),
                ((+1, 8, 0xD3A751C9F5C2C109uL, 0x12E10A3AD4AB2D59uL), (-1, 13, 0xE3C5C292478B5418uL, 0x7F5AFCA84C0986E4uL)),
                ((+1, 7, 0xB93A5BD1714605CAuL, 0x5F1D57AB21D5CB80uL), Zero),
            }));

            public static ddouble Value(ddouble x) {
                x = Abs(x);

                ddouble y;
                if (x <= 1d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                }
                else if (x <= 3d) {
                    y = ApproxUtil.Pade(x - 2d, pade_plus_2_3);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 3d, pade_plus_3_4);
                }
                else if (x <= 6d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_6);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 6d, pade_plus_6_8);
                }
                else if (x <= 10d) {
                    y = ApproxUtil.Pade(x - 8d, pade_plus_8_10);
                }
                else if (x <= 12d) {
                    y = ApproxUtil.Pade(x - 10d, pade_plus_10_12);
                }
                else if (x <= 16d) {
                    y = ApproxUtil.Pade(x - 12d, pade_plus_12_16);
                }
                else if (x <= 32d) {
                    y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                }
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    ddouble u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u / x;
                }

                return y;
            }
        }

        private static class CDFPade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x8000000000000000uL, 0x0000000000000000uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -3, 0xFDC640DE5489B2EAuL, 0x909B82CFF67C57ABuL), (+1, -4, 0xA1E62F651F27C09CuL, 0xA5871FF9286EB9ECuL)),
                ((+1, -4, 0xEF8565AE305F7102uL, 0xC7E92049C30B2CBFuL), (+1, -2, 0x8F058948832AB6B9uL, 0x0441AD1BCDBC6F4FuL)),
                ((-1, -6, 0xF6F443C19FD5699AuL, 0xD2FC42FC20B832FEuL), (+1, -6, 0xF1C053008537B3B5uL, 0x6520BD8F534696B9uL)),
                ((+1, -7, 0x9F83D72B8D8413E4uL, 0xD7777E79C5118832uL), (+1, -6, 0xFCA4878BA9278651uL, 0x5BD2F3208BF47C9DuL)),
                ((-1, -10, 0xE8150D00D3A17256uL, 0x74A302B63A5BFA1CuL), (+1, -8, 0x8B3D47838441F3BEuL, 0x4B0803C29D4A6D24uL)),
                ((+1, -13, 0xF55A996F61769909uL, 0xFF811DFDB3F2119EuL), (+1, -10, 0xD1A0D573B30F42F3uL, 0x605AA419467EEBD5uL)),
                ((+1, -16, 0xAF14F828AAE0FC87uL, 0x2593CAB361646D7BuL), (+1, -12, 0x959086830214F032uL, 0x42F2EAF55B6B41CAuL)),
                ((-1, -17, 0xBC900AA40BB2F478uL, 0x63FB2471F30FDFBAuL), (+1, -15, 0x8A1B1759B71FC803uL, 0xA7B879FFD4CFA240uL)),
                ((+1, -19, 0x984D58626505C184uL, 0x03D8F54440B8745FuL), (+1, -17, 0x806735A4BC6C0B93uL, 0x8D9BC3D10D723BFCuL)),
                ((-1, -22, 0x8847620734F969E2uL, 0xE5E32FC78A7F5FB1uL), Zero),
                ((+1, -27, 0xEB9AC5A866B42954uL, 0x468FC63727F16871uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p5_1 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xB8A0025CAD721F6DuL, 0x8459A0279F4773FCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -3, 0xB11A7785DBB26275uL, 0x3ACF9DFC8B0191E9uL), (+1, -3, 0xFDB70E2901649F99uL, 0xA133005305984D1DuL)),
                ((+1, -4, 0xC1E61C11FB677E40uL, 0x415EE1E04EDC04E2uL), (+1, -2, 0x9F8C5B7871AAD51CuL, 0xD24A9F22828817F7uL)),
                ((-1, -6, 0xBB899FD6452C02F1uL, 0x128D36B11EB7605CuL), (+1, -4, 0x8154D48C490E0E38uL, 0x88839447ACEFD483uL)),
                ((+1, -7, 0x851FB7182AE5BCFFuL, 0x3D425481A17CA116uL), (+1, -5, 0x9C27E3AB1FDD90DBuL, 0xC0847E8FD4F535CFuL)),
                ((-1, -10, 0xEC117304F8C5728EuL, 0x1075D86E38BDA342uL), (+1, -8, 0xCA10314B9A9FAE2EuL, 0x9F91127445237D37uL)),
                ((+1, -12, 0xE98468D83C599AF8uL, 0x33AD076CAA07BF3CuL), (+1, -9, 0x907B0BD1342A04DCuL, 0x411854596C5CB3A3uL)),
                ((-1, -15, 0xF4AE699CD4384609uL, 0x53FACC55A801E4D3uL), (+1, -12, 0x8D718AAAF9AEB920uL, 0x8D3F45120479186FuL)),
                ((+1, -18, 0xAA6A2D79EE7B4B53uL, 0x18C2407E36D98F23uL), (+1, -15, 0xDF1311E23AAF207EuL, 0xF90B887D109D589BuL)),
                ((-1, -24, 0x9A959C19BCE5615AuL, 0xCF54F1ABFCAB25ECuL), (+1, -18, 0x908832EDB687B559uL, 0xE0DAEE4FD5EF8FDDuL)),
                ((-1, -27, 0xC8DA4CC7B6F64335uL, 0xBDA3F3CAEBD88B56uL), (+1, -23, 0xA06C9005E15FEA77uL, 0xDD8F2DB7C3F8A3ABuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_1_2 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xF98179F27573BC4FuL, 0x30ACA8D6CECA9DA8uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -5, 0xC3294BA026F37EF8uL, 0xD54574726FA60F1BuL), (+1, -1, 0xA236355775AA0336uL, 0x0D31214E488BFB3EuL)),
                ((+1, -4, 0x829A7E7D8532927AuL, 0x265224A2AB7740B5uL), (+1, -1, 0x824402BC582B798BuL, 0x54374A4D0EF1D5A8uL)),
                ((-1, -9, 0xE57F8049533ABE31uL, 0x0FED4D4212D3A795uL), (+1, -3, 0xDFFDC621C7F8A068uL, 0xA75FDCB03163DBECuL)),
                ((+1, -8, 0xC8466B7482908616uL, 0xA70571D2E59D89E8uL), (+1, -4, 0xC6AED00015D3EF90uL, 0x2E34F1EE8F79C5BDuL)),
                ((-1, -14, 0xAB820788BCA51A11uL, 0x873A4CBA798A56C9uL), (+1, -6, 0xFAECBC447087DA52uL, 0x78EC8FCDC5CF19D5uL)),
                ((+1, -12, 0xC512BF61537AE819uL, 0xABB773AE74DAA725uL), (+1, -7, 0x97CAB55E46CA43F0uL, 0x4D18FD010C3E61C8uL)),
                ((-1, -18, 0x837EB04DDD597339uL, 0xD09D541E80122DA0uL), (+1, -9, 0x8FB57A7B298FB8A2uL, 0x2DAF069AB91F5A5AuL)),
                ((+1, -17, 0xB81FDBB82EE1F083uL, 0x23AE388A28F90863uL), (+1, -12, 0xF740F3FA3442536CuL, 0x991E26B72EA3A680uL)),
                ((+1, -22, 0x99D044AE58C6ABA7uL, 0xA570D2BAE8E8ECCEuL), (+1, -14, 0xAAF801BA9059E72AuL, 0x1DF5D0BEABAD44ACuL)),
                ((+1, -24, 0xBE1EACC78F55559AuL, 0xEE62C21B81847B32uL), (+1, -17, 0xC79A9D09B778A6D8uL, 0xFFA3C4D574E2BFDCuL)),
                ((+1, -28, 0xCD48E4A0083C72A3uL, 0xB2BB184F8202D354uL), (+1, -20, 0xB49E7CC1CD403000uL, 0xF40DAFB78455EC56uL)),
                (Zero, (+1, -24, 0xE5AF798008B1B77BuL, 0xB5DA46D4E365456BuL)),
                (Zero, (+1, -28, 0xA0498E6F8FB5A3CFuL, 0x95469BCAC304F80BuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_2_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -4, 0xD71F1F4928C283B6uL, 0x6E7CB23276CE927FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -4, 0x9921A05E8B2B4233uL, 0xCA24AAE828868F24uL), (+1, 0, 0xC22240C0DD12FF7DuL, 0x6802D0DE53CA8148uL)),
                ((+1, -5, 0xDE51B484417976ACuL, 0xCBF11B64837A1C1EuL), (+1, 0, 0xAB102F3ED41A23CEuL, 0xD46518EDE62603D3uL)),
                ((+1, -6, 0xF47EDB892A4DBA18uL, 0x340F46BE93EE6206uL), (+1, -1, 0xDF110FC4CAF3F268uL, 0x5CB9411E5E41F2D9uL)),
                ((+1, -7, 0xB367FADCA8610F45uL, 0x11D242FD8DC014DDuL), (+1, -2, 0xE28637897CA49E42uL, 0x7D91ACFB1202DE6FuL)),
                ((+1, -8, 0x96BFE977DDC910ACuL, 0x9897A49F23D6A86AuL), (+1, -3, 0xBFA680646D15DE82uL, 0xFB596EB1094E68EDuL)),
                ((+1, -10, 0x95AF9004BD9360E8uL, 0xDEF28B12E285AC79uL), (+1, -4, 0x87A6E35BFB581048uL, 0x5E097BE544721EEDuL)),
                ((+1, -12, 0xB9798A36883F2CB6uL, 0xE84FADB36814FF03uL), (+1, -6, 0xA549376CF25D7EF6uL, 0xB1B4290FEB69920FuL)),
                ((+1, -14, 0x88560CB6BBA0A7DEuL, 0x632FAA0DA8237F4AuL), (+1, -8, 0xAC30AA7782F88FF9uL, 0x378DC725252E5C8DuL)),
                ((+1, -17, 0xEF627F6F1B47D5FAuL, 0x4576E17C8E8EE83BuL), (+1, -10, 0x9AF6D6C17AA47F21uL, 0x55680DBA063E7800uL)),
                ((+1, -20, 0xFFD159FAFFA58B5DuL, 0x8BED9508DE6830E5uL), (+1, -13, 0xED593994EBBB7AB4uL, 0x5BDF422DD1CE5398uL)),
                ((+1, -22, 0x9491D56680996FBAuL, 0xCC4AF2D5ABCDA016uL), (+1, -15, 0x99B459DFCE429824uL, 0xAE45BD8EB3270BD8uL)),
                ((+1, -26, 0xC46533729C2051D2uL, 0xAE43E5B5D2EFF9BBuL), (+1, -18, 0xA338664E81C9F757uL, 0x8AC228C83374F9D3uL)),
                ((+1, -30, 0xF3C623510FD8DB34uL, 0x39CDE81F90682B6AuL), (+1, -21, 0x883DC79F16699FF1uL, 0xFB42601F9809D4FFuL)),
                ((+1, -36, 0xF6C99AFB550D05CBuL, 0x13A2A1B18E0DFB40uL), (+1, -25, 0xA37C01A80246258DuL, 0x2A3E51787D8E6A99uL)),
                ((-1, -43, 0x975212BFDAE3ED63uL, 0xC5B79AA06EDADB1BuL), (+1, -30, 0xDA842D644F538F0FuL, 0x9A962EB17170F42DuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_4_5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0xFA79610FA0059D97uL, 0xCAC50CED8ABE19C7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -5, 0x9980894A1F01F016uL, 0x8F3BC33B11019E08uL), (+1, 0, 0xD620E650FF8DBAB2uL, 0x6885E288E9E08576uL)),
                ((+1, -6, 0xD0B64D8569C84492uL, 0x209DF9C32514ABE2uL), (+1, 0, 0xB54E44023CAB54E8uL, 0x0482587535E6ECF7uL)),
                ((+1, -7, 0xB7A12E6A6B3E72FEuL, 0x1529F4C0F3483569uL), (+1, -1, 0xC74B444F81DA6B05uL, 0x77871AAE1DD5F0E2uL)),
                ((+1, -9, 0xE64B900D568B95FAuL, 0x0528BB4C42ABF243uL), (+1, -2, 0x9BDD3793A35088A5uL, 0xD8BE2E185F4C1913uL)),
                ((+1, -11, 0xCF59F1B2D4B0E7A2uL, 0x5C597043AD46385EuL), (+1, -4, 0xB41C520FFE5F3137uL, 0x4550C67E9F7B1114uL)),
                ((+1, -13, 0x86328897D3D2FCD0uL, 0xF6CCF85A50FFD9BCuL), (+1, -6, 0x9B4C136A5B8F96D1uL, 0xF87584546ADBB9C1uL)),
                ((+1, -17, 0xEB234AEB04CCDEB6uL, 0x7D0F341390A87815uL), (+1, -9, 0xC563D5ED075CF492uL, 0x00768AE6BBB5800EuL)),
                ((+1, -21, 0xFE30891EB80DAE2DuL, 0x583E72D3F679DEBDuL), (+1, -12, 0xB21C4435950A564BuL, 0x4A92899E668EE431uL)),
                ((+1, -26, 0xF08E3660D2177D8DuL, 0x7517AD591549FEF4uL), (+1, -16, 0xD2C03529996CD69BuL, 0x4B4D34B5CC72776BuL)),
                ((+1, -33, 0x8A82F6DFF6C12836uL, 0xB6816E897C6BAD6FuL), (+1, -20, 0x896179F92D6A93F1uL, 0x0FC3B9149B3006E9uL)),
                (Zero, (+1, -27, 0xE97B70ECC307F33BuL, 0x4EE5A02ED74CB52FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_5_6 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -6, 0xA95237A9C2572CE5uL, 0xB7432DF035BBC062uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -6, 0xDBC0A5E05C67D95DuL, 0x97E8DE474A7FD3A7uL), (+1, 0, 0xD22A72C01A9DC41AuL, 0x318A43645E4BB270uL)),
                ((+1, -6, 0x89E5964F4925C4F0uL, 0xF5FCE72A94EB833FuL), (+1, 0, 0xA3CE6E5828569F65uL, 0x16A1574B2D51AF75uL)),
                ((+1, -8, 0xD282351FFC1B90E8uL, 0xECC932D1BFDBB70BuL), (+1, -1, 0x9D44B398D122EA26uL, 0xBC8FF24981488FFDuL)),
                ((+1, -10, 0xD274DB2549835E9FuL, 0xA125EDBB22AFD342uL), (+1, -3, 0xCB7EB0BD9A13E74BuL, 0x596F3C7038360A36uL)),
                ((+1, -12, 0x8A6A8FC2D35FE5FAuL, 0x67EB04BE2AB26253uL), (+1, -5, 0xB6FD54CE5253BAFDuL, 0xA305DBD03477162BuL)),
                ((+1, -16, 0xE51E5C7A33ACDDDCuL, 0x5C326FF0B84B7ABBuL), (+1, -8, 0xE3D170FA1BA87350uL, 0x1A4E215EB749E770uL)),
                ((+1, -20, 0xD20292BBBFA96CC2uL, 0xE4910C2098D4C647uL), (+1, -11, 0xBD41528991000B47uL, 0x302056E869571064uL)),
                ((+1, -25, 0x9F12DBF3DEB9FDB6uL, 0x63399A25AC616BD2uL), (+1, -15, 0xC1072D61A480BF19uL, 0x362D306DBF4E7B58uL)),
                ((+1, -33, 0xB8D65976662F89F0uL, 0x213C8BE7FB7DF2C2uL), (+1, -20, 0xCAFD1E17279F666AuL, 0x0D76DF1EAC8FCE13uL)),
                ((-1, -42, 0xB0AC4A2E410E4408uL, 0x15D44EC224EA52DAuL), (+1, -26, 0x926D96348B4374C8uL, 0x06F82AF9CC18ED89uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -7, 0xF8A5B186A9627AA8uL, 0x483470D522CE0D0AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -6, 0xA76E5CC31270F994uL, 0xC5B78362C7A8FE55uL), (+1, 0, 0xD001027FAC417914uL, 0x71D53101DCD13D6AuL)),
                ((+1, -7, 0xD0BFFDCF21B61796uL, 0x678CB47D0A1B47AFuL), (+1, 0, 0x9CF9E4A312182E71uL, 0x83554284196031C5uL)),
                ((+1, -8, 0x9B77AB04EC0B68DDuL, 0x6C6A1B67674D30C3uL), (+1, -1, 0x8FE2810F3AABD39BuL, 0x19BF83AA6EEF4E09uL)),
                ((+1, -10, 0x966186B970FDEB41uL, 0xF8A4EEF16B99AB1FuL), (+1, -3, 0xB04EA5B830671376uL, 0xC2884C4E3918DA91uL)),
                ((+1, -13, 0xC0E94C907DF74A07uL, 0xB616DC25DA4EFE05uL), (+1, -5, 0x95FD443B900AFDD2uL, 0xED52D8320E951B45uL)),
                ((+1, -16, 0xA0E9E8CB1DCEE9EAuL, 0x0097259BC43BA313uL), (+1, -8, 0xB25DFB27CA7D146CuL, 0x50FCB431F4036666uL)),
                ((+1, -20, 0xA3CC849CA3773798uL, 0x34FC84587E94952DuL), (+1, -11, 0x9194FCC9978F8AE8uL, 0x752FC2D465BC65F5uL)),
                ((+1, -25, 0xB2BDBAFE9122918CuL, 0x1F41CF30D5FC208FuL), (+1, -15, 0x9BBD4CD99B8CCED9uL, 0xF9F9580C103A8418uL)),
                ((+1, -31, 0xA281DFD20F480B78uL, 0xE3C45CACDBB76E86uL), (+1, -20, 0xC7F3B33DAAA86702uL, 0x69E273C62BB277DCuL)),
                ((+1, -40, 0xF07018E83B0D233DuL, 0x32097599B818DB1DuL), (+1, -25, 0x825A99BFECDE20D5uL, 0x7CF6100195107177uL)),
                ((-1, -49, 0xA20E2570CE7B4929uL, 0xBE0D48DC4B8C4BEEuL), (+1, -33, 0xEC4B71E2284CEF83uL, 0x420D90D6E3E24187uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_8_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -7, 0x9B3930FD3A803B45uL, 0x4BF483743555B3DDuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -7, 0xC6A4C5051720E8AEuL, 0xE327A1104BDD6E1EuL), (+1, 0, 0xBD901AF99D0A6CBFuL, 0x7391C61C6E993151uL)),
                ((+1, -8, 0xE55864F0ACBE34C9uL, 0x5F55D78F173C2C38uL), (+1, 0, 0x8056738E705A4445uL, 0x9BE64C2360A28E26uL)),
                ((+1, -9, 0x9CBE82D2892F7B51uL, 0x6E16A7EAF29EFAE7uL), (+1, -2, 0xD175FB87CB6729B6uL, 0x61CB78B4E6D9BAB6uL)),
                ((+1, -11, 0x8C119E148E66FF15uL, 0x8DF32CDD9A77A4D1uL), (+1, -4, 0xE4CDA6DF55E3E10DuL, 0x14BC3B4675876FFDuL)),
                ((+1, -14, 0xAAC5B74DA8805CD7uL, 0x871735E732C7FC36uL), (+1, -6, 0xAFBE9B9CBF90756EuL, 0xDA4977F3EB877B90uL)),
                ((+1, -17, 0x900C3F40660E661CuL, 0x173191D28434B3BAuL), (+1, -9, 0xC237E814D5A4D079uL, 0x544FD5DA50512690uL)),
                ((+1, -21, 0xA74B384C7B75992AuL, 0xD7D36B628E73B37EuL), (+1, -12, 0x9B61412336CD7BB6uL, 0x9602A1A076A7D69FuL)),
                ((+1, -25, 0x82EB7C103E6D4822uL, 0xAC884EE4EF2E3FC0uL), (+1, -16, 0xB2D13FCE6ED7FA7FuL, 0x5DA14206156BA080uL)),
                ((+1, -30, 0x84BA95F2FFADC674uL, 0xADD5AB25304356B7uL), (+1, -20, 0x91448D7709166E60uL, 0x6C66B5AF22535A4EuL)),
                ((+1, -36, 0xA36C01DC571043EDuL, 0xD84B08A0D54DB571uL), (+1, -25, 0xA13516DE47BE34AFuL, 0xA40AA059A2D512CFuL)),
                ((+1, -43, 0xDB355D0027FEE51DuL, 0x1C36F1E09040A260uL), (+1, -31, 0xE7D8E714315F7E29uL, 0xCD8D57A250988332uL)),
                ((+1, -50, 0x819B279ABC2D9BCBuL, 0x38486899BCB34DD9uL), (+1, -37, 0xC6825690CD4A9963uL, 0xBB635BAF167DCB57uL)),
                ((+1, -60, 0x8D1CB7B49E18B178uL, 0x586FCFADF438D99EuL), (+1, -44, 0xAE414C70E2FD8864uL, 0x73F8C8991082A5DAuL)),
                ((-1, -71, 0xBD00FDDEAB470F17uL, 0x8D40D84020FE3968uL), (+1, -53, 0xDFB2544A28A44BBEuL, 0xC60A7F2B5A17BB3BuL)),
                ((+1, -81, 0xAF8143A01DEDBBA0uL, 0xA428C52660BEE733uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -9, 0xD175D391E97C0CCCuL, 0xA2664757647D6D1DuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -10, 0xA5C66A0A4188D21FuL, 0xCF34C3576756E18DuL), (+1, -2, 0xFBD347BAF43A5291uL, 0xF764EF226038610BuL)),
                ((+1, -13, 0xE0F42770CB3D6009uL, 0xF0D1AA766994C7F3uL), (+1, -4, 0xDA7360429F82AFC7uL, 0xA6105AC65B5E6511uL)),
                ((+1, -16, 0xAB2A96F94BC51B5EuL, 0xCFF2BAACB2C167AAuL), (+1, -7, 0xDBCD205B2A328D47uL, 0x7A764B9A22850AB6uL)),
                ((+1, -20, 0xA08DC6F32941F6E5uL, 0xC09AC734C750F40DuL), (+1, -10, 0x8DF3C028A1712F32uL, 0x5D1D811B1BC63ACFuL)),
                ((+1, -25, 0xC0676D01A9CEEF24uL, 0x92B7839317A6E0FBuL), (+1, -15, 0xF61496854E3BA404uL, 0x65A277478DF1E4B6uL)),
                ((+1, -30, 0x93D49966BF30D66AuL, 0x4BA639B5148981D0uL), (+1, -19, 0x91698B627A2A2F0BuL, 0x1D2D59DB60F8240EuL)),
                ((+1, -36, 0x8E7A7EBE2A3C500EuL, 0x948A85DD8E34BD9BuL), (+1, -25, 0xE96383A5DB5515E3uL, 0x4E90C997776E6E4AuL)),
                ((+1, -43, 0xA3994A23A83E5AA3uL, 0x58CA4125A1706C3BuL), (+1, -31, 0xF8CB9D597E5D60C6uL, 0x49045D73EC61CDFBuL)),
                ((+1, -51, 0xCAC3FC991BBC9288uL, 0x70120919DA59D577uL), (+1, -37, 0xA8929A6D39EF89ECuL, 0xC7A2B23C1B00CEE5uL)),
                ((+1, -60, 0xDD2C6F0C7C0DF7E9uL, 0x47DC14D9EF3D28A2uL), (+1, -44, 0x865F7DBE9598A141uL, 0x08F325B5D2778FE7uL)),
                ((+1, -71, 0xDF04D11BE5CA8CA4uL, 0x963C759093AD07AAuL), (+1, -53, 0xDA673A78CE5D7715uL, 0x19F432644AF2B791uL)),
                ((-1, -82, 0x8B72296409381FDBuL, 0x980934DCEE98B4B6uL), (+1, -62, 0x81C7C1DECAE07AF0uL, 0x030D27E706C5DDB5uL)),
                ((+1, -94, 0xF466FB3BFF1C4AFAuL, 0x74D2396A7F68881CuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -10, 0x91B725448A329675uL, 0xEA0081A9319681A8uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -13, 0xD0D6947DC3236446uL, 0x041835428BD23EDFuL), (+1, -3, 0xE7DFD72B6AC977BCuL, 0x296983D7798D8EB4uL)),
                ((+1, -17, 0xFFF61D601B760A3FuL, 0xD0AAF95BCBF4020AuL), (+1, -6, 0xB8DA06371F91DD4BuL, 0xF8F97247FA2A2006uL)),
                ((+1, -21, 0xAFBBDB72D4740D15uL, 0x3CF7C5074C6D8E66uL), (+1, -10, 0xAAAC1088A3A77C1EuL, 0x1F3979859832CEF5uL)),
                ((+1, -26, 0x94C9DDDA8EB92E7FuL, 0x6CECE60BD69F3A85uL), (+1, -15, 0xCA11B008DB89687BuL, 0x75C9E3A2F39BFA7CuL)),
                ((+1, -32, 0xA12FFFF2B90CD45AuL, 0x0957C5A8941EDD30uL), (+1, -20, 0xA071CCD02932877AuL, 0x77F9BD4CFD47095DuL)),
                ((+1, -39, 0xE0758B18F8F499FDuL, 0xB1BE64FBA6DC444CuL), (+1, -26, 0xADAAEE659B5BA92DuL, 0x3246994199690408uL)),
                ((+1, -46, 0xC4A96A9CCDEA5F19uL, 0x58D69EC3E4A2D246uL), (+1, -33, 0xFF5D2AF6497E9518uL, 0x7A27FEAC5B760207uL)),
                ((+1, -54, 0xCE0B711D023A72E9uL, 0xD54005E2D55CE48AuL), (+1, -40, 0xF9952F79FBB0F99AuL, 0xD84695B736A82042uL)),
                ((+1, -63, 0xE9FC05536F259C83uL, 0x4263B7A32EF94E4BuL), (+1, -47, 0x9B3D48959EB698FDuL, 0x2053D246C324728DuL)),
                ((+1, -73, 0xEAE665C3D269D483uL, 0x93851150E95C75CDuL), (+1, -56, 0xE39A77B7E2DD7041uL, 0x23C4C8E1A53455C6uL)),
                ((+1, -85, 0xDB07AA9BB0E90502uL, 0x1163BFF16FDDC673uL), (+1, -65, 0xAA811D61B8C783D8uL, 0xED855789A321682AuL)),
                ((-1, -98, 0xFE8A0BF4FA220282uL, 0x9F321F7736512234uL), (+1, -76, 0xBB554E7C91A0B48AuL, 0x0A59AB6051E47464uL)),
                ((+1, -110, 0xD055F1D745C8C27BuL, 0x286D5DA5432FE0E9uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_limit = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -3, 0xCC42299EA1B28468uL, 0x7E59E2805D5C717FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 1, 0xC21CD31A673ACB1BuL, 0xC34E86FD057E1F55uL), (-1, 4, 0x86687EB67831EEC5uL, 0x399F0F936F368806uL)),
                ((+1, 4, 0x97F9CE2D16F829D2uL, 0xAAC29B4146F12919uL), (+1, 6, 0xEFB8330899EAD3ABuL, 0xF9AD243B4415CCD9uL)),
                ((-1, 5, 0xD2EC2E421BE53A3AuL, 0x2590BC9452A0B0C0uL), (-1, 8, 0xD16F2BD950B37702uL, 0x7C53FDEF0F455785uL)),
                ((+1, 5, 0xDD5C43BFF82483A5uL, 0x9279B2ABB4F9C045uL), (+1, 9, 0xAF50427AED649AE1uL, 0x2D44497D15F19A7CuL)),
                ((-1, 1, 0x9F715D0912A1BA5BuL, 0xE6237C65490196C6uL), (-1, 8, 0xD2C22A1615C584E0uL, 0xE1CF791EF4FAAA56uL)),
            }));

            public static ddouble Value(ddouble x) {
                if (IsNegative(x)) {
                    return 1d - Value(-x);
                }

                ddouble y;
                if (x <= 0.5d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_0p5);
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
                else if (x <= 5d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_5);
                }
                else if (x <= 6d) {
                    y = ApproxUtil.Pade(x - 5d, pade_plus_5_6);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 6d, pade_plus_6_8);
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
                    ddouble u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u;
                }

                return y;
            }
        }

        private static class QuantilePade {
            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p125_0p1328125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xDB90088D17A85C87uL, 0x554CFEAFDE412BB9uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 5, 0x94462234BDA3EDABuL, 0xBB9C42ADB9127A90uL), (-1, 4, 0xD6A9AFC19720ABF5uL, 0x241F195FE9B3EC88uL)),
                ((+1, 8, 0x88589D274AAF620EuL, 0x36782F7751843B10uL), (+1, 8, 0x8A71EE36AD2899D0uL, 0x3722EA9AC71D8199uL)),
                ((-1, 9, 0x8CAB74C585555021uL, 0x3370F6E69C25E730uL), (-1, 10, 0xA4F9AA499791885DuL, 0x8B9D46287BB600DAuL)),
                ((-1, 10, 0xD8B10C965F66C023uL, 0xDA2BB7B7C872713EuL), (+1, 11, 0xA0988E003A8F755BuL, 0x435CA1E79E16C39DuL)),
                ((+1, 12, 0xD7D549D8431FF0F9uL, 0x7E065540BD49BE5DuL), (-1, 8, 0x8277F37A3FCD12F1uL, 0xA65A8BA821A29B16uL)),
                ((-1, 11, 0x90D96B06D9356D00uL, 0xF7191CC5D8162C4CuL), (-1, 11, 0xC2BC6ED426641D1BuL, 0x24E07A634BBCF91BuL)),
                ((-1, 11, 0xE9234D5615D92338uL, 0xA94F5324671CE8E5uL), (+1, 9, 0x8BE63EC55CFC1766uL, 0xC174D73E50B275CEuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p1328125_0p140625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xD2E33744E0A44F7BuL, 0x3C8C63CB9872AF54uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 4, 0xF34B335F46B8C157uL, 0x79A287A1FBCA7706uL), (-1, 4, 0xBC7AD97D933F8483uL, 0xA25CDBD69DFF9D3BuL)),
                ((+1, 7, 0xB6A08A1D7993E071uL, 0x351D8DCCEBA7CEF3uL), (+1, 7, 0xD305A5E9960DF0EAuL, 0x0DD1A1BFE6F5F67FuL)),
                ((-1, 7, 0xCDFF54748D948CC4uL, 0x7224C3CA4EE2A53CuL), (-1, 9, 0xD4915D5AEAE32A4AuL, 0x38EEE00838790243uL)),
                ((-1, 10, 0xBDE6D6F1BE95CCEDuL, 0xEDD587A1DF5DE29EuL), (+1, 10, 0x9CC85DE78F42D641uL, 0x3568BC047B832FC7uL)),
                ((+1, 11, 0xE935D1987C7DA981uL, 0x70D31BAF99BCD0F2uL), (+1, 9, 0x802A0978A11DE83CuL, 0xD26702D0F4E513F8uL)),
                ((+1, 7, 0xB8EA39B960D37095uL, 0x0183FB00862C545BuL), (-1, 10, 0xD846A7DE18BF37C9uL, 0x5981A7BCC1ECE942uL)),
                ((-1, 11, 0x80110295E9AFA01FuL, 0x92F4CFBB9C52213CuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p140625_0p1484375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xCABA9871574D7163uL, 0x333408B411FBD423uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 4, 0xDAA662E59BFF8A6EuL, 0x120D76441E4573B6uL), (-1, 4, 0xB2132C5E656E679EuL, 0x4D5933078FB59CC6uL)),
                ((+1, 7, 0x95BA2843B2534C59uL, 0x89C09FBA5BF10A6AuL), (+1, 7, 0xBBBD26A1E4EFF057uL, 0xD283578ADD310368uL)),
                ((-1, 6, 0xDFDD26E2C163920CuL, 0x3D22D8E267D41BC9uL), (-1, 9, 0xB12E2E9A47433D57uL, 0xEFBA9B1C34042175uL)),
                ((-1, 10, 0xA7EDFD95D8740599uL, 0x3436C018CACBE43FuL), (+1, 9, 0xEFB4DA95BA3B41ECuL, 0x30B0A7D31B58FF06uL)),
                ((+1, 11, 0xB546B15A591D11D8uL, 0xFBB07AA22731BBE4uL), (+1, 8, 0xF256654F00BE7C1BuL, 0x6E99C53B824F0068uL)),
                ((+1, 8, 0x83A5026278A461DBuL, 0x6F126E8685CE94DCuL), (-1, 10, 0xA732D8F9C3E2DDA2uL, 0xFF0D99B6EB3E8554uL)),
                ((-1, 10, 0xC43507BA8FBC3959uL, 0x01525CEA5CA956B2uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p1484375_0p15625 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xC3053DBD8E45F59AuL, 0xD2C43346F8C0448BuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 4, 0xC88F349FB1E6B53FuL, 0xE04EA19751FF9F4AuL), (-1, 4, 0xAB0738C38D6F29DDuL, 0x7D6DDE8CF57D2738uL)),
                ((+1, 7, 0x807A3D19DCA39245uL, 0x5D1782E38385B8F2uL), (+1, 7, 0xAD02E20EFCDE3574uL, 0xC9055F3918864B8DuL)),
                ((-1, 5, 0xEA9C49712AE1349DuL, 0x923D3B39934078E5uL), (-1, 9, 0x9C3EE9DB7D910BE7uL, 0x3D99C6DD50BA86DBuL)),
                ((-1, 10, 0x99F72F37A3ACEC7BuL, 0x06C193E29D40F2B0uL), (+1, 9, 0xC70AA93EA145724AuL, 0x0C0F6CA488BDFF83uL)),
                ((+1, 11, 0x98F2522971AF9EE3uL, 0x5856CF3F87AF2425uL), (+1, 8, 0xE743C03930B5A0D2uL, 0xEC1F30DE1BB10323uL)),
                ((+1, 8, 0x92F6A318DB926E9CuL, 0xAFF862C813CF5D7DuL), (-1, 10, 0x8DCD72293565CAF7uL, 0x7EF515E1D32B0A5CuL)),
                ((-1, 10, 0xA670BB104FA2BE1DuL, 0x741C85764836F633uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p15625_0p1875 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0xA7B80AF21ED0232DuL, 0x4383461F92CC69B7uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 0, 0x9C9C6BAE7F570037uL, 0x6D519486380AC176uL), (-1, 1, 0xF53B318537383CB1uL, 0x0CDBE976A454CC6DuL)),
                ((+1, 4, 0xE03F6DEA1E79B0C0uL, 0xCD51190A6A2170A4uL), (+1, 4, 0xDBA32F3CA635E8C9uL, 0x5BC2645CA6229230uL)),
                ((-1, 13, 0x882FEED23BE7E703uL, 0xE82868656163C5FCuL), (-1, 12, 0xD4334EE3B2660CC0uL, 0xF70435850F9F5A83uL)),
                ((+1, 16, 0xFF90F01812B51609uL, 0xA9FD3079CEA915AAuL), (+1, 17, 0x80C09BF341EAD2A3uL, 0x5E01E3D0342141B6uL)),
                ((-1, 19, 0xB02A2F989D70A7D5uL, 0x2265519345E0E20CuL), (-1, 20, 0x8606449E0F660A5FuL, 0x32F54645AECA4FA4uL)),
                ((+1, 19, 0xE008CFCBFCB9D832uL, 0xBA12D6A2C71CB251uL), (+1, 22, 0x8EFD339B8A3E6EF8uL, 0x435F61613EC85EC4uL)),
                ((+1, 22, 0xC7565181CC5FB491uL, 0x9C5B7AF774063C72uL), (-1, 23, 0x90365D693D7F27D2uL, 0xCCE1AD0F3DA68020uL)),
                ((-1, 24, 0xCC23B5E29FA7DAE2uL, 0x3274E347899F2B15uL), (+1, 21, 0x80DE3C5B1FE0D88BuL, 0xAFD91DC64C73E83FuL)),
                ((+1, 24, 0xB4D8FE6925AC3195uL, 0x1E0AB1A6885CBA31uL), (+1, 24, 0xB1F3382A4F6168C8uL, 0x489111A73A19B158uL)),
                ((+1, 25, 0x8F8E2D30202AB31EuL, 0x749BEDBD69190D8AuL), (-1, 24, 0xCAC44D315B6A2DD4uL, 0x646398FCBADE7AB1uL)),
                ((-1, 25, 0xE38B8867E8E0AB06uL, 0x3DEC998D4C08CFCCuL), (-1, 23, 0xA7A6FAC1CB2C34F8uL, 0x6AC7DBE374B44E75uL)),
                ((-1, 23, 0x84339F9CA4B6F5CEuL, 0x64FB49EC23ABEDC3uL), (+1, 24, 0x9E6F622334AEB5FCuL, 0x35198B0DDE8AAB33uL)),
                ((+1, 24, 0xBA8B3EABC6066627uL, 0x2574A8C77218CFA1uL), (+1, 18, 0xE1C8B15DFE2D6F7EuL, 0xA452A2779F783070uL)),
                (Zero, (-1, 20, 0xF737A3C35216338DuL, 0xF958EBDEFE8D779CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p1875_0p21875 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, 0, 0x907F0A2CA11CA7D5uL, 0xDE283D13A1288E32uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, 0, 0xF9CFA22722DDC6D9uL, 0xB49F72CBB83DCD51uL), (-1, 2, 0xD0F8D063528694C1uL, 0x39D9D464C6B2801BuL)),
                ((-1, 7, 0xA5C7D0651947940CuL, 0x5634747CA79A492EuL), (-1, 6, 0xFA327EFC9FAC82B2uL, 0x9090E0B41B851A12uL)),
                ((+1, 10, 0xA04968028B95F048uL, 0xC8F234B061E61647uL), (+1, 10, 0xDC4432BFC4BFC904uL, 0x25FD384DF8BB7F16uL)),
                ((-1, 10, 0xC760E6CEAED359DFuL, 0xC8F69529FA2E82FDuL), (-1, 13, 0x85FBFE3A87BAA935uL, 0xBB56BF27FB1596DFuL)),
                ((-1, 13, 0xE0A3BDF2EFA367F2uL, 0x92A68EBD0DBF461CuL), (+1, 14, 0x80EEC28B8114AC36uL, 0x09E7173DACA7D423uL)),
                ((+1, 15, 0xCC61A1EFDC6D64C2uL, 0x0BACAF469F227DF4uL), (+1, 11, 0xD738F40EEDACBDC4uL, 0x3481163E9DB8A7F4uL)),
                ((-1, 14, 0xD3787D9EFBFA5E4EuL, 0xC28E31803FAF4863uL), (-1, 15, 0xD24F0761267FD345uL, 0x97448A47D842EDBCuL)),
                ((-1, 16, 0xCD4CF2489E57599AuL, 0xF3C6823B9941BE77uL), (+1, 15, 0x9CF65D113C33825BuL, 0x0A55ACBA81CD08EFuL)),
                ((+1, 16, 0xC7FFBF9D8349EA12uL, 0xBF0CA5A37FE082B0uL), (+1, 15, 0xAC03778FF3905CF2uL, 0x5E4408A78C159FC5uL)),
                ((+1, 15, 0xEAA634CAEBDA6D78uL, 0x638BE37943728A2EuL), (-1, 15, 0x92306631CC738F45uL, 0xE145390441FC5F1BuL)),
                ((-1, 15, 0xA9D35AF15BC9F6B9uL, 0x73B1F1266732096CuL), (-1, 13, 0xBD37A8E3A487C858uL, 0x59119EB175FDA696uL)),
                ((-1, 13, 0xBD685AE81DA9C204uL, 0x17C3A4EBA244343DuL), (-1, 9, 0xF9F716D823459BA3uL, 0x6A208A8E0D2E01EFuL)),
                ((-1, 13, 0x9FEA5A3BF4C7573DuL, 0xC1C8FC1BF45CA47EuL), (+1, 9, 0xEFF44C8C9F518780uL, 0x9C7756250CCB4FC1uL)),
                (Zero, (+1, 10, 0xC7C17A709C48B482uL, 0x7FC79A7E9D2F856FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p21875_0p25 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xF80C01477814982DuL, 0xCF420B8BC042CDF3uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0xFE382E2EA29E9B28uL, 0xAA33C4069E554656uL), (-1, -1, 0xE7950240AADDFE91uL, 0xFF89FD3AA311AE81uL)),
                ((-1, 7, 0x85CF889AC3E0FAC0uL, 0x7FA2844314F4B3EAuL), (-1, 7, 0x8D82DC5B02AABF74uL, 0x983677B4EAA74B60uL)),
                ((+1, 8, 0xE427B460A9DB41C3uL, 0x36F5EB33CAC7BE10uL), (+1, 10, 0x90D4C8AE76EB8624uL, 0x34C071644389E30CuL)),
                ((+1, 10, 0xA7217BA711753651uL, 0xE0556A1EB749FEBBuL), (-1, 11, 0xD2164B25E717732CuL, 0x784E3F03AE9ECB3AuL)),
                ((-1, 13, 0x9035620C33CDBD31uL, 0x3A782937090030B1uL), (+1, 10, 0xDB7B11F906A7F43CuL, 0x2D0CF569D09A1D6BuL)),
                ((+1, 13, 0xA4D64075F85F990DuL, 0x7695B1A6E8088D26uL), (+1, 13, 0x892A24AFE9CB0C02uL, 0x499A3F04794AC5AFuL)),
                ((+1, 13, 0xE85E2E2DDAE316ABuL, 0x1DE62B6AB1F74277uL), (-1, 13, 0xB9FABD1521966D72uL, 0x765ABE1843F4D65AuL)),
                ((-1, 14, 0xD411884A403BF65EuL, 0x17C2F4D4ACE5763FuL), (-1, 12, 0x94AF87C1B408CC6BuL, 0x5EBC4BCB1D1BCEE7uL)),
                ((-1, 12, 0x80DCB094C2DD92DFuL, 0xD433B9A7D633EABAuL), (+1, 13, 0x9D9BC7B473C48BADuL, 0x31F62B31E8BD23A5uL)),
                ((+1, 13, 0xC1D9A1C160EBD874uL, 0xFD145ACDA9E98DB4uL), (+1, 7, 0xFE43932903A8C510uL, 0xFD1C5199742CE546uL)),
                (Zero, (-1, 10, 0x887D94F0CAC006F2uL, 0xFE7051F483109E7CuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p25_0p3125 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0xB0C481C57F835A9FuL, 0xB633E969709A1187uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0x98004C46A78637CBuL, 0x7C61BD3E2410339FuL), (-1, 1, 0xA303B2598C44863EuL, 0x75B0A1F3FC54CCDAuL)),
                ((-1, 5, 0x88C3E45A07993C81uL, 0x6C2FC3F66887331DuL), (-1, 5, 0xA132A760746CEC82uL, 0x4C2F46F89CE98004uL)),
                ((-1, 4, 0xBF173A1B44035E68uL, 0x0F6E8E0A0AA0BB36uL), (+1, 7, 0xCAFCD1F355059846uL, 0xF33E22C702D4B326uL)),
                ((+1, 9, 0x8ACF5C47FE79938FuL, 0xFD26F70DEF09E221uL), (-1, 7, 0xA589FFC99C75371DuL, 0x7017E616B689A688uL)),
                ((-1, 9, 0xBC56E4B6F4AD1CBFuL, 0x5F6F5B7F39655AE6uL), (-1, 9, 0x904412B5B1009CDEuL, 0xD5F9F1B360EBEA35uL)),
                ((-1, 10, 0x8539ED0F80C8013FuL, 0x8BBF88E3E27EE790uL), (+1, 9, 0xD4E02914A0421BFAuL, 0xED0728F746C954C4uL)),
                ((+1, 10, 0xF5A2B69E02A2CE74uL, 0x0B3343A1FED9084FuL), (+1, 8, 0xC0498480A88A920FuL, 0xBA7DDFD5B702EB65uL)),
                ((+1, 8, 0xC639908A0D4B04C8uL, 0xB59F770D87012BEBuL), (-1, 9, 0xBEDBFD77798BF3E6uL, 0x1F1ABC117DCC1D21uL)),
                ((-1, 9, 0xF13CB302F046DB34uL, 0xB4AA0BD9A86D1A92uL), (-1, 4, 0xF900ABBEF659A968uL, 0x148850E00927DE04uL)),
                (Zero, (+1, 6, 0xAFBE51503EE4587DuL, 0x4607F30B3009EDB6uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p3125_0p34375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -1, 0x90A4902D23E8EE0FuL, 0x86C5C68A9644E2D6uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 4, 0x81003AD8B31871FFuL, 0x8897A3E856833764uL), (+1, 4, 0xACF9770700C1BA7DuL, 0x60EB921D276BB457uL)),
                ((+1, 4, 0xAA8A222E1125D295uL, 0x0B2A4BB03AF54E34uL), (-1, 6, 0xEAA1FF924C26939BuL, 0x6401D8D61A0BABE4uL)),
                ((-1, 8, 0xB68A5A10A78A8C81uL, 0x2F2FD992E962839EuL), (+1, 4, 0xEB71014E20E66A03uL, 0x969F2AE70DDC2806uL)),
                ((+1, 7, 0xE414C9A0ED98BAA4uL, 0x66B591DBFFCE4A1EuL), (+1, 8, 0xD49FC5E4174569EDuL, 0xF459534B4849C90BuL)),
                ((+1, 9, 0xECCDB6613431C712uL, 0xC9AE6693FCD6F200uL), (-1, 7, 0xE83F3F63EB359B1FuL, 0xF132047671141561uL)),
                ((-1, 9, 0x82784FE2B0857DF2uL, 0xD8F0DCF0451CEB62uL), (-1, 8, 0xCAF7A1582FF2DBA3uL, 0x315F4F245EA52A48uL)),
                ((-1, 9, 0x88A9055B5B3037DDuL, 0x1249707FE28024D3uL), (+1, 7, 0x8C3E9E7466284648uL, 0x054E2D098CD2BE43uL)),
                ((+1, 6, 0xE081F948546EC2C9uL, 0xED568D9CE07C80A3uL), (+1, 5, 0xE15BEC2E94CE466BuL, 0x89809FBE37621D4BuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p34375_0p375 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xE426054CB845E48BuL, 0xDE602BE4EF6BFF1FuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0x8ACC2D9852B2ABFBuL, 0xBCEF03B5123CE1D5uL), (-1, 1, 0xE1FF7A419A63C33AuL, 0x3100AA253D107048uL)),
                ((-1, 3, 0xC5ABDE2B3CE7108DuL, 0x03C0A97FCE287BD4uL), (-1, 1, 0xC9455A7D5BAB9565uL, 0xC7AED03A812779AFuL)),
                ((-1, 2, 0x86C14A529BA779EFuL, 0x33F01360470988B6uL), (+1, 4, 0x88F31B3216134F0DuL, 0xE717DF956D62F5D6uL)),
                ((+1, 5, 0xAAFBDBC1BC56D233uL, 0x825C576C286A1E67uL), (+1, 0, 0xC71B5D1A145EFBE7uL, 0x0CE099094D76FE5EuL)),
                ((-1, -2, 0xF1421D8077A39565uL, 0x72BD4557377C7474uL), (-1, 4, 0xBDC53AD899FD5043uL, 0x8B024036D52A434DuL)),
                ((-1, 5, 0x9D9D764C0A5CBAA0uL, 0x36CEF9CABB42A4E1uL), (+1, -2, 0xFCF59CF622B87627uL, 0x77ECADFC8A1B225EuL)),
                ((+1, -1, 0xB57709BDE14CD806uL, 0xA419A07CC44E5B39uL), (+1, 2, 0xFE1D487911CD632CuL, 0xD3918D2248A56B97uL)),
                ((+1, 2, 0xAB29552A406353B7uL, 0xA5DAAC5D094FD566uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_0p375_0p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                (Zero, (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, 1, 0xDEB907476B686EE5uL, 0x8639B974191ECE31uL), (-1, -1, 0xF36EB70D36023FF0uL, 0x42F0B3EC04E30C44uL)),
                ((-1, 1, 0xD3C9F4ABAE3AB127uL, 0x50518EEA99191096uL), (-1, 3, 0xFA86E10E8A0A8F5AuL, 0xDBB9F6780A6DBEFAuL)),
                ((-1, 5, 0xC5362F1EF8799D76uL, 0xB78D3618C3FEC973uL), (+1, 4, 0x8576B4754D6CD3ECuL, 0xD27CF2BA573241A7uL)),
                ((+1, 5, 0xD47FB853866FCBFDuL, 0xDB55151B5611413BuL), (+1, 6, 0xB691D841E27C5291uL, 0x78147FC5A36870C5uL)),
                ((+1, 7, 0xFAEC1916058A921EuL, 0xD9AF01E8CE1443C6uL), (-1, 6, 0xCF99C011A67D356AuL, 0xDABCDBDB3BC41CCDuL)),
                ((-1, 8, 0x903977373B2FDF8AuL, 0xA5F0AD64548A306CuL), (-1, 7, 0xF7A28071F84C7B88uL, 0xF6956C6BFDA545E9uL)),
                ((-1, 9, 0x8DD0648FDD8042CCuL, 0x1186C3970401CA04uL), (+1, 8, 0x919BFBAF884FB067uL, 0x5EC030BFBFDBCF9FuL)),
                ((+1, 9, 0xA7CB3D7C822131D2uL, 0xC0BE4794B0A0485CuL), (+1, 8, 0x9E1D010BF01A14E2uL, 0x8D2DB9370BBB3417uL)),
                ((+1, 9, 0x8C0D6CD00A35015FuL, 0xE4C3567A883936C6uL), (-1, 8, 0xBBA8970C966A7D5AuL, 0xD9CCA352AD6BD6E6uL)),
                ((-1, 9, 0xA5F5F02560D26DF7uL, 0x8A7547C061118F58uL), (-1, 7, 0xA6CB384B0B597B77uL, 0xB0AE33E25919D51CuL)),
                ((-1, 7, 0xC40E809556A1A4A0uL, 0x9F05673569022403uL), (+1, 7, 0xC32E179E7E6124A1uL, 0xB70341B58BAA4864uL)),
                ((+1, 7, 0xE278864EDD20B612uL, 0x524E631CD53C8A25uL), (+1, 4, 0xBF369EE6A3F64BA0uL, 0xD16E107BF5DCFEC7uL)),
                ((+1, 3, 0xB17214D3B47E55EDuL, 0x74B22EEF250E7C1FuL), (-1, 4, 0xD691691F9611D1A3uL, 0x1D1CD7298DE664E1uL)),
                ((-1, 3, 0xC0EF6CAC3327C1C2uL, 0x9B4D6C30028CFF6DuL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm3_3p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xE4D5C09E32CAE547uL, 0x155F21CCC573F8F4uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x85BAE68E7A5B7C0CuL, 0xB687E9E7B51BBBC4uL), (+1, -2, 0x93ECED2120EA5EC0uL, 0x83A6D6A8C6581E9EuL)),
                ((-1, -4, 0xA9B0DFC3E7ADDA24uL, 0xC75587721BB197F5uL), (-1, -4, 0xD43F14B3C18DFC3CuL, 0xCC1A9B59474755CEuL)),
                ((+1, -6, 0x9BC2A138BA03746BuL, 0x7B0A00A6D1BA9DB0uL), (+1, -6, 0xDE55035379231B9DuL, 0xA8194FF372BA9424uL)),
                ((+1, -7, 0xD51E68F3C190C037uL, 0x2E970120D35D5917uL), (+1, -6, 0xA0597A9E58571C95uL, 0x91C8704377C62BE1uL)),
                ((-1, -9, 0xDF0172407188F3D6uL, 0x337A6F60BE16D759uL), (-1, -11, 0xDC6D3DD8D783E712uL, 0x20AEDF5842F44225uL)),
                ((+1, -12, 0x9011B6CD8DE847EDuL, 0xC7AD69A6C5076D0CuL), (-1, -17, 0xF0D9DA72F4BC1065uL, 0x5338532F52195AC7uL)),
                ((+1, -13, 0xD0C4C87A71E409A5uL, 0x1A295FF7DE3CD7E5uL), (+1, -13, 0xBBB518525CFF6018uL, 0xF9E7C64AB50BC0C5uL)),
                ((-1, -15, 0x93BB2C7E102596E1uL, 0xF23209B14B68BD48uL), (+1, -16, 0xAC4302DCA09A744CuL, 0x0C80AF13AA87B4A3uL)),
                ((+1, -18, 0x97A5F413E5B1EF21uL, 0x94A44DC56726264AuL), (+1, -20, 0xAF376A68A265F341uL, 0x1E2BFBE0D09EEFC8uL)),
                (Zero, (+1, -24, 0xC26127934FD625D1uL, 0xAAB6809CF490A2CEuL)),
                (Zero, (+1, -25, 0x832A5995030738C8uL, 0x7E906ACCB72307A0uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm3p5_4 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xE17A5C974382201CuL, 0x2F2B89A84B0DE1A5uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -3, 0xB79FAB5530B716F0uL, 0xBC9DB43C7796C6EFuL), (-1, -2, 0xB4457FA24CA47A5BuL, 0xB0D00AB14CE89D8CuL)),
                ((+1, -5, 0x973FE50111AC62C1uL, 0x0B4DFF2DAE5CE138uL), (+1, -4, 0xCFB8B43242E8B8DEuL, 0x817E2BE837EC8586uL)),
                ((+1, -5, 0x83447865097C6C9DuL, 0xE4538A17164C1781uL), (+1, -5, 0xA9C99EA13FA9DF04uL, 0x7C7681C025F6D5D4uL)),
                ((-1, -7, 0xD5769FA5CA1742E5uL, 0x6250569FBF4511E8uL), (-1, -7, 0x930239164BFB5024uL, 0xBF3BD3BD6A754BFCuL)),
                ((+1, -9, 0xD566546FC863029FuL, 0xFDBD2F054DE551F6uL), (+1, -9, 0xAB405520B075589CuL, 0x93B5AAE7C100B8C6uL)),
                ((+1, -13, 0xC1C8BE554DB49F38uL, 0x47DBF7AA6740ACA3uL), (+1, -11, 0x9F72E92113397BA7uL, 0xBB455F7E150A1F5DuL)),
                ((-1, -13, 0x8FB31851DF5955F9uL, 0x28ED651B45B62F8DuL), (-1, -18, 0xE3DBCD51A16C020AuL, 0x86C72F56D867E701uL)),
                ((+1, -15, 0xAC1BE5C237C625E0uL, 0xFB4285B9E8FD34ABuL), (+1, -17, 0xCDDBAD0BF54693C9uL, 0x14296FC7E17348C9uL)),
                ((-1, -19, 0xF3BB68C39E073FABuL, 0xC2F10C4D72753574uL), (+1, -19, 0xA97D4CFE3930D91DuL, 0x9D52E64576A54A43uL)),
                ((+1, -22, 0xA87EBB8AEC5375C1uL, 0x85AF04BA40AE1F4EuL), (+1, -22, 0x90257E9634356AB3uL, 0x4E54F89ABAE9F58FuL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm4_4p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xD9C6C01833F2F4D9uL, 0x2A3A95F640266D5CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -4, 0x8DE63E7DF273A7FBuL, 0x41936D476F66EEBCuL), (-1, -4, 0xAA894235C3BBB81CuL, 0xBAABB4EEF042BFA9uL)),
                ((+1, -4, 0x8B4FA97952A17670uL, 0xFB39482F41A572C1uL), (+1, -3, 0xA6A317364B68705DuL, 0x4027282445E3CE41uL)),
                ((+1, -6, 0xF845FB5CE312F20FuL, 0x2DFCFA4184AE77D8uL), (+1, -4, 0x8D3E95467BB8CF93uL, 0xC1966FEF5EBB5450uL)),
                ((-1, -10, 0xEF10A6192381B439uL, 0xACD08204514B7BBEuL), (+1, -7, 0x84FD3511E17A7DC3uL, 0x43886579B5865104uL)),
                ((+1, -8, 0x81AD6FD38C60E9BAuL, 0xC38AAA9F887F7642uL), (+1, -8, 0xEFA6624B5F2865FBuL, 0x76D2358DFED385CFuL)),
                ((+1, -11, 0x85BD9CDDB0A9DA3BuL, 0x6D809E4F05ABC693uL), (+1, -10, 0xE333B9ADCA5A94EEuL, 0xEDCDFFD013FBFC35uL)),
                ((+1, -16, 0xB03A34AE4A3BB179uL, 0x015706BA27CECB50uL), (+1, -12, 0x8FFDC0E9514B9786uL, 0xBE51001057E1F720uL)),
                ((+1, -15, 0xC7D2F070041BFE45uL, 0x4404A81612247AADuL), (+1, -14, 0x8922BDB3BB2584F4uL, 0x68E641423F9510F7uL)),
                ((-1, -21, 0xBDC775554C69BC11uL, 0x1D984AD97A79E7AAuL), (+1, -17, 0x96B06B7E9F23E7C2uL, 0xF8C1BB4F9460AB0AuL)),
                ((+1, -21, 0x83A31444E83608B4uL, 0x94AF26B26180FF71uL), (+1, -21, 0xCD36BD82B89473B9uL, 0x47B173859E238414uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm4p5_5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xD0EEBAF10D0A5E6EuL, 0xB3ABFA48CD5DC452uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((-1, -4, 0x8E3193DB8E41D894uL, 0xC4ADD1BEDC5C2760uL), (-1, -4, 0xB15DAC0077D4E1C7uL, 0xC9882D314DCD2FB9uL)),
                ((+1, -3, 0x8FAE0BAFE2A3479DuL, 0xF3194CCEBF70D0AAuL), (+1, -2, 0xA8CFF9DAAA215D5AuL, 0x685EF5292E2A835FuL)),
                ((+1, -8, 0xF9FB0E82F3A9689CuL, 0xBBBAE308FCE2792AuL), (+1, -5, 0xA09BB8DC80D1F985uL, 0x873F37C66F6355DBuL)),
                ((+1, -7, 0xB1AF183F34B18709uL, 0xCAB7AFD7F2A735E6uL), (+1, -5, 0x80B2B9FB26CD2929uL, 0x0C8FDAE41E98C9F7uL)),
                ((+1, -8, 0x8D4BD6872204BBE1uL, 0xD035C513C6FF5FA3uL), (+1, -7, 0x9F78013346805967uL, 0x349CA4D87CF728A6uL)),
                ((+1, -13, 0xA0E0FF0A8EF88D57uL, 0x0FC236E0F35C0167uL), (+1, -10, 0xE659AD32DE54331EuL, 0xA9F72F37806D86E4uL)),
                ((+1, -12, 0x87B7EB29BF43A3B7uL, 0x65163F898AD57EB3uL), (+1, -12, 0xE48AB62CB9A47ECEuL, 0x4957FD90F96AC2CEuL)),
                ((+1, -18, 0xA91B13791F06811CuL, 0x6BEA934F086DA65BuL), (+1, -15, 0xFB82AB477DC1E324uL, 0x2730B611A95FD8A8uL)),
                ((+1, -19, 0x94BDF4ADF4AD2F56uL, 0xC6526A7DB264A903uL), (+1, -18, 0x99A0B497E12F6CDDuL, 0xD82DC08F0AE80459uL)),
                ((+1, -23, 0x8702362C711AEB2AuL, 0x7A38046DB093A371uL), Zero),
                ((-1, -26, 0x818B76EE5A641DEFuL, 0x805D1BD16FE51773uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm5_5p5 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xC8B9DB7F28D4B6AFuL, 0xA7D3DA6692EFC105uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x96FB1C9BB7EE9AB7uL, 0x19A8C307F5643FDFuL), (+1, -2, 0xE72265CA4FABCB38uL, 0xB1F7CAE948D1F1FCuL)),
                ((+1, -3, 0xF16A596E0926A8DDuL, 0x7CDA965582761F0CuL), (+1, -1, 0x9EE79354B22A031BuL, 0x26EC5EC7F435B406uL)),
                ((+1, -4, 0xB696DEC9E898E66BuL, 0x4921151647E3E7A0uL), (+1, -2, 0x87C449E56115825EuL, 0x1C013F4F7CF70EA4uL)),
                ((+1, -5, 0xDE0B7CA2896E4799uL, 0x4261724A722E3F07uL), (+1, -3, 0x99E01B9AE16D028CuL, 0x72684396AF98E404uL)),
                ((+1, -6, 0x96C6DC01AAF5C94FuL, 0xDC7DDB4179DAD0C7uL), (+1, -5, 0xDAB023E2B95BBBE9uL, 0xB4C2C0490529B19AuL)),
                ((+1, -8, 0xB9D373C4015BBC10uL, 0x343F5B3728DD2B4DuL), (+1, -6, 0x8959FAA77F147C1BuL, 0xCCBF06E1B6252B10uL)),
                ((+1, -10, 0xBFF94326F798FFEAuL, 0xC9D1EBCFBBE619B2uL), (+1, -8, 0x87BDC94AF7FB8E6DuL, 0xE341F33C9455325FuL)),
                ((+1, -13, 0xF3488699D4AF584FuL, 0xD2AF109D71126180uL), (+1, -11, 0xC08F4990A4499959uL, 0xE2879D85437D8C96uL)),
                ((+1, -15, 0x8F93DE7413A09269uL, 0x0C440A759A12F7BBuL), (+1, -14, 0xBF8FEF557B114482uL, 0xE7EF653173611110uL)),
                ((+1, -20, 0xF391994B464B3DD2uL, 0xF9FE1F3E096B3D9EuL), (+1, -18, 0xCBBC7BBE86A1AEE6uL, 0xB80B19359EA0D740uL)),
                ((+1, -27, 0xAEEAB8ECF75FC2ACuL, 0x7207E77107E8CB42uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm5p5_6 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xC1EF53E5C4A19BD7uL, 0x18B450B545B228BFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -2, 0xC932829990618F75uL, 0xC75625C935CFD0CDuL), (+1, 0, 0x8CB88A580BFDDA65uL, 0x5BDE60472A9E2F45uL)),
                ((+1, -2, 0xD8C7354BFD94AF6AuL, 0x6BFA12A010F3CC1DuL), (+1, 0, 0x95BA5311AE13FA9FuL, 0xF501FB48C1751592uL)),
                ((+1, -2, 0x847B52A6A36A0449uL, 0x00A0B4C33C636CEFuL), (+1, -1, 0xBD15F2E78B2C1B14uL, 0x40FDFE2970FCAC32uL)),
                ((+1, -3, 0x8E1FBBBDBD2F4109uL, 0x8EBBDD4B04AC0582uL), (+1, -2, 0xCA54E72B118A09D0uL, 0x5A7CAE3941801D51uL)),
                ((+1, -5, 0xD5E06A20A7B00AA6uL, 0x358699D34F3DD764uL), (+1, -3, 0x9BA9E9F04D87F84AuL, 0x8994440BF914D8D2uL)),
                ((+1, -7, 0xFFB987561FF98A83uL, 0xC0B41D4FDEFF4E53uL), (+1, -5, 0xBAD01F6849EE3DB6uL, 0x57B9F39F80DE7EBEuL)),
                ((+1, -9, 0xD7D90023C30AD315uL, 0x9E8F22269226A377uL), (+1, -7, 0x9F1D08482DB5F476uL, 0xB4F91E453FB487F3uL)),
                ((+1, -12, 0xEA948320B78510CBuL, 0x49D57EAD7D785A41uL), (+1, -10, 0xAEEAF757CE6732CBuL, 0x061DB8A330B8E074uL)),
                ((+1, -15, 0x99C361BD38899979uL, 0x67238A60434DFE20uL), (+1, -14, 0xCF0309C18B8AB59BuL, 0x01D4C83ADAD1E37CuL)),
                ((-1, -22, 0xBB3A327612D16E69uL, 0xCAF59E5FBAD07E9DuL), Zero),
                ((+1, -26, 0x876AD1442511E8ABuL, 0xF0DF65E80F632D97uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm6_8 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xBCAEB5E59DF3D830uL, 0x5E3ECCF617AF5504uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0xDEACA1921DB9962FuL, 0xADF3C151AB479EC6uL), (+1, 1, 0x9A230493F4233E99uL, 0xAEE39B70E1B8DA62uL)),
                ((+1, 0, 0xAC61EB398FF71129uL, 0x715D1DFD61CC7ABFuL), (+1, 1, 0xF061835007A36EF5uL, 0xF38EED5F1C90286DuL)),
                ((+1, 0, 0xB9AB0832A0CF2898uL, 0xD9F082AB09212BA8uL), (+1, 2, 0x82B221A8FF9D00C6uL, 0x02B7642FC3172790uL)),
                ((+1, 0, 0x9C62710DFF02A24CuL, 0x78DA7627928720F2uL), (+1, 1, 0xDDA1775E0BE87FBAuL, 0xE865155D71A0CEB2uL)),
                ((+1, -1, 0xD0AED6FD35A3A19AuL, 0x23F04C18F311949CuL), (+1, 1, 0x94FC1184B36470E4uL, 0x0B7F529260B49C46uL)),
                ((+1, -2, 0xE2C45A40EC297A82uL, 0xE5D3AA3FD9EE0D79uL), (+1, 0, 0xA2E92DEE4EE0BAC6uL, 0xDE3532A909EADE59uL)),
                ((+1, -3, 0xC807FB0CEA0AC69DuL, 0xA0DA6AD1D78DE72CuL), (+1, -1, 0x90AA23A43800C244uL, 0x580FC0E2F9634693uL)),
                ((+1, -4, 0x8F8035D3E05F037FuL, 0xB80DFFFFF01DD7E2uL), (+1, -3, 0xD0B9D04F15664C79uL, 0xB335BB7627C3C351uL)),
                ((+1, -6, 0xA4DC3CDCA95F4D55uL, 0x482DC45FA7A16528uL), (+1, -5, 0xF1111A13FD18CDB4uL, 0xE89A743348F8A380uL)),
                ((+1, -8, 0x95AFD811D8D2FE8BuL, 0x1619D8A6AF84FEFEuL), (+1, -7, 0xDB8844D84873A65EuL, 0x9ACB140DC21E0350uL)),
                ((+1, -11, 0xD0FEF840645D2A7CuL, 0x8868999D225CF858uL), (+1, -9, 0x9962B98A4C7C48DCuL, 0x4E7E71B8B65718F3uL)),
                ((+1, -14, 0xDC1F2A48F57F1CAAuL, 0xE4911975A048B525uL), (+1, -12, 0xA0F92132C087DB36uL, 0x77C42C931817E4B7uL)),
                ((+1, -17, 0xA993581845D0C5FAuL, 0x5806E5C34920C741uL), (+1, -16, 0xF7A0DDC6CEC3BEE2uL, 0xDE0407B75AA079C1uL)),
                ((+1, -21, 0xBC9CF039AD4AB594uL, 0x21E132E27B13E8F9uL), (+1, -19, 0x8B366DCEB934D777uL, 0x713551302C5C1818uL)),
                ((+1, -25, 0x908E5515F236E075uL, 0x10AE80A13074680FuL), (+1, -24, 0xD136008A8E231958uL, 0x9DF0648C21624BC7uL)),
                ((+1, -30, 0x95F0595A2F465C81uL, 0x5E8B22E6628743F0uL), (+1, -29, 0xDE7A7D5969DEEA4FuL, 0xFA778D9A5FB1E397uL)),
                ((+1, -43, 0xF50A0FD06F659A86uL, 0x9F1424A6D45177A3uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm8_12 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xB265C5F2C8B85237uL, 0x5C819DAD8F75B4CCuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -1, 0x99C831C6B8EA292EuL, 0x477B313F6C18A58CuL), (+1, 0, 0xDE73514FFA38FAA8uL, 0x08CA446D8D7A9DC9uL)),
                ((+1, -2, 0xF19DD6314B1FC913uL, 0x015BA2DF8964F0A4uL), (+1, 0, 0xAFD7577B85692209uL, 0xE172247A21CA7803uL)),
                ((+1, -3, 0xE4A97270A623FE03uL, 0xB2EC3D27D75DA862uL), (+1, -1, 0xA71D0628EE0C4F52uL, 0x4C37D03784AF3CD6uL)),
                ((+1, -4, 0x91E9B25BEC0CC082uL, 0x54486BB86A65112AuL), (+1, -3, 0xD5B77DC00A9F82D1uL, 0x11164026E590D4D9uL)),
                ((+1, -6, 0x859E2EEA83E2A8BBuL, 0x5EA5720ACD566594uL), (+1, -5, 0xC3C642B6F98AC97FuL, 0x06E839734E66168DuL)),
                ((+1, -9, 0xB6B7DE22CC9051A8uL, 0x28BBECEB2166FE87uL), (+1, -7, 0x85CE53CF413A915BuL, 0xDED42A34FE6D8F71uL)),
                ((+1, -12, 0xBF8E65B8EBF166D8uL, 0xB24E170EA26C5B0AuL), (+1, -10, 0x8C4245A93976F8A6uL, 0x44C265C48E4F58DFuL)),
                ((+1, -15, 0x9C96B692F510AEA8uL, 0x30C795834BA4724EuL), (+1, -14, 0xE556AEB03B942389uL, 0x3799ADEBDEBB53BFuL)),
                ((+1, -19, 0xC968A5ED2A967DB5uL, 0x0AD62BBB672335A5uL), (+1, -17, 0x937E3DCBE8225CBDuL, 0x03E3AF74914E1D53uL)),
                ((+1, -23, 0xCBD5E3C24BDF7539uL, 0x461425DE654A764AuL), (+1, -21, 0x954BCE75A326670AuL, 0xE1ADCE1C8092D464uL)),
                ((+1, -27, 0xA0C102C4B67094DDuL, 0x046FEE747EA71558uL), (+1, -26, 0xEB4FAAF5EA4B6EADuL, 0x2AC05B4837E5B387uL)),
                ((+1, -32, 0xBF67496CA58423CCuL, 0x9A3A804F1E1BBC57uL), (+1, -30, 0x8C4B1344EF57E098uL, 0xE43583F6C12D6453uL)),
                ((+1, -37, 0xA17078430C41EDA5uL, 0x6FF1A47119593518uL), (+1, -36, 0xEC25DC3369150F94uL, 0x3528F2805F0F0AB2uL)),
                ((+1, -43, 0xA1121433B8E09E64uL, 0xD332162496A807D2uL), (+1, -42, 0xEC471095F8541D3FuL, 0x9ED75C812B8A15E7uL)),
                ((+1, -59, 0x935B889C96AEC17DuL, 0xF3EAFB75455B41C1uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm12_16 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xAF052A5E2FCD244BuL, 0x3F9AC983661C623AuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xE913B1B84EE56174uL, 0xB58EB9700A045646uL), (+1, -1, 0xAAB0A0110789E99CuL, 0xE37EA95988F1E5C1uL)),
                ((+1, -4, 0x8F79C67B2ABC22EFuL, 0x8800922A2682F30BuL), (+1, -3, 0xD227E2059B5351B2uL, 0x74AF8B0A8F6FCA55uL)),
                ((+1, -7, 0xDA232E11C0A6290BuL, 0x35B5CA78017FE5B6uL), (+1, -5, 0x9FBBC419DB787299uL, 0x788EBF4D01E262FDuL)),
                ((+1, -10, 0xE6C6F3971F92112EuL, 0x6AA3BF1D92CEF640uL), (+1, -8, 0xA8FEE63CF479ED6DuL, 0x37B9813253D9CFAEuL)),
                ((+1, -13, 0xB52877C2D854CF04uL, 0x26FDB8EB8BFE9DFCuL), (+1, -11, 0x84AA5327286AC9C9uL, 0xA4D6C7903FDA150AuL)),
                ((+1, -17, 0xDAF5D28C3AAB6696uL, 0x51CE3D92F4A58F9AuL), (+1, -15, 0xA056DE6188A87331uL, 0xD3EEDEC8246C4D74uL)),
                ((+1, -21, 0xCFBBEFB36F0CC6BFuL, 0x6205F4F2EB5F68B2uL), (+1, -19, 0x982009F66586C531uL, 0xB331955315301C7BuL)),
                ((+1, -25, 0x9BAA29F554635E46uL, 0x51221D53769E5281uL), (+1, -24, 0xE3FC03F3D9CC7AFFuL, 0x8EE2BD74913A0AADuL)),
                ((+1, -30, 0xB6B303539CF61BD8uL, 0x59C1B68F709C0932uL), (+1, -28, 0x85CA1E471A0AE683uL, 0x475E4031192F50CAuL)),
                ((+1, -35, 0xA2F1A3CBA50B4B4FuL, 0x5A028E6D90F8E016uL), (+1, -34, 0xEEA5C1D54BCCE7F0uL, 0xD4A9980DAEC4EBC8uL)),
                ((+1, -41, 0xCD1F5AB5301B9EA9uL, 0x857430F5C693AAA2uL), (+1, -39, 0x9635B572D4AEE71DuL, 0x4BA2E52C0DD8A68AuL)),
                ((+1, -47, 0x95894397D96AFC59uL, 0x50ECFE2D6A23163AuL), (+1, -46, 0xDB0285CDA37A738CuL, 0x93484D4F029A9A17uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm16_32 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xAECE93E7466925F5uL, 0x68DCA68F2D7C286CuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xB756F7F9ACF3FB85uL, 0x3ED14ED4CACA71E8uL), (+1, -1, 0x86432FC7ED62C649uL, 0xE3DB88983E4B67A2uL)),
                ((+1, -5, 0xB73DAA7ED29810E9uL, 0xDB9829EE8862430FuL), (+1, -3, 0x862F97E64C5546D4uL, 0xA010821CCA385298uL)),
                ((+1, -8, 0xE903889E8329103CuL, 0xF7534CEBC8654E57uL), (+1, -6, 0xAAA27E28A6A8E328uL, 0xC5D2AEE45797C107uL)),
                ((+1, -11, 0xD407274F8D22DF13uL, 0x3176E400F1868E11uL), (+1, -9, 0x9B4483F424D3A5C5uL, 0xCE37437F101989BDuL)),
                ((+1, -14, 0x93058948B7D8EF26uL, 0xABF85650778C7A0DuL), (+1, -13, 0xD75386DC4D64EC0CuL, 0xBEA7ECB46A66F9CBuL)),
                ((+1, -18, 0xA15B562F8082A358uL, 0x01AC689750F74BC2uL), (+1, -17, 0xEC526FF6BBC2A641uL, 0x48CAF1D57856D611uL)),
                ((+1, -22, 0x8F80B58D1F44696CuL, 0x524CD9CE0029A6AFuL), (+1, -21, 0xD22C4857F6528A14uL, 0xDB9DFEE0C627BE29uL)),
                ((+1, -27, 0xD1CE936B02B392D5uL, 0x1619777C50C12054uL), (+1, -25, 0x99A40787BEE36D10uL, 0xF9B6C54551D88758uL)),
                ((+1, -32, 0xFDE248C9357B1421uL, 0xB763B7277AF38BBFuL), (+1, -30, 0xB9EB0E4284882D9EuL, 0xA2CFD871D3C72DB7uL)),
                ((+1, -37, 0xFE4888F88DAFBC00uL, 0x5AA845461CB5EFCDuL), (+1, -35, 0xBA35EFCAB7BF3F3AuL, 0x8A836A3F154CC3E9uL)),
                ((+1, -42, 0xD122DB537A379469uL, 0x4E4ECCE67428780DuL), (+1, -40, 0x99264490A2F67346uL, 0x9DCF6DA3492AAFCEuL)),
                ((+1, -47, 0x8A918EA5CF9937D7uL, 0x26BC590987886D7CuL), (+1, -46, 0xCAF253A8492C2ABFuL, 0x624306A4A9DB367AuL)),
                ((+1, -53, 0x8E5D705673B7DEA5uL, 0x4C5F985496A57BCBuL), (+1, -52, 0xD081A2AC192F01FFuL, 0x72E1DFABD1B0737FuL)),
                ((+1, -60, 0xD0C6956C04F214D3uL, 0xB293D5A0B0EA9B7AuL), (+1, -58, 0x98E2C5669889C032uL, 0xCD8AE77FDBB0CFF4uL)),
                ((+1, -67, 0xB1DD5D76D6A4BF6BuL, 0x81A0EABE1CEE7D87uL), (+1, -65, 0x823FD499BF1AA572uL, 0x9BB4477D5598195BuL)),
                ((-1, -93, 0x90B0CD8BE67D0506uL, 0xD8AA5D3163961786uL), Zero),
                ((+1, -102, 0xCCB32674FF74D132uL, 0xEC96647261D42080uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm32_64 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xAECAEFB98D91C092uL, 0xA8993952B50D24BFuL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x8CC130B9D31A5F35uL, 0xD6235951A59A28E6uL), (+1, -2, 0xCE25EB156AEAA1DFuL, 0xB01CB4CB39502648uL)),
                ((+1, -6, 0xDDA168662FFEAE57uL, 0x3BCD3CD687CE9264uL), (+1, -4, 0xA24C92D5AF172B82uL, 0x9718210E884BFF8DuL)),
                ((+1, -9, 0xE2F2F90DC928BF67uL, 0xC4701FFBD6D1BDB7uL), (+1, -7, 0xA631A46657B3B1FDuL, 0x0E0E58DA256C52F8uL)),
                ((+1, -12, 0xA993F620BBE40E49uL, 0x007F7A7FB1589ACEuL), (+1, -11, 0xF85CCE786D28E2A6uL, 0x8E110194CF5D7C0FuL)),
                ((+1, -16, 0xC4A1E58983CB1C9AuL, 0x580A660CB3EBFB92uL), (+1, -14, 0x8FFE3CCE97CA7A2DuL, 0xCC600B1C637F14ABuL)),
                ((+1, -20, 0xB797931A9AF3AA06uL, 0x8C77DDDD351136F8uL), (+1, -18, 0x8671990CE382C1B3uL, 0x32F56448EB1E80D7uL)),
                ((+1, -24, 0x8D527752B3DDFADEuL, 0x6BE9EF94122ED668uL), (+1, -23, 0xCEFAB01E96C9AA14uL, 0x7B7436FE16373B93uL)),
                ((+1, -29, 0xB6117270F3F93B86uL, 0xA011412C36AB77A8uL), (+1, -27, 0x8553E8B48AD6E522uL, 0x3EC45141872C209FuL)),
                ((+1, -34, 0xC5F722AF14B62E41uL, 0x8F8113309ECDAD1BuL), (+1, -32, 0x90F82016F2F236FEuL, 0x8BF01B0BDB888E67uL)),
                ((+1, -39, 0xB634A36C7C6BEB76uL, 0x7C895FB87306C371uL), (+1, -37, 0x856DADF26762E59BuL, 0xB97C590C4F876818uL)),
                ((+1, -44, 0x8D90E302BB03E158uL, 0xB83C77C5B3B05BD4uL), (+1, -43, 0xCF561BD6E36EFB85uL, 0xF1425199A286683EuL)),
                ((+1, -50, 0xB7E2C6DB93F0B054uL, 0xFEF1375E5A7958E8uL), (+1, -48, 0x86A8AB0859039729uL, 0xBDFAF745F9B8D3BDuL)),
                ((+1, -56, 0xC3A53814773FA863uL, 0xE42E1DE714880AD6uL), (+1, -54, 0x8F4533E8421CB649uL, 0x1C07665BEE4EE3BDuL)),
                ((+1, -62, 0xA3F6171E41836522uL, 0xBC6047CD8986C0BBuL), (+1, -61, 0xF022E98299CA697DuL, 0x034934656DBA5AF7uL)),
                ((+1, -69, 0xC76DCF65485ECFA2uL, 0x8306406A552D8C6AuL), (+1, -67, 0x920A7F82234D66FEuL, 0x0CB55391C07DDD9BuL)),
                ((+1, -76, 0x8E76BACBC4D498E0uL, 0x093644EE253C6A8AuL), (+1, -75, 0xD0A6BC1C8E17EA60uL, 0xDE1A1763B3CB49E2uL)),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm64_128 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xAECAEFB5E9576CD5uL, 0x114F76C8A28C3260uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0x9B09F9C848FA0218uL, 0x5B407A149C2AD303uL), (+1, -2, 0xE3119DCF64E37DFDuL, 0xBC594FFA3D957578uL)),
                ((+1, -5, 0x8715C569F28E9DC7uL, 0x2BB3F6EC67CD351CuL), (+1, -4, 0xC5D82E7584D0D105uL, 0x7B80D37C7E498D89uL)),
                ((+1, -8, 0x99DFF8A93D1F8ED9uL, 0x398B49D17149AC2DuL), (+1, -7, 0xE15D2961CB1787F7uL, 0x638670C814517928uL)),
                ((+1, -11, 0x809E5BC293A5F95AuL, 0xA018B24D7404A1D3uL), (+1, -10, 0xBC5FAD6DAED92015uL, 0x630C3266B8A140BBuL)),
                ((+1, -15, 0xA7DB9968FE4D056EuL, 0xEEE5B3D1DC7FED84uL), (+1, -14, 0xF5D7DB0D2460BD3BuL, 0x08D1A9AEE8AF2520uL)),
                ((+1, -19, 0xB1958C34116D6F67uL, 0xC9161E382703BC7DuL), (+1, -17, 0x820B4E94B147BA97uL, 0xBA395AE49EF1C8BCuL)),
                ((+1, -23, 0x9C027C0E99FD29D2uL, 0x91D146F0B87BE84BuL), (+1, -22, 0xE47D949E7EB4AED4uL, 0xE67BDC1E167FFBA9uL)),
                ((+1, -28, 0xE75DCEE513174A7AuL, 0x34573460E074C93FuL), (+1, -26, 0xA96DBF8F7DB9CE2BuL, 0x54F23AA37B6B0EDBuL)),
                ((+1, -32, 0x918B917F3519638DuL, 0x35ECD48E3ABBB63FuL), (+1, -31, 0xD52A0FC0C9E84404uL, 0x68F4C410FA0D9FA8uL)),
                ((+1, -37, 0x9EBDBD07B79FCA84uL, 0x46139CD2BD249DD8uL), (+1, -36, 0xE87DB3B57CD83F20uL, 0x1D356FB5BB063A66uL)),
                ((+1, -42, 0x8B4429D87D31BC34uL, 0xE5C67746CC0C92A2uL), (+1, -41, 0xCBF7DE8F7D18357FuL, 0xF5BC36D860CA43FFuL)),
                ((+1, -48, 0xF6356F60F71875DFuL, 0x0BB4429941F511ECuL), (+1, -46, 0xB44C331BCD818D73uL, 0x43E435F563250F33uL)),
                ((+1, -54, 0xCE98467DDB46B487uL, 0xA2D57D7436456E92uL), (+1, -52, 0x9749DE98BE0C85D0uL, 0x4A9FC47FB42DFD16uL)),
                ((+1, -59, 0xD7F47DCFF6A2F287uL, 0x9F77DCA8F5229A10uL), (+1, -57, 0x9E249BDD56FC9147uL, 0x8FA81D90C8CA3A03uL)),
                ((-1, -151, 0xC56A6744D3AAF2ADuL, 0x53885206871505D9uL), Zero),
            }));

            private static readonly ReadOnlyCollection<(ddouble c, ddouble d)> pade_plus_expm128_192 = new(Array.AsReadOnly(new (ddouble c, ddouble)[]{
                ((+1, -2, 0xAECAEFB5E9576CD1uL, 0x6D1522FE71B6C326uL), (+1, 0, 0x8000000000000000uL, 0x0000000000000000uL)),
                ((+1, -3, 0xA03903CE7321F9A2uL, 0x0E74DB1AF62304C9uL), (+1, -2, 0xEAA92FEE3DA03E4CuL, 0xE1055A5EC96BF1A6uL)),
                ((+1, -5, 0x8F8B5931D5105742uL, 0x2FF6CE2DDC59B7EAuL), (+1, -4, 0xD23BDE4E2D6E690EuL, 0xECB2DCCDD73C751FuL)),
                ((+1, -8, 0xA740820377578211uL, 0x3B3EFEA66E5FB9C0uL), (+1, -7, 0xF4F4B5C75E33B5E2uL, 0x8884ABCB0C9A9561uL)),
                ((+1, -11, 0x8C4D5FFBE49A5D1CuL, 0x8C9D56322DFA3804uL), (+1, -10, 0xCD7C4BC3C1123B9FuL, 0xA8F1B2E691DD94E3uL)),
                ((+1, -15, 0xC031E1B6C1066D42uL, 0x4E40B373DFA5DF49uL), (+1, -13, 0x8CBE56A0D3ACDC06uL, 0x3E4CF16DC45637FDuL)),
                ((+1, -19, 0xA17F40A2A8418672uL, 0x70B23916E7A4B9E9uL), (+1, -18, 0xEC86FE0031EC77EEuL, 0x0D4F807B9AE3B8DAuL)),
                ((+1, -23, 0xEC82D9CE62DF5702uL, 0xC4A83640B826F293uL), (+1, -21, 0xAD323697C894F5B9uL, 0x6ECF23B764DEB4A7uL)),
                ((+1, -34, 0xEE62F8523ED4C8A1uL, 0x8C11A81636829251uL), (+1, -32, 0xAE91CD76EC370244uL, 0x05A87889BACF331CuL)),
                ((+1, -31, 0x90DF80AEFF3E10E3uL, 0xA74B28D73028E1F1uL), (+1, -30, 0xD42E0E37EB3A95ADuL, 0x52C78576C2303F1EuL)),
                ((+1, -180, 0xFC06BB0047A8661FuL, 0xB306200B412CDECCuL), Zero),
                ((-1, -189, 0xFB986F5E5D955C69uL, 0x18D9B59FDF513AEEuL), Zero),
            }));


            public static ddouble Value(ddouble x) {
                if (x > 0.5) {
                    return -Value(1d - x);
                }

                if (x >= 0.125d) {
                    ddouble y;
                    if (x <= 0.1328125d) {
                        y = ApproxUtil.Pade(0.1328125d - x, pade_plus_0p125_0p1328125);
                    }
                    else if (x <= 0.140625d) {
                        y = ApproxUtil.Pade(0.140625d - x, pade_plus_0p1328125_0p140625);
                    }
                    else if (x <= 0.1484375d) {
                        y = ApproxUtil.Pade(0.1484375d - x, pade_plus_0p140625_0p1484375);
                    }
                    else if (x <= 0.15625d) {
                        y = ApproxUtil.Pade(0.15625d - x, pade_plus_0p1484375_0p15625);
                    }
                    else if (x <= 0.1875d) {
                        y = ApproxUtil.Pade(0.1875d - x, pade_plus_0p15625_0p1875);
                    }
                    else if (x <= 0.21875d) {
                        y = ApproxUtil.Pade(0.21875d - x, pade_plus_0p1875_0p21875);
                    }
                    else if (x <= 0.25d) {
                        y = ApproxUtil.Pade(0.25d - x, pade_plus_0p21875_0p25);
                    }
                    else if (x <= 0.3125d) {
                        y = ApproxUtil.Pade(0.3125d - x, pade_plus_0p25_0p3125);
                    }
                    else if (x <= 0.34375d) {
                        y = ApproxUtil.Pade(0.34375d - x, pade_plus_0p3125_0p34375);
                    }
                    else if (x <= 0.375d) {
                        y = ApproxUtil.Pade(0.375d - x, pade_plus_0p34375_0p375);
                    }
                    else {
                        y = ApproxUtil.Pade(0.5d - x, pade_plus_0p375_0p5);
                    }

                    return y;
                }
                else {
                    ddouble v;
                    int exponent = double.ILogB((double)x);

                    if (exponent >= -4) {
                        ddouble u = -Log2(Ldexp(x, 3));

                        if (u <= 0.5d) {
                            v = ApproxUtil.Pade(u, pade_plus_expm3_3p5);
                        }
                        else {
                            v = ApproxUtil.Pade(u - 0.5d, pade_plus_expm3p5_4);
                        }
                    }
                    else if (exponent >= -5) {
                        ddouble u = -Log2(Ldexp(x, 4));

                        if (u <= 0.5d) {
                            v = ApproxUtil.Pade(u, pade_plus_expm4_4p5);
                        }
                        else {
                            v = ApproxUtil.Pade(u - 0.5d, pade_plus_expm4p5_5);
                        }
                    }
                    else if (exponent >= -6) {
                        ddouble u = -Log2(Ldexp(x, 5));

                        if (u <= 0.5d) {
                            v = ApproxUtil.Pade(u, pade_plus_expm5_5p5);
                        }
                        else {
                            v = ApproxUtil.Pade(u - 0.5d, pade_plus_expm5p5_6);
                        }
                    }
                    else if (exponent >= -8) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 6)), pade_plus_expm6_8);
                    }
                    else if (exponent >= -12) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 8)), pade_plus_expm8_12);
                    }
                    else if (exponent >= -16) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 12)), pade_plus_expm12_16);
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
                    else if (exponent >= -192) {
                        v = ApproxUtil.Pade(-Log2(Ldexp(x, 128)), pade_plus_expm128_192);
                    }
                    else {
                        v = 1d / Ldexp(Cbrt(PI), 1);
                    }

                    ddouble y = v / ExMath.Pow2d3(x);

                    return y;
                }
            }
        }
    }
}
