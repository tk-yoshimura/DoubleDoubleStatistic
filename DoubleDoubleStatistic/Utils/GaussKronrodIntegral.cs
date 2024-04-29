using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DoubleDouble.ddouble;

namespace DoubleDoubleStatistic {
    internal static class GaussKronrodIntegral {
        const long gk_points = 63;

        public static (ddouble value, ddouble error) Integrate(Func<ddouble, ddouble> f, ddouble a, ddouble b) {
            ReadOnlyCollection<(ddouble x, ddouble wg, ddouble wk)> ps = gauss_kronrod_table;

            ddouble sg = Zero, sk = Zero;
            ddouble r = b - a;

            if (!IsFinite(r)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            for (int i = 0; i < ps.Count; i++) {
                ddouble x = ps[i].x;
                ddouble x_shifted = x * r + a;

                ddouble y = f(x_shifted);

                sk += ps[i].wk * y;

                if ((i & 1) == 1) {
                    sg += ps[i].wg * y;
                }
            }

            sk *= r;
            sg *= r;

            ddouble error = Abs(sk - sg);

            return (sk, error);
        }

        private static (ddouble value, ddouble error, long eval_points) IntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, long discontinue_eval_points) {
            Debug.Assert(discontinue_eval_points >= 0);

            PriorityQueue<(ddouble a, ddouble b, ddouble eps), long> queue = new();
            queue.Enqueue((a, b, eps), 0);

            long eval_points_sum = 0;
            ddouble value_sum = 0d, error_sum = 0d;

            while (queue.Count > 0) {
                (a, b, eps) = queue.Dequeue();

                (ddouble value, ddouble error) = Integrate(f, a, b);

                long eval_points = gk_points;
                eval_points_sum += eval_points;

                if (!(error > eps) || eval_points_sum > discontinue_eval_points) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                ddouble c = Ldexp(a + b, -1), eps_half = Ldexp(eps, -1);
                long priority = double.ILogB((double)error);
                queue.Enqueue((a, c, eps_half), -priority);
                queue.Enqueue((c, b, eps_half), -priority);
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (ddouble value, ddouble error, long eval_points) AdaptiveIntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, long discontinue_eval_points) {
            return IntegrateFiniteInterval(f, a, b, eps, discontinue_eval_points);
        }

        private static (ddouble value, ddouble error, long eval_points) AdaptiveIntegrateInfiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, long discontinue_eval_points) {
            if (IsNaN(a) || IsNaN(b)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            if (a > b) {
                (ddouble value, ddouble error, long eval_points) = AdaptiveIntegrateInfiniteInterval(f, b, a, eps, discontinue_eval_points);

                return (-value, error, eval_points);
            }

            if (IsInfinity(a) && IsInfinity(b)) {
                if (Sign(a) == Sign(b)) {
                    throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
                }

                ddouble g(ddouble t) {
                    if (IsZero(t)) {
                        return Zero;
                    }

                    ddouble u = (1d - t) / t;

                    return (f(u) + f(-u)) / (t * t);
                }

                return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, discontinue_eval_points);
            }

            if (IsFinite(a) && IsInfinity(b)) {
                if (IsZero(a)) {
                    ddouble g(ddouble t) {
                        if (IsZero(t)) {
                            return Zero;
                        }

                        ddouble u = (1d - t) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, discontinue_eval_points);
                }
                else {
                    ddouble g(ddouble t) {
                        if (IsZero(t)) {
                            return Zero;
                        }

                        ddouble u = (1d - t) / t;

                        return f(a + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, discontinue_eval_points);
                }
            }

            if (IsInfinity(a) && IsFinite(b)) {
                if (IsZero(b)) {
                    ddouble g(ddouble t) {
                        if (IsZero(t)) {
                            return Zero;
                        }

                        ddouble u = (t - 1d) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, discontinue_eval_points);
                }
                else {
                    ddouble g(ddouble t) {
                        if (IsZero(t)) {
                            return Zero;
                        }

                        ddouble u = (t - 1d) / t;

                        return f(b + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, discontinue_eval_points);
                }
            }

            throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
        }

        public static (ddouble value, ddouble error, long eval_points) AdaptiveIntegrate(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, long discontinue_eval_points = -1) {
            if (discontinue_eval_points < -1) {
                throw new ArgumentOutOfRangeException(nameof(discontinue_eval_points), "Invalid param. discontinue_eval_points=-1: infinite, discontinue_eval_points>=0: finite");
            }
            if (!(eps >= 0d)) {
                throw new ArgumentOutOfRangeException(nameof(eps), "Invalid param. eps must be nonnegative value");
            }

            if (IsFinite(a) && IsFinite(b)) {
                return AdaptiveIntegrateFiniteInterval(f, a, b, eps, discontinue_eval_points);
            }
            else {
                return AdaptiveIntegrateInfiniteInterval(f, a, b, eps, discontinue_eval_points);
            }
        }

        private static readonly ReadOnlyCollection<(ddouble x, ddouble wg, ddouble wk)> gauss_kronrod_table =
            new([
                ((+1, -13, 0xFD8340F5FF8608A4uL, 0x5D0C3EC2D8485580uL), Zero, (+1, -11, 0xAAB76582A85A393FuL, 0x4CD98C758D4C6984uL)),
                ((+1, -10, 0xBEDFF255980D822FuL, 0xC58E40E2F554990EuL), (+1, -9, 0xF4CDE0A745242B3DuL, 0xEE87BF6947587077uL), (+1, -10, 0xEF04BBC9225B48B5uL, 0x0A286FFC3F906D8CuL)),
                ((+1, -8, 0x808E3995728DDE05uL, 0xE10B625EA85AD919uL), Zero, (+1, -9, 0xCBBBAE7C75F193C0uL, 0x4C2E641D8A2F26D7uL)),
                ((+1, -8, 0xFAE7F349FE858CA8uL, 0x61A38D937967D797uL), (+1, -7, 0x8DDFC7BCFFAB5408uL, 0x00C43E4A0BAFEEB3uL), (+1, -8, 0x8E8057273DF5A6FCuL, 0xB0E3B35A3B1B2C22uL)),
                ((+1, -7, 0xCE88A2CF1CA7EC33uL, 0x9A9B364269EEB324uL), Zero, (+1, -8, 0xB5B58E2C73FDD585uL, 0xE30D0543BD9AD03FuL)),
                ((+1, -6, 0x99957BF744301BCDuL, 0x9721F191D947676EuL), (+1, -7, 0xDD4204C316B7CD9DuL, 0xF8D06E802A89279AuL), (+1, -8, 0xDCE515C3F665AAD6uL, 0x7E322459455CBF31uL)),
                ((+1, -6, 0xD5BABD4C4A0B7F63uL, 0xAC486CA486BFEE88uL), Zero, (+1, -7, 0x821BAD280529B1E6uL, 0x65A70CACDC5C3B6EuL)),
                ((+1, -5, 0x8DCF466C1606E5C4uL, 0x6D1328D53B2D2696uL), (+1, -6, 0x953A020F8FCF04A9uL, 0xCC5F88778796412FuL), (+1, -7, 0x9558F0F476FCC29CuL, 0x4E104521FFBB88D6uL)),
                ((+1, -5, 0xB57C96652F29D61EuL, 0x54918880FD479083uL), Zero, (+1, -7, 0xA7FBAD2280BA659FuL, 0x28774041883A01FBuL)),
                ((+1, -5, 0xE1C57E50AC6E9577uL, 0x8031084C6F8FE6CBuL), (+1, -6, 0xBA579C200CC5374EuL, 0xFA53C7F4C6A8295EuL), (+1, -7, 0xBA40EDCCDBFE7A47uL, 0x66101D5E51DA7274uL)),
                ((+1, -4, 0x894B956529E67E35uL, 0x6DC14E65DF6787A4uL), Zero, (+1, -7, 0xCC3BF603CB587378uL, 0x428C8C78E6560548uL)),
                ((+1, -4, 0xA3EBED76E75AACFAuL, 0xD45FE7B96E82D659uL), (+1, -6, 0xDD9B319A0B82C604uL, 0x9E62EEAD81757094uL), (+1, -7, 0xDDACE10BF6194661uL, 0x09CD2FB4432BAF79uL)),
                ((+1, -4, 0xC0AF5ACA39CAA887uL, 0x83A3424459018FCBuL), Zero, (+1, -7, 0xEE6DF39B20D25E78uL, 0x4058EA88C49338BCuL)),
                ((+1, -4, 0xDF81587567EAD17CuL, 0x9EF0F2A898710063uL), (+1, -6, 0xFEAAFD2EA577BEF2uL, 0xBF640B8CB7FC8C06uL), (+1, -7, 0xFE9C855AA7962747uL, 0x3105CB7071AB14CCuL)),
                ((+1, -3, 0x8028834E10858DDCuL, 0x61614421B55F4DA9uL), Zero, (+1, -6, 0x8722D0A18230EB7EuL, 0xECCD18D5BA2375C4uL)),
                ((+1, -3, 0x9185B48DED467882uL, 0x604966D19CA7D8FEuL), (+1, -5, 0x8E996E3EE738D84CuL, 0xA8882C4BBB10DAF1uL), (+1, -6, 0x8E9F8FA1EEBF8CD8uL, 0x0AA959736129B236uL)),
                ((+1, -3, 0xA3CC122446A1BE59uL, 0xF6D1E5609222E0A9uL), Zero, (+1, -6, 0x95B564815AC27BD8uL, 0x5BF5985ACD8384F3uL)),
                ((+1, -3, 0xB6EF24F83DEF8F8AuL, 0x7F5344F94FB12990uL), (+1, -5, 0x9C7290215D454022uL, 0xF9B8D5414C324FF8uL), (+1, -6, 0x9C6D35AF3D4003BAuL, 0x2B559041F9201065uL)),
                ((+1, -3, 0xCAE3B4D2A6535A09uL, 0x16917F24F826CCF9uL), Zero, (+1, -6, 0xA2CC186B4109676AuL, 0x81F425DF46EF754AuL)),
                ((+1, -3, 0xDF9DCE508384472AuL, 0xC72B3DEC1F16D133uL), (+1, -5, 0xA8BDA82848F4E0BBuL, 0x70395AA23939BDD9uL), (+1, -6, 0xA8C2720262AC4566uL, 0x30F6161C9B48432CuL)),
                ((+1, -3, 0xF50F691912A065B5uL, 0xF8663ABEDD6DABBFuL), Zero, (+1, -6, 0xAE43DF194875369FuL, 0x344EA898C86E36A1uL)),
                ((+1, -2, 0x859517A0F7729022uL, 0x4223F4972617FBD9uL), (+1, -5, 0xB35B6F266DF3C626uL, 0x71513E8B63AB3BC3uL), (+1, -6, 0xB3570E86BAEDF839uL, 0x92C3251B419C5D75uL)),
                ((+1, -2, 0x90F068205DBE73AEuL, 0xAC49D75110E2A377uL), Zero, (+1, -6, 0xB800CF381E76697BuL, 0x949BC5A0B9CB9171uL)),
                ((+1, -2, 0x9C92BDD65D7B7FE9uL, 0x76410791CFB61F9BuL), (+1, -5, 0xBC30E24FFDC2196BuL, 0xB5FCC5DEDB046A3EuL), (+1, -6, 0xBC34F5F01177A88CuL, 0xFF916202C697CAD8uL)),
                ((+1, -2, 0xA87458D64D0B1B18uL, 0xDD71FAFC756E86C1uL), Zero, (+1, -6, 0xBFE9362C44442D88uL, 0xAEFFA5B2AFA9D1B4uL)),
                ((+1, -2, 0xB48D5B47E36B0D62uL, 0x61163BA0A8CD0CA1uL), (+1, -5, 0xC32787E77D8F1A4BuL, 0x3F826BF1F049627EuL), (+1, -6, 0xC323AAB5760D5E0CuL, 0x05B5CC315214F0E3uL)),
                ((+1, -2, 0xC0D6647107E28578uL, 0xE571861DDC049B0AuL), Zero, (+1, -6, 0xC5E99F04DD11F3CAuL, 0x62DACAC3C6559882uL)),
                ((+1, -2, 0xCD47EDD8359D05B4uL, 0x779A17FE25808B5CuL), (+1, -5, 0xC82DA8671705BB23uL, 0x64F0E798D1E7CDCAuL), (+1, -6, 0xC83160EFFCE98DE8uL, 0x97E6649913333FF8uL)),
                ((+1, -2, 0xD9D9BDE741F14BF0uL, 0x5DA7BAA0550BE76FuL), Zero, (+1, -6, 0xC9F24F3F934A2FC9uL, 0x371456393DC6F258uL)),
                ((+1, -2, 0xE6838B0629A19512uL, 0x402D8482B5499185uL), (+1, -5, 0xCB367B914F71580EuL, 0xF0D57001090A26FCuL), (+1, -6, 0xCB32D7E28095C48AuL, 0x87DCC296F4C41DE7uL)),
                ((+1, -2, 0xF33D9139622A0024uL, 0x7A0B69C9A3A46B93uL), Zero, (+1, -6, 0xCBF94D81D0B8B273uL, 0x152C0D545D41CB04uL)),
                ((+1, -1, 0x8000000000000000uL, 0x0000000000000000uL), (+1, -5, 0xCC3A48F5043D1C5FuL, 0x82C81EBC9F783387uL), (+1, -6, 0xCC3DE5A9AB0D3362uL, 0xACDC08E394621031uL)),
                ((+1, -1, 0x866137634EEAFFEDuL, 0xC2FA4B1B2E2DCA36uL), Zero, (+1, -6, 0xCBF94D81D0B8B273uL, 0x152C0D545D41CB04uL)),
                ((+1, -1, 0x8CBE3A7CEB2F3576uL, 0xDFE93DBEA55B373DuL), (+1, -5, 0xCB367B914F71580EuL, 0xF0D57001090A26FCuL), (+1, -6, 0xCB32D7E28095C48AuL, 0x87DCC296F4C41DE7uL)),
                ((+1, -1, 0x9313210C5F075A07uL, 0xD12C22AFD57A0C48uL), Zero, (+1, -6, 0xC9F24F3F934A2FC9uL, 0x371456393DC6F258uL)),
                ((+1, -1, 0x995C0913E5317D25uL, 0xC432F400ED3FBA51uL), (+1, -5, 0xC82DA8671705BB23uL, 0x64F0E798D1E7CDCAuL), (+1, -6, 0xC83160EFFCE98DE8uL, 0x97E6649913333FF8uL)),
                ((+1, -1, 0x9F94CDC77C0EBD43uL, 0x8D473CF111FDB27AuL), Zero, (+1, -6, 0xC5E99F04DD11F3CAuL, 0x62DACAC3C6559882uL)),
                ((+1, -1, 0xA5B9525C0E4A794EuL, 0xCF74E22FAB9979AFuL), (+1, -5, 0xC32787E77D8F1A4BuL, 0x3F826BF1F049627EuL), (+1, -6, 0xC323AAB5760D5E0CuL, 0x05B5CC315214F0E3uL)),
                ((+1, -1, 0xABC5D394D97A7273uL, 0x91470281C548BC9FuL), Zero, (+1, -6, 0xBFE9362C44442D88uL, 0xAEFFA5B2AFA9D1B4uL)),
                ((+1, -1, 0xB1B6A114D142400BuL, 0x44DF7C371824F032uL), (+1, -5, 0xBC30E24FFDC2196BuL, 0xB5FCC5DEDB046A3EuL), (+1, -6, 0xBC34F5F01177A88CuL, 0xFF916202C697CAD8uL)),
                ((+1, -1, 0xB787CBEFD120C628uL, 0xA9DB1457778EAE44uL), Zero, (+1, -6, 0xB800CF381E76697BuL, 0x949BC5A0B9CB9171uL)),
                ((+1, -1, 0xBD35742F8446B7EEuL, 0xDEEE05B46CF40213uL), (+1, -5, 0xB35B6F266DF3C626uL, 0x71513E8B63AB3BC3uL), (+1, -6, 0xB3570E86BAEDF839uL, 0x92C3251B419C5D75uL)),
                ((+1, -1, 0xC2BC25B9BB57E692uL, 0x81E6715048A49510uL), Zero, (+1, -6, 0xAE43DF194875369FuL, 0x344EA898C86E36A1uL)),
                ((+1, -1, 0xC8188C6BDF1EEE35uL, 0x4E353084F83A4BB3uL), (+1, -5, 0xA8BDA82848F4E0BBuL, 0x70395AA23939BDD9uL), (+1, -6, 0xA8C2720262AC4566uL, 0x30F6161C9B48432CuL)),
                ((+1, -1, 0xCD4712CB566B297DuL, 0xBA5BA036C1F64CC1uL), Zero, (+1, -6, 0xA2CC186B4109676AuL, 0x81F425DF46EF754AuL)),
                ((+1, -1, 0xD24436C1F0841C1DuL, 0x602B2EC1AC13B59BuL), (+1, -5, 0x9C7290215D454022uL, 0xF9B8D5414C324FF8uL), (+1, -6, 0x9C6D35AF3D4003BAuL, 0x2B559041F9201065uL)),
                ((+1, -1, 0xD70CFB76EE579069uL, 0x824B86A7DB7747D5uL), Zero, (+1, -6, 0x95B564815AC27BD8uL, 0x5BF5985ACD8384F3uL)),
                ((+1, -1, 0xDB9E92DC84AE61DFuL, 0x67EDA64B98D609C0uL), (+1, -5, 0x8E996E3EE738D84CuL, 0xA8882C4BBB10DAF1uL), (+1, -6, 0x8E9F8FA1EEBF8CD8uL, 0x0AA959736129B236uL)),
                ((+1, -1, 0xDFF5DF2C7BDE9C88uL, 0xE7A7AEF792A82C95uL), Zero, (+1, -6, 0x8722D0A18230EB7EuL, 0xECCD18D5BA2375C4uL)),
                ((+1, -1, 0xE40FD4F15302A5D0uL, 0x6C21E1AAECF1DFF3uL), (+1, -6, 0xFEAAFD2EA577BEF2uL, 0xBF640B8CB7FC8C06uL), (+1, -7, 0xFE9C855AA7962747uL, 0x3105CB7071AB14CCuL)),
                ((+1, -1, 0xE7EA14A6B8C6AAEFuL, 0x0F8B97B774DFCE06uL), Zero, (+1, -7, 0xEE6DF39B20D25E78uL, 0x4058EA88C49338BCuL)),
                ((+1, -1, 0xEB8282512314AA60uL, 0xA5740308D22FA534uL), (+1, -6, 0xDD9B319A0B82C604uL, 0x9E62EEAD81757094uL), (+1, -7, 0xDDACE10BF6194661uL, 0x09CD2FB4432BAF79uL)),
                ((+1, -1, 0xEED68D535AC33039uL, 0x5247D63344130F0BuL), Zero, (+1, -7, 0xCC3BF603CB587378uL, 0x428C8C78E6560548uL)),
                ((+1, -1, 0xF1E3A81AF53916A8uL, 0x87FCEF7B39070193uL), (+1, -6, 0xBA579C200CC5374EuL, 0xFA53C7F4C6A8295EuL), (+1, -7, 0xBA40EDCCDBFE7A47uL, 0x66101D5E51DA7274uL)),
                ((+1, -1, 0xF4A83699AD0D629EuL, 0x1AB6E777F02B86F7uL), Zero, (+1, -7, 0xA7FBAD2280BA659FuL, 0x28774041883A01FBuL)),
                ((+1, -1, 0xF7230B993E9F91A3uL, 0xB92ECD72AC4D2D96uL), (+1, -6, 0x953A020F8FCF04A9uL, 0xCC5F88778796412FuL), (+1, -7, 0x9558F0F476FCC29CuL, 0x4E104521FFBB88D6uL)),
                ((+1, -1, 0xF9522A159DAFA404uL, 0xE29DBC9ADBCA008BuL), Zero, (+1, -7, 0x821BAD280529B1E6uL, 0x65A70CACDC5C3B6EuL)),
                ((+1, -1, 0xFB33542045DE7F21uL, 0x9346F0737135C4C4uL), (+1, -7, 0xDD4204C316B7CD9DuL, 0xF8D06E802A89279AuL), (+1, -8, 0xDCE515C3F665AAD6uL, 0x7E322459455CBF31uL)),
                ((+1, -1, 0xFCC5DD74C38D604FuL, 0x31959326F6584533uL), Zero, (+1, -8, 0xB5B58E2C73FDD585uL, 0xE30D0543BD9AD03FuL)),
                ((+1, -1, 0xFE0A30196C02F4E6uL, 0xAF3CB8E4D90D3050uL), (+1, -7, 0x8DDFC7BCFFAB5408uL, 0x00C43E4A0BAFEEB3uL), (+1, -8, 0x8E8057273DF5A6FCuL, 0xB0E3B35A3B1B2C22uL)),
                ((+1, -1, 0xFEFEE38CD51AE443uL, 0xF43DE93B42AF4A4DuL), Zero, (+1, -9, 0xCBBBAE7C75F193C0uL, 0x4C2E641D8A2F26D7uL)),
                ((+1, -1, 0xFFA09006D533F93EuL, 0xE81D38DF8E8555B3uL), (+1, -9, 0xF4CDE0A745242B3DuL, 0xEE87BF6947587077uL), (+1, -10, 0xEF04BBC9225B48B5uL, 0x0A286FFC3F906D8CuL)),
                ((+1, -1, 0xFFF027CBF0A0079FuL, 0x75BA2F3C13D27B7AuL), Zero, (+1, -11, 0xAAB76582A85A393FuL, 0x4CD98C758D4C6984uL)),
            ]);
    }
}
