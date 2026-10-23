/* static/pfdvideo.js — dynamic start-up film of the last run
   NavNevAI_mAbWasteWater · CarbonTesseract 2026

   The PFD starts greyed out. Each unit and utility lights up, with an iOS-style
   zoom and highlight, at the hour it first carries flow, in flow order. Flowing
   streams animate; vessel levels, column states and values follow the run hour
   by hour. Quiet stretches where no new unit joins (e.g. the N-1 seed phase) are
   fast-forwarded to at most MAX_GAP_S seconds of film.
   When the NevAI watch flags titre drifting from the forecast, the top pill says so during playback.
   Recording uses MediaRecorder: MP4 (H.264) where the browser supports
   it (Safari, recent Chrome/Edge), otherwise WebM. Recording runs in real time.
   Page ids: vidStart vidHour vidLen vidRate vidTour vidRes vidPreview vidRecord vidStop vidInfo vidCanvas vidOut
   API: PFDVideo.init()  PFDVideo.load(simResult)  PFDVideo.stop()
*/
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const num = v => (typeof v === 'number' && isFinite(v)) ? v : 0;
  const clamp = (x, a = 0, b = 1) => x < a ? a : x > b ? b : x;
  const easeIO = x => x < 0.5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2;
  const easeOutBack = x => { const c1 = 1.4, c3 = c1 + 1; return 1 + c3 * Math.pow(x - 1, 3) + c1 * Math.pow(x - 1, 2); };
  const hexA = (h, a) => `rgba(${parseInt(h.slice(1, 3), 16)},${parseInt(h.slice(3, 5), 16)},${parseInt(h.slice(5, 7), 16)},${a})`;
  const mmss = t => `${Math.floor(t / 60)}:${String(Math.floor(t % 60)).padStart(2, '0')}`;
  function fmt(v, d) {
    if (v === null || v === undefined || !isFinite(v)) return '—';
    const a = Math.abs(v);
    if (a < 5e-4) return '0';
    return a >= 1000 ? Math.round(v).toLocaleString() : a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : v.toFixed(d === undefined ? 2 : d);
  }

  const FONT = '-apple-system, BlinkMacSystemFont, "SF Pro Display", "SF Pro Text", "Helvetica Neue", "Segoe UI", Roboto, Arial, sans-serif';
  const VW = 1920, VH = 1080;                       // virtual frame, scaled for 720p
  const VP = {x: 24, y: 24, w: 1376, h: 1032};      // diagram card
  const HUD = {x: 1424, y: 24, w: 472, h: 1032};    // side panel
  const CAP_H = 176;
  const TEAL = '#0e7c7b', TEAL2 = '#14b8a6', INK = '#16202c', GREY = '#6b7280';
  const KCOL = {pp: '#1d2733', pb: '#7c3aed', pw: '#2563eb', pa: '#dc2626', pd: '#a16207', px: '#64748b', ps: '#16a34a'};
  const FCOL = Object.assign({}, KCOL, {pp: TEAL});
  const KDASH = {pp: [], pb: [6, 3], pw: [6, 3], pa: [6, 3], pd: [6, 3], px: [2, 3], ps: [9, 4]};
  const TITLE_S = 3.2, FOCUS_S = 2.8, FLY_S = 0.9, OUT_S = 1.0, END_S = 5;
  const MAX_GAP_S = 5;   // longest quiet stretch (no new unit) in film seconds
  const WATERMARK = 'CarbonTesseract 2026';

  function U(key, short, title, tag, util, box, act, desc) { return {key, short, title, tag, util: !!util, box, act, desc}; }
  // flow order; box = [x0, y0, x1, y1] in diagram coordinates; act = streams that put the unit in the circuit
  const UNITS = [
    U('n1', 'N-1 seed', 'N-1 perfusion seed', 'R-101', 0, [38, 28, 190, 240], ['U1', 'U2'],
      'Media in = cell-free permeate out to waste; cells stay and build to seed density.'),
    U('seed', 'Seed transfer', 'Seed transfer to the N reactor', 'S1', 0, [95, 236, 265, 290], ['S1'],
      'The N-1 culture moves into the N reactor at inoculation.'),
    U('n', 'N reactor', 'N production reactor (IDB)', 'R-102', 0, [190, 28, 325, 252], ['U3'],
      'Continuous media feed; outflow = cells + broth + product (in = out).'),
    U('tff', 'TFF cell removal', 'TFF cell removal', 'F-101', 0, [318, 132, 510, 272], ['U4', 'U5', '1'],
      'Removes all cells: concentrate to bio-kill, permeate (broth + product) forward.'),
    U('divert', 'Pre-capture divert', 'Divert valve before capture', 'XV-101', 0, [470, 100, 565, 242], ['2'],
      'Permeate goes to drain until capture starts.'),
    U('surge', 'Surge tank', 'Conditioned media surge tank', 'T-101', 0, [548, 66, 690, 244], ['3', '4'],
      'Collects permeate; PAT titre × volume decides when a Protein A load is ready.'),
    U('det', 'Brij-9 dosing', 'Brij-9 detergent dosing', 'Brij-9', 1, [636, 0, 744, 150], ['8'],
      '10 %wt concentrate dosed in-line into each load: kills enveloped viruses before the column.'),
    U('pabuf', 'Prot A buffers', 'Protein A buffers + WFI', 'Buf A', 1, [28, 296, 712, 434], ['7', '5', '6'],
      'Post-load wash, wash 1 (arginine pH 9), wash 2 (pH 6.5), elution pH 4.45, strip, 0.5 N NaOH, re-equilibration.'),
    U('pa', 'Protein A', 'Protein A capture: 2 trains', 'C-101', 0, [722, 50, 968, 304], ['9', '10', '11', '12', '13', '14', 'W1'],
      'Columns 1&2 and 3&4 alternate: one loads while the other washes, elutes, strips and cleans.'),
    U('vi', 'Prot A vessel SUM', 'Protein A vessel SUM', 'T-102', 0, [978, 78, 1104, 248], ['14', '15'],
      'Pools eluates; depth filtration starts after three elution volumes.'),
    U('vibuf', 'Depth solutions', 'Depth filtration solutions + WFI', 'Buf DF', 1, [938, 292, 1428, 434], ['17', '16', '18'],
      'Acetate pH 5.0, 0.5 M NaCl (2 M in-line), arginine, 1 N NaOH and Tris flushes (order assumed).'),
    U('depth', 'Depth filtration', 'Depth filtration + surge', 'FH-101', 0, [1112, 78, 1438, 300], ['19', '20', '21'],
      'New filter each batch: WFI flush, equilibrate, filter, chase, post-use flushes.'),
    U('base', 'Base addition', 'On-line base addition', 'SM-201', 1, [18, 470, 182, 640], ['25'],
      '2 M Tris base and 2 M NaCl dosed into the AEX load to set pH and conductivity.'),
    U('aexbuf', 'AEX buffers', 'AEX buffers + WFI', 'Buf AEX', 1, [28, 692, 432, 842], ['23', '22', '24'],
      'Equilibration, chase, strip, 1 N NaOH CIP and benzyl alcohol storage.'),
    U('aex', 'AEX', 'AEX chromatography: 2 columns', 'C-201', 0, [198, 455, 345, 650], ['26', '27', '28', '29', '30'],
      'Flow-through: product passes and is chased out; then strip, CIP, re-equilibrate.'),
    U('aexsum', 'AEX SUM', 'AEX vessel SUM', 'T-201', 0, [350, 462, 488, 628], ['31', '32'],
      'AEX pool hold before acid conditioning.'),
    U('acid', 'Acid addition', 'On-line acid addition', 'SM-202', 1, [482, 470, 588, 640], ['36'],
      '1 M acetic acid dosed into the CEX load to bring it to pH 5.5.'),
    U('cexbuf', 'CEX buffers', 'CEX buffers + WFI', 'Buf CEX', 1, [428, 692, 912, 842], ['34', '33', '35'],
      'Equilibration/wash pH 5.5, elution 280 mM NaCl, strip, 1 N NaOH CIP, storage.'),
    U('cex', 'CEX', 'CEX chromatography: 2 columns', 'C-203', 0, [600, 455, 800, 674], ['37', '38', '39', '40', '41', 'W2'],
      'Bind and elute: load, wash, elute the peak, strip, CIP, re-equilibrate.'),
    U('cexsum', 'CEX SUM', 'CEX vessel SUM', 'T-202', 0, [802, 462, 930, 628], ['42', '43'],
      'Pools CEX eluates; virus filtration starts after three elution volumes.'),
    U('vfbuf', 'VF buffer', 'Virus-filter buffer + WFI', 'Buf VF', 1, [906, 692, 1074, 842], ['45', '44', '46'],
      'WFI flush and equilibration per new filter; chase per batch.'),
    U('vrf', 'Virus filtration', 'Virus retentive filtration', 'F-301', 0, [26, 905, 232, 1072], ['47', '48', '49', '50', '51', '52'],
      'Holders A and B alternate at each filter change.'),
    U('vrfsum', 'VRF SUM', 'VRF vessel SUM', 'T-301', 0, [238, 858, 372, 1018], ['53', '54'],
      'Collects virus filtrate for one lot window.'),
    U('dfbuf', 'UF/DF buffers', 'UF/DF buffers + WFI', 'Buf DF', 1, [388, 1000, 632, 1124], ['56', '55', '57'],
      'Diafiltration pH 5.1 proline, recovery pH 5.3 proline, NaOH CIP and storage.'),
    U('ufdf', 'UF/DF', 'Ultrafiltration / diafiltration', 'UF-301', 0, [470, 904, 650, 1052], ['58', 'W3', '59'],
      'Per lot: transfer in, concentrate, diafilter, concentrate, recover, spike, transfer out, clean.'),
    U('ufsum', 'UF/DF SUM', 'UF/DF vessel SUM', 'T-302', 0, [640, 858, 775, 1018], ['60'],
      'Recirculation vessel for the lot being processed.'),
    U('form', 'Formulation spike', 'Formulation spike', 'Form.', 1, [758, 1000, 842, 1124], ['61'],
      'Acetate, proline, 1.17 M sucrose and poloxamer 188 added to reach the final formulation.'),
    U('bds', 'BDS hold', 'BDS hold mixer', 'M-301', 0, [848, 856, 985, 1018], ['62'],
      'Each lot is mixed, then filtered through 0.2 µm.'),
    U('cryo', 'Cryovault', 'Cryovault containers', 'Cryo', 0, [1004, 926, 1126, 1062], ['63'],
      'Bulk drug substance filled into containers and frozen.'),
  ];

  let J = null, G = null, BASE = null, BASE_K = 0, PLAN = null, RUN = null, LASTURL = null;
  let OWNER = {}, VTU = {}, TKU = {}, COLU = {};

  // ── helpers ──
  function rr(c, x, y, w, h, r) {
    r = Math.min(r, w / 2, h / 2);
    c.moveTo(x + r, y); c.arcTo(x + w, y, x + w, y + h, r); c.arcTo(x + w, y + h, x, y + h, r);
    c.arcTo(x, y + h, x, y, r); c.arcTo(x, y, x + w, y, r); c.closePath();
  }
  function card(c, x, y, w, h, r, fill, blur) {
    c.save();
    c.shadowColor = 'rgba(15,23,42,0.14)'; c.shadowBlur = blur === undefined ? 30 : blur; c.shadowOffsetY = 8;
    c.beginPath(); rr(c, x, y, w, h, r); c.fillStyle = fill || '#fff'; c.fill();
    c.restore();
  }
  function fit(c, s, w) {
    s = String(s);
    if (c.measureText(s).width <= w) return s;
    while (s.length > 1 && c.measureText(s + '…').width > w) s = s.slice(0, -1);
    return s + '…';
  }
  function path(c, pts) { c.beginPath(); c.moveTo(pts[0][0], pts[0][1]); for (let i = 1; i < pts.length; i++) c.lineTo(pts[i][0], pts[i][1]); }
  function tankPath(c, g) {
    c.moveTo(g.x0, g.yTop);
    c.ellipse(g.cx, g.yTop, g.w / 2, g.ry, 0, Math.PI, 2 * Math.PI);
    c.lineTo(g.x1, g.yb);
    c.ellipse(g.cx, g.yb, g.w / 2, g.ry, 0, 0, Math.PI);
    c.closePath();
  }
  const isOn = (st, id) => { const v = st && st.vals[id]; return !!v && num(v.L) > 1e-6; };
  const lit = (u, t) => { const L = PLAN && PLAN.litAt[u.key]; return L === undefined ? 0 : clamp((t - L) / 0.6); };
  const lerpCam = (a, b, q) => ({x: a.x + (b.x - a.x) * q, y: a.y + (b.y - a.y) * q, s: a.s * Math.pow(b.s / a.s, q)});
  const lerpBox = (a, b, q) => a.map((v, i) => v + (b[i] - v) * q);
  function info(t, err) { const e = $('vidInfo'); if (e) { e.textContent = t; e.className = err ? 'err' : 'muted'; } }

  // ── diagram bitmap and geometry ──
  function indexGeo() {
    OWNER = {};
    UNITS.forEach(u => u.act.forEach(id => { if (!(id in OWNER)) OWNER[id] = u; }));
    const inBox = (x, y) => {
      let best = null, ba = Infinity;
      for (const u of UNITS) {
        const [a, b, c, d] = u.box;
        if (x >= a && x <= c && y >= b && y <= d) { const ar = (c - a) * (d - b); if (ar < ba) { ba = ar; best = u; } }
      }
      return best;
    };
    VTU = {}; for (const id in G.vt) VTU[id] = inBox(G.vt[id].x, G.vt[id].y);
    TKU = {}; for (const id in G.tanks) { const g = G.tanks[id]; TKU[id] = inBox(g.cx, (g.yTop + g.yb) / 2); }
    COLU = {}; for (const id in G.cols) { const g = G.cols[id]; COLU[id] = inBox(g.x0 + g.w / 2, g.yTop + g.h / 2); }
    UNITS.forEach(u => { u.vts = Object.keys(VTU).filter(id => VTU[id] === u); });
    G.pipes.forEach(p => {
      p.seg = []; p.len = 0;
      for (let k = 1; k < p.pts.length; k++) {
        const l = Math.hypot(p.pts[k][0] - p.pts[k - 1][0], p.pts[k][1] - p.pts[k - 1][1]);
        p.seg.push(l); p.len += l;
      }
    });
  }

  async function ensureBase(K) {
    if (BASE && BASE_K === K) return;
    const svg = PFD.svg({video: true});
    G = PFD.geo();
    const sized = svg.replace('<svg ', `<svg width="${Math.round(G.W * K)}" height="${Math.round(G.H * K)}" `);
    const img = new Image();
    await new Promise((res, rej) => {
      img.onload = res;
      img.onerror = () => rej(new Error('Could not draw the diagram for the video'));
      img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(sized);
    });
    const cv = document.createElement('canvas');
    cv.width = Math.round(G.W * K); cv.height = Math.round(G.H * K);
    const x = cv.getContext('2d');
    x.fillStyle = '#fff'; x.fillRect(0, 0, cv.width, cv.height);
    x.drawImage(img, 0, 0, cv.width, cv.height);
    BASE = cv; BASE_K = K;
    indexGeo();
  }

  // ── the script of the film ──
  function opts() {
    return {start: $('vidStart').value, hour: +$('vidHour').value || 0, len: $('vidLen').value,
            rate: $('vidRate').value, tour: $('vidTour').checked, res: +$('vidRes').value || 1080};
  }
  function OV() { const s = Math.min((VP.w - 40) / G.W, (VP.h - 30) / G.H); return {x: G.W / 2, y: G.H / 2, s}; }
  function camFor(u) {
    const pad = 26, [a, b, c, d] = u.box, bw = c - a + 2 * pad, bh = d - b + 2 * pad;
    const s = clamp(Math.min((VP.w - 80) / bw, (VP.h - CAP_H - 60) / bh), OV().s * 1.35, 3.1);
    return {x: (a + c) / 2, y: (b + d) / 2 + (CAP_H + 30) / 2 / s, s};
  }

  function buildPlan(o) {
    G = G || PFD.geo();
    const S = PFD.series(), n = S.length;
    const find = f => { const k = S.findIndex(f); return k < 0 ? null : k; };
    const paLoads = Object.keys(LST).filter(k => k.indexOf('pa') === 0)
      .map(k => (LST[k].find(s => s.kind === 'load') || {}).t0).filter(v => v !== undefined);
    const tCap = paLoads.length ? Math.min(...paLoads) : null;   // first Protein A load starts
    const iInoc = find(r => r.phase === 1), iDs = find(r => r.ds_g_h > 0);
    const iCap = tCap !== null ? find(r => r.h >= Math.floor(tCap)) : find(r => r.capture_switch > 0);
    let i0 = 0;
    if (o.start === 'inoc' && iInoc !== null) i0 = iInoc;
    else if (o.start === 'cap' && iCap !== null) i0 = Math.max(0, iCap - 1);
    else if (o.start === 'ds' && iDs !== null) i0 = Math.max(0, iDs - 1);
    else if (o.start === 'custom') i0 = clamp(Math.round(o.hour), 0, n - 1);
    let i1;
    if (o.len === 'startup') i1 = iDs === null ? n - 1 : Math.max(i0, iDs + 6);
    else if (o.len === 'all') i1 = n - 1;
    else i1 = i0 + Math.max(1, +o.len) - 1;
    i1 = Math.min(n - 1, i1);

    const ST = [];
    for (let i = i0; i <= i1; i++) ST.push(PFD.state(i));
    const firstS = {}, first = {};
    ST.forEach((st, k) => { for (const id in st.vals) if (firstS[id] === undefined && isOn(st, id)) firstS[id] = k; });
    UNITS.forEach(u => {
      let f = null;
      u.act.forEach(id => { if (firstS[id] !== undefined && (f === null || firstS[id] < f)) f = firstS[id]; });
      first[u.key] = f;
    });

    const H0 = S[i0].h, H1 = S[i1].h + 1, hours = H1 - H0;
    const rate = o.rate === 'auto' ? hours / clamp(8 + 4 * Math.log2(Math.max(hours, 1)), 8, 45) : Math.max(+o.rate, 0.01);

    const segs = [], litAt = {};
    let T = 0, cam = OV(), last = null, clock = H0;
    const push = s => { s.t0 = T; T += s.dur; s.t1 = T; segs.push(s); last = s; return s; };
    push({type: 'title', dur: TITLE_S, clock: H0, from: cam, to: cam});
    const groups = {};
    UNITS.forEach(u => { const f = first[u.key]; if (f !== null) (groups[f] = groups[f] || []).push(u); });
    Object.keys(groups).map(Number).sort((a, b) => a - b).forEach(k => {
      const hk = H0 + k;
      if (hk > clock) {
        const full = (hk - clock) / rate;
        push({type: 'play', dur: Math.min(full, MAX_GAP_S), c0: clock, c1: hk, from: cam, to: cam,
              ff: full > MAX_GAP_S});
        clock = hk;
      }
      if (o.tour) {
        groups[k].forEach(u => {
          const to = camFor(u);
          const s = push({type: 'focus', dur: FOCUS_S, clock, unit: u, from: cam, to,
                          prevUnit: last.type === 'focus' ? last.unit : null});
          litAt[u.key] = s.t0 + FLY_S * 0.7;
          cam = to;
        });
        const fu = last.unit;
        push({type: 'out', dur: OUT_S, clock, from: cam, to: OV(), prevUnit: fu});
        cam = OV();
      } else groups[k].forEach(u => { litAt[u.key] = T; });
    });
    if (H1 > clock) push({type: 'play', dur: (H1 - clock) / rate, c0: clock, c1: H1, from: cam, to: cam});
    push({type: 'end', dur: END_S, clock: H1 - 1e-6, from: cam, to: cam});
    segs.forEach((s, i) => { s.nextFocus = !!(segs[i + 1] && segs[i + 1].type === 'focus'); });

    let g = 0, ds = 0, buf = 0, ww = 0;
    ST.forEach(st => { g += st.r.harvest_g_h; ds += st.r.ds_g_h; buf += st.r.buffer_L_h; ww += st.r.ww_total_Lh; });
    const cyc = num(S[i1].cycles) - (i0 > 0 ? num(S[i0 - 1].cycles) : 0);
    const p = J.params;
    return {
      S, ST, i0, i1, H0, H1, rate, segs, T, litAt, first, firstS, tour: o.tour,
      hN: iInoc === null ? null : S[iInoc].h,
      fcUsed: !!J.forecaster_used,
      title: (J.design && J.design.molecule_label) || 'Legacy ATF perfusion model',
      subtitle: [p.upstream.mode === 'mechanistic' ? 'Legacy ATF perfusion' : 'IDB titre curve',
                 p.capture.mode === 'pcc' ? '3-column PCC' : 'Protein A 2 trains',
                 p.vi && p.vi.mode === 'low_ph' ? 'low-pH VI' : 'Brij-9 VI',
                 p.polish.order === 'cex_aex' ? 'CEX → AEX' : 'AEX → CEX',
                 J.forecaster_used ? 'NevAI on' : 'NevAI off'].join(' · '),
      sum: {g, ds, buf, ww, cyc},
    };
  }

  function segAt(t) {
    const s = PLAN.segs;
    for (let i = 0; i < s.length; i++) if (t < s[i].t1) return s[i];
    return s[s.length - 1];
  }

  // ── drawing ──
  function dots(c, p, t) {
    const gap = 12, ph = (t * 30) % gap;
    c.fillStyle = '#fff';
    let si = 0, acc = 0;
    for (let d = ph; d < p.len; d += gap) {
      while (si < p.seg.length - 1 && d > acc + p.seg[si]) { acc += p.seg[si]; si++; }
      const l = p.seg[si] || 1, f = (d - acc) / l, a = p.pts[si], b = p.pts[si + 1];
      c.beginPath(); c.arc(a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, 1.25, 0, 6.2832); c.fill();
    }
  }

  function drawPipes(c, t, st, k) {
    c.lineJoin = 'round'; c.lineCap = 'round';
    for (const p of G.pipes) {
      if (!p.id) continue;
      const f = PLAN.firstS[p.id];
      if (f === undefined || f > k) continue;
      const own = OWNER[p.id], a = own ? lit(own, t) : 1;
      if (a <= 0) continue;
      c.save(); c.globalAlpha = a;
      if (isOn(st, p.id)) {
        const col = FCOL[p.k] || TEAL;
        c.setLineDash([]);
        c.strokeStyle = hexA(p.k === 'pp' ? TEAL2 : col, 0.22); c.lineWidth = 8; path(c, p.pts); c.stroke();
        c.strokeStyle = col; c.lineWidth = 3; path(c, p.pts); c.stroke();
        dots(c, p, t);
      } else {
        c.strokeStyle = KCOL[p.k] || INK; c.lineWidth = 1.6; c.setLineDash(KDASH[p.k] || []);
        path(c, p.pts); c.stroke();
      }
      c.restore();
    }
  }

  function drawEquipment(c, t, st) {
    for (const id in st.eq.lv) {
      const f = st.eq.lv[id], g = G.tanks[id], u = TKU[id];
      if (!g || f <= 0 || !u) continue;
      const a = lit(u, t); if (a <= 0) continue;
      c.save(); c.globalAlpha = a;
      c.beginPath(); tankPath(c, g); c.clip();
      const y = g.base - f * g.span, grd = c.createLinearGradient(0, y, 0, g.base);
      grd.addColorStop(0, 'rgba(20,184,166,0.45)'); grd.addColorStop(1, 'rgba(14,124,123,0.30)');
      c.fillStyle = grd; c.fillRect(g.x0, y, g.w, g.base - y);
      c.strokeStyle = 'rgba(14,124,123,0.85)'; c.lineWidth = 1;
      c.beginPath(); c.moveTo(g.x0, y); c.lineTo(g.x1, y); c.stroke();
      c.restore();
    }
    for (const id in st.eq.col) {
      const cls = st.eq.col[id], g = G.cols[id], u = COLU[id];
      if (!cls || !g || !u) continue;
      const a = lit(u, t); if (a <= 0) continue;
      c.save();
      c.beginPath(); rr(c, g.x0, g.yTop, g.w, g.h, 4);
      c.globalAlpha = a * 0.85; c.fillStyle = cls === 'load' ? '#99e0d8' : '#fde68a'; c.fill();
      c.globalAlpha = a; c.strokeStyle = '#1d2733'; c.lineWidth = 1.3; c.stroke();
      if (g.label) { c.fillStyle = INK; c.font = `700 8px ${FONT}`; c.textAlign = 'center'; c.fillText(g.label, g.x0 + g.w / 2, g.yTop + g.h / 2 + 3); }
      c.restore();
    }
    c.font = `700 9px ${FONT}`; c.lineJoin = 'round';
    for (const id in st.eq.vt) {
      const [l1, l2] = st.eq.vt[id], g = G.vt[id], u = VTU[id];
      if (!g || !u || !(l1 || l2)) continue;
      const a = lit(u, t); if (a <= 0) continue;
      c.save(); c.globalAlpha = a;
      c.textAlign = g.anchor === 'start' ? 'left' : g.anchor === 'end' ? 'right' : 'center';
      [l1, l2].forEach((s, j) => {
        if (!s) return;
        c.strokeStyle = '#fff'; c.lineWidth = 3; c.strokeText(s, g.x, g.y + j * 10);
        c.fillStyle = '#0b5e5d'; c.fillText(s, g.x, g.y + j * 10);
      });
      c.restore();
    }
  }

  function drawTags(c, t, st) {
    c.font = `700 7.5px ${FONT}`; c.textAlign = 'center';
    for (const id in G.tags) {
      if (!isOn(st, id)) continue;
      const own = OWNER[id], a = own ? lit(own, t) : 1;
      if (a <= 0) continue;
      const w = 7 + String(id).length * 3.2;
      for (const [x, y] of G.tags[id]) {
        c.save(); c.globalAlpha = a;
        c.beginPath(); c.moveTo(x - w, y); c.lineTo(x, y - 8); c.lineTo(x + w, y); c.lineTo(x, y + 8); c.closePath();
        c.fillStyle = '#ccefeb'; c.fill(); c.strokeStyle = TEAL; c.lineWidth = 0.9; c.stroke();
        c.fillStyle = INK; c.fillText(id, x, y + 2.8);
        c.restore();
      }
    }
  }

  function drawSpot(c, sp, cam) {
    const pad = 6, [x0, y0, x1, y1] = sp.box;
    c.save();
    c.beginPath(); c.rect(-3000, -3000, G.W + 6000, G.H + 6000);
    rr(c, x0 - pad, y0 - pad, x1 - x0 + 2 * pad, y1 - y0 + 2 * pad, 14);
    c.fillStyle = `rgba(15,23,42,${0.42 * sp.dim})`; c.fill('evenodd');
    if (sp.ringA > 0) {
      const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
      c.translate(cx, cy); c.scale(sp.ringS, sp.ringS); c.translate(-cx, -cy);
      c.beginPath(); rr(c, x0 - pad, y0 - pad, x1 - x0 + 2 * pad, y1 - y0 + 2 * pad, 14);
      c.globalAlpha = sp.ringA; c.lineWidth = 3.2 / cam.s; c.strokeStyle = TEAL2;
      c.shadowColor = 'rgba(20,184,166,0.9)'; c.shadowBlur = 26; c.stroke();
    }
    c.restore();
  }

  function drawPFD(c, t, cam, st, k, spot) {
    card(c, VP.x, VP.y, VP.w, VP.h, 30);
    c.save();
    c.beginPath(); rr(c, VP.x, VP.y, VP.w, VP.h, 30); c.clip();
    c.translate(VP.x + VP.w / 2, VP.y + VP.h / 2); c.scale(cam.s, cam.s); c.translate(-cam.x, -cam.y);
    c.imageSmoothingEnabled = true; c.imageSmoothingQuality = 'high';
    c.drawImage(BASE, 0, 0, G.W, G.H);
    c.fillStyle = 'rgba(243,246,250,0.82)'; c.fillRect(-4, -4, G.W + 8, G.H + 8);   // not yet in the circuit
    for (const u of UNITS) {
      const a = lit(u, t); if (a <= 0) continue;
      const [x0, y0, x1, y1] = u.box;
      c.save(); c.globalAlpha = a;
      c.beginPath(); rr(c, x0, y0, x1 - x0, y1 - y0, 10); c.clip();
      c.drawImage(BASE, 0, 0, G.W, G.H);
      c.restore();
    }
    drawPipes(c, t, st, k);
    drawEquipment(c, t, st);
    drawTags(c, t, st);
    if (spot && spot.dim > 0) drawSpot(c, spot, cam);
    c.restore();
  }

  function live(u, st) {
    const parts = [];
    (u.vts || []).forEach(id => { const v = st.eq.vt[id]; if (v && (v[0] || v[1])) parts.push([v[0], v[1]].filter(Boolean).join(' · ')); });
    u.act.forEach(id => {
      const v = st.vals[id];
      if (v && num(v.L) > 1e-6 && parts.length < 4)
        parts.push(`stream ${id}: ${fmt(num(v.L))} L/h` + (v.g !== undefined && v.g !== null && num(v.g) > 0 ? ` · ${fmt(num(v.g))} g/h` : ''));
    });
    return parts.length ? parts.join('   |   ') : 'joins the circuit';
  }

  function island(c, o) {
    c.save(); c.globalAlpha = o.a;
    c.font = `600 19px ${FONT}`;
    const tw = c.measureText(o.text).width, fw = tw + 84, w = 140 + (fw - 140) * o.w, h = 46;
    const x = VP.x + VP.w / 2 - w / 2, y = VP.y + 18;
    c.beginPath(); rr(c, x, y, w, h, 23); c.fillStyle = '#0b0f14'; c.fill();
    if (o.w > 0.6) {
      c.globalAlpha = o.a * clamp((o.w - 0.6) / 0.4);
      c.beginPath(); c.arc(x + 26, y + h / 2, 6, 0, 6.2832); c.fillStyle = '#34d399'; c.fill();
      c.fillStyle = '#fff'; c.textAlign = 'left'; c.fillText(o.text, x + 44, y + h / 2 + 7);
    }
    c.restore();
  }

  function caption(c, u, st, a, dy) {
    const x = VP.x + 28, w = VP.w - 56, h = CAP_H - 24, y = VP.y + VP.h - CAP_H + dy;
    c.save(); c.globalAlpha = a;
    c.shadowColor = 'rgba(15,23,42,0.22)'; c.shadowBlur = 40; c.shadowOffsetY = 12;
    c.beginPath(); rr(c, x, y, w, h, 30); c.fillStyle = 'rgba(255,255,255,0.95)'; c.fill();
    c.shadowColor = 'transparent'; c.shadowBlur = 0; c.shadowOffsetY = 0;
    const isz = 104, ix = x + 24, iy = y + (h - isz) / 2, grd = c.createLinearGradient(ix, iy, ix + isz, iy + isz);
    if (u.util) { grd.addColorStop(0, '#a78bfa'); grd.addColorStop(1, '#6d28d9'); }
    else { grd.addColorStop(0, '#2dd4bf'); grd.addColorStop(1, '#0e7c7b'); }
    c.beginPath(); rr(c, ix, iy, isz, isz, 24); c.fillStyle = grd; c.fill();
    c.fillStyle = '#fff'; c.textAlign = 'center';
    c.font = `700 ${u.tag.length > 5 ? 19 : 24}px ${FONT}`; c.fillText(u.tag, ix + isz / 2, iy + isz / 2 + 6);
    c.font = `600 11px ${FONT}`; c.globalAlpha = a * 0.85;
    c.fillText(u.util ? 'UTILITY' : 'EQUIPMENT', ix + isz / 2, iy + isz - 14);
    c.globalAlpha = a;
    const k = PLAN.first[u.key], hj = k === null ? null : PLAN.H0 + k;
    const chip = hj === null ? '' : `Day ${(hj / 24).toFixed(2)} · hour ${hj}`;
    c.font = `600 17px ${FONT}`;
    const cw = c.measureText(chip).width + 36, cx = x + w - 24 - cw, cy = y + 22;
    if (chip) {
      c.beginPath(); rr(c, cx, cy, cw, 36, 18); c.fillStyle = '#e6f5f3'; c.fill();
      c.fillStyle = TEAL; c.textAlign = 'center'; c.fillText(chip, cx + cw / 2, cy + 24);
    }
    const tx = ix + isz + 26, tw = x + w - 24 - tx;
    c.textAlign = 'left';
    c.fillStyle = INK; c.font = `700 32px ${FONT}`; c.fillText(fit(c, u.title, tw - cw - 20), tx, y + 52);
    c.fillStyle = GREY; c.font = `400 21px ${FONT}`; c.fillText(fit(c, u.desc, tw), tx, y + 88);
    c.fillStyle = TEAL; c.font = `600 21px ${FONT}`; c.fillText(fit(c, live(u, st), tw), tx, y + 122);
    c.restore();
  }

  function hud(c, t, clock, st) {
    const P = PLAN;
    card(c, HUD.x, HUD.y, HUD.w, HUD.h, 30);
    const x = HUD.x + 32, xr = HUD.x + HUD.w - 32, W = HUD.w - 64;
    let y = HUD.y + 50;
    c.textAlign = 'left';
    c.fillStyle = TEAL; c.font = `700 14px ${FONT}`; c.fillText('IDB mAb DIGITAL TWIN', x, y);
    y += 36; c.fillStyle = INK; c.font = `700 28px ${FONT}`; c.fillText(fit(c, P.title, W), x, y);
    y += 26; c.fillStyle = GREY; c.font = `400 16px ${FONT}`; c.fillText(fit(c, P.subtitle, W), x, y);
    y += 74; c.fillStyle = INK; c.font = `700 56px ${FONT}`; c.fillText('Day ' + (clock / 24).toFixed(2), x, y);
    y += 32; c.fillStyle = GREY; c.font = `500 18px ${FONT}`;
    const nd = P.hN === null ? null : (clock - P.hN) / 24;
    c.fillText((nd === null || nd < 0 ? 'N-1 seed phase' : 'N day ' + nd.toFixed(2)) + ' · hour ' + Math.floor(clock), x, y);
    y += 22;
    const pr = clamp((clock - P.H0) / Math.max(P.H1 - P.H0, 1e-9));
    c.beginPath(); rr(c, x, y, W, 8, 4); c.fillStyle = '#e5e7eb'; c.fill();
    if (pr > 0) { c.beginPath(); rr(c, x, y, Math.max(8, W * pr), 8, 4); c.fillStyle = TEAL2; c.fill(); }
    y += 30; c.font = `500 13px ${FONT}`; c.fillStyle = GREY;
    c.fillText(`Window: hour ${P.H0} → ${P.H1} (${P.H1 - P.H0} h)`, x, y);
    y += 20;
    const r = st.r, rows = [
      ['Permeate titre', st.r0.titre_harvest === null || st.r0.titre_harvest === undefined ? '—' : fmt(r.titre_harvest) + ' g/L'],
      ['NevAI +24 h', (st.fc === null || st.fc === undefined ? (P.fcUsed ? '—' : 'off') : fmt(st.fc) + ' g/L') +
                      (st.watch ? '  ⚠' : '')],
      ['Surge tank', `${fmt(r.sump_L)} L · ${fmt(r.sump_g)} g`],
      ['Protein A cycles', String(Math.round(r.cycles))],
      ['Buffer use', fmt(r.buffer_L_h) + ' L/h'],
      ['Wastewater', fmt(r.ww_total_Lh) + ' L/h'],
      ['Drug substance', fmt(r.ds_g_cum / 1000, 3) + ' kg'],
    ];
    c.beginPath(); rr(c, x - 12, y, W + 24, rows.length * 34 + 10, 16); c.fillStyle = '#f4f6f9'; c.fill();
    rows.forEach(([l, v], j) => {
      const yy = y + 29 + j * 34;
      c.textAlign = 'left'; c.fillStyle = GREY; c.font = `400 17px ${FONT}`; c.fillText(l, x + 4, yy);
      c.textAlign = 'right'; c.fillStyle = INK; c.font = `600 17px ${FONT}`; c.fillText(v, xr - 4, yy);
      if (j < rows.length - 1) { c.fillStyle = '#e3e7ee'; c.fillRect(x + 4, yy + 11, W - 8, 1); }
    });
    y += rows.length * 34 + 46;
    const joined = UNITS.filter(u => lit(u, t) > 0.5).length;
    c.textAlign = 'left'; c.fillStyle = TEAL; c.font = `700 14px ${FONT}`; c.fillText('IN CIRCUIT', x, y);
    c.textAlign = 'right'; c.fillStyle = GREY; c.font = `500 14px ${FONT}`; c.fillText(`${joined} of ${UNITS.length} units`, xr, y);
    y += 10;
    const colW = W / 2, half = Math.ceil(UNITS.length / 2);
    UNITS.forEach((u, j) => {
      const ux = x + (j < half ? 0 : colW), uy = y + 24 + (j % half) * 26;
      const a = lit(u, t), now = a > 0 && u.act.some(id => isOn(st, id));
      c.beginPath(); c.arc(ux + 7, uy - 5, 6, 0, 6.2832);
      if (a > 0) {
        if (now) { c.fillStyle = TEAL2; c.fill(); } else { c.strokeStyle = TEAL2; c.lineWidth = 2; c.stroke(); }
      } else { c.fillStyle = '#d1d5db'; c.fill(); }
      c.textAlign = 'left'; c.fillStyle = a > 0 ? INK : '#9ca3af';
      c.font = `${a > 0 ? 600 : 400} 15px ${FONT}`; c.fillText(fit(c, u.short, colW - 30), ux + 22, uy);
    });
    const ly = HUD.y + HUD.h - 28;
    c.font = `500 13px ${FONT}`; c.textAlign = 'left';
    let lx = x;
    [['fill', TEAL2, 'flowing now'], ['ring', TEAL2, 'in circuit'], ['fill', '#d1d5db', 'not yet']].forEach(([m, col, s]) => {
      c.beginPath(); c.arc(lx + 6, ly - 4, 5, 0, 6.2832);
      if (m === 'fill') { c.fillStyle = col; c.fill(); } else { c.strokeStyle = col; c.lineWidth = 2; c.stroke(); }
      c.fillStyle = GREY; c.fillText(s, lx + 16, ly); lx += c.measureText(s).width + 40;
    });
  }

  function titleCard(c, lt, dur) {
    const P = PLAN, a = 1 - clamp((lt - (dur - 0.6)) / 0.6);
    if (a <= 0) return;
    c.save(); c.globalAlpha = a;
    const w = 900, h = 300, x = VP.x + (VP.w - w) / 2;
    const y = VP.y + (VP.h - h) / 2 - 40 + (1 - easeOutBack(clamp(lt / 0.6))) * 30;
    card(c, x, y, w, h, 36, 'rgba(255,255,255,0.97)', 50);
    c.textAlign = 'center';
    c.fillStyle = TEAL; c.font = `700 16px ${FONT}`; c.fillText('INTEGRATED DYNAMIC BIOPROCESSING', x + w / 2, y + 62);
    c.fillStyle = INK; c.font = `700 44px ${FONT}`; c.fillText('Dynamic start-up of the mAb line', x + w / 2, y + 120);
    c.fillStyle = GREY; c.font = `400 22px ${FONT}`;
    c.fillText(fit(c, P.title + ' · ' + P.subtitle, w - 80), x + w / 2, y + 166);
    c.fillText(`Campaign day ${(P.H0 / 24).toFixed(2)} → ${(P.H1 / 24).toFixed(2)}  ·  ${P.H1 - P.H0} h of operation`, x + w / 2, y + 202);
    c.fillStyle = TEAL; c.font = `600 19px ${FONT}`;
    c.fillText(P.tour ? 'Each unit and utility zooms into view as it joins the circuit'
                      : 'Units and utilities light up as they join the circuit', x + w / 2, y + 252);
    c.restore();
  }

  function endCard(c, t, lt) {
    const P = PLAN, a = clamp(lt / 0.5);
    c.save(); c.globalAlpha = a;
    const w = 940, h = 400, x = VP.x + (VP.w - w) / 2, y = VP.y + (VP.h - h) / 2 + (1 - easeOutBack(a)) * 30;
    card(c, x, y, w, h, 36, 'rgba(255,255,255,0.97)', 50);
    c.textAlign = 'center';
    c.fillStyle = TEAL; c.font = `700 16px ${FONT}`; c.fillText('WINDOW SUMMARY', x + w / 2, y + 56);
    c.fillStyle = INK; c.font = `700 36px ${FONT}`;
    c.fillText(`Campaign day ${(P.H0 / 24).toFixed(2)} → ${(P.H1 / 24).toFixed(2)}`, x + w / 2, y + 102);
    const joined = UNITS.filter(u => lit(u, t) > 0.5).length;
    const rows = [['Product harvested', fmt(P.sum.g) + ' g'], ['Drug substance', fmt(P.sum.ds) + ' g'],
                  ['Protein A cycles', String(Math.round(P.sum.cyc))], ['Buffer used', fmt(P.sum.buf) + ' L'],
                  ['Wastewater', fmt(P.sum.ww) + ' L'], ['Units in circuit', `${joined} of ${UNITS.length}`]];
    const tw = (w - 120) / 3, th = 100;
    rows.forEach(([l, v], j) => {
      const tx = x + 40 + (j % 3) * (tw + 20), ty = y + 140 + Math.floor(j / 3) * (th + 20);
      c.beginPath(); rr(c, tx, ty, tw, th, 20); c.fillStyle = '#f4f6f9'; c.fill();
      c.textAlign = 'left'; c.fillStyle = GREY; c.font = `500 16px ${FONT}`; c.fillText(l, tx + 20, ty + 34);
      c.fillStyle = INK; c.font = `700 30px ${FONT}`; c.fillText(v, tx + 20, ty + 76);
    });
    c.restore();
  }

  // ── steps running now, from the run's step log ──
  const SOLC = {wfi: '#7dd3fc', brij_10: '#a855f7', tris_74: '#0e7c7b', pa_wash1: '#f59e0b', pa_wash2: '#fbbf24',
    pa_elute: '#dc2626', pa_strip: '#f43f5e', naoh_05n: '#475569', naoh_1n: '#334155', naoh_01n: '#64748b',
    storage: '#16a34a', df_ac50: '#fb923c', nacl_05m: '#60a5fa', nacl_2m: '#3b82f6', df_arg: '#c084fc',
    tris_base_2m: '#14b8a6', tris_nacl_74: '#2563eb', acetic_1m: '#e11d48', cex_eq: '#eab308',
    cex_elute: '#be123c', uf_df1: '#8b5cf6', uf_df2: '#6366f1', uf_spike: '#db2777', vi_acid: '#b91c1c', vi_base: '#0891b2'};
  const KINDC = {load: '#64748b', hold: '#94a3b8', perm: '#a5b4fc', out: '#0ea5e9'};
  let LST = {}, LLAB = {};
  function indexSteps(j) {
    LST = {}; LLAB = {};
    (j.lanes || []).forEach(l => { LLAB[l.key] = l.label; LST[l.key] = []; });
    (j.steps || []).forEach(s => { (LST[s.lane] = LST[s.lane] || []).push(s); });
    for (const k in LST) LST[k].sort((a, b) => a.t0 - b.t0);
  }
  function stepsAt(t) {
    const out = [];
    for (const k of Object.keys(LLAB)) {
      const a = LST[k] || [];
      let lo = 0, hi = a.length - 1, f = -1;
      while (lo <= hi) { const m = (lo + hi) >> 1; if (a[m].t0 <= t) { f = m; lo = m + 1; } else hi = m - 1; }
      for (let i = f; i >= 0 && i > f - 2; i--) if (a[i].t1 > t) { out.push({label: LLAB[k], s: a[i]}); break; }
    }
    return out;
  }
  function stepStrip(c, t, a) {
    if (a <= 0.01 || !Object.keys(LLAB).length) return;
    const now = stepsAt(t), cols = 3, rowH = 34, rows = Math.max(1, Math.ceil(now.length / cols));
    const w = VP.w - 56, h = 52 + rows * rowH, x = VP.x + 28, y = VP.y + VP.h - h - 20;
    c.save(); c.globalAlpha = a;
    c.shadowColor = 'rgba(15,23,42,0.18)'; c.shadowBlur = 30; c.shadowOffsetY = 10;
    c.beginPath(); rr(c, x, y, w, h, 24); c.fillStyle = 'rgba(255,255,255,0.94)'; c.fill();
    c.shadowColor = 'transparent'; c.shadowBlur = 0; c.shadowOffsetY = 0;
    c.textAlign = 'left'; c.fillStyle = TEAL; c.font = `700 14px ${FONT}`; c.fillText('STEPS RUNNING NOW', x + 24, y + 30);
    c.textAlign = 'right'; c.fillStyle = GREY; c.font = `500 14px ${FONT}`;
    c.fillText(now.length ? `${now.length} unit${now.length > 1 ? 's' : ''} busy` : 'all units idle', x + w - 24, y + 30);
    const cw = (w - 48) / cols;
    now.forEach((o, i) => {
      const cx = x + 24 + (i % cols) * cw, cy = y + 46 + Math.floor(i / cols) * rowH, bw = cw - 12, bh = rowH - 8;
      const col = (o.s.sol && SOLC[o.s.sol]) || KINDC[o.s.kind] || '#94a3b8';
      const p = clamp((t - o.s.t0) / Math.max(o.s.t1 - o.s.t0, 1e-6));
      c.beginPath(); rr(c, cx, cy, bw, bh, 8); c.fillStyle = '#f4f6f9'; c.fill();
      c.beginPath(); rr(c, cx, cy, Math.max(10, bw * p), bh, 8); c.fillStyle = hexA(col, 0.30); c.fill();
      c.beginPath(); c.arc(cx + 14, cy + bh / 2, 6, 0, 6.2832); c.fillStyle = col; c.fill();
      c.textAlign = 'left'; c.fillStyle = INK; c.font = `600 15px ${FONT}`;
      c.fillText(fit(c, o.label + ' · ' + o.s.step, bw - 36), cx + 28, cy + bh / 2 + 5);
    });
    c.restore();
  }

  // drawn last in frame(), right-aligned in the HUD header row: never cropped, never covered
  function watermark(c) {
    c.save();
    c.font = `700 15px ${FONT}`;
    c.textAlign = 'right';
    c.fillStyle = INK;
    c.globalAlpha = 0.8;
    c.fillText(WATERMARK, HUD.x + HUD.w - 32, HUD.y + 50);   // same baseline as "IDB mAb DIGITAL TWIN"
    c.restore();
  }

  function frame(c, t) {
    const P = PLAN, SC = c.canvas.width / VW;
    const s = segAt(t), lt = t - s.t0, p = s.dur > 0 ? clamp(lt / s.dur) : 1;
    const clock = s.type === 'play' ? s.c0 + (s.c1 - s.c0) * p : s.clock;
    const k = clamp(Math.floor(clock - P.H0 + 1e-6), 0, P.ST.length - 1);
    const st = P.ST[k];
    let cam = s.to, spot = null, cap = null, isl = null;
    if (s.type === 'focus') {
      const q = easeIO(clamp(lt / FLY_S)), rt = lt - FLY_S;
      cam = lerpCam(s.from, s.to, q);
      spot = {box: s.prevUnit ? lerpBox(s.prevUnit.box, s.unit.box, q) : s.unit.box, dim: s.prevUnit ? 1 : q,
              ringA: rt < 0 ? 0 : clamp(rt / 0.15), ringS: rt < 0 ? 1 : 1 + 0.07 * Math.exp(-5 * rt) * Math.cos(11 * rt)};
      const ain = clamp((lt - FLY_S * 0.6) / 0.35), aout = s.nextFocus ? clamp((s.dur - lt) / 0.2) : 1;
      cap = {u: s.unit, a: Math.min(ain, aout), dy: (1 - easeOutBack(ain)) * 40};
      isl = {text: 'Joining the circuit · ' + s.unit.tag, a: 1, w: easeOutBack(clamp(lt / 0.5))};
    } else if (s.type === 'out') {
      cam = lerpCam(s.from, s.to, easeIO(p));
      const f = 1 - clamp(lt / 0.45);
      spot = {box: s.prevUnit.box, dim: f, ringA: f, ringS: 1};
      cap = {u: s.prevUnit, a: 1 - clamp(lt / 0.3), dy: 0};
      isl = {text: 'Joining the circuit · ' + s.prevUnit.tag, a: 1 - clamp(lt / 0.4), w: 1};
    }
    if (s.type === 'play' && s.ff)
      isl = {text: `Fast-forward · next unit joins at hour ${s.c1}`, a: 1, w: easeOutBack(clamp(lt / 0.5))};
    const wThr = num(J.params.control && J.params.control.nevai_watch_threshold) || 0.12;
    if (s.type === 'play' && !s.ff && st.r0.nevai_watch !== null && st.r0.nevai_watch !== undefined &&
        st.r0.nevai_watch > wThr)
      isl = {text: 'NevAI watch · permeate titre off forecast', a: 1, w: 1};
    c.setTransform(SC, 0, 0, SC, 0, 0);
    const bg = c.createLinearGradient(0, 0, 0, VH);
    bg.addColorStop(0, '#eef1f6'); bg.addColorStop(1, '#dde3ec');
    c.fillStyle = bg; c.fillRect(0, 0, VW, VH);
    drawPFD(c, t, cam, st, k, spot);
    if (s.type !== 'title' && s.type !== 'end') stepStrip(c, clock, 1 - (cap ? clamp(cap.a) : 0));
    if (isl && isl.a > 0) island(c, isl);
    if (cap && cap.a > 0) caption(c, cap.u, st, cap.a, cap.dy);
    hud(c, t, clock, st);
    if (s.type === 'title') titleCard(c, lt, s.dur);
    if (s.type === 'end') endCard(c, t, lt);
    watermark(c);
  }

  // ── preview and recording ──
  function pickMime() {
    if (!window.MediaRecorder) return '';
    const c = ['video/mp4;codecs=avc1.640028', 'video/mp4;codecs=avc1.42E028', 'video/mp4;codecs=avc1', 'video/mp4',
               'video/webm;codecs=vp9', 'video/webm;codecs=vp8', 'video/webm'];
    return c.find(m => { try { return MediaRecorder.isTypeSupported(m); } catch (e) { return false; } }) || '';
  }
  function setBusy(b) {
    $('vidPreview').disabled = b; $('vidRecord').disabled = b; $('vidStop').disabled = !b;
  }
  function stop() { if (RUN) RUN.stop = true; }

  function loop(c, record) {
    return new Promise(resolve => {
      const t0 = performance.now(), T = PLAN.T;
      const step = now => {
        if (!RUN || RUN.stop) return resolve(false);
        const t = Math.min((now - t0) / 1000, T);
        frame(c, t);
        info(`${record ? '● Recording' : 'Preview'} ${mmss(t)} / ${mmss(T)}` + (record ? ' — keep this tab visible' : ''));
        if (t >= T) return resolve(true);
        requestAnimationFrame(step);
      };
      requestAnimationFrame(step);
    });
  }

  function offer(blob, mime, complete) {
    const ext = mime.indexOf('video/mp4') === 0 ? 'mp4' : 'webm';
    if (LASTURL) URL.revokeObjectURL(LASTURL);
    LASTURL = URL.createObjectURL(blob);
    const mol = (J.params.upstream && J.params.upstream.molecule) || 'run';
    const name = `IDB_${mol}_h${PLAN.H0}-${PLAN.H1}.${ext}`;
    const notes = [`${(blob.size / 1e6).toFixed(1)} MB · ${ext.toUpperCase()}`];
    if (ext === 'webm') notes.push(`for iPhone/Keynote convert with: ffmpeg -i ${name} -c:v libx264 -pix_fmt yuv420p ${name.replace('.webm', '.mp4')}`);
    if (RUN && RUN.hidden) notes.push('the tab was hidden during recording, so some frames may be frozen');
    if (!complete) notes.push('stopped early');
    $('vidOut').innerHTML = `<video controls playsinline src="${LASTURL}"></video>` +
      `<div class="vidlinks"><a class="btnlink" href="${LASTURL}" download="${name}">⬇ Download ${name}</a>` +
      `<span class="muted">${notes.join(' · ')}</span></div>`;
    info('Video ready.');
  }

  async function start(record) {
    if (RUN) return;
    if (!J || !window.PFD || !PFD.length()) { info('Run a simulation first.', true); return; }
    const o = opts(), cv = $('vidCanvas');
    cv.width = o.res === 720 ? 1280 : 1920; cv.height = o.res === 720 ? 720 : 1080;
    RUN = {stop: false, hidden: false};
    setBusy(true);
    try {
      info('Preparing the diagram…');
      PFD.pause();
      await ensureBase(o.res === 720 ? 2 : 2.6);
      PLAN = buildPlan(o);
      const c = cv.getContext('2d');
      frame(c, 0);
      let rec = null, chunks = [], mime = '';
      if (record) {
        mime = pickMime();
        if (!mime) throw new Error('This browser cannot record video (no MediaRecorder). Use Safari, Chrome, Edge or Firefox.');
        const stream = cv.captureStream(30);
        rec = new MediaRecorder(stream, {mimeType: mime, videoBitsPerSecond: o.res === 720 ? 6e6 : 12e6});
        rec.ondataavailable = e => { if (e.data && e.data.size) chunks.push(e.data); };
        rec.start(1000);
      }
      const done = await loop(c, record);
      if (rec) {
        await new Promise(r => setTimeout(r, 300));
        await new Promise(r => { rec.onstop = r; rec.stop(); });
        offer(new Blob(chunks, {type: mime.split(';')[0]}), mime, done);
      } else info(done ? 'Preview finished. Press Record video to save it.' : 'Stopped.');
    } catch (e) {
      info(e && e.name === 'SecurityError' ? 'The browser blocked recording of this canvas; try Chrome or Safari.' : (e.message || String(e)), true);
    } finally {
      RUN = null; setBusy(false);
    }
  }

  function estimate() {
    if (RUN || !J || !window.PFD || !PFD.length()) return;
    try {
      const p = buildPlan(opts());
      const n = Object.values(p.first).filter(v => v !== null).length;
      info(`Film ≈ ${mmss(p.T)} for ${p.H1 - p.H0} h of operation (campaign day ${(p.H0 / 24).toFixed(2)} → ` +
           `${(p.H1 / 24).toFixed(2)}); ${n} of ${UNITS.length} units join.` +
           (p.T > 600 ? ' Long film: pick a faster speed or a shorter window.' : '') +
           (n <= 2 ? ' Few units run in this window: start at N inoculation or pick a longer length.' : ''));
    } catch (e) { info(e.message, true); }
  }

  async function poster() {
    if (RUN || !J || !window.PFD || !PFD.length()) return;
    try {
      const o = opts(), cv = $('vidCanvas');
      cv.width = o.res === 720 ? 1280 : 1920; cv.height = o.res === 720 ? 720 : 1080;
      await ensureBase(o.res === 720 ? 2 : 2.6);
      if (RUN) return;
      PLAN = buildPlan(o);
      frame(cv.getContext('2d'), Math.min(1.2, PLAN.T));
    } catch (e) { info(e.message, true); }
  }

  function init() {
    if (!$('vidCanvas')) return;
    $('vidPreview').onclick = () => start(false);
    $('vidRecord').onclick = () => start(true);
    $('vidStop').onclick = stop;
    ['vidStart', 'vidHour', 'vidLen', 'vidRate', 'vidTour', 'vidRes'].forEach(id =>
      $(id).addEventListener('change', () => {
        $('vidHour').hidden = $('vidStart').value !== 'custom';
        estimate(); poster();
      }));
    document.addEventListener('visibilitychange', () => { if (RUN && document.hidden) RUN.hidden = true; });
    const c = $('vidCanvas').getContext('2d');
    c.fillStyle = '#e9edf3'; c.fillRect(0, 0, 1920, 1080);
    c.fillStyle = GREY; c.font = `500 40px ${FONT}`; c.textAlign = 'center';
    c.fillText('Run a simulation, then Preview or Record', 960, 550);
  }

  function load(j) { J = j; PLAN = null; indexSteps(j); estimate(); poster().then(estimate); }

  window.PFDVideo = {init, load, stop};
})();
