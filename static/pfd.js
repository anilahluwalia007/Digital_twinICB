/* static/pfd.js — live process flow diagram for the IDB mAb digital twin
   NavNevAI_mAbWasteWater · CarbonTesseract 2026

   Draws the whole line as a PFD (equipment symbols, stream numbers as on the
   site flowsheet) and plays the last simulation hour by hour. Flowing streams
   turn teal and animate; hover a line or stream number for L/h, g/h and g/L.
   Page ids used: pfdBox pfdPlay pfdSlider pfdSpeed pfdSvg pfdClock pfdNote pfdTable
   API: PFD.init()  PFD.load(simResult)  PFD.pause()  PFD.clock()
   For the video exporter (static/pfdvideo.js):
        PFD.svg({video: true})  clean diagram (no live text) as an SVG string
        PFD.geo()               pipe, vessel, column, value-label and tag coordinates
        PFD.state(i)            stream values, equipment state and readout at series row i
        PFD.series()  PFD.params()  PFD.length()
*/
(function () {
  'use strict';
  const W = 1440, H = 1130;
  const $ = id => document.getElementById(id);
  const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const num = v => (typeof v === 'number' && isFinite(v)) ? v : 0;
  function fmt(v) {
    if (v === null || v === undefined || !isFinite(v)) return '—';
    const a = Math.abs(v);
    if (a < 5e-4) return '0';
    return a >= 1000 ? Math.round(v).toLocaleString() : a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : v.toFixed(2);
  }

  const CSS = `
.pfd text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;fill:#16202c}
.pfd .lb{font-size:9px}.pfd .en{font-size:9px;font-weight:600}.pfd .cn{font-size:8px;font-weight:700}
.pfd .tt{font-size:7.6px;font-weight:600;fill:#3b2a66}.pfd .tt.wf{font-size:11px;font-weight:700;fill:#1e3a8a}
.pfd .tg{font-size:8px;font-weight:700;fill:#5b6675}.pfd .cap{font-size:10px;font-weight:700;fill:#0e7c7b}
.pfd .val{font-size:9px;font-weight:700;fill:#0b5e5d}.pfd .ttl{font-size:16px;font-weight:700}.pfd .sub{font-size:10px;fill:#5b6675}
.pfd .halo{paint-order:stroke;stroke:#fff;stroke-width:3px;stroke-linejoin:round}
.pfd .pl{fill:none;stroke-width:1.5;stroke-linejoin:round}
.pfd .pp .pl{stroke:#1d2733}.pfd .pb .pl{stroke:#7c3aed;stroke-dasharray:6 3}.pfd .pw .pl{stroke:#2563eb;stroke-dasharray:6 3}
.pfd .pa .pl{stroke:#dc2626;stroke-dasharray:6 3}.pfd .pd .pl{stroke:#a16207;stroke-dasharray:6 3}
.pfd .px .pl{stroke:#64748b;stroke-dasharray:2 3}.pfd .ps .pl{stroke:#16a34a;stroke-dasharray:9 4}
.pfd .st.on .pl{stroke-width:2.8}.pfd .pp.on .pl{stroke:#0e7c7b}.pfd .st.nm .pl{stroke:#94a3b8;opacity:.6}
.pfd .fl{fill:none;stroke:#fff;stroke-width:1.6;stroke-linecap:round;stroke-dasharray:1 9;opacity:0}
.pfd .st.on .fl{opacity:1;animation:pfdflow .9s linear infinite}
@keyframes pfdflow{to{stroke-dashoffset:-10}}
.pfd .vb{fill:#f1f8f8}.pfd .vo{fill:none;stroke:#1d2733;stroke-width:1.5}.pfd .lvl{fill:#0e7c7b;fill-opacity:.28}
.pfd .sym{fill:#fff;stroke:#1d2733;stroke-width:1.3}.pfd .pt{fill:#1d2733}
.pfd .mem{stroke:#1d2733;stroke-width:1;stroke-dasharray:3 2}.pfd .zz{fill:none;stroke:#1d2733;stroke-width:1}
.pfd .jk{fill:none;stroke:#94a3b8;stroke-width:2}
.pfd .ag rect{fill:#fff;stroke:#1d2733;stroke-width:1.2}.pfd .ag line,.pfd .ag path{stroke:#1d2733;stroke-width:1.4;fill:none}
.pfd .col .cb{fill:#fff;stroke:#1d2733;stroke-width:1.3}.pfd .col .bed{fill:url(#pfd-hatch)}.pfd .col .cl{stroke:#1d2733;stroke-width:1}
.pfd .col.load .cb{fill:#99e0d8}.pfd .col.regen .cb{fill:#fde68a}
.pfd .tote rect{fill:#efeafd;stroke:#6d28d9;stroke-width:1}.pfd .tote.wfi rect{fill:#e0ecff;stroke:#2563eb}
.pfd .tote .cage{fill:none;stroke:#6d28d9;stroke-opacity:.18}.pfd .tote.wfi .cage{stroke:#2563eb}
.pfd .tag path{fill:#fff;stroke:#1d2733;stroke-width:.9}.pfd .tag text{font-size:7.5px;font-weight:700;text-anchor:middle}
.pfd .tag.on path{fill:#ccefeb;stroke:#0e7c7b}.pfd .tag.nm path{stroke:#94a3b8;stroke-dasharray:2 2}
.pfd .box{fill:#f8fafc;stroke:#cbd5e1}.pfd .cryo{fill:#ecfdf5;stroke:#15803d;stroke-width:1.3}.pfd .cap2{fill:#1d2733}
@media print{.pfd .fl{display:none}}`;

  // ── drawing primitives (three layers: pipes, equipment, labels) + geometry for the video ──
  let LP, LE, LT, VIDEO = false;
  let GEO = {pipes: [], tanks: {}, cols: {}, vt: {}, tags: {}, W, H};
  const LV = {};
  const KIND = {p: 'pp', b: 'pb', w: 'pw', a: 'pa', d: 'pd', x: 'px', s: 'ps'};
  const COLOR = {pp: '#1d2733', pb: '#7c3aed', pw: '#2563eb', pa: '#dc2626', pd: '#a16207', px: '#64748b', ps: '#16a34a'};

  function txt(x, y, s, cls, anchor) {
    LT.push(`<text x="${x}" y="${y}" class="${cls || 'lb'} halo" text-anchor="${anchor || 'middle'}">${esc(s)}</text>`);
  }
  function mtxt(x, yc, arr, cls, dy) {
    dy = dy || 9;
    const y0 = yc - (arr.length - 1) * dy / 2 + 3;
    LT.push(`<text x="${x}" y="${y0}" class="${cls}" text-anchor="middle">` +
      arr.map((t, i) => `<tspan x="${x}" dy="${i ? dy : 0}">${esc(t)}</tspan>`).join('') + '</text>');
  }
  function vt(id, x, y, anchor) {
    GEO.vt[id] = {x, y, anchor: anchor || 'middle'};
    LT.push(`<text id="vt-${id}" x="${x}" y="${y}" class="val halo" text-anchor="${anchor || 'middle'}">` +
      `<tspan x="${x}"></tspan><tspan x="${x}" dy="10"></tspan></text>`);
  }
  function tag(id, x, y) {
    (GEO.tags[id] = GEO.tags[id] || []).push([x, y]);
    const w = 7 + String(id).length * 3.2;
    LT.push(`<g class="tag" data-s="${id}"><path d="M${x - w},${y} L${x},${y - 8} L${x + w},${y} L${x},${y + 8} Z"/>` +
      `<text x="${x}" y="${y + 2.8}">${esc(id)}</text><title></title></g>`);
  }
  function pipe(id, pts, kind, badge, noArrow) {
    const k = KIND[kind || 'p'];
    GEO.pipes.push({id: id || null, pts, k});
    const d = 'M' + pts.map(q => q[0] + ',' + q[1]).join(' L');
    LP.push(`<g class="st ${k}"${id ? ` data-s="${id}"` : ''}><path class="pl" d="${d}"` +
      (noArrow ? '' : ` marker-end="url(#ar-${k})"`) + '/>' +
      (id ? `<path class="fl" d="${d}"/><title></title>` : '') + '</g>');
    if (badge) tag(id, badge[0], badge[1]);
  }
  function tee(id, cx, y, half, yCol) {            // manifold to/from a column pair
    const inlet = yCol > y;
    for (const s of [-1, 1]) {
      const xs = cx + s * half;
      pipe(id, inlet ? [[cx, y], [xs, y]] : [[xs, y], [cx, y]], 'p', null, true);
      pipe(id, inlet ? [[xs, y], [xs, yCol]] : [[xs, yCol], [xs, y]], 'p', null, true);
    }
  }
  function tank(id, cx, yTop, h, name, tagNo, o) {
    o = o || {};
    const w = o.w || 70, ry = w * 0.18, x0 = cx - w / 2, x1 = cx + w / 2, yb = yTop + h;
    const d = `M${x0},${yTop} A${w / 2},${ry} 0 0 1 ${x1},${yTop} L${x1},${yb} A${w / 2},${ry} 0 0 1 ${x0},${yb} Z`;
    LV[id] = {base: yb + ry, span: h + 2 * ry};
    GEO.tanks[id] = {cx, x0, x1, yTop, yb, w, ry, base: yb + ry, span: h + 2 * ry};
    LE.push(`<g id="eq-${id}"><clipPath id="cp-${id}"><path d="${d}"/></clipPath><path class="vb" d="${d}"/>` +
      `<rect class="lvl" id="lv-${id}" clip-path="url(#cp-${id})" x="${x0}" y="${yb + ry}" width="${w}" height="0"/>` +
      `<path class="vo" d="${d}"/></g>`);
    mtxt(cx, yTop + h * 0.3, name, 'en halo');
    vt(id, cx, yTop + h * 0.72);
    if (o.tagLeft) txt(x0 - 5, yTop + 2, tagNo, 'tg', 'end'); else txt(x1 + 5, yTop + 2, tagNo, 'tg', 'start');
  }
  function agitator(cx, yTop, h, w) {
    const ym = yTop - w * 0.18 - 12;
    LE.push(`<g class="ag"><rect x="${cx - 8}" y="${ym - 10}" width="16" height="10" rx="2"/>` +
      `<line x1="${cx}" y1="${ym}" x2="${cx}" y2="${yTop + h - 6}"/>` +
      `<path d="M${cx - 15},${yTop + h * 0.5} h30 M${cx - 15},${yTop + h - 6} h30"/></g>`);
  }
  function jacket(cx, yTop, h, w) {
    LE.push(`<path class="jk" d="M${cx - w / 2 - 4},${yTop + h * 0.2} V${yTop + h + 4} M${cx + w / 2 + 4},${yTop + h * 0.2} V${yTop + h + 4}"/>`);
  }
  function column(id, cx, yTop, h, label) {
    const w = 20, x0 = cx - w / 2;
    GEO.cols[id] = {x0, yTop, w, h, label: label || ''};
    LE.push(`<g class="col" id="eq-${id}"><rect class="cb" x="${x0}" y="${yTop}" width="${w}" height="${h}" rx="4"/>` +
      `<rect class="bed" x="${x0 + 1.5}" y="${yTop + 8}" width="${w - 3}" height="${h - 16}"/>` +
      `<path class="cl" d="M${x0},${yTop + 8} h${w} M${x0},${yTop + h - 8} h${w}"/></g>`);
    if (label) txt(cx, yTop + h / 2 + 3, label, 'cn');
  }
  function filt(cx, cy, tagNo, big) {
    const w = big ? 26 : 20, h = big ? 36 : 30;
    LE.push(`<g><rect class="sym" x="${cx - w / 2}" y="${cy - h / 2}" width="${w}" height="${h}" rx="${w / 2}"/>` +
      `<line class="mem" x1="${cx}" y1="${cy - h / 2 + 5}" x2="${cx}" y2="${cy + h / 2 - 5}"/></g>`);
    if (tagNo) txt(cx, cy - h / 2 - 4, tagNo, 'tg');
  }
  function pump(cx, cy, tagNo) {
    LE.push(`<g><circle class="sym" cx="${cx}" cy="${cy}" r="10"/>` +
      `<path class="pt" d="M${cx - 4},${cy - 6} L${cx + 7},${cy} L${cx - 4},${cy + 6} Z"/></g>`);
    txt(cx, cy + 21, tagNo, 'tg');
  }
  function valve(cx, cy, vertical, branch) {
    let d = vertical ? `M${cx - 5},${cy - 7} L${cx + 5},${cy + 7} L${cx - 5},${cy + 7} L${cx + 5},${cy - 7} Z`
                     : `M${cx - 7},${cy - 5} L${cx + 7},${cy + 5} L${cx + 7},${cy - 5} L${cx - 7},${cy + 5} Z`;
    if (branch === 'up') d += ` M${cx},${cy} L${cx - 5},${cy - 9} L${cx + 5},${cy - 9} Z`;
    if (branch === 'down') d += ` M${cx},${cy} L${cx - 5},${cy + 9} L${cx + 5},${cy + 9} Z`;
    LE.push(`<path class="sym" d="${d}"/>`);
  }
  function smix(cx, cy) {
    LE.push(`<g><rect class="sym" x="${cx - 20}" y="${cy - 8}" width="40" height="16" rx="2"/>` +
      `<path class="zz" d="M${cx - 16},${cy} l4,-5 l4,10 l4,-10 l4,10 l4,-10 l4,10 l4,-10 l4,5"/></g>`);
  }
  function tff(cx, cy, w) {                         // membrane module, membrane 6 px above the centre line
    const x0 = cx - w / 2;
    LE.push(`<g><rect class="sym" x="${x0}" y="${cy - 18}" width="${w}" height="36" rx="3"/>` +
      `<line class="mem" x1="${x0 + 3}" y1="${cy - 6}" x2="${x0 + w - 3}" y2="${cy - 6}"/></g>`);
  }
  function holder(cx, cy) {
    let s = `<g><rect class="sym" x="${cx - 26}" y="${cy - 28}" width="52" height="56" rx="2"/>`;
    for (let i = -18; i <= 18; i += 6) s += `<line class="mem" x1="${cx + i}" y1="${cy - 24}" x2="${cx + i}" y2="${cy + 24}"/>`;
    LE.push(s + '</g>');
  }
  function tote(cx, yTop, lines, wfi, capBottom) {
    const w = 72, h = 58, x0 = cx - w / 2, yc = capBottom ? yTop + h - 1 : yTop - 4;
    LE.push(`<g class="tote${wfi ? ' wfi' : ''}"><rect x="${cx - 7}" y="${yc}" width="14" height="5" rx="1"/>` +
      `<rect x="${x0}" y="${yTop}" width="${w}" height="${h}" rx="4"/>` +
      `<path class="cage" d="M${x0},${yTop + h / 3} h${w} M${x0},${yTop + 2 * h / 3} h${w} M${x0 + w / 3},${yTop} v${h} M${x0 + 2 * w / 3},${yTop} v${h}"/></g>`);
    mtxt(cx, yTop + h / 2, lines, wfi ? 'tt wf' : 'tt', 8.6);
  }
  function totes(list, yTop) { list.forEach(([x, l, w]) => tote(x, yTop, l, !!w)); }
  function header(x0, x1, y, drops, yDrop, kind) {
    pipe(null, [[x0, y], [x1, y]], kind || 'b', null, true);
    drops.forEach(x => pipe(null, [[x, yDrop], [x, y]], kind || 'b', null, true));
  }
  function cryo(cx, cy) {
    let s = `<g><rect class="cryo" x="${cx - 48}" y="${cy - 42}" width="96" height="84" rx="6"/>`;
    for (const dx of [-26, 0, 26]) {
      const bx = cx + dx, by = cy + 6;
      s += `<path class="sym" d="M${bx - 8},${by + 18} V${by - 4} Q${bx - 8},${by - 11} ${bx - 3},${by - 12} V${by - 17} H${bx + 3} V${by - 12} Q${bx + 8},${by - 11} ${bx + 8},${by - 4} V${by + 18} Z"/>` +
           `<rect class="cap2" x="${bx - 4}" y="${by - 21}" width="8" height="4" rx="1"/>`;
    }
    LE.push(s + '</g>');
  }
  const WFI = ['WFI'];

  // ── the flowsheet ──
  function drawRow1() {
    LT.push(`<text x="${W - 10}" y="24" class="ttl" text-anchor="end">IDB mAb line — process flow diagram</text>`);
    if (!VIDEO) LT.push(`<text id="pfd-clk" x="${W - 10}" y="40" class="sub" text-anchor="end">Run a simulation to start the live view</text>`);
    // USP: N-1 seed, N reactor, TFF
    txt(90, 38, 'Media');
    pipe('U1', [[90, 45], [90, 111]], 'p', [90, 74]);
    tank('R101', 110, 120, 100, ['N-1', 'perfusion', 'seed'], 'R-101', {w: 64, tagLeft: true});
    agitator(110, 120, 100, 64); jacket(110, 120, 100, 64);
    pipe('U2', [[142, 140], [168, 140], [168, 62]], 'x', [168, 100]); txt(168, 55, 'Waste');
    pipe('S1', [[110, 233], [110, 262], [250, 262], [250, 244]], 's', [180, 262]);
    txt(180, 280, 'Seed transfer at N inoculation');
    txt(245, 38, 'Media');
    pipe('U3', [[245, 45], [245, 100]], 'p', [245, 74]);
    tank('R102', 270, 110, 120, ['N reactor', '(IDB)'], 'R-102', {w: 80, tagLeft: true});
    agitator(270, 110, 120, 80); jacket(270, 110, 120, 80);
    pipe('U4', [[310, 215], [380, 215], [380, 203], [399, 203]], 'p', [362, 215]);
    pump(336, 215, 'P-101');
    tff(435, 215, 70); txt(415, 191, 'F-101', 'tg'); txt(435, 246, 'TFF cell removal'); vt('F101', 435, 259);
    pipe('U5', [[445, 197], [445, 152]], 'x', [445, 174]); txt(445, 144, 'Cell concentrate → bio-kill');
    pipe('1', [[470, 215], [513, 215]], 'p', [492, 215]);
    valve(515, 215, false, 'up'); txt(515, 236, 'XV-101', 'tg');
    pipe('2', [[515, 215], [515, 120]], 'x', [515, 165]); txt(515, 112, 'Waste (before capture)');
    pipe('3', [[515, 215], [545, 215], [545, 75], [585, 75], [585, 84]], 'p', [545, 150]);
    tank('T101', 600, 95, 90, ['Conditioned', 'media surge', 'tank'], 'T-101');
    pipe('4', [[600, 198], [600, 215], [698, 215]], 'p', [676, 215]);
    pump(646, 215, 'P-102');
    // detergent, buffers, WFI to the load line
    tote(700, 6, ['10 %wt Brij-9', 'concentrate'], false, true);
    pipe('8', [[700, 68], [700, 205]], 'd', [700, 165]);
    valve(700, 118, true); txt(690, 104, 'Detergent', 'lb', 'end');
    pipe('6', [[70, 348], [70, 305], [698, 305]], 'w', [600, 305]);
    header(150, 630, 330, [150, 230, 310, 390, 470, 550, 630], 348);
    pipe('7', [[630, 330], [700, 330], [700, 307]], 'b', [665, 330]);
    pipe('5', [[700, 305], [700, 225]], 'b', [700, 265]);
    totes([[70, WFI, 1], [150, ['50 mM Tris', 'pH 7.4']],
      [230, ['25 mM Tris', '0.42 M Arg', '0.5 M NaCl', '0.1% Virodex', 'pH 9.0']],
      [310, ['25 mM NaOAc', 'pH 6.5']], [390, ['25 mM NaOAc', 'pH 4.45']], [470, ['100 mM', 'acetic acid']],
      [550, ['0.5 N NaOH', '(ILD w/ WFI)']], [630, ['2% benzyl alc.', '100 mM NaOAc', 'pH 5.0']]], 352);
    txt(350, 428, 'Protein A buffers and WFI', 'tg');
    // Protein A capture: columns 1 & 2 and 3 & 4
    pipe('9', [[700, 215], [735, 215], [735, 92], [760, 92]], 'p', [735, 165]);
    txt(860, 62, 'Capture chromatography — Protein A', 'cap');
    pipe('10', [[760, 92], [800, 92], [800, 106]], 'p', [780, 92], true); tee('10', 800, 106, 15, 112);
    pipe('11', [[760, 92], [760, 80], [880, 80], [880, 106]], 'p', [830, 80], true); tee('11', 880, 106, 15, 112);
    column('C101', 785, 112, 82, '1'); column('C102', 815, 112, 82, '2');
    column('C103', 865, 112, 82, '3'); column('C104', 895, 112, 82, '4');
    tee('12', 800, 203, 15, 194); pipe('12', [[800, 203], [800, 235], [878, 235]], 'p', [838, 235]);
    tee('13', 880, 203, 15, 194); pipe('13', [[880, 203], [880, 233]], 'p', [880, 219]);
    txt(800, 252, 'Columns 1 & 2'); txt(880, 252, 'Columns 3 & 4'); vt('CAP', 840, 268);
    valve(912, 235, false, 'down'); txt(912, 224, 'XV-102', 'tg');
    pipe('W1', [[912, 235], [912, 285]], 'x', [912, 262]); txt(912, 298, 'Column waste → drain');
    filt(945, 235, 'FL-101');
    pipe('14', [[880, 235], [985, 235], [985, 90], [1015, 90], [1015, 99]], 'p', [970, 235]);
    tank('T102', 1030, 110, 95, ['Protein A', 'vessel SUM'], 'T-102');
    pipe('15', [[1030, 218], [1030, 235], [1098, 235]], 'p', [1084, 235]);
    pump(1056, 235, 'P-103');
    pipe('16', [[1100, 305], [1100, 237]], 'b', [1100, 272]);
    pipe('17', [[1100, 332], [1100, 307]], 'b', [1100, 319]);
    header(980, 1300, 332, [980, 1060, 1140, 1220, 1300], 348);
    pipe('18', [[1385, 348], [1385, 305], [1102, 305]], 'w', [1250, 305]);
    totes([[980, ['50 mM Tris', 'pH 7.4']], [1060, ['2 M NaCl', '(0.5 M w/ WFI)']],
      [1140, ['1 M arginine', '50 mM NaPi', 'pH 7']], [1220, ['1 N NaOH']],
      [1300, ['50 mM NaOAc', 'pH 5.0']], [1385, WFI, 1]], 352);
    txt(1140, 428, 'Post-capture buffers and WFI', 'tg');
    // depth filtration
    pipe('19', [[1100, 235], [1159, 235]], 'p', [1130, 235]);
    holder(1185, 235); txt(1185, 198, 'Depth filtration', 'cap'); txt(1185, 277, 'FH-101', 'tg'); vt('FH', 1185, 287);
    pipe('20', [[1211, 235], [1250, 235], [1250, 90], [1285, 90], [1285, 99]], 'p', [1230, 235]);
    tank('T103', 1300, 110, 95, ['Depth filt.', 'surge tank'], 'T-103');
    pipe('21', [[1300, 218], [1300, 235], [1432, 235], [1432, 455], [20, 455], [20, 600], [79, 600]], 'p', [1392, 235]);
    tag('21', 20, 530);
    pump(1342, 235, 'P-104');
  }

  function drawRow2() {
    // on-line base addition + AEX
    tote(60, 478, ['2 M Tris', 'base'], false, true); tote(140, 478, ['2 M NaCl'], false, true);
    header(60, 140, 556, [60, 140], 540, 'a');
    pipe('25', [[100, 556], [100, 591]], 'a', [100, 574]);
    smix(100, 600); txt(100, 621, 'SM-201', 'tg'); txt(100, 634, 'On-line base addition');
    pipe('26', [[121, 600], [205, 600], [205, 478], [229, 478]], 'p', [205, 540]);
    pipe('22', [[165, 700], [165, 603]], 'b', [165, 652]);
    pipe('23', [[165, 740], [165, 702]], 'b', [165, 721]);
    pipe('24', [[390, 756], [390, 700], [167, 700]], 'w', [300, 700]);
    header(70, 310, 740, [70, 150, 230, 310], 756);
    totes([[70, ['50 mM Tris', '0.5 M NaCl', 'pH 7.4']], [150, ['50 mM Tris', 'pH 7.4']], [230, ['1.0 N NaOH']],
      [310, ['2% benzyl alc.', '100 mM NaOAc', 'pH 5.0']], [390, WFI, 1]], 760);
    txt(230, 836, 'AEX buffers and WFI', 'tg');
    pipe('27', [[229, 478], [250, 478], [250, 505]], 'p', [250, 492], true);
    pipe('28', [[229, 478], [229, 466], [290, 466], [290, 505]], 'p', [262, 466], true);
    column('C201', 250, 505, 74, '1'); column('C202', 290, 505, 74, '2');
    pipe('29', [[250, 579], [250, 600], [288, 600]], 'p', [269, 600]);
    pipe('30', [[290, 579], [290, 598]], 'p', [290, 589]);
    txt(262, 622, 'AEX chromatography', 'cap'); vt('AEX', 262, 636);
    filt(325, 600, 'FL-201');
    pipe('31', [[290, 600], [360, 600], [360, 472], [395, 472], [395, 484]], 'p', [360, 540]);
    tank('T201', 410, 495, 80, ['AEX', 'vessel SUM'], 'T-201');
    pipe('32', [[410, 588], [410, 600], [513, 600]], 'p', [490, 600]);
    pump(460, 600, 'P-201');
    // on-line acid addition + CEX
    tote(535, 478, ['1 M', 'acetic acid'], false, true);
    pipe('36', [[535, 540], [535, 591]], 'a', [535, 566]);
    smix(535, 600); txt(535, 621, 'SM-202', 'tg'); txt(535, 634, 'On-line acid addition');
    pipe('33', [[590, 700], [590, 603]], 'b', [590, 655]);
    pipe('34', [[590, 740], [590, 702]], 'b', [590, 721]);
    pipe('35', [[470, 756], [470, 700], [588, 700]], 'w', [530, 700]);
    header(550, 870, 740, [550, 630, 710, 790, 870], 756);
    totes([[470, WFI, 1], [550, ['50 mM NaOAc', '120 mM NaCl', 'pH 5.5']], [630, ['50 mM NaOAc', '280 mM NaCl', 'pH 6.0']],
      [710, ['50 mM Tris', '0.5 M NaCl', 'pH 7.4']], [790, ['1.0 N NaOH']],
      [870, ['2% benzyl alc.', '100 mM NaOAc', 'pH 5.0']]], 760);
    txt(670, 836, 'CEX buffers and WFI', 'tg');
    pipe('37', [[555, 600], [620, 600], [620, 478], [649, 478]], 'p', [620, 540]);
    pipe('38', [[649, 478], [670, 478], [670, 505]], 'p', [670, 492], true);
    pipe('39', [[649, 478], [649, 466], [710, 466], [710, 505]], 'p', [682, 466], true);
    column('C203', 670, 505, 74, '1'); column('C204', 710, 505, 74, '2');
    pipe('40', [[670, 579], [670, 600], [708, 600]], 'p', [689, 600]);
    pipe('41', [[710, 579], [710, 598]], 'p', [710, 589]);
    txt(680, 622, 'CEX chromatography', 'cap'); vt('CEX', 676, 636);
    valve(748, 600, false, 'down'); txt(748, 589, 'XV-201', 'tg');
    pipe('W2', [[748, 600], [748, 652]], 'x', [748, 628]); txt(748, 666, 'Flow-through → drain');
    filt(780, 600, 'FL-202');
    pipe('42', [[710, 600], [810, 600], [810, 472], [840, 472], [840, 484]], 'p', [810, 540]);
    tank('T202', 855, 495, 80, ['CEX', 'vessel SUM'], 'T-202');
    pipe('43', [[855, 588], [855, 600], [948, 600]], 'p', [926, 600]);
    pump(898, 600, 'P-202');
    pipe('44', [[950, 730], [950, 603]], 'b', [950, 662]);
    pipe('45', [[1030, 756], [1030, 700], [952, 700]], 'b', [992, 700]);
    pipe('46', [[950, 756], [950, 732]], 'w', [950, 744]);
    totes([[950, WFI, 1], [1030, ['50 mM NaOAc', '280 mM NaCl', 'pH 6.0']]], 760);
    txt(990, 836, 'VF flush buffer and WFI', 'tg');
    pipe('47', [[950, 600], [1432, 600], [1432, 850]], 'p', [1060, 600]);
    pipe('48', [[1432, 850], [20, 850], [20, 990], [54, 990]], 'p', [20, 920]);
    // live readout box (the video draws its own side panel instead)
    if (!VIDEO) {
      LE.push('<rect class="box" x="1090" y="470" width="330" height="116" rx="6"/>');
      txt(1102, 488, 'Live readout', 'cap', 'start');
      for (let i = 0; i < 5; i++) LT.push(`<text id="ro${i}" x="1102" y="${506 + i * 15}" class="lb" text-anchor="start"></text>`);
    }
  }

  function drawRow3() {
    // virus retentive filtration
    txt(110, 918, 'Virus retentive filtration', 'cap'); txt(110, 932, 'F-301 A / B', 'tg');
    pipe('49', [[54, 990], [54, 955], [110, 955]], 'p', [78, 955], true);
    pipe('50', [[54, 990], [54, 1025], [110, 1025]], 'p', [78, 1025], true);
    filt(110, 955, '', true); filt(110, 1025, '', true);
    pipe('51', [[110, 955], [165, 955], [165, 987]], 'p', [140, 955]);
    pipe('52', [[110, 1025], [165, 1025], [165, 993]], 'p', [140, 1025]);
    vt('VRF', 110, 1058);
    filt(205, 990, 'FL-301');
    pipe('53', [[165, 990], [245, 990], [245, 868], [285, 868], [285, 880]], 'p', [245, 930]);
    tank('T301', 300, 890, 75, ['VRF', 'vessel SUM'], 'T-301');
    pipe('54', [[300, 978], [300, 990], [428, 990]], 'p', [392, 990]);
    pump(345, 990, 'P-301');
    // UF/DF
    pipe('55', [[430, 1030], [430, 993]], 'b', [430, 1011]);
    pipe('57', [[430, 1056], [430, 1032]], 'w', [430, 1044]);
    pipe('56', [[550, 1045], [550, 1030], [432, 1030]], 'b', [492, 1030]);
    header(510, 590, 1045, [510, 590], 1056);
    totes([[430, WFI, 1], [510, ['15 mM acetate', '130 mM proline', 'pH 5.1']],
      [590, ['15 mM acetate', '113 mM proline', 'pH 5.3']]], 1060);
    pipe('58', [[430, 990], [519, 990]], 'p', [475, 990]);
    tff(560, 990, 80); txt(606, 968, 'UF-301', 'tg', 'start');
    txt(560, 1022, 'Ultrafiltration / diafiltration', 'cap'); vt('UF', 606, 1040, 'start');
    pipe('W3', [[560, 972], [560, 925]], 'x', [560, 950]); txt(560, 916, 'UF/DF permeate → drain');
    pipe('59', [[600, 990], [650, 990], [650, 868], [685, 868], [685, 880]], 'p', [650, 930]);
    tank('T302', 700, 890, 75, ['UF/DF', 'vessel SUM'], 'T-302');
    pipe('60', [[700, 978], [700, 990], [798, 990]], 'p', [775, 990]);
    pump(742, 990, 'P-302');
    // formulation, BDS hold, cryovault
    pipe('61', [[800, 1056], [800, 993]], 'b', [800, 1025]);
    tote(800, 1060, ['25 mM acetate', '113 mM proline', '1.17 M sucrose', '6 mg/mL P188', 'pH 5.3']);
    pipe('62', [[800, 990], [860, 990], [860, 880], [890, 880], [890, 891]], 'p', [860, 935]);
    tank('M301', 910, 900, 65, ['BDS hold', 'mixer'], 'M-301'); agitator(910, 900, 65, 70);
    pipe('63', [[910, 978], [910, 990], [1016, 990]], 'p', [988, 990]);
    pump(952, 990, 'P-303');
    cryo(1065, 990); txt(1065, 940, 'Cryovault containers', 'cap'); vt('CV', 1065, 1046);
    // legend
    LE.push('<rect class="box" x="1135" y="868" width="295" height="258" rx="6"/>');
    txt(1147, 886, 'Legend', 'cap', 'start');
    [['pp', 'Process stream'], ['pp on', 'Process stream, flowing this hour'], ['pb', 'Buffer'], ['pw', 'WFI'],
     ['pa', 'On-line acid / base addition'], ['pd', 'Detergent'], ['px', 'Waste, drain or bio-kill'],
     ['ps', 'Cell seed transfer']].forEach(([k, t], i) => {
      const y = 902 + i * 15, d = `M1147,${y} H1185`;
      LT.push(`<g class="st ${k}"><path class="pl" d="${d}" marker-end="url(#ar-${k.split(' ')[0]})"/>` +
        (k.indexOf('on') > 0 ? `<path class="fl" d="${d}"/>` : '') + '</g>');
      txt(1192, y + 3, t, 'lb', 'start');
    });
    let y = 1026;
    LT.push(`<g class="tag"><path d="M1156,${y} L1166,${y - 8} L1176,${y} L1166,${y + 8} Z"/><text x="1166" y="${y + 2.8}">n</text></g>`);
    txt(1192, y + 3, 'Stream number: hover for live values', 'lb', 'start');
    y += 18;
    LT.push(`<rect x="1150" y="${y - 7}" width="9" height="14" rx="2" fill="#99e0d8" stroke="#1d2733"/>` +
            `<rect x="1166" y="${y - 7}" width="9" height="14" rx="2" fill="#fde68a" stroke="#1d2733"/>`);
    txt(1192, y + 3, 'Column loading / eluting', 'lb', 'start');
    y += 18;
    LT.push(`<rect x="1152" y="${y - 7}" width="26" height="14" fill="#0e7c7b" fill-opacity=".28" stroke="#1d2733"/>`);
    txt(1192, y + 3, 'Vessel shading = liquid level', 'lb', 'start');
    y += 18; txt(1147, y + 3, 'Capsule = 0.2 µm filter · circle = pump', 'lb', 'start');
    y += 14; txt(1147, y + 3, 'Bow-tie = valve · zig-zag box = static mixer', 'lb', 'start');
    y += 14; txt(1147, y + 3, 'Grey dashed = on the flowsheet, not simulated', 'lb', 'start');
    y += 14; txt(1147, y + 3, 'Values: hourly, from the last simulation run', 'lb', 'start');
  }

  function defs() {
    let s = '<defs>';
    for (const [k, c] of Object.entries(COLOR))
      s += `<marker id="ar-${k}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" ` +
           `markerUnits="userSpaceOnUse" orient="auto"><path d="M0,0 L10,5 L0,10 Z" fill="${c}"/></marker>`;
    s += '<pattern id="pfd-hatch" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">' +
         '<line x1="0" y1="0" x2="0" y2="5" stroke="#64748b" stroke-width="0.8"/></pattern>';
    return s + '</defs>';
  }

  function build(opts) {
    VIDEO = !!(opts && opts.video);
    LP = []; LE = []; LT = [];
    GEO = {pipes: [], tanks: {}, cols: {}, vt: {}, tags: {}, W, H};
    drawRow1(); drawRow2(); drawRow3();
    return `<svg xmlns="http://www.w3.org/2000/svg" class="pfd" viewBox="0 0 ${W} ${H}" role="img" ` +
      `aria-label="Process flow diagram of the IDB mAb line"><style>${CSS}</style>${defs()}` +
      `<rect x="0" y="0" width="${W}" height="${H}" fill="#fff"/>` + LP.join('') + LE.join('') + LT.join('') + '</svg>';
  }

  // ── stream values: [id, description, fn(row, ctx) -> {L: L/h, g: g/h} | null = not modelled] ──
  const Z = {L: 0, g: 0};
  const duty = (c, which, v) => (c.duty === which ? v : Z);
  const STREAMS = [
    ['U1', 'Media feed to the N-1 perfusion seed', r => ({L: r.phase ? 0 : r.ww_n1_permeate})],
    ['U2', 'N-1 cell-free permeate to waste', r => ({L: r.ww_n1_permeate})],
    ['S1', 'N-1 cells to the N reactor at inoculation', (r, c) => ({L: c.inoc ? num(c.P.n1.transfer_L) : 0})],
    ['U3', 'Media feed to the N reactor', r => ({L: r.Fin_Lh})],
    ['U4', 'N outflow (cells + broth + product) to TFF', r => ({L: r.Fin_Lh, g: r.harvest_g_h + r.Fb_Lh * r.titre_harvest})],
    ['U5', 'TFF cell concentrate to bio-kill', r => ({L: r.Fb_Lh, g: r.Fb_Lh * r.titre_harvest})],
    ['1', 'TFF permeate (perfusion harvest from USP)', r => ({L: r.Fh_Lh, g: r.harvest_g_h})],
    ['2', 'Permeate to waste before capture starts', r => ({L: r.ww_ramp_harvest, g: r.ww_ramp_harvest * r.titre_harvest})],
    ['3', 'Permeate to the conditioned media surge tank', r => { const L = Math.max(0, r.Fh_Lh - r.ww_ramp_harvest); return {L, g: L * r.titre_harvest}; }],
    ['4', 'Surge tank to Protein A load (g ≈ sump concentration)', (r, c) => ({L: r.Fload_Lh, g: r.Fload_Lh * c.cs})],
    ['8', 'Detergent (10 %wt Brij-9) — not modelled', null],
    ['6', 'WFI — included in buffer volumes in the model', null],
    ['7', 'Protein A buffers (equil., wash, elution, strip, CIP, re-equil.)', r => ({L: r.buf_capture})],
    ['5', 'Protein A buffers + WFI', r => ({L: r.buf_capture})],
    ['9', 'Column feed (load + buffers, hourly totals)', (r, c) => c.feed],
    ['10', 'Feed to Protein A columns 1 & 2', (r, c) => (c.loadA ? c.feed : Z)],
    ['11', 'Feed to Protein A columns 3 & 4', (r, c) => (c.loadA ? Z : c.feed)],
    ['12', 'Protein A eluate, columns 1 & 2', (r, c) => ({L: c.elutA ? r.capture_switch * c.elu : 0})],
    ['13', 'Protein A eluate, columns 3 & 4', (r, c) => ({L: c.elutA ? 0 : r.capture_switch * c.elu})],
    ['14', 'Protein A eluate to Protein A vessel SUM', (r, c) => ({L: r.capture_switch * c.elu})],
    ['W1', 'Protein A flow-through, buffers, strip and CIP to drain',
      r => ({L: r.ww_capture_flowthrough + r.ww_capture_buffers + r.ww_capture_strip_cip})],
    ['15', 'Low-pH VI pool leaving the Protein A vessel SUM', (r, c) => ({L: r.F_vi_Lh, g: r.g_vi / c.viy})],
    ['18', 'WFI — not modelled', null],
    ['17', 'VI acid and base titrants', r => ({L: r.buf_vi})],
    ['16', 'VI titrants + WFI', r => ({L: r.buf_vi})],
    ['19', 'Neutralised VI pool to depth filtration', (r, c) => c.vi],
    ['20', 'Depth filtrate (depth filtration not modelled: no loss)', (r, c) => c.vi],
    ['21', 'Depth filtration surge tank to AEX', (r, c) => c.vi],
    ['25', 'On-line base addition (2 M Tris base / 2 M NaCl) — not modelled', null],
    ['24', 'WFI — not modelled', null],
    ['23', 'AEX buffers', r => ({L: r.buf_aex})],
    ['22', 'AEX buffers + WFI', r => ({L: r.buf_aex})],
    ['26', 'AEX load', (r, c) => c.aexIn],
    ['27', 'AEX column 1 load', (r, c) => duty(c, 0, c.aexIn)],
    ['28', 'AEX column 2 load', (r, c) => duty(c, 1, c.aexIn)],
    ['29', 'AEX column 1 flow-through pool', (r, c) => duty(c, 0, c.aex)],
    ['30', 'AEX column 2 flow-through pool', (r, c) => duty(c, 1, c.aex)],
    ['31', 'AEX pool to AEX vessel SUM', (r, c) => c.aex],
    ['32', 'AEX vessel SUM to on-line acid addition', (r, c) => c.aex],
    ['36', 'On-line acid addition (1 M acetic acid) — not modelled', null],
    ['35', 'WFI — not modelled', null],
    ['34', 'CEX buffers', r => ({L: r.buf_cex})],
    ['33', 'CEX buffers + WFI', r => ({L: r.buf_cex})],
    ['37', 'CEX load', (r, c) => c.cexIn],
    ['38', 'CEX column 1 load', (r, c) => duty(c, 0, c.cexIn)],
    ['39', 'CEX column 2 load', (r, c) => duty(c, 1, c.cexIn)],
    ['W2', 'CEX load flow-through to drain', r => ({L: r.ww_cex_flowthrough})],
    ['40', 'CEX column 1 eluate', (r, c) => duty(c, 0, c.cex)],
    ['41', 'CEX column 2 eluate', (r, c) => duty(c, 1, c.cex)],
    ['42', 'CEX eluate pool to CEX vessel SUM', (r, c) => c.cex],
    ['43', 'CEX vessel SUM out', (r, c) => c.cex],
    ['46', 'WFI — not modelled', null],
    ['45', 'Virus filter pre-use flush buffer', r => ({L: r.buf_vf})],
    ['44', 'Flush buffer + WFI', r => ({L: r.buf_vf})],
    ['47', 'Feed to virus retentive filtration', (r, c) => c.vfIn],
    ['48', 'Feed to virus retentive filtration', (r, c) => c.vfIn],
    ['49', 'Virus filter A feed', (r, c) => duty(c, 0, c.vfIn)],
    ['50', 'Virus filter B feed', (r, c) => duty(c, 1, c.vfIn)],
    ['51', 'Virus filter A filtrate', (r, c) => duty(c, 0, c.vf)],
    ['52', 'Virus filter B filtrate', (r, c) => duty(c, 1, c.vf)],
    ['53', 'Virus filtrate to VRF vessel SUM', (r, c) => c.vf],
    ['54', 'VRF vessel SUM to UF/DF', (r, c) => c.vf],
    ['57', 'WFI — not modelled', null],
    ['56', 'UF/DF diafiltration buffer', r => ({L: r.L_df})],
    ['55', 'DF buffer + WFI', r => ({L: r.L_df})],
    ['58', 'UF/DF feed (VRF pool + DF buffer)', r => ({L: r.L_vf + r.L_df, g: r.g_vf})],
    ['W3', 'UF/DF permeate to drain', r => ({L: r.ww_sptff_permeate})],
    ['59', 'UF/DF retentate at target concentration', r => ({L: r.L_ret, g: r.ds_g_h})],
    ['60', 'UF/DF vessel SUM to formulation', r => ({L: r.L_ret, g: r.ds_g_h})],
    ['61', 'Formulation buffer (excipients)', r => ({L: r.L_exc})],
    ['62', 'Formulated drug substance to BDS hold mixer', r => ({L: r.L_ds, g: r.ds_g_h})],
    ['63', 'Bulk drug substance to cryovault containers', r => ({L: r.L_ds, g: r.ds_g_h})],
  ];

  // ── playback state ──
  let S = [], P = null, ELU = 0, idx = 0, pos = 0, timer = null, clockTxt = '—';
  const ELS = {}, ROWS = {};

  function ctx(r0, r) {
    const aexFirst = !(P.polish && P.polish.order === 'cex_aex');
    const vi = {L: r.L_vi, g: r.g_vi}, aex = {L: r.L_aex, g: r.g_aex}, cex = {L: r.L_cex, g: r.g_cex};
    const cs = r.sump_L > 1 ? r.sump_g / r.sump_L : r.titre_harvest;
    const cyc = Math.round(r.cycles);
    return {
      P, aexFirst, vi, aex, cex, cs, elu: ELU,
      aexIn: aexFirst ? vi : cex, cexIn: aexFirst ? aex : vi, vfIn: aexFirst ? cex : aex,
      vf: {L: r.L_vf, g: r.g_vf},
      feed: {L: r.Fload_Lh + r.buf_capture, g: r.Fload_Lh * cs},
      duty: Math.floor(Math.max(r.n_day, 0) * 24 / Math.max(num(P.lot_window_h) || 24, 1)) % 2,
      loadA: cyc % 2 === 0, elutA: cyc % 2 === 1,
      inoc: !!(P.n1 && P.n1.enabled) && r0.phase === 1 && Math.abs(r.n_day) < 1e-6,
      viy: Math.max(num(P.vi.vi_yield), 1e-9),
    };
  }

  // DOM writers for the live diagram
  function domLvl(id, f) {
    const g = LV[id], el = $('lv-' + id);
    if (!g || !el) return;
    f = Math.max(0, Math.min(1, num(f)));
    el.setAttribute('y', g.base - f * g.span);
    el.setAttribute('height', f * g.span);
  }
  function domVT(id, a, b) {
    const t = $('vt-' + id);
    if (!t) return;
    t.children[0].textContent = a || '';
    t.children[1].textContent = b || '';
  }
  function domCol(id, cls) {
    const e = $('eq-' + id);
    if (!e) return;
    e.classList.remove('load', 'regen');
    if (cls) e.classList.add(cls);
  }

  // equipment state for one hour, written into E = {lv, vt, col}
  function equipment(r0, r, c, v, E) {
    const lvl = (id, f) => { E.lv[id] = Math.max(0, Math.min(1, num(f))); };
    const setVT = (id, a, b) => { E.vt[id] = [a || '', b || '']; };
    const colCls = (id, cls) => { E.col[id] = cls || ''; };
    const inN = r0.phase === 1, profile = P.upstream.mode !== 'mechanistic';
    lvl('R101', inN ? 0 : 0.8);
    setVT('R101', inN ? 'transferred' : `${fmt(r.X_n1)} ×10⁶/mL`, '');
    lvl('R102', inN ? 0.8 : 0);
    if (!inN) setVT('R102', 'awaiting seed', '');
    else if (profile) setVT('R102', `${fmt(r.titre_harvest)} g/L`, `${fmt(r.VVD)} VVD`);
    else setVT('R102', `VCD ${fmt(r.Xv)}`, `${fmt(100 * r.viability)} % viable`);
    setVT('F101', r.Fh_Lh > 0 ? `${fmt(r.Fh_Lh)} L/h permeate` : '', '');
    const sumMax = P.capture.mode === 'pcc' ? P.surge.V_max_L : P.capture.sump_max_L;
    lvl('T101', r.sump_L / Math.max(num(sumMax), 1e-9));
    setVT('T101', `${fmt(r.sump_L)} L`, `${fmt(r.sump_g)} g`);
    const st = Math.round(r.col_state);
    const clsA = st === 1 && c.loadA ? 'load' : st === 2 && c.elutA ? 'regen' : '';
    const clsB = st === 1 && !c.loadA ? 'load' : st === 2 && !c.elutA ? 'regen' : '';
    colCls('C101', clsA); colCls('C102', clsA); colCls('C103', clsB); colCls('C104', clsB);
    setVT('CAP', `cycle ${Math.round(r.cycles)} · ` +
      (st === 1 ? 'loading' : st === 2 ? 'eluting / regenerating' : 'waiting for a full load'), '');
    lvl('T102', r.eluate_surge_L / Math.max(num(P.vi.eluate_surge_max_L), 1e-9));
    setVT('T102', `${fmt(r.eluate_surge_L)} L`, r0.vi_pH !== null && r0.vi_pH !== undefined ? `pH ${num(r0.vi_pH).toFixed(2)}` : '');
    setVT('FH', 'not modelled', '(no loss)');
    const pass = (id, s) => {
      const on = s && s.L > 1e-6;
      lvl(id, on ? 0.12 : 0);
      setVT(id, on ? `${fmt(s.L)} L/h` : 'idle', on && s.g ? `${fmt(s.g)} g/h` : '');
    };
    pass('T103', c.vi); pass('T201', c.aex); pass('T202', c.cex); pass('T301', c.vf); pass('T302', v['59']);
    const aOn = c.aexIn.L > 1e-6, cOn = c.cexIn.L > 1e-6;
    colCls('C201', aOn && c.duty === 0 ? 'load' : ''); colCls('C202', aOn && c.duty === 1 ? 'load' : '');
    colCls('C203', cOn && c.duty === 0 ? 'load' : ''); colCls('C204', cOn && c.duty === 1 ? 'load' : '');
    setVT('AEX', `col. ${c.duty + 1} on duty` + (aOn ? ` · ${fmt(c.aexIn.g)} g/h` : ''), '');
    setVT('CEX', `col. ${c.duty + 1} on duty` + (cOn ? ` · ${fmt(c.cexIn.g)} g/h` : ''), '');
    setVT('VRF', `filter ${c.duty ? 'B' : 'A'} on duty` + (c.vfIn.L > 1e-6 ? ` · ${fmt(c.vfIn.L)} L/h` : ''), '');
    setVT('UF', v['59'] && v['59'].L > 1e-6 ? `${fmt(v['59'].L)} L/h retentate` : '', '');
    const lotH = Math.max(num(P.lot_window_h) || 24, 1);
    lvl('M301', inN && r.ds_g_cum > 0 ? ((Math.max(r.n_day, 0) * 24) % lotH) / lotH : 0);
    setVT('M301', r.ds_g_h > 0 ? `${fmt(r.ds_g_h)} g/h` : '', '');
    setVT('CV', `DS ${fmt(r.ds_g_cum / 1000)} kg`, '');
  }

  function lastFc(i) {
    for (let k = i; k >= 0 && k > i - 12; k--) {
      const v = S[k].titre_fc_24h;
      if (v !== null && v !== undefined) return v;
    }
    return null;
  }

  function readout(i, r0, r) {
    const fc = lastFc(i);
    const lines = [
      `Campaign day ${fmt(r.day)} · ` + (r0.phase ? `N day ${fmt(r.n_day)}` : 'N-1 seed phase'),
      `Permeate titre ${r0.titre_harvest != null ? fmt(r.titre_harvest) + ' g/L' : '—'} · NevAI +24 h ` +
        (fc != null ? fmt(fc) + ' g/L' : '—'),
      `Surge tank ${fmt(r.sump_L)} L, ${fmt(r.sump_g)} g · Protein A cycles ${Math.round(r.cycles)}`,
      `Buffer use ${fmt(r.buffer_L_h)} L/h · wastewater ${fmt(r.ww_total_Lh)} L/h`,
      `Drug substance so far ${fmt(r.ds_g_cum / 1000)} kg`,
    ];
    const clock = `hour ${r.h} of ${S[S.length - 1].h} · campaign day ${fmt(r.day)}` +
      (r0.phase ? ` · N day ${fmt(r.n_day)}` : ' · N-1 phase');
    return {lines, clock, fc};
  }

  // everything the diagram shows at series row i (used by the live view and the video)
  function state(i) {
    if (!S.length || !P) return null;
    i = Math.max(0, Math.min(S.length - 1, i | 0));
    const r0 = S[i], r = {};
    for (const k in r0) r[k] = num(r0[k]);
    const c = ctx(r0, r), vals = {};
    for (const [id, , fn] of STREAMS) vals[id] = fn ? fn(r, c) : null;
    const eq = {lv: {}, vt: {}, col: {}};
    equipment(r0, r, c, vals, eq);
    const ro = readout(i, r0, r);
    return {i, r0, r, c, vals, eq, lines: ro.lines, clock: ro.clock, fc: ro.fc};
  }

  function apply(st) {
    for (const [id, desc] of STREAMS) {
      const v = st.vals[id];
      const L = v ? num(v.L) : 0;
      const g = v && v.g !== undefined && v.g !== null ? num(v.g) : null;
      const on = L > 1e-6;
      const tip = `Stream ${id}: ${desc}\n` + (v ? `${fmt(L)} L/h` +
        (g !== null ? ` · ${fmt(g)} g/h` + (on ? ` · ${fmt(g / L)} g/L` : '') : '') : 'not modelled in the simulator');
      for (const el of ELS[id] || []) {
        el.classList.toggle('on', on);
        el.classList.toggle('nm', !v);
        const t = el.querySelector('title');
        if (t) t.textContent = tip;
      }
      const tr = ROWS[id];
      if (tr) {
        const td = tr.children;
        td[2].textContent = v ? fmt(L) : 'not modelled';
        td[3].textContent = g !== null ? fmt(g) : '—';
        td[4].textContent = g !== null && on ? fmt(g / L) : '—';
        tr.classList.toggle('on', on);
      }
    }
    for (const id in st.eq.lv) domLvl(id, st.eq.lv[id]);
    for (const id in st.eq.vt) domVT(id, st.eq.vt[id][0], st.eq.vt[id][1]);
    for (const id in st.eq.col) domCol(id, st.eq.col[id]);
    st.lines.forEach((s, k) => { const e = $('ro' + k); if (e) e.textContent = s; });
    clockTxt = st.clock;
    const k1 = $('pfdClock'); if (k1) k1.textContent = clockTxt;
    const k2 = $('pfd-clk'); if (k2) k2.textContent = 'Live playback · ' + clockTxt;
    const sl = $('pfdSlider');
    if (sl && +sl.value !== st.i) sl.value = st.i;
  }

  function show(i) {
    const st = state(i);
    if (!st) return;
    idx = st.i;
    apply(st);
  }

  function speed() { const s = $('pfdSpeed'); return s ? (+s.value || 6) : 6; }
  function pause() {
    clearInterval(timer); timer = null;
    const b = $('pfdPlay'); if (b) b.textContent = '▶ Play';
  }
  function play() {
    if (!S.length) return;
    if (pos >= S.length - 1) pos = 0;
    clearInterval(timer);
    timer = setInterval(() => {
      pos = Math.min(pos + speed() / 10, S.length - 1);
      show(Math.floor(pos));
      if (pos >= S.length - 1) pause();
    }, 100);
    const b = $('pfdPlay'); if (b) b.textContent = '⏸ Pause';
  }

  function note() {
    const e = $('pfdNote');
    if (!e) return;
    const legacy = P.polish && P.polish.order === 'cex_aex';
    e.innerHTML = (legacy ? '<b>Note:</b> polish.order = cex_aex, so the simulator runs CEX before AEX; ' +
      'each unit still shows its own values, but the drawn order differs. ' : '') +
      'Grey dashed items are on the site flowsheet but not in the simulator (detergent, separate WFI, depth filtration, ' +
      'on-line acid/base, holds in the intermediate SUMs). Protein A pairs alternate each cycle and AEX/CEX/virus-filter ' +
      'duty swaps each lot; the model itself runs one column train.' +
      (P.upstream.mode === 'mechanistic' ? ' Legacy ATF mode: U5 is the cell bleed.' : '');
  }

  function buildTable() {
    const t = $('pfdTable');
    if (!t) return;
    t.innerHTML = '<thead><tr><th>Stream</th><th>Description</th><th>L/h</th><th>g/h</th><th>g/L</th></tr></thead><tbody>' +
      STREAMS.map(([id, d]) => `<tr data-row="${id}"><td><b>${esc(id)}</b></td><td>${esc(d)}</td><td></td><td></td><td></td></tr>`).join('') +
      '</tbody>';
    t.querySelectorAll('tr[data-row]').forEach(tr => { ROWS[tr.dataset.row] = tr; });
  }

  function downloadSvg() {
    const s = document.querySelector('#pfdBox svg');
    if (!s) return;
    const blob = new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(s)],
                          {type: 'image/svg+xml'});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'IDB_PFD_hour_' + (S[idx] ? S[idx].h : 0) + '.svg';
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }

  function init() {
    const box = $('pfdBox');
    if (!box) return;
    box.innerHTML = build();
    box.querySelectorAll('[data-s]').forEach(el => (ELS[el.dataset.s] = ELS[el.dataset.s] || []).push(el));
    buildTable();
    const pb = $('pfdPlay'); if (pb) pb.onclick = () => (timer ? pause() : play());
    const sl = $('pfdSlider'); if (sl) sl.oninput = e => { pause(); pos = +e.target.value; show(pos); };
    const sv = $('pfdSvg'); if (sv) sv.onclick = downloadSvg;
  }

  function load(j) {
    if (!j || !j.series || !j.params) return;
    S = j.series; P = j.params;
    const cp = P.capture, pc = P.pcc;
    ELU = cp.mode === 'pcc' ? num(pc.eluate_CV) * num(pc.CV_L)
                            : num(cp.eluate_CV) * num(cp.column_capacity_g) / Math.max(num(cp.DBC_gL), 1e-9);
    const sl = $('pfdSlider'); if (sl) sl.max = Math.max(0, S.length - 1);
    note();
    pos = 0;
    show(0);
    play();
  }

  window.PFD = {
    init, load, pause, clock: () => clockTxt,
    svg: opts => build(opts), geo: () => GEO, state,
    series: () => S, params: () => P, length: () => S.length,
  };
})();
