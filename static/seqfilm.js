/* static/seqfilm.js — step-sequence film for the IDB mAb digital twin
   NavNevAI_mAbWasteWater · CarbonTesseract 2026

   Moving Gantt chart of every downstream step of the last run: one lane per
   column, filter or in-line addition, coloured by the solution in use, with
   vessel levels underneath. Play, scrub, jump to a unit, record MP4 / WebM.
   Uses result.steps, result.lanes, result.series, result.solution_labels.
   API: SeqFilm.init()  SeqFilm.load(simResult)  SeqFilm.stop()
*/
(function () {
  'use strict';
  const W = 1920, HH = 1080;
  const FONT = '-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Arial,sans-serif';
  const SOL = {wfi: '#7dd3fc', brij_10: '#a855f7', tris_74: '#0e7c7b', pa_wash1: '#f59e0b', pa_wash2: '#fbbf24',
    pa_elute: '#dc2626', pa_strip: '#f43f5e', naoh_05n: '#475569', naoh_1n: '#334155', naoh_01n: '#64748b',
    storage: '#16a34a', df_ac50: '#fb923c', nacl_05m: '#60a5fa', nacl_2m: '#3b82f6', df_arg: '#c084fc',
    tris_base_2m: '#14b8a6', tris_nacl_74: '#2563eb', acetic_1m: '#e11d48', cex_eq: '#eab308',
    cex_elute: '#be123c', uf_df1: '#8b5cf6', uf_df2: '#6366f1', uf_spike: '#db2777', vi_acid: '#b91c1c', vi_base: '#0891b2'};
  const KIND = {load: '#94a3b8', hold: '#cbd5e1', perm: '#c7d2fe', out: '#0ea5e9'};
  const VES = [['sump_L', 'Surge tank'], ['eluate_surge_L', 'Protein A vessel'], ['vip_L', 'VI pool'],
    ['dfs_L', 'Depth filtr. surge'], ['aexv_L', 'AEX vessel'], ['cexv_L', 'CEX vessel'], ['vrfv_L', 'VRF vessel'],
    ['ufdf_L', 'UF/DF vessel'], ['bdsm_L', 'BDS mixer']];
  const $ = id => document.getElementById(id);
  let D = null, now = 0, playing = false, raf = 0, last = 0, rec = null, cv = null, ctx = null;

  function info(t) { const e = $('seqInfo'); if (e) e.textContent = t; }
  function col(s) { return (s.sol && SOL[s.sol]) || KIND[s.kind] || '#94a3b8'; }
  function ink(hex) {
    const n = parseInt(hex.slice(1), 16);
    return (0.299 * (n >> 16) + 0.587 * ((n >> 8) & 255) + 0.114 * (n & 255)) / 255 > 0.6 ? '#0f172a' : '#ffffff';
  }
  function rr(x, y, w, h, r) {
    r = Math.max(0, Math.min(r, w / 2, h / 2));
    ctx.beginPath(); ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
    ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
  }
  function txt(s, x, y, size, color, align, weight, maxW) {
    ctx.font = (weight || '500') + ' ' + size + 'px ' + FONT;
    s = String(s);
    if (maxW) {
      if (maxW < 12) return;
      while (s.length > 1 && ctx.measureText(s).width > maxW) s = s.slice(0, -2) + '…';
    }
    ctx.fillStyle = color; ctx.textAlign = align || 'left'; ctx.textBaseline = 'middle'; ctx.fillText(s, x, y);
  }
  function val(k, t) {
    const S = D.S; if (!S.length) return null;
    const x = Math.max(0, t / D.dh - 1), i = Math.min(S.length - 1, Math.floor(x)), j = Math.min(S.length - 1, i + 1);
    const a = S[i][k], b = S[j][k];
    if (typeof a !== 'number') return typeof b === 'number' ? b : null;
    if (typeof b !== 'number') return a;
    return a + (b - a) * Math.min(1, x - i);
  }
  function speed() { return parseFloat($('seqSpeed').value) || 0.5; }
  function sync() { const s = $('seqSlider'); if (s) s.value = String(now); }

  function load(j) {
    stop();
    const steps = (j.steps || []).slice().sort((a, b) => a.t0 - b.t0);
    const S = j.series || [];
    const dh = S.length > 1 ? (S[1].h - S[0].h) || 1 : 1;
    const labels = {}; (j.lanes || []).forEach(l => { labels[l.key] = l.label; });
    const lanes = (j.lanes || []).map(l => l.key);
    const byLane = {}; lanes.forEach(k => { byLane[k] = steps.filter(s => s.lane === k); });
    const vess = VES.filter(([k]) => S.some(r => typeof r[k] === 'number' && r[k] > 1e-6));
    const peaks = {}; vess.forEach(([k]) => { peaks[k] = Math.max(1e-6, ...S.map(r => (typeof r[k] === 'number' ? r[k] : 0))); });
    D = {steps, S, dh, labels, lanes, byLane, vess, peaks, H: S.length ? S[S.length - 1].h + dh : 0,
         solLab: j.solution_labels || {}};
    const first = steps.find(s => s.lane.indexOf('pa') === 0 && s.kind === 'load');
    now = first ? Math.max(0, first.t0 - 1) : 0;
    const sl = $('seqSlider'); if (sl) sl.max = String(D.H);
    sync();
    info(`${steps.length} steps logged in ${lanes.length} lanes. Play, scrub or jump to a unit.`);
    draw();
  }

  function draw() {
    if (!ctx) return;
    ctx.fillStyle = '#f8fafc'; ctx.fillRect(0, 0, W, HH);
    if (!D || !D.lanes.length) { txt('Run a simulation to see the step sequence', W / 2, HH / 2, 40, '#64748b', 'center'); return; }
    const win = parseFloat($('seqWin').value) || 12, X0 = 360, X1 = W - 40, TOP = 150;
    const nL = D.lanes.length, LH = Math.min(50, 600 / nL);
    const tA = now - win * 0.3, tB = tA + win, sx = t => X0 + (t - tA) / (tB - tA) * (X1 - X0);
    txt('IDB downstream — every step in sequence (Case 5 flowsheet)', 40, 44, 34, '#16202c', 'left', '700');
    const d = Math.floor(now / 24), hh = now - d * 24, hr = Math.floor(hh), mn = Math.floor((hh - hr) * 60);
    const nd = val('n_day', now);
    txt(`Campaign day ${d + 1}  ${String(hr).padStart(2, '0')}:${String(mn).padStart(2, '0')}` +
        (nd !== null && nd >= 0 ? `  ·  N day ${nd.toFixed(2)}` : '  ·  N-1 seed'), W - 40, 44, 28, '#0e7c7b', 'right', '600');
    const active = [];
    for (const k of D.lanes) for (const s of D.byLane[k]) {
      if (s.t0 > now) break;
      if (s.t1 > now && s.sol && active.indexOf(s.sol) < 0) active.push(s.sol);
    }
    txt('Solutions in use:', 40, 98, 20, '#5b6675', 'left', '600');
    let cx = 220;
    for (const k of active) {
      const lab = D.solLab[k] || k;
      ctx.font = '500 18px ' + FONT;
      const w = Math.min(ctx.measureText(lab).width + 46, 520);
      if (cx + w > W - 40) break;
      rr(cx, 82, w, 32, 16); ctx.fillStyle = '#ffffff'; ctx.fill();
      ctx.strokeStyle = SOL[k] || '#94a3b8'; ctx.lineWidth = 2; ctx.stroke();
      ctx.beginPath(); ctx.arc(cx + 18, 98, 8, 0, 6.2832); ctx.fillStyle = SOL[k] || '#94a3b8'; ctx.fill();
      txt(lab, cx + 32, 98, 18, '#16202c', 'left', '500', w - 40);
      cx += w + 10;
    }
    if (!active.length) txt('none (idle, loading or holding)', cx, 98, 20, '#94a3b8');
    const tick = win <= 3 ? 0.5 : win <= 12 ? 1 : win <= 24 ? 2 : 6;
    ctx.lineWidth = 1;
    for (let t = Math.max(0, Math.ceil(tA / tick) * tick); t <= tB; t += tick) {
      const x = sx(t);
      ctx.strokeStyle = '#e2e8f0'; ctx.beginPath(); ctx.moveTo(x, TOP - 8); ctx.lineTo(x, TOP + nL * LH); ctx.stroke();
      const dd = Math.floor(t / 24), h2 = t - dd * 24;
      txt(h2 < 1e-6 ? 'Day ' + (dd + 1) : String(Math.floor(h2)).padStart(2, '0') + ':' + String(Math.round((h2 % 1) * 60)).padStart(2, '0'),
          x, TOP - 22, 16, '#5b6675', 'center');
    }
    D.lanes.forEach((k, i) => {
      const y = TOP + i * LH;
      ctx.fillStyle = i % 2 ? '#eef2f7' : '#ffffff'; ctx.fillRect(X0, y, X1 - X0, LH - 2);
      txt(D.labels[k] || k, X0 - 16, y + LH / 2, Math.min(20, LH * 0.45), '#16202c', 'right', '600', X0 - 40);
      for (const s of D.byLane[k]) {
        if (s.t0 > tB) break;
        const t1 = Math.max(s.t1, s.t0 + 0.02);
        if (t1 < tA) continue;
        const a = sx(Math.max(s.t0, tA)), b = sx(Math.min(t1, tB)), on = s.t0 <= now && now < t1, done = t1 <= now;
        const c = col(s);
        ctx.globalAlpha = on ? 1 : done ? 0.6 : 0.3;
        rr(a, y + 5, Math.max(b - a, 3), LH - 12, 6); ctx.fillStyle = c; ctx.fill();
        ctx.globalAlpha = 1;
        if (on) { ctx.strokeStyle = '#0f172a'; ctx.lineWidth = 3; ctx.stroke(); }
        if (b - a > 60) txt(s.step, a + 8, y + LH / 2, Math.min(17, LH * 0.38), on || done ? ink(c) : '#334155', 'left', '600', b - a - 14);
      }
    });
    const xn = sx(now);
    ctx.strokeStyle = '#0f172a'; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.moveTo(xn, TOP - 10); ctx.lineTo(xn, TOP + nL * LH + 4); ctx.stroke();
    const vy = TOP + nL * LH + 30, vh = HH - vy - 70, nV = D.vess.length;
    if (nV) {
      txt('Vessel levels (bar = share of the run peak)', 40, vy + 4, 18, '#5b6675', 'left', '600');
      const gap = 24, vw = (W - 80 - gap * (nV - 1)) / nV;
      D.vess.forEach(([k, lab], i) => {
        const x = 40 + i * (vw + gap), top = vy + 22, hgt = Math.max(vh - 24, 30);
        const v = Math.max(0, val(k, now) || 0), f = Math.min(1, v / D.peaks[k]);
        rr(x, top, vw, hgt, 12); ctx.fillStyle = '#ffffff'; ctx.fill(); ctx.strokeStyle = '#cbd5e1'; ctx.lineWidth = 2; ctx.stroke();
        ctx.save(); rr(x, top, vw, hgt, 12); ctx.clip();
        ctx.fillStyle = '#7dd3c8'; ctx.fillRect(x, top + hgt * (1 - f), vw, hgt * f); ctx.restore();
        txt((v >= 100 ? v.toFixed(0) : v.toFixed(1)) + ' L', x + vw / 2, top + hgt / 2, 24, '#0f172a', 'center', '700');
        txt(lab, x + vw / 2, top + hgt + 18, 17, '#334155', 'center', '600', vw);
      });
    }
    const dsg = val('ds_g_cum', now), cyc = val('cycles', now);
    txt(`Drug substance so far ${dsg === null ? '—' : (dsg / 1000).toFixed(3) + ' kg'} · Protein A cycles ${cyc === null ? '—' : Math.round(cyc)}` +
        ' · grey = loading/filtering · pale = concentrating · outline = running now', W - 40, HH - 18, 16, '#5b6675', 'right');
  }

  function tick(ts) {
    if (!playing) return;
    const dt = last ? (ts - last) / 1000 : 0; last = ts;
    now += dt * speed();
    const end = rec ? rec.endAt : D.H;
    if (now >= end) { now = end; sync(); draw(); pause(); if (rec) finishRec(); return; }
    sync(); draw();
    raf = requestAnimationFrame(tick);
  }
  function play() { if (!D || playing) return; playing = true; last = 0; $('seqPlay').textContent = '❚❚ Pause'; raf = requestAnimationFrame(tick); }
  function pause() { playing = false; cancelAnimationFrame(raf); const b = $('seqPlay'); if (b) b.textContent = '▶ Play'; }
  function finishRec() {
    const r = rec; rec = null;
    if (r && r.mr.state !== 'inactive') setTimeout(() => { try { r.mr.stop(); } catch (e) { /* ignore */ } }, 300);
  }
  function stop() { pause(); if (rec) finishRec(); }
  function pickMime() {
    const c = ['video/mp4;codecs=avc1.640028', 'video/mp4;codecs=avc1', 'video/mp4', 'video/webm;codecs=vp9', 'video/webm'];
    return c.find(m => { try { return MediaRecorder.isTypeSupported(m); } catch (e) { return false; } }) || '';
  }

  function record() {
    if (!D) { info('Run a simulation first.'); return; }
    if (!window.MediaRecorder || !cv.captureStream) { info('This browser cannot record canvas video.'); return; }
    if (rec) return;
    pause();
    const len = $('seqLen').value;
    const endAt = len === 'all' ? D.H : Math.min(D.H, now + parseFloat(len));
    const mime = pickMime();
    let mr;
    try { mr = new MediaRecorder(cv.captureStream(30), mime ? {mimeType: mime, videoBitsPerSecond: 8e6} : {videoBitsPerSecond: 8e6}); }
    catch (e) { info('Recording failed: ' + e.message); return; }
    const chunks = [];
    mr.ondataavailable = e => { if (e.data && e.data.size) chunks.push(e.data); };
    mr.onstop = () => {
      const type = (mr.mimeType || mime || 'video/webm').split(';')[0];
      const blob = new Blob(chunks, {type}), url = URL.createObjectURL(blob);
      const ext = type.indexOf('mp4') >= 0 ? 'mp4' : 'webm';
      $('seqOut').innerHTML = `<video controls playsinline src="${url}"></video><div class="vidlinks">` +
        `<a class="btnlink" href="${url}" download="IDB_step_sequence.${ext}">⬇ Download ${ext.toUpperCase()}</a>` +
        `<span class="muted">${(blob.size / 1e6).toFixed(1)} MB</span></div>`;
      $('seqRecord').disabled = false; $('seqStopRec').disabled = true;
      info('Recording finished.');
    };
    rec = {mr, endAt};
    $('seqRecord').disabled = true; $('seqStopRec').disabled = false;
    info(`Recording ${(endAt - now).toFixed(1)} h of run time ≈ ${((endAt - now) / speed()).toFixed(0)} s. Keep this tab visible.`);
    draw(); mr.start(500); play();
  }

  function jump(v) {
    if (!D || !v) return;
    const s = D.steps.find(x => v === 'end' ? x.step.indexOf('Store') === 0 : (x.lane.indexOf(v) === 0 && x.kind === 'load'));
    if (!s) { info('That unit has no step in this run.'); return; }
    now = Math.max(0, s.t0 - 0.5); sync(); draw();
  }

  function init() {
    cv = $('seqCanvas'); if (!cv) return;
    ctx = cv.getContext('2d');
    $('seqPlay').onclick = () => (playing ? pause() : play());
    $('seqSlider').oninput = e => { now = +e.target.value; if (!playing) draw(); };
    $('seqWin').onchange = draw;
    $('seqJump').onchange = e => { jump(e.target.value); e.target.value = ''; };
    $('seqRecord').onclick = record;
    $('seqStopRec').onclick = () => { pause(); if (rec) finishRec(); };
    const ps = $('pfdSlider');
    if (ps) ps.addEventListener('input', () => {
      if (D && !playing && $('seqFollow').checked) { now = Math.min(+ps.value * D.dh, D.H); sync(); draw(); }
    });
    draw();
  }

  window.SeqFilm = {init, load, stop};
})();
