export async function api(path, options = {}) {
  const response = await fetch(path, options);
  const data = await response.json();
  if (!response.ok) throw new Error(typeof data.detail === 'string' ? data.detail : JSON.stringify(data.detail || data));
  return data;
}
export const $ = id => document.getElementById(id);
export const escapeHTML = value => String(value ?? '').replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
export const fmt = (value, digits = 2) => value === null || value === undefined ? '—' : Number(value).toLocaleString(undefined, {maximumFractionDigits:digits});
export function filteredRows(rows, query, tier, sorts) {
  const q = query.trim().toLowerCase();
  return rows.filter(r => (!tier || r.tier === tier) && ['Gene','ProteinChange','Neopeptide_9mer','WT_9mer'].map(k => r[k] || '').join(' ').toLowerCase().includes(q)).sort((a,b) => {
    for (const [key,direction] of sorts) {
      const av=a[key], bv=b[key];
      if (av == null && bv == null) continue;
      if (av == null) return 1;
      if (bv == null) return -1;
      const cmp=typeof av === 'number' ? av-bv : (av < bv ? -1 : av > bv ? 1 : 0);
      if(cmp) return direction === 'asc' ? cmp : -cmp;
    }
    return 0;
  });
}
export const tierNames = {0:'Unclassified',1:'High Priority / Vaccine Payload',2:'Secondary / Backup',3:'Deprioritized'};
export const badge = r => `<span class="tier-badge tier-${r.tier}">Tier ${r.tier || '?'} · ${tierNames[r.tier]}</span>`;
export function peptide(r) {
  return [...String(r.Neopeptide_9mer || '')].map((ch,i) => i+1 === r.MutPos ? `<mark>${escapeHTML(ch)}</mark>` : escapeHTML(ch)).join('');
}
