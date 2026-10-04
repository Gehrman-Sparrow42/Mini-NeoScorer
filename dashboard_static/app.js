import {api,$,escapeHTML as esc,fmt,filteredRows,badge,peptide} from './data.js';
import {renderCharts,clearCharts} from './charts.js';
const columns=[['Rank','Rank'],['Gene','Gene'],['ProteinChange','Protein change'],['Neopeptide_9mer','Neopeptide (9-mer)'],['WT_9mer','WT 9-mer'],['MutPos','Mut. pos.'],['AF','VAF'],['Reads','Reads'],['TPM','RNA TPM'],['MT_IC50_nM','IC50 (nM)'],['Percentile_Rank','Rank (%)'],['Agretopicity','Agretopicity (×)'],['Score','Score'],['ClinicalTier','Clinical tier'],['Action','Action'],['Predictor','Predictor']];
let report=null, visible=[], selected=null, page=0, generation=0, jobTimer=null;
let sorts=[['Score','desc'],['MT_IC50_nM','asc'],['Agretopicity','desc']];
const pageSize=20;
function notice(message){$('notice').textContent=message;$('notice').hidden=!message;}
function resetState(){ $('search').value='';$('tier').value='0';$('sort').value='Score:desc';sorts=[['Score','desc'],['MT_IC50_nM','asc'],['Agretopicity','desc']];page=0;selected=null; }
async function discover(preferred){
  const patients=await api('/api/patients');
  const current=preferred || $('patient').value;
  $('patient').replaceChildren(...patients.map(p=>new Option(`Patient: ${p.patient}`,p.id)));
  if(!patients.length){generation++;report=null;clearCharts();$('dashboard').hidden=true;$('loading').hidden=false;$('loading').textContent='No patient reports yet. Run a new patient analysis to get started.';$('patient').append(new Option('No reports available',''));return;}
  $('patient').value=patients.some(p=>p.id===current)?current:patients[0].id;
  await loadPatient();
}
async function loadPatient(){
  const token=++generation;
  report=null;clearCharts();$('dashboard').hidden=true;$('loading').hidden=false;$('loading').textContent='Loading patient report…';notice('');resetState();
  try{
    const next=await api(`/api/patients/${encodeURIComponent($('patient').value)}`);
    if(token!==generation)return;
    report=next;$('patient-name').textContent=next.patient;$('hla').textContent=next.hla;$('hla-source').textContent=next.hla_source;
    $('updated').textContent=new Date(next.modified*1000).toLocaleString();$('source').textContent=next.source;
    $('total').textContent=fmt(next.total,0);$('strong').textContent=fmt(next.strong,0);$('tier1').textContent=fmt(next.tier1,0);$('tier2').textContent=fmt(next.rows.filter(r=>r.tier===2).length,0);
    $('strong-note').textContent=`IC50 ≤ 50 nM or rank ≤ 0.5% · ${next.strong_affinity} affinity-only`;
    const stats=next.metadata.stats || {};
    $('provenance').innerHTML=[['Source',next.source],['HLA evidence',next.hla_source],['Somatic variants in source',stats.total_somatic_mutations ?? 'Not recorded'],['9-mers evaluated',stats.total_9mers_evaluated ?? 'Not recorded'],['RNA gate (TPM)',next.metadata.tpm_threshold ?? 'Not recorded']].map(([k,v])=>`<div class="metric-row"><span>${esc(k)}</span><strong>${esc(v)}</strong></div>`).join('');
    $('dashboard').hidden=false;$('loading').hidden=true;render();
  }catch(error){if(token===generation){report=null;clearCharts();$('dashboard').hidden=true;$('loading').textContent='This report could not be loaded. Select another patient or refresh.';notice(error.message);}}
}
function render(){
  if(!report)return;
  visible=filteredRows(report.rows,$('search').value,Number($('tier').value),sorts);
  page=Math.min(page,Math.max(0,Math.ceil(visible.length/pageSize)-1));
  $('result-count').textContent=`${visible.length} / ${report.total}`;
  $('sort-description').textContent='Sort: '+sorts.map(([k,d],i)=>`${i+1}. ${columns.find(c=>c[0]===k)?.[1] || k} ${d==='asc'?'↑':'↓'}`).join(' · ');
  if(!visible.some(r=>r.uid===selected?.uid))selected=visible.find(r=>r.tier===1)||visible[0]||null;
  renderTable();renderConsensus();renderCharts(visible,selectRow);
  $('export-tsv').disabled=$('export-csv').disabled=!visible.length;
}
function renderTable(){
  $('columns').innerHTML=columns.map(([key,label])=>{const pos=sorts.findIndex(s=>s[0]===key),dir=pos>=0?sorts[pos][1]:'';return `<th scope="col" ${pos===0?`aria-sort="${dir==='asc'?'ascending':'descending'}"`:''}><button data-sort="${key}">${label}${pos>=0?` ${dir==='asc'?'▲':'▼'}${pos+1}`:''}</button></th>`;}).join('');
  $('rows').innerHTML=visible.slice(page*pageSize,(page+1)*pageSize).map(r=>`<tr class="${r.uid===selected?.uid?'selected':''}">${columns.map(([key])=>{
    let value=esc(r[key] ?? '—'),cls='';
    if(key==='Gene')value=`<button class="gene-button" data-uid="${r.uid}" aria-label="Inspect ${esc(r.Gene)} ${esc(r.Neopeptide_9mer)}">${esc(r.Gene)} ↗</button>`;
    else if(key==='ClinicalTier')value=badge(r);
    else if(key==='Neopeptide_9mer'){value=peptide(r);cls='peptide';}
    else if(key==='WT_9mer')cls='peptide';
    else if(typeof r[key]==='number')value=fmt(r[key],key==='AF'?3:key==='Rank'||key==='MutPos'?0:2);
    if(key==='Score')cls='score';if(key==='Action')cls='action-cell';if(key==='Predictor')cls='predictor-cell';
    return `<td class="${cls}">${value}</td>`;
  }).join('')}</tr>`).join('');
  $('empty-results').hidden=visible.length>0;
  $('empty-results').textContent=report.total?'No candidates match these filters. Try another gene or reset the filters.':'This report contains no evaluated candidates. Review the pipeline log and input data.';
  $('page-info').textContent=visible.length?`Showing ${page*pageSize+1}–${Math.min((page+1)*pageSize,visible.length)} of ${visible.length} candidates`:'0 candidates';
  $('page-number').textContent=`${page+1} / ${Math.max(1,Math.ceil(visible.length/pageSize))}`;$('previous').disabled=page===0;$('next').disabled=(page+1)*pageSize>=visible.length;
}
function selectRow(row){selected=row;renderTable();renderConsensus();$('consensus').scrollIntoView({behavior:'smooth',block:'center'});}
function metric(label,value,suffix=''){return `<div class="metric-row"><span>${label}</span><strong>${value==null?'Not recorded':fmt(value,4)+suffix}</strong></div>`;}
function renderConsensus(){
  if(!selected){$('consensus-content').innerHTML='<p class="muted">Select a candidate from the table or chart to inspect its model evidence.</p>';return;}
  const r=selected, isFlurry=/MHCflurry|Ensemble/i.test(r.Predictor),isNet=/NetMHCpan/i.test(r.Predictor),isSingleNet=/^NetMHCpan_4/.test(r.Predictor);
  const mf=r.MHCflurry_IC50 ?? (/^MHCflurry_2/.test(r.Predictor)?r.MT_IC50_nM:null);
  const net=/Filter_NonBinder|Unavailable_Fallback/.test(r.Predictor)?null:(r.NetMHCpan_IC50 ?? (isSingleNet?r.MT_IC50_nM:null));
  $('consensus-content').innerHTML=`<div class="selected-title"><h3>${esc(r.Gene)}</h3><span class="muted">${esc(r.ProteinChange)}</span><span class="peptide">${peptide(r)}</span></div><div class="model-grid"><div class="model-box"><h3>MHCflurry 2.0</h3><small>Binding affinity & antigen presentation</small>${metric('Predicted affinity',mf,' nM')}${metric('Presentation score',isFlurry?r.Presentation_Score:null)}</div><div class="model-box"><h3>NetMHCpan-4.1</h3><small>Pan-allele binding · NIH IEDB</small>${metric('Predicted affinity',net,' nM')}${metric('Model percentile rank',isSingleNet?r.Percentile_Rank:null,' %')}</div></div><div class="metric-row"><span>Report IC50 / percentile</span><strong>${fmt(r.MT_IC50_nM)} nM / ${fmt(r.Percentile_Rank)}%</strong></div><p class="model-note">${esc(r.Predictor)}<br>${isFlurry&&isNet?'Ensemble percentile is a combined metric; individual model ranks are not stored in this report. ':''}${/Filter_NonBinder/.test(r.Predictor)?'NetMHCpan was skipped by the screening gate; its placeholder is not a model prediction. ':''}${!isFlurry?'This predictor does not provide an MHCflurry processing score. ':''}Presentation is an integrated model score, not a separately measured cleavage probability.</p><p class="model-note">${esc(r.Action)}</p>`;
}
function exportReport(format){if(!report)return;const q=new URLSearchParams({q:$('search').value.trim(),tier:$('tier').value,sorts:sorts.map(s=>s.join(':')).join(','),format});window.location.href=`/api/patients/${encodeURIComponent(report.id)}/export?${q}`;}
async function openRun(){
  $('run-dialog').showModal();
  try{const inputs=await api('/api/inputs');for(const [kind,id]of [['mutations','mutation-file'],['expression','expression-file']])$(''+id).replaceChildren(new Option('Select a local file…',''),...inputs[kind].map(name=>new Option(name,name)));await pollJob();}catch(e){$('run-status').textContent=e.message;}
}
async function upload(id){const file=$(id).files[0];if(!file)return undefined;if(file.size>18_000_000)throw new Error('JSON files must be smaller than 18 MB.');try{return JSON.parse((await file.text()).replace(/^\uFEFF/,''));}catch{throw new Error(`Invalid JSON in ${file.name}`);}}
async function submitRun(event){
  event.preventDefault();$('run-submit').disabled=true;$('run-status').textContent='Validating inputs…';
  try{
    const [mutations,expression]=await Promise.all([upload('mutation-upload'),upload('expression-upload')]);
    const body={patient_id:$('run-patient').value.trim()||null,mutations_file:$('mutation-file').value||null,expression_file:$('expression-file').value||null,mutations,expression,predictor:$('predictor').value,hla:$('run-hla').value.trim(),tpm_threshold:Number($('run-tpm').value),top:Number($('run-top').value)};
    await api('/api/jobs',{method:'POST',headers:{'Content-Type':'application/json','X-NeoScorer':'dashboard'},body:JSON.stringify(body)});await pollJob();
  }catch(e){$('run-status').textContent=e.message;$('run-submit').disabled=false;}
}
let completedJob=null;
async function pollJob(){
  clearTimeout(jobTimer);
  try{const job=await api('/api/jobs/current');$('run-submit').disabled=job.status==='running';$('run-log').textContent=job.log||'No active run.';
    if(job.status!=='idle')$('run-status').textContent=`${job.status.toUpperCase()} · ${job.patient}\n${job.message}`;
    if(job.status==='running')jobTimer=setTimeout(pollJob,1500);
    if(job.status==='completed'&&completedJob!==job.id){completedJob=job.id;await discover(job.report_id);notice(`Analysis complete for ${job.patient}. The new report is selected.`);}
  }catch(e){$('run-status').textContent=`Status unavailable: ${e.message}. Retrying…`;jobTimer=setTimeout(pollJob,4000);}
}
$('patient').addEventListener('change',loadPatient);$('refresh').addEventListener('click',()=>discover().catch(e=>notice(e.message)));
$('search').addEventListener('input',()=>{page=0;render();});$('tier').addEventListener('change',()=>{page=0;render();});
$('sort').addEventListener('change',()=>{const primary=$('sort').value.split(':');sorts=[primary,...[['Score','desc'],['MT_IC50_nM','asc'],['Agretopicity','desc']].filter(s=>s[0]!==primary[0])];page=0;render();});
$('reset').addEventListener('click',()=>{resetState();render();});
$('columns').addEventListener('click',event=>{const button=event.target.closest('[data-sort]');if(!button)return;const key=button.dataset.sort,existing=sorts.find(s=>s[0]===key);const spec=[key,existing?.[1]==='asc'?'desc':'asc'];sorts=event.shiftKey?[...sorts.filter(s=>s[0]!==key),spec]:[spec];page=0;render();});
$('rows').addEventListener('click',event=>{const button=event.target.closest('[data-uid]');if(button)selectRow(visible.find(r=>r.uid===Number(button.dataset.uid)));});
$('previous').addEventListener('click',()=>{page--;renderTable();});$('next').addEventListener('click',()=>{page++;renderTable();});
$('export-tsv').addEventListener('click',()=>exportReport('tsv'));$('export-csv').addEventListener('click',()=>exportReport('csv'));
$('open-run').addEventListener('click',openRun);$('new-run').addEventListener('click',openRun);$('close-run').addEventListener('click',()=>$('run-dialog').close());$('run-form').addEventListener('submit',submitRun);
await discover().catch(e=>{notice(e.message);$('loading').textContent='Unable to discover patient reports. Use refresh to retry.';});
await pollJob();
