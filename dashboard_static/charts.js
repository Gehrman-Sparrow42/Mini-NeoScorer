import {$, fmt} from './data.js';
const colors = {0:'#b1a0c7',1:'#34d399',2:'#e9b65c',3:'#73839c'};
let scatter, foreign;
const quadrantPlugin = {
  id:'quadrants', beforeDatasetsDraw(chart) {
    const {ctx, chartArea:a, scales:{x,y}} = chart;
    if(!a) return;
    const gx = Math.max(a.left, Math.min(a.right, x.getPixelForValue(Math.log2(5))))
    const gy = Math.max(a.top, Math.min(a.bottom, y.getPixelForValue(-Math.log10(50))));
    ctx.save(); ctx.fillStyle='#34d39908'; ctx.fillRect(gx,a.top,a.right-gx,gy-a.top);
    ctx.strokeStyle='#57716b88';ctx.setLineDash([4,5]);ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(gx,a.top);ctx.lineTo(gx,a.bottom);ctx.moveTo(a.left,gy);ctx.lineTo(a.right,gy);ctx.stroke();
    ctx.setLineDash([]);ctx.font='9px Segoe UI';ctx.fillStyle='#6e9d90';ctx.textAlign='right';ctx.fillText('HIGH EXPRESSION / HIGH AFFINITY',a.right-8,a.top+14);ctx.restore();
  }
};
function baseOptions(){return {responsive:true,maintainAspectRatio:false,animation:false,plugins:{legend:{display:false},tooltip:{backgroundColor:'#0b1626',borderColor:'#3a5365',borderWidth:1,padding:12,titleColor:'#eef6ff',bodyColor:'#b7c7db'}},scales:{x:{grid:{color:'#26354966'},ticks:{color:'#8e9fb5',font:{size:10}},border:{display:false}},y:{grid:{color:'#26354966'},ticks:{color:'#8e9fb5',font:{size:10}},border:{display:false}}}};}
export function clearCharts(){scatter?.destroy(); foreign?.destroy(); scatter=foreign=null;}
export function renderCharts(rows, select) {
  clearCharts();
  if(!window.Chart) throw new Error('The local chart library could not be loaded. Refresh the page.');
  window.Chart.defaults.font.family="'Segoe UI', Arial, sans-serif";
  const valid=rows.filter(r=>r.TPM>0 && r.MT_IC50_nM>0);
  $('plot-count').textContent=`${valid.length} candidates`;
  const omitted=rows.length-valid.length;
  $('scatter-note').textContent=`Reference gates: TPM 5 and IC50 50 nM. Upper-right: higher expression and affinity.${omitted ? ` ${omitted} rows with non-positive/missing values omitted.` : ''}`;
  const options=baseOptions();
  options.scales.x.title={display:true,text:'RNA expression · log₂(TPM)',color:'#9eafc4',font:{size:10}};
  options.scales.y.title={display:true,text:'HLA affinity · −log₁₀(IC50)',color:'#9eafc4',font:{size:10}};
  options.scales.x.suggestedMin=0;options.scales.x.suggestedMax=8;
  options.scales.y.suggestedMin=-5;options.scales.y.suggestedMax=0;
  options.plugins.tooltip.callbacks={title:items=>items[0]?.raw.row.Gene || '',label:context=>{const r=context.raw.row;return [r.Neopeptide_9mer,`Score ${fmt(r.Score,1)} · VAF ${fmt(r.AF,3)}`,`TPM ${fmt(r.TPM)} · IC50 ${fmt(r.MT_IC50_nM)} nM`];}};
  options.onClick=(_,items)=>{if(items.length){const p=items[0];select(scatter.data.datasets[p.datasetIndex].data[p.index].row);}};
  scatter=new window.Chart($('scatter'),{type:'bubble',plugins:[quadrantPlugin],data:{datasets:[3,2,1,0].map(t=>({label:`Tier ${t}`,data:valid.filter(r=>r.tier===t).map(row=>({x:Math.log2(row.TPM),y:-Math.log10(row.MT_IC50_nM),r:3+Math.sqrt(Math.min(1,Math.max(0,row.AF || 0)))*10,row})),backgroundColor:colors[t]+'80',borderColor:colors[t],borderWidth:1,hoverBorderWidth:2}))},options});
  const top=rows.filter(r=>r.tier===1).sort((a,b)=>(b.Score||0)-(a.Score||0)).slice(0,8);
  $('foreign-empty').hidden=top.length>0;
  const foreignOptions=baseOptions();
  foreignOptions.scales.y.type='logarithmic';foreignOptions.scales.y.title={display:true,text:'IC50 (nM) · logarithmic scale',color:'#9eafc4',font:{size:10}};
  foreignOptions.scales.x.grid.display=false;foreignOptions.scales.x.ticks.maxRotation=35;foreignOptions.scales.x.ticks.minRotation=0;
  foreignOptions.plugins.tooltip.callbacks={afterTitle:items=>top[items[0].dataIndex]?.Neopeptide_9mer,label:ctx=>`${ctx.dataset.label}: ${fmt(ctx.raw)} nM`,afterBody:items=>{const r=top[items[0].dataIndex];return `Agretopicity: ${fmt(r.Agretopicity)}× · Mutation position ${fmt(r.MutPos,0)}`;}};
  foreignOptions.onClick=(_,items)=>{if(items.length)select(top[items[0].index]);};
  foreign=new window.Chart($('foreign'),{type:'bar',data:{labels:top.map(r=>r.Gene),datasets:[{label:'Mutant',data:top.map(r=>r.MT_IC50_nM>0?r.MT_IC50_nM:null),backgroundColor:'#34d399',borderRadius:3,maxBarThickness:19},{label:'Wild type',data:top.map(r=>r.WT_IC50_nM>0?r.WT_IC50_nM:null),backgroundColor:'#729fe3aa',borderRadius:3,maxBarThickness:19}]},options:foreignOptions});
}
