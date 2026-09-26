#!/usr/bin/env python3

from pathlib import Path
import json, hashlib, subprocess, time, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
# ==================== COMMON CONFIGURATION ====================
DATA = Path('/Users/L/Desktop/ESCC/shuffle')
TARGET_FILE = DATA / 'overlap_SNP/SNP.csv'
BACKGROUND_FILE = DATA / 'background SNP/phs000361.pha003107.txt'
PEAK_DIR = DATA / 'peak'
CHROMOSOME_SIZES = PEAK_DIR / 'hg19.chrom.sizes'
ROOT = Path('/Users/L/Desktop/ESCC/shuffle/fig1fg_complete_pipeline/output')
B = 10000
SEED = 20260924
# Re-running replaces files only inside the output directory above.
(ROOT / 'scripts').mkdir(parents=True, exist_ok=True)
for d in ['inputs','results','figures']: (ROOT/d).mkdir(exist_ok=True)
# Compile the included Figure 1F engine using the installed Apple compiler.
CPP_SOURCE = r"""
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <string>
using namespace std;
int main(int argc,char**argv){
 if(argc!=7){cerr<<"sizes snps peaks B seed output\n";return 1;}
 map<string,long long> sizes; string c;long long a,b;ifstream sz(argv[1]);while(sz>>c>>a)sizes[c]=a;
 map<string,vector<long long>> snps;ifstream ss(argv[2]);while(ss>>c>>a>>b)snps[c].push_back(a);
 for(auto &kv:snps)sort(kv.second.begin(),kv.second.end());
 map<string,vector<long long>> lens;ifstream pp(argv[3]);while(pp>>c>>a>>b){if(a<0||b<=a||!sizes.count(c)||b>sizes[c])return 2;lens[c].push_back(b-a);}
 mt19937_64 rng(stoull(argv[5]));ofstream out(argv[6]);int B=stoi(argv[4]);
 for(int k=0;k<B;k++){int count=0;
  for(auto &kv:snps){auto &pos=kv.second;vector<unsigned char> hit(pos.size(),0);long long L=sizes[kv.first];
   for(long long len:lens[kv.first]){long long start=uniform_int_distribution<long long>(0,L-len)(rng);long long end=start+len;
    if(start>pos.back()||end<=pos.front())continue;
    auto it=lower_bound(pos.begin(),pos.end(),start);for(;it!=pos.end()&&*it<end;++it)hit[it-pos.begin()]=1;
   }
   for(auto h:hit)count+=h;
  }out<<count<<'\n';
 }
}
"""
cpp_file = ROOT / 'scripts/shuffle_peaks.cpp'
cpp_file.write_text(CPP_SOURCE)
subprocess.run(['/usr/bin/clang++', '-O3', '-std=c++17', str(cpp_file),
                '-o', str(ROOT / 'scripts/shuffle_peaks')], check=True)

# ==================== SHARED INPUT PREPARATION ====================
# SNP.csv pvalue < 1e-5 selects 167 targets.
# Coordinates are already hg19 BED 0-based starts.
sizes=dict(pd.read_csv(CHROMOSOME_SIZES,sep='\t',header=None).itertuples(index=False,name=None))
b=pd.read_csv(BACKGROUND_FILE,sep='\t',header=None,names=['id','chr','pos','pvalue','unknown5','unknown6','unknown7','unknown8','unknown9','unknown10'])
b['chrom']='chr'+b.chr.astype(str)
raw_s=pd.read_csv(TARGET_FILE)
s=raw_s.loc[raw_s.pvalue<1e-5].copy().reset_index(drop=True)
s['chrom']='chr'+s.chr.astype(str)
assert len(s)==167 and (s.end==s.pos+1).all()
assert s.id.is_unique and not s.duplicated(['chrom','pos']).any()
s.to_csv(ROOT/'results/target_snps_167.csv',index=False)
s[['chrom','pos','end']].to_csv(ROOT/'inputs/target.bed',sep='\t',header=False,index=False)


fixed={'KYSE140_HiChIP_anchors':23,'KYSE180_HiChIP_anchors':5,'KYSE70_HiChIP_anchors':21,'TT_HiChIP_anchors':23,'KYSE140_H3K27ac':7,'KYSE180_H3K27ac':2,'KYSE70_H3K27ac':11,'TT_H3K27ac':2,'TE5_H3K27ac':13,'KYSE70_H3K4me1':4,'TE5_H3K4me1':2,'KYSE70_H3K4me3':4,'TE5_H3K4me3':4,'KYSE140_SOX2':3,'TE5_SOX2':0,'promoter_new1':5,'genebody_new1':22,'other_region':3,'KYSE180_H3K27ac_deduplicated_sensitivity':2}
audit={'seed':SEED,'B':B,'target_used':len(s),'target_source':'overlap_SNP/SNP.csv, pvalue < 1e-5','comparison':'Fixed historical Figure 1C count thresholds against a 167-SNP null, NOT an enrichment test of the actual 167-SNP overlaps','background_raw':len(b),'background_p_gt_threshold':int((b.pvalue>1e-5).sum()),'background_p_equal_threshold':int((b.pvalue==1e-5).sum()),'coordinate_convention':'Supplied hg19 BED 0-based starts; no additional subtraction.'}
b=b.loc[(b.pvalue>1e-5)&b.chrom.isin(sizes)&(b.pos>=0)].copy();b=b.loc[b.pos< b.chrom.map(sizes)].copy(); before=len(b)
keys=set(zip(s.chrom,s.pos)); keep=~b.id.isin(s.id)&~pd.Series(list(zip(b.chrom,b.pos)),index=b.index).isin(keys);b=b.loc[keep];audit['background_target_removed']=before-len(b)
b=b.drop_duplicates(['chrom','pos']).drop_duplicates('id').reset_index(drop=True);audit['background_eligible']=len(b)
b[['id','chrom','pos','pvalue']].to_csv(ROOT/'results/background_eligible.tsv.gz',sep='\t',index=False)
# ==================== FIGURE 1G: RANDOMIZE SNPs (PART 1) ====================
# Real peaks stay fixed. For each of 10,000 replicates, draw 167 SNPs from
# P > 1e-5 background, excluding all target rsIDs/positions. Sampling is
# WITHOUT replacement WITH exact chromosome matching. The same replicate
# is reused across all features. No MAF/TSS/LD matching is performed.
rng=np.random.default_rng(SEED); samples=[]; chromosome=[]
for c,sc in s.groupby('chrom',sort=True):
 ix=np.flatnonzero(b.chrom.to_numpy()==c); n=len(sc); assert len(ix)>=n
 samples.append(np.stack([rng.choice(ix,n,replace=False) for _ in range(B)]))
 chromosome.append({'chrom':c,'target_N':n,'background_N':len(ix)})
draws=np.concatenate(samples,axis=1); assert draws.shape==(B,len(s));np.savez_compressed(ROOT/'results/background_draw_indices.npz',indices=draws)
pd.DataFrame(chromosome).to_csv(ROOT/'results/chromosome_matching.csv',index=False)

def intervals_mask(points,peaks,offset=0):
 ans=np.zeros(len(points),bool)
 for c,g in points.groupby('chrom',sort=False):
  p=peaks.loc[peaks.chrom==c].sort_values('start');
  if len(p)==0:continue
  st=p.start.to_numpy(); ends=np.maximum.accumulate(p.end.to_numpy());q=g.pos.to_numpy()+offset
  i=np.searchsorted(st,q,side='right')-1;ans[g.index.to_numpy()]=(i>=0)&(q<ends[np.maximum(i,0)])
 return ans
# ==================== SHARED FEATURE PREPARATION ====================
# BED: retain rows/lengths, including duplicates. HiChIP BEDPE: use both
# anchors and remove identical anchors; never use the whole loop span.
# Keep this feature ordering: the Figure 1F seeds depend on its indices.
# 15 chromatin tracks + 3 annotations form the 18-test BH family.
# Extra KYSE180 deduplication sensitivity is excluded from that family.
features=[]; peak_audits=[]
for path in sorted((PEAK_DIR).glob('*.bed'))+sorted((PEAK_DIR).glob('*.bedpe')):
 if '_chr' in path.stem or path.name=='gene.bed':continue # chromosome-only examples are subsets, not additional features
 d=pd.read_csv(path,sep='\t',header=None);raw=len(d)
 if path.suffix=='.bedpe':
  a=d.iloc[:,:3].copy();z=d.iloc[:,3:6].copy();a.columns=z.columns=['chrom','start','end'];p=pd.concat([a,z],ignore_index=True);raw_intervals=len(p);p=p.drop_duplicates().reset_index(drop=True);name=path.stem+'_HiChIP_anchors';kind='HiChIP'
 else:
  p=d.iloc[:,:3].copy();p.columns=['chrom','start','end'];raw_intervals=len(p);name=path.stem;kind='annotation' if name in ['gene','genebody_new1','promoter_new1'] else 'ChIP'
 valid=p.chrom.isin(sizes)&(p.start>=0)&(p.end>p.start); valid &=p.end<=p.chrom.map(sizes)
 if not valid.all():
  p.loc[~valid].to_csv(ROOT/'results'/f'{name}_invalid_intervals.csv',index=False)
  raise ValueError(f'{name}: invalid intervals found; inspect before proceeding')
 peak_audits.append({'feature':name,'kind':kind,'input_file':path.name,'raw_rows':raw,'raw_intervals':raw_intervals,'tested_intervals':len(p),'exact_duplicate_intervals_in_test':int(p.duplicated().sum())})
 features.append((name,kind,p.reset_index(drop=True)))
 if name=='KYSE180_H3K27ac':features.append((name+'_deduplicated_sensitivity','sensitivity',p.drop_duplicates().reset_index(drop=True)))
 # Confirm KYSE70 chromosome-only BEDPE files contain only records already in full file.
 if path.name=='KYSE70.bedpe':
  full=pd.MultiIndex.from_frame(d.iloc[:,:6])
  for sub in (PEAK_DIR).glob('KYSE70_H3K27ac_chr*.bedpe'):
   t=pd.read_csv(sub,sep='\t',header=None);assert pd.MultiIndex.from_frame(t.iloc[:,:6]).isin(full).all()
# Other regions: complement of promoter/gene-body interval union.
parts=[p for name,kind,p in features if name in ['promoter_new1','genebody_new1']]
other=pd.concat(parts,ignore_index=True)
features.append(('other_region','annotation',other))
peak_audits.append({'feature':'other_region','kind':'annotation','input_file':'complement of promoter_new1 + genebody_new1','tested_intervals':len(other)})
pd.DataFrame(peak_audits).to_csv(ROOT/'results/feature_input_audit.csv',index=False)
obs=pd.DataFrame({'id':s.id,'chrom':s.chrom,'start':s.pos,'end':s.end,'GWAS_P':s.pvalue})
rows=[];tasks=[];sensitivity=[]
for i,(name,kind,p) in enumerate(features):
 path=ROOT/'inputs'/f'{name}.bed';p.to_csv(path,sep='\t',index=False,header=False)
 mask=intervals_mask(s,p);bm=intervals_mask(b,p)
 if name=='other_region':mask=~mask;bm=~bm
 obs[name]=mask.astype(int)
 # FIGURE 1G (PART 2): count distinct sampled SNPs hitting >=1 real interval.
 nullg=bm[draws].sum(axis=1);np.savetxt(ROOT/'results'/f'{name}.G.null.tsv',nullg,fmt='%d')
 # boundary +/-1 sensitivity for observed counts and background overlap fractions
 sensitivity.append({'feature':name,'observed_pos_minus1':(len(s)-int(intervals_mask(s,p,-1).sum()) if name=="other_region" else int(intervals_mask(s,p,-1).sum())),'observed_as_supplied':int(mask.sum()),'observed_pos_plus1':(len(s)-int(intervals_mask(s,p,1).sum()) if name=="other_region" else int(intervals_mask(s,p,1).sum())),'background_hits_minus1':(len(b)-int(intervals_mask(b,p,-1).sum()) if name=="other_region" else int(intervals_mask(b,p,-1).sum())),'background_hits_as_supplied':int(bm.sum()),'background_hits_plus1':(len(b)-int(intervals_mask(b,p,1).sum()) if name=="other_region" else int(intervals_mask(b,p,1).sum()))})
 rows.append({'feature':name,'kind':kind,'N':len(s),'observed':fixed[name],'actual_overlap_167_audit_only':int(mask.sum()),'fixed_threshold_matches_actual':fixed[name]==int(mask.sum()),'intervals_shuffled':len(p),'background_overlap':int(bm.sum())})
 tasks.append((name,path,SEED+1000+i))
obs.to_csv(ROOT/'results/target_feature_overlap_matrix.csv',index=False)
pd.DataFrame(sensitivity).to_csv(ROOT/'results/coordinate_sensitivity.csv',index=False)
# Validate interval lookup independently using brute-force for all targets in every feature.
for name,kind,p in features:
 brute=[]
 for q in s.itertuples():brute.append(bool(((p.chrom==q.chrom)&(p.start<=q.pos)&(p.end>q.pos)).any()))
 if name=='other_region':brute=np.logical_not(brute)
 assert np.array_equal(brute,obs[name].to_numpy().astype(bool)),name
# Confirm no overlapping target IDs or coordinates in eligible background, fixed chromosome counts.
assert not b.id.isin(s.id).any();assert not set(zip(b.chrom,b.pos))&keys
assert all(len(set(row))==len(s) for row in draws)
assert all(np.all((b.chrom.to_numpy()[draws]==r['chrom']).sum(axis=1)==r['target_N']) for r in chromosome)
checks={'brute_force_overlap_all_targets_all_features':'passed','sampling_no_replacement_all_10000':'passed','sampling_chromosome_counts_all_10000':'passed','target_exclusion':'passed','KYSE70_chr_subsets':'passed'}
(ROOT/'results/validation.json').write_text(json.dumps(checks,indent=2))
(ROOT/'results/input_audit.json').write_text(json.dumps(audit,indent=2))
# fingerprint inputs, including all source features
manifest={str(p.relative_to(DATA)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(list((DATA/'overlap_SNP').glob('*'))+list((PEAK_DIR).glob('*'))+list((DATA/'background SNP').glob('*'))) if p.is_file() and ROOT not in p.parents and p.name!='.DS_Store'}
(ROOT/'results/source_sha256.json').write_text(json.dumps(manifest,indent=2))
# ==================== FIGURE 1F: RANDOMIZE PEAKS ====================
# Keep target SNP positions fixed. Relocate each peak independently within
# its original chromosome. Uniform integer starts range from 0 to L-length,
# inclusive. Preserve every peak length and count; shuffled peaks may overlap.
# No blacklist/gap mask is used. Each SNP contributes at most 1 overlap.
# C++ accelerates this exact simulation; bedtools is NOT required.
# Chromosomes without target SNPs can be skipped because they contribute 0.
# For 'other_region', count the complement of the shuffled promoter/body union.
def run(t):
 name,path,seed=t;out=ROOT/'results'/f'{name}.F.null.tsv';start=time.time()
 subprocess.run([str(ROOT/'scripts/shuffle_peaks'),str(CHROMOSOME_SIZES),str(ROOT/'inputs/target.bed'),str(path),str(B),str(seed),str(out)],check=True)
 if name=='other_region':
  counts=np.loadtxt(out,dtype=int);np.savetxt(out,len(s)-counts,fmt='%d')
 print(name,'finished',round(time.time()-start,1),'seconds',flush=True)
with ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(run,tasks))
# ==================== FIGURES 1F AND 1G: STATISTICS ====================
# Tail P = (1 + number of null counts >= fixed Fig1C count) / (B + 1).
# Difference = fixed count - median(null), signed, NOT absolute.
# FE = fixed count / median(null); zero denominators remain inf/NaN.
# BH correction is separate for F and G across 18 primary features.
def bh(p):
 p=np.asarray(p);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
for row in rows:
 for method in ['F','G']:
  null=np.loadtxt(ROOT/'results'/f'{row["feature"]}.{method}.null.tsv',dtype=int);o=row['observed'];med=float(np.median(null));k=int((null>=o).sum());assert len(null)==B and np.all((null>=0)&(null<=len(s)))
  for key,val in {'null_mean':float(null.mean()),'null_median':med,'null_q025':float(np.quantile(null,.025)),'null_q975':float(np.quantile(null,.975)),'exceedances':k,'P':(1+k)/(B+1),'FE':o/med if med else (float('inf') if o else float('nan')),'difference':o-med}.items():row[f'{method}_{key}']=val
res=pd.DataFrame(rows)
for m in ['F','G']:
 main=res.kind!='sensitivity';res.loc[main,f'{m}_BH_q']=bh(res.loc[main,f'{m}_P']);res[f'{m}_stars_P']=res[f'{m}_P'].map(lambda p:'***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns')
res['inference_status']='historical_count_threshold_vs_167SNP_null_not_actual167_enrichment'
res.to_csv(ROOT/'results/enrichment_summary.csv',index=False)
print(res[['feature','observed','F_P','F_BH_q','F_FE','G_P','G_BH_q','G_FE']].to_string(index=False),flush=True)
print(json.dumps(audit,indent=2),flush=True)

# ==================== FIGURE 1F / FIGURE 1G: DRAW BOTH PANELS ====================
# R code below uses the original Fig1.R style, without reading that Chinese-
# directory file. No titles; 5 x 4.3 inches; PNG 300 dpi; PDF also exported.
# Effect Size = signed difference, not FE. P > 0.05 maps to white for display.
# The R loop method="F" draws Figure 1F; method="G" draws Figure 1G.
R_SOURCE = r"""
# Historical Fig1C count thresholds against the SNP.csv-derived 167-SNP null.
# This is NOT a test of actual 167-SNP observed enrichment.
# Run from any working directory: Rscript /absolute/path/to/this/script.R
library(ggplot2)
args <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1]))
root <- dirname(dirname(script))
out <- file.path(root, "figures", "Fig1R_style")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
s <- read.csv(file.path(root, "results", "enrichment_summary.csv"))
s <- s[s$kind %in% c("ChIP", "HiChIP"), ]
s$cellline <- sub("_.*", "", s$feature)
s$factor <- sub("^[^_]+_", "", s$feature)
s$factor[s$factor == "HiChIP_anchors"] <- "H3K27ac anchors"
grid <- expand.grid(cellline=c("KYSE140","KYSE180","KYSE70","TT","TE5"),
                    factor=c("SOX2","H3K4me3","H3K4me1","H3K27ac","H3K27ac anchors"),
                    stringsAsFactors=FALSE)
key <- function(d) paste(d$cellline,d$factor,sep="|")
data3 <- grid[!key(grid) %in% key(s), ]
write.csv(data3,file.path(out,"heat3.csv"),row.names=FALSE)
for (method in c("F","G")) {
  data2 <- data.frame(cellline=s$cellline,factor=s$factor,
                      P_raw=s[[paste0(method,"_P")]],
                      FE_raw=s[[paste0(method,"_FE")]],
                      difference=s[[paste0(method,"_difference")]],
                      observed=s$observed,
                      null_median=s[[paste0(method,"_null_median")]])
  # Map nonsignificant P to the white endpoint; preserve original P separately.
  data2$P <- pmin(data2$P_raw,0.05)
  # Original manuscript effect-size definition: fixed count minus null median.
  # All 15 plotted chromatin features have nonnegative differences in this run.
  data2$effectsize <- data2$difference
  stopifnot(all(is.finite(data2$effectsize)),all(data2$effectsize>=0),
            all(data2$effectsize<=50),nrow(data2)==15,nrow(data3)==10)
  write.csv(data2,file.path(out,paste0("plot_data_",method,".csv")),row.names=FALSE)
p1 <- ggplot(data2, aes(cellline,factor))+
  geom_point(aes(size = effectsize,fill=P),shape=21,stroke=0.6) +
  scale_fill_gradient(name = 'P',
                      limit = c(-0.00001,0.050001),
                      breaks = c(0.001,0.01,0.05),
                      low = 'deeppink',
                      high = 'white')+
  scale_size_continuous(name = 'Effect Size',
                        limit = c(-0.001,50),
                        breaks = c(0,5,10,15),
                        range = c(3,17))+
  
  scale_y_discrete(limits=as.character(c("SOX2","H3K4me3","H3K4me1","H3K27ac","H3K27ac anchors")))+ 
  scale_x_discrete(limits=as.character(c("KYSE140","KYSE180","KYSE70","TT","TE5")))+ 
  theme_bw()+
  xlab(NULL) + 
  ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1,colour = "black"),  
        axis.text.y = element_text(size = 14, angle = 45,colour = "black"),  
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12))  
p1
p2<-p1+geom_point(data=data3,
                  mapping =aes(cellline,factor),
                  shape=4,
                  stroke=1,
                  size=8,
                  color='black')
p2

  stem <- if (method=="F") "Fig1F-shuffle" else "Fig1G-shuffle_snp"
  ggsave(filename=file.path(out,paste0(stem,".png")),plot=p2,
         width=5,height=4.3,dpi=300,bg="white")
  ggsave(filename=file.path(out,paste0(stem,".pdf")),plot=p2,
         width=5,height=4.3,bg="white")
}
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
"""
r_file = ROOT / 'scripts/plot_fig1FG.R'
r_file.write_text(R_SOURCE)
subprocess.run(['/usr/local/bin/Rscript', str(r_file)], check=True)
print('Finished. Results:', ROOT / 'results', flush=True)
print('Figure 1F:', ROOT / 'figures/Fig1R_style/Fig1F-shuffle.pdf', flush=True)
print('Figure 1G:', ROOT / 'figures/Fig1R_style/Fig1G-shuffle_snp.pdf', flush=True)
