#!/usr/bin/env python3
"""Проверка трёх генов на случай для учебного набора, Python 3."""
import argparse
import hashlib
import io
import json
import os
import re
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

SESSION = requests.Session()
SESSION.mount('https://', HTTPAdapter(max_retries=Retry(
    total=3, backoff_factor=1, status_forcelist=[429,500,502,503,504],
    allowed_methods=['GET','POST'], respect_retry_after_header=True)))
SESSION.headers['User-Agent'] = 'EducationalVariantReview/1.0'


def read(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')


def cached_json(cache, label, url, params=None, payload=None):
    key = hashlib.sha256(json.dumps([label, url, params, payload], sort_keys=True).encode()).hexdigest()
    path = cache / (label + '_' + key[:20] + '.json')
    if path.exists():
        return json.loads(path.read_text())
    try:
        if payload is None:
            response = SESSION.get(url, params=params, timeout=(20,90))
        else:
            response = SESSION.post(url, json=payload, timeout=(20,90))
        if not response.ok:
            return {'status':'request_failed', 'http_status':response.status_code}
        data = response.json()
        result = {'status':'retrieved', 'data':data}
        dump(path, result)
        return result
    except (requests.RequestException, ValueError):
        return {'status':'request_failed'}


def download(cache, name, url):
    path = cache / name
    if path.exists():
        return path
    try:
        with SESSION.get(url, stream=True, timeout=(20,120)) as response:
            response.raise_for_status()
            tmp = path.with_suffix('.tmp')
            with tmp.open('wb') as handle:
                for chunk in response.iter_content(1024*1024):
                    if chunk:
                        handle.write(chunk)
            tmp.replace(path)
        return path
    except requests.RequestException:
        return None


def normalize(db, number):
    if db.upper() == 'OMIM' and re.fullmatch(r'\d{6}', number):
        return 'OMIM:' + number
    if db.upper() in ('ORPHA','ORPHANET'):
        return 'ORPHA:' + number.split(':')[-1]
    return number if number.startswith(db+':') else db+':'+number


def read_orphanet(cache):
    links, disorders, gene_mim = defaultdict(set), {}, defaultdict(set)
    status = {}
    for product in (1,6):
        name = f'en_product{product}.xml'
        path = download(cache, name, 'https://www.orphadata.com/data/xml/'+name)
        if path is None:
            status[str(product)] = 'download_failed'
            continue
        try:
            count = 0
            for _, node in ET.iterparse(path, events=['end']):
                if node.tag != 'Disorder':
                    continue
                code = node.findtext('OrphaCode','')
                did = 'ORPHA:' + code
                ids = {did}
                exact_ids = {did}
                for x in node.findall('./ExternalReferenceList/ExternalReference'):
                    source, ref = x.findtext('Source',''), x.findtext('Reference','')
                    if source and ref:
                        identifier = normalize(source, ref)
                        ids.add(identifier)
                        relation = x.findtext('DisorderMappingRelation/Name','')
                        if relation.startswith('E ('):
                            exact_ids.add(identifier)
                if product == 1:
                    disorders[did] = {'name':node.findtext('Name',''), 'ids':sorted(ids), 'exact_ids':sorted(exact_ids)}
                    for identifier in ids:
                        links[identifier].add(did)
                else:
                    item = disorders.setdefault(did, {'name':node.findtext('Name',''), 'ids':[did]})
                    associations = item.setdefault('genes', [])
                    for assoc in node.findall('./DisorderGeneAssociationList/DisorderGeneAssociation'):
                        gene = assoc.find('Gene')
                        if gene is None:
                            continue
                        symbol = gene.findtext('Symbol','')
                        associations.append({
                            'gene':symbol,
                            'type':assoc.findtext('DisorderGeneAssociationType/Name',''),
                            'status':assoc.findtext('DisorderGeneAssociationStatus/Name',''),
                        })
                        for x in gene.findall('./ExternalReferenceList/ExternalReference'):
                            if x.findtext('Source','').upper()=='OMIM':
                                gene_mim[symbol].add(x.findtext('Reference',''))
                count += 1
                node.clear()
            status[str(product)] = 'retrieved' if count else 'unexpected_schema'
        except ET.ParseError:
            status[str(product)] = 'parse_failed'
    return links, disorders, gene_mim, status


def read_clingen(cache):
    url = 'https://search.clinicalgenome.org/kb/gene-validity/download'
    path = download(cache, 'clingen_validity.csv', url)
    if path is None:
        return None
    text = path.read_text(encoding='utf-8-sig')
    lines = text.splitlines()
    for i,line in enumerate(lines):
        if line.startswith('"GENE SYMBOL"') or line.startswith('GENE SYMBOL,'):
            table = pd.read_csv(io.StringIO('\n'.join(lines[i:])), dtype=str).fillna('')
            if 'GENE SYMBOL' in table:
                return table[~table['GENE SYMBOL'].str.startswith('+')]
    return None


def variant_notations(tx):
    terms = set()
    for value in (tx.get('hgvsc',''), tx.get('hgvsp','')):
        if ':' in value:
            short = value.split(':',1)[1]
            terms.add(short)
            if short.startswith(('p.','c.')):
                terms.add(short[2:])
    return sorted(terms)


def publications(cache, gene, tx):
    terms = variant_notations(tx)
    if not terms:
        return {'status':'no_variant_notation', 'hits':[]}
    query = gene + ' AND (' + ' OR '.join('"'+x+'"' for x in terms) + ')'
    response = cached_json(cache, 'literature',
        'https://www.ebi.ac.uk/europepmc/webservices/rest/search',
        params={'query':query,'format':'json','resultType':'core','pageSize':5})
    if response['status']!='retrieved':
        return response
    hits = []
    for hit in response['data'].get('resultList',{}).get('result',[]):
        text = hit.get('title','')+' '+hit.get('abstractText','')
        mentions = [term for term in terms if term.lower() in text.lower()]
        # Совпадение поискового запроса не равно доказательству патогенности.
        hits.append({'title':hit.get('title',''), 'pmid':hit.get('pmid',''),
            'pmcid':hit.get('pmcid',''), 'doi':hit.get('doi',''),
            'exact_notation_in_title_or_abstract':mentions,
            'url':'https://europepmc.org/article/'+hit.get('source','MED')+'/'+hit.get('id','')})
    return {'status':'retrieved', 'query':query, 'hits':hits,
            'interpretation':'search_hits_require_reading; not_causality_confirmation'}


def gnomad(cache, variant):
    chrom,pos,ref,alt = variant.split(':')
    query = '''query($id: String!) { variant(variantId: $id, dataset: gnomad_r4) {
      variant_id exome { ac an af populations { id ac an } }
      genome { ac an af populations { id ac an } }
    } }'''
    result = cached_json(cache, 'gnomad_r4',
        'https://gnomad.broadinstitute.org/api',
        payload={'query':query,'variables':{'id':f'{chrom.removeprefix("chr")}-{pos}-{ref}-{alt}'}})
    if result['status']=='retrieved':
        data = result['data']
        if data.get('errors'):
            return {'status':'variant_not_returned' if all('not found' in x.get('message','').lower() for x in data['errors']) else 'api_error', 'details':data['errors']}
        if not (data.get('data') or {}).get('variant'):
            return {'status':'no_record_returned_not_zero_frequency'}
    return result


def omim(cache, did, gene, tx, gene_mim):
    key = os.environ.get('OMIM_API_KEY','').strip()
    if not key:
        return {'status':'not_checked_API_key_required'}
    numbers = set(gene_mim.get(gene,set()))
    if did.startswith('OMIM:'):
        numbers.add(did.split(':')[1])
    numbers = sorted(x for x in numbers if re.fullmatch(r'\d{6}',x))
    if not numbers:
        return {'status':'not_checked_no_OMIM_entry_mapping'}
    result = cached_json(cache, 'omim', 'https://api.omim.org/api/entry', params={
        'mimNumber':','.join(numbers), 'include':'text,allelicVariantList,geneMap',
        'format':'json', 'apiKey':key})
    if result['status'] != 'retrieved':
        return result
    data = result['data']
    entries = data.get('omim',{}).get('entryList',[])
    if not entries:
        return {'status':'unexpected_response_or_access_error'}
    rawtext = json.dumps(entries, ensure_ascii=False)
    return {'status':'retrieved', 'entries':entries,
        'lexical_variant_mentions':[t for t in variant_notations(tx) if t in rawtext],
        'interpretation':'lexical_match_requires_manual_context_review'}



def live_clinvar(cache, variation_id, variant_id):
    if not re.fullmatch(r"\d+", variation_id):
        return {"status":"no_VariationID_in_previous_VCF", "records":[]}
    path=cache/("clinvar_live_"+variation_id+".xml")
    try:
        if path.exists():
            root=ET.fromstring(path.read_bytes())
        else:
            time.sleep(0.4)
            response=SESSION.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={"db":"clinvar","rettype":"vcv","is_variationid":"1","id":variation_id},timeout=(20,90))
            response.raise_for_status(); root=ET.fromstring(response.content)
            for node in root.iter():node.tag=node.tag.split("}")[-1]
            path.write_bytes(ET.tostring(root,encoding="utf-8"))
        archives=[root] if root.tag=="VariationArchive" else root.findall(".//VariationArchive")
        archive=next((x for x in archives if x.get("VariationID")==variation_id),None)
        if archive is None:return {"status":"unexpected_XML", "records":[]}
        keys=set()
        for x in archive.findall("./ClassifiedRecord/SimpleAllele/.//SequenceLocation"):
            if x.get("Assembly")=="GRCh38":
                keys.add(":".join(["chr"+x.get("Chr","").removeprefix("chr"),x.get("positionVCF",""),x.get("referenceAlleleVCF",""),x.get("alternateAlleleVCF","")]))
        if variant_id not in keys:return {"status":"allele_not_verified", "records":[]}
        maps=defaultdict(lambda:defaultdict(set))
        for ts in archive.findall(".//TraitSet"):
            if not ts.get("ID"):continue
            for trait in ts.findall("./Trait"):
                ids={normalize(x.get("DB",""),x.get("ID","")) for x in trait.findall(".//XRef") if x.get("DB") and x.get("ID")}
                ids={x for x in ids if not (x.startswith("OMIM:") and not re.fullmatch(r"OMIM:\d{6}",x))}
                for x in ids:maps[ts.get("ID")][x].update(ids)
        records=[]
        for r in archive.findall("./ClassifiedRecord/RCVList/RCVAccession"):
            cc=r.find("./ClassifiedConditionList")
            if cc is None:continue
            conditions=[]
            for x in cc.findall("./ClassifiedCondition"):
                ident=normalize(x.get("DB",""),x.get("ID",""))
                conditions.append({"name":x.text or "", "ids":sorted({ident}|maps[cc.get("TraitSetID","")].get(ident,set()))})
            gl=r.find("./RCVClassifications/GermlineClassification")
            records.append({"rcv":r.get("Accession",""),"version":r.get("Version",""),"conditions":conditions,
                "multiple_conditions":len(conditions)>1,
                "classification":gl.findtext("Description","") if gl is not None else "",
                "review_status":gl.findtext("ReviewStatus","") if gl is not None else ""})
        return {"status":"retrieved_verified_allele", "records":records}
    except (requests.RequestException,ET.ParseError):
        return {"status":"request_or_parse_failed", "records":[]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/content'))
    parser.add_argument('--local-only', action='store_true', help='Проверка отбора без сети')
    args=parser.parse_args(); root=args.root
    out=root/'stage5_results'; out.mkdir(exist_ok=True)
    cache=out/'cache'; cache.mkdir(exist_ok=True)
    candidates=read(root/'stage4_results/phenotype_candidates_with_ClinVar_RCV.tsv')
    transcripts=read(root/'stage2_results/cases_annotated.tsv')
    evidence=read(root/'stage3_results/variant_evidence.tsv').set_index('variant_id').to_dict('index')
    rcv=json.loads((root/'stage4_results/clinvar_rcv_records.json').read_text())
    txgroups={key:group for key,group in transcripts.groupby(['case_id','variant_id'])}
    # Не приписываем заболевание всем соседним генам в объединённом поле.
    gp=read(root/'stage4_results/genes_to_phenotype.txt')
    gp.columns=gp.columns.str.lstrip('#')
    allowed_pairs=set(zip(gp['gene_symbol'],gp['disease_id']))
    # Важно: кандидат должен действительно затрагивать ЭТОТ ген,
    # а не только соседний ген из объединённого списка.
    selected_rows=[]
    for row in candidates.to_dict('records'):
        cid,vid=row['case_id'],row['variant_id']
        if row['genotype_incomplete'].lower()=='true' or row['allelic_state']=='not_carried':
            continue
        ev=evidence[vid]
        freq=json.loads(ev['gnomad_frequencies_json'])
        overall=[float(v) for k,vs in freq.items() if k in ('gnomad','gnomade','gnomadg') for v in vs]
        # 1% — только исследовательский флаг для отбора, не порог ACMG.
        common=bool(overall and max(overall)>=0.01)
        for gene in row['gene_symbols'].split(';'):
            if (gene,row['disease_id']) not in allowed_pairs:
                continue
            group=txgroups[(cid,vid)]
            group=group[group['gene_symbol']==gene].copy()
            if group.empty: continue
            group['_mane']=group['mane_select'].ne('').astype(int)
            group['_canonical']=group['canonical'].eq('1').astype(int)
            group['_coding']=group['impact'].isin(['HIGH','MODERATE']).astype(int)
            tx=group.sort_values(['_mane','_canonical','_coding'],ascending=False).iloc[0].to_dict()
            relevant=group['impact'].isin(['HIGH','MODERATE']).any()
            records=json.loads(row['clinvar_matching_rcv']) if row['clinvar_matching_rcv'] else []
            single=[x for x in records if not x.get('multiple_conditions')]
            classes=[x.get('classification','').lower() for x in single]
            lp=any(x in ('pathogenic','likely pathogenic','pathogenic/likely pathogenic') for x in classes)
            benign=any(x in ('benign','likely benign','benign/likely benign') for x in classes)
            # Прозрачная эвристика отбора; коэффициенты не обучены и не являются вероятностями.
            tier=0 if (lp and not benign and not common) else 1 if (relevant and not common and not benign) else 2
            selected_rows.append({**row,'gene':gene,'_tx':tx,'_tier':tier,
                'common_AF_flag':common,'coding_or_splice_candidate':bool(relevant),
                'condition_specific_P_LP':lp,'condition_specific_benign':benign,
                'selection_score':float(row['hpo_match_fraction'])+0.15*lp+0.05*bool(relevant)-0.8*common-0.8*benign})
    if not selected_rows: raise ValueError('Нет кандидатов с генами.')
    # Берём лучшую пару «вариант–заболевание» каждого гена, затем три разных гена.
    selected_rows.sort(key=lambda x:(x['case_id'],-x['selection_score'],-float(x['hpo_match_fraction']),int(x['excluded_overlap_count']),x['variant_id']))
    top=[]; seen=defaultdict(set)
    for row in selected_rows:
        cid,gene=row['case_id'],row['gene']
        if gene in seen[cid] or len(seen[cid])>=3:continue
        seen[cid].add(gene); row['gene_rank']=len(seen[cid]);top.append(row)
    missing=[c for c in sorted(set(candidates['case_id'])) if len(seen[c])<3]
    print('Выбрано:',len(top),'пар ген–вариант; случаев:',len(seen),'Меньше 3 генов:',missing,flush=True)
    if args.local_only:
        dump(out/'selection_local_check.json',[{k:v for k,v in r.items() if k!='_tx'} for r in top])
        for r in top:print(r['case_id'],r['gene_rank'],r['gene'],r['variant_id'],r['disease_id'],r['hpo_match_fraction'])
        return
    print('Загружаем официальные справочники Orphanet и ClinGen…',flush=True)
    links,disorders,gene_mim,orpha_status=read_orphanet(cache)
    clingen=read_clingen(cache)
    hpoa=root/'stage4_results/phenotype.hpoa'
    lines=hpoa.read_text().splitlines()
    index=next(i for i,x in enumerate(lines) if x.lstrip('#').startswith(('database_id\t','DatabaseID\t')))
    hpotable=pd.read_csv(io.StringIO('\n'.join([lines[index].lstrip('#')]+[x for x in lines[index+1:] if x and not x.startswith('#')])),sep='\t',dtype=str).fillna('')
    hpotable.columns=hpotable.columns.str.lower()
    hpotable=hpotable.rename(columns={'databaseid':'database_id'})
    modes={'HP:0000006':'AD','HP:0000007':'AR','HP:0001417':'XL'}
    output_rows=[]
    for n,row in enumerate(top,1):
        cid,gene,vid,did=row['case_id'],row['gene'],row['variant_id'],row['disease_id']
        tx=row['_tx']; print(f'[{n}/{len(top)}] {cid}: {gene} {vid}',flush=True)
        orpha_ids=links.get(did,set()) | ({did} if did.startswith('ORPHA:') else set())
        orpha=[]
        disease_ids={did}
        for oid in sorted(orpha_ids):
            data=disorders.get(oid,{})
            disease_ids.update(data.get('exact_ids',[]))
            assoc=[x for x in data.get('genes',[]) if x['gene']==gene]
            orpha.append({'orpha_id':oid,'name':data.get('name',''),
                'cross_reference_only_not_equivalence_proof':True,
                'gene_associations':assoc,'url':'https://www.orpha.net/en/disease/detail/'+oid.split(':')[1]})
        cg=[]
        if clingen is not None:
            for x in clingen[clingen['GENE SYMBOL']==gene].to_dict('records'):
                x['exact_disease_ID_match']=x.get('DISEASE ID (MONDO)','') in disease_ids
                cg.append(x)
        livecv=live_clinvar(cache,evidence[vid].get('clinvar_id',''),vid)
        records=livecv.get('records',[]) if livecv['status']=='retrieved_verified_allele' else rcv.get(vid,[])
        matching=[x for x in records if any(did in c.get('ids',[]) for c in x.get('conditions',[]))]
        inh=set(hpotable.loc[(hpotable['database_id']==did)&(~hpotable['qualifier'].str.contains('NOT',regex=False)),'hpo_id']) & modes.keys()
        mode_names=sorted(modes[x] for x in inh)
        state=row['allelic_state']
        compatibility='not_assessed_no_disease_specific_model'
        if mode_names==['AD']:
            compatibility='heterozygous_compatible_not_proof' if state=='heterozygous' else 'requires_manual_review'
        elif mode_names==['AR']:
            compatibility='biallelic_genotype_compatible_not_proof' if state=='homozygous_alt' else 'second_allele_and_phase_not_established'
        gn=gnomad(cache,vid)
        pub=publications(cache,gene,tx)
        om=omim(cache,did,gene,tx,gene_mim)
        source_data={'orphanet_download_status':orpha_status,'orphanet':orpha,
            'clingen_status':'retrieved' if clingen is not None else 'download_or_parse_failed',
            'clingen_gene_curations':cg,'gnomad_live':gn,'publications':pub,'omim':om,
            'clinvar_RCV_from_previous_stage':matching,
            'clinvar_live_status':livecv['status'],
            'clinvar_scope':'live_or_cached_online_response' if livecv['status']=='retrieved_verified_allele' else 'fallback_stage4_snapshot_or_no_data',
            'parents_and_ancestry':'not_supplied',
            'mechanism_and_literature_causality_review':'manual_review_required'}
        detail=out/(cid+'_'+gene+'_evidence.json');dump(detail,source_data)
        ev=evidence[vid]
        output_rows.append({
            'case_id':cid,'gene_rank':row['gene_rank'],'gene':gene,'variant_id':vid,
            'GT':row['GT'],'selection_score':row['selection_score'],'disease_id':did,'disease_name':row['disease_name'],
            'hpo_match_fraction':row['hpo_match_fraction'],'matched_count':row['matched_count'],
            'patient_present_count':row['patient_present_count'],
            'excluded_overlap_count':row['excluded_overlap_count'],
            'consequence':tx['consequence'],'impact':tx['impact'],'transcript_id':tx['transcript_id'],
            'hgvsc':tx['hgvsc'],'hgvsp':tx['hgvsp'],
            'gnomad_from_VEP':ev['gnomad_frequencies_json'],'gnomad_live_status':gn['status'],
            'common_AF_flag':row['common_AF_flag'],'ancestry_frequency_status':'unknown_ancestry',
            'clinvar_live_status':livecv['status'],'clinvar_RCV':json.dumps(matching,ensure_ascii=False),
            'clinvar_condition_specific_P_LP':row['condition_specific_P_LP'],
            'orphanet_status':json.dumps(orpha_status),
            'orphanet_gene_disease_association_found':any(x['gene_associations'] for x in orpha),
            'orphanet_urls':';'.join(x['url'] for x in orpha),
            'clingen_status':'retrieved' if clingen is not None else 'download_or_parse_failed',
            'clingen_exact_disease_curations':json.dumps([x for x in cg if x['exact_disease_ID_match']],ensure_ascii=False),
            'clingen_other_gene_curations':json.dumps([x for x in cg if not x['exact_disease_ID_match']],ensure_ascii=False),
            'HPO_inheritance':';'.join(mode_names),'genotype_inheritance_check':compatibility,
            'omim_status':om['status'],'omim_url':'https://omim.org/entry/'+did.split(':')[1] if did.startswith('OMIM:') else '',
            'literature_status':pub['status'],'literature_hits':json.dumps(pub.get('hits',[]),ensure_ascii=False),
            'other_HPO_diseases_for_this_gene':json.dumps(sorted({(x['disease_id'],x['disease_name'],x['hpo_match_fraction']) for x in selected_rows if x['case_id']==cid and x['gene']==gene}),ensure_ascii=False),
            'phenotype_support':'none' if float(row['hpo_match_fraction'])==0 else 'present',
            'overall_status':'candidate_requires_manual_causality_review',
            'evidence_file':detail.name})
        pd.DataFrame(output_rows).to_csv(out/'top3_genes_checks.tsv',sep='\t',index=False)
    dump(out/'run_info.json',{'completed_at_utc':pd.Timestamp.now(tz='UTC').isoformat(),
        'selection':'heuristic HPO + 0.15 condition-specific P/LP + 0.05 coding/splice - 0.8 common_AF - 0.8 condition-specific benign; not calibrated',
        'not_probability_model':True,'uses_development_labels':False,
        'missing_parents_and_ancestry':True,'clinical_confirmation_automated':False})
    print('\nГотово:',out/'top3_genes_checks.tsv',flush=True)
    print(pd.DataFrame(output_rows)[['case_id','gene_rank','gene','disease_name','hpo_match_fraction','omim_status']].to_string(index=False))

if __name__=='__main__':
    main()
