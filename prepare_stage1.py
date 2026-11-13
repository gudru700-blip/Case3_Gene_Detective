#!/usr/bin/env python3
"""Стадия 1: подготовка данных."""
import argparse
import csv
import gzip
import json
import re
from pathlib import Path


def write_tsv(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def case_id(path):
    match = re.search(r'case_\d+', str(path))
    if not match:
        raise ValueError(f'No case_XX identifier: {path}')
    return match.group()


def read_vcf(path, sample):
    opener = gzip.open if path.name.endswith('.gz') else open
    rows, samples = [], None
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                samples = line.rstrip('\r\n').split('\t')[9:]
                if sample and sample not in samples:
                    raise ValueError(f'{path}: sample {sample!r} not found')
                if len(samples) > 1 and not sample:
                    raise ValueError(f'{path}: multiple samples; specify --sample')
                continue
            if line.startswith('#') or not line.strip():
                continue
            if samples is None:
                raise ValueError(f'{path}: missing #CHROM header')
            cols = line.rstrip('\r\n').split('\t')
            if len(cols) < 8:
                raise ValueError(f'{path}: invalid VCF row')
            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
            if int(pos) < 1 or len(ref) != 1 or ref.upper() not in 'ACGT':
                raise ValueError(f'{path}: expected positive position and SNV: {line.strip()}')
            selected = sample or (samples[0] if samples else '')
            fmt = cols[8].split(':') if len(cols) > 8 else []
            values = cols[9 + samples.index(selected)].split(':') if selected else []
            data = dict(zip(fmt, values))
            gt = data.get('GT', '.')
            alleles = re.split(r'[/|]', gt)
            # Для каждого ALT создаём отдельную строку; GT остаётся исходным.
            for index, allele in enumerate(alt.split(','), 1):
                if len(allele) != 1 or allele.upper() not in 'ACGT':
                    raise ValueError(f'{path}: only SNVs supported: {alt}')
                canonical_chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
                missing = any(a == '.' for a in alleles)
                count = alleles.count(str(index))
                state = ('missing_or_partial' if missing else
                         'not_carried' if count == 0 else
                         'haploid_alt' if len(alleles) == 1 else
                         'homozygous_alt' if count == len(alleles) else 'heterozygous')
                rows.append(dict(case_id=case_id(path), assembly='GRCh38',
                    variant_id=f'{canonical_chrom}:{pos}:{ref}:{allele}',
                    chrom=canonical_chrom, pos=int(pos), ref=ref, alt=allele,
                    alt_index=index, sample=selected, GT=gt, allelic_state=state,
                    phase_separator_present='|' in gt, QUAL=qual, FILTER=filt,
                    DP=data.get('DP', '.'), GQ=data.get('GQ', '.'),
                    FORMAT=cols[8] if len(cols)>8 else '.',
                    sample_data=cols[9+samples.index(selected)] if selected else '.',
                    INFO=info, source=str(path)))
    if samples is None:
        raise ValueError(f'{path}: missing #CHROM header')
    return rows


def read_hpo(path):
    obj = json.loads(path.read_text(encoding='utf-8'))
    if isinstance(obj, list):
        features = obj
    elif isinstance(obj, dict):
        features = next((obj[k] for k in ('phenotypicFeatures', 'phenotypic_features',
                        'phenotypes', 'features', 'hpo_terms') if k in obj), None)
    else:
        features = None
    if not isinstance(features, list):
        raise ValueError(f'{path}: expected a feature list or phenotypicFeatures/phenotypes/features/hpo_terms list')
    rows = []
    for feature in features:
        if isinstance(feature, str):
            feature = {'id': feature}
        term = feature.get('type', feature)
        identifier = term.get('id', term.get('hpo_id', ''))
        if not re.fullmatch(r'HP:\d{7}', identifier):
            raise ValueError(f'{path}: invalid HPO identifier: {identifier!r}')
        # excluded=true означает отсутствие признака.
        excluded = feature.get('excluded', False)
        if not isinstance(excluded, bool):
            raise ValueError(f'{path}: excluded must be JSON true/false')
        rows.append(dict(case_id=case_id(path), hpo_id=identifier,
                         label=term.get('label', term.get('name', '')),
                         excluded=excluded, source=str(path)))
    states = {}
    for row in rows:
        key = row['hpo_id']
        if key in states and states[key] != row['excluded']:
            raise ValueError(f'{path}: conflicting present/excluded states for {key}')
        states[key] = row['excluded']
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--sample', help='Sample name for multisample VCF')
    args = parser.parse_args()
    root, out = args.input.resolve(), args.output.resolve()
    if root == out or root in out.parents:
        parser.error('--output must be outside the input directory')
    # Рекурсивно ищем исходные файлы во вложенных папках.
    vcfs = sorted(p for p in root.rglob('*') if p.name.endswith(('.vcf', '.vcf.gz')))
    phenos = sorted(p for p in root.rglob('*.json') if 'phenotype' in p.name.lower())
    if not vcfs or not phenos:
        parser.error('Need both case_XX*.vcf[.gz] and case_XX*phenotype*.json files')
    vmap, pmap = {}, {}
    for paths, mapping in ((vcfs, vmap), (phenos, pmap)):
        for path in paths:
            cid = case_id(path)
            if cid in mapping:
                raise ValueError(f'Multiple input files for {cid}: {mapping[cid]} / {path}')
            mapping[cid] = path
    # У каждого случая должны быть и VCF, и файл фенотипов.
    if vmap.keys() != pmap.keys():
        raise ValueError(f'Unpaired cases: {sorted(vmap.keys() ^ pmap.keys())}')
    # Сначала читаем и проверяем все входные данные, затем записываем результат.
    variants, hpo, summary = [], [], []
    for cid in sorted(vmap):
        vr = read_vcf(vmap[cid], args.sample)
        hr = read_hpo(pmap[cid])
        variants.extend(vr)
        hpo.extend(hr)
        summary.append(dict(case_id=cid, variant_allele_rows=len(vr),
                            hpo_present=len({r['hpo_id'] for r in hr if not r['excluded']}),
                            hpo_excluded=len({r['hpo_id'] for r in hr if r['excluded']})))
    out.mkdir(parents=True, exist_ok=True)
    fields = ['case_id','assembly','variant_id','chrom','pos','ref','alt','alt_index',
              'sample','GT','allelic_state','phase_separator_present','QUAL','FILTER',
              'DP','GQ','FORMAT','sample_data','INFO','source']
    write_tsv(out/'variants.tsv', fields, variants)
    hf = ['case_id','hpo_id','label','excluded','source']
    write_tsv(out/'hpo_present.tsv', hf, [r for r in hpo if not r['excluded']])
    write_tsv(out/'hpo_excluded.tsv', hf, [r for r in hpo if r['excluded']])
    write_tsv(out/'summary.tsv', list(summary[0]), summary)
    unique = {(r['chrom'],r['pos'],r['ref'],r['alt']) for r in variants}
    with (out/'vep_input.vcf').open('w', encoding='utf-8') as handle:
        handle.write('##fileformat=VCFv4.2\n##reference=GRCh38\n')
        handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref, alt in sorted(unique):
            handle.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n')
    print(f'Prepared {len(summary)} cases; {len(variants)} allele rows. Output: {out}')
    print('Assembly GRCh38 is an input assumption. REF has NOT been checked against FASTA.')
    print('Phase separators are preserved but do not establish phase in this synthetic dataset.')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, IndexError, OSError) as exc:
        raise SystemExit(f'ERROR: {exc}')
