"""Explain the original-China cohort intersection without changing any data."""
import collections
import json
import pandas as pd
from global_common import ROOT, OUT, CFG, path, ogr
from global_statistics import CHINA_IDS
from spatial_rules import exclusion_reason


def main():
    original = pd.read_csv(CHINA_IDS)
    screened = pd.read_csv(ROOT / 'tables/global_species_screened.csv')
    bridge = pd.read_csv(OUT / 'tables/china_matched_species_comparison.csv')
    old_ids = set(original.id_no.astype(int))
    kept_ids = set(screened.species_id.astype(int))
    matched_ids = set(bridge.species_id.astype(int))
    missing = old_ids - matched_ids
    evidence = collections.defaultdict(collections.Counter)
    raw_present = set()
    for source in sorted(path(CFG['species_root']).glob('*/*.shp')):
        dataset = ogr.Open(str(source), 0)
        layer = dataset.GetLayer(0)
        layer.SetIgnoredFields(['OGR_GEOMETRY'])
        for feature in layer:
            key = feature.GetField('id_no')
            if key is None or int(key) not in missing:
                continue
            key = int(key)
            raw_present.add(key)
            evidence[key][exclusion_reason(dict(feature.items())) or 'retained_record'] += 1
        dataset = None
    rows = []
    for row in original.to_dict('records'):
        key = int(row['id_no'])
        if key in matched_ids:
            status = 'matched_China_common_land'
        elif key in kept_ids:
            status = 'globally_screened_but_no_retained_China_common_land'
        elif key not in raw_present:
            status = 'ID_absent_from_current_global_source_files'
        else:
            status = 'excluded_by_global_attribute_rules'
        rows.append({'species_id': key, 'original_taxon': row['class'],
                     'original_category': row['category'], 'match_status': status,
                     'source_record_screening': json.dumps(dict(evidence[key]), sort_keys=True) if key in missing else ''})
    assert len(rows) == 1820 and len(matched_ids) == sum(r['match_status'] == 'matched_China_common_land' for r in rows)
    pd.DataFrame(rows).to_csv(OUT / 'tables/china_cohort_match_status.csv', index=False, encoding='utf-8-sig')
    summary = {'original_species': len(old_ids), 'matched_species': len(matched_ids),
               'status_counts': dict(collections.Counter(r['match_status'] for r in rows)),
               'note': 'No records, filters, habitat metrics or original-China results were changed.'}
    (OUT / 'progress/china_cohort_match_review.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
