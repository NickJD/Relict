"""Cumulative project overviews and decision deltas."""

from __future__ import annotations

import csv
import html
import json
from pathlib import Path

from branchmanager.pipeline.workflow_helpers import build_selection_decision


def write_performance_review_dashboard(selection_summary: str | Path, outdir: str | Path) -> str:
    """Write a compact, static decision dashboard linked to detailed evidence."""
    with open(selection_summary, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    counts = {}
    for row in rows:
        decision = row.get('Recommendation', 'UNKNOWN')
        counts[decision] = counts.get(decision, 0) + 1
    order = {
        'PRIORITISE - SET PRIMARY': 0, 'RESERVE - SET BACKUP': 1,
        'STRONG CANDIDATE': 2, 'SECONDARY CANDIDATE': 3,
        'REVIEW BEFORE SELECTION': 4, 'ALREADY SEQUENCED': 9,
    }
    rows.sort(key=lambda row: (order.get(row.get('Recommendation', ''), 8), row.get('SequenceID', '')))
    columns = [
        'SequenceID', 'PartnerID', 'Recommendation', 'SequencingSetRole',
        'EvidenceQuality', 'MarkerQC', 'GTDBTaxonomy', 'CulturedGap',
        'ProjectCoverage', 'GenomeCoverage', 'LocalTreeFigure', 'RecommendationReason',
    ]
    available = list(rows[0]) if rows else []
    columns = [column for column in columns if column in available]
    cards = ''.join(
        f'<div class="card"><b>{count}</b>{html.escape(decision)}</div>'
        for decision, count in sorted(counts.items(), key=lambda item: (order.get(item[0], 8), item[0]))
    )
    head = ''.join(f'<th>{html.escape(column)}</th>' for column in columns)
    body_rows = []
    for row in rows:
        cells = []
        for column in columns:
            value = str(row.get(column, ''))
            if column == 'LocalTreeFigure' and value not in {'', 'NA', 'None'}:
                target = value if value.startswith('assessment/') else f'assessment/{value}'
                value = f'<a href="{html.escape(target)}">view</a>'
            else:
                value = html.escape(value)
            cells.append(f'<td>{value}</td>')
        body_rows.append('<tr>' + ''.join(cells) + '</tr>')
    output = Path(outdir) / 'performance_review_dashboard.html'
    output.write_text(
        '<!doctype html><html><head><meta charset="utf-8"><title>BranchManager Performance Review</title>'
        '<style>body{font:14px system-ui;margin:28px;color:#17202a}h1{color:#123b2d}.cards{display:flex;flex-wrap:wrap;gap:10px}'
        '.card{border:1px solid #cbd7d1;border-radius:6px;padding:11px;min-width:150px}.card b{display:block;font-size:22px}'
        'table{border-collapse:collapse;width:100%;margin-top:20px}th,td{border-bottom:1px solid #dce3df;padding:7px;text-align:left;vertical-align:top}'
        'th{position:sticky;top:0;background:#edf4f0}</style></head><body><h1>Performance Review</h1>'
        '<p>BranchManager Hiring Panel summary. Detailed scientific evidence remains in <a href="assessment/sequence_assessment.tsv">sequence_assessment.tsv</a>.</p>'
        f'<div class="cards">{cards}</div><table><thead><tr>{head}</tr></thead><tbody>'
        + ''.join(body_rows) + '</tbody></table></body></html>\n'
    )
    return str(output)


def _latest_two_snapshots(db):
    with db.connect() as conn:
        ids = [row[0] for row in conn.execute(
            'SELECT snapshot_id FROM assessment_snapshots GROUP BY snapshot_id '
            'ORDER BY MAX(created_at) DESC, MAX(rowid) DESC LIMIT 2'
        ).fetchall()]
        snapshots = []
        for snapshot_id in ids:
            rows = conn.execute(
                'SELECT sequence_id, assessment_json FROM assessment_snapshots WHERE snapshot_id = ?',
                (snapshot_id,),
            ).fetchall()
            snapshots.append((snapshot_id, {
                str(sequence_id): json.loads(payload) for sequence_id, payload in rows
            }))
    return snapshots


def write_decision_changes(db, path: str | Path) -> str:
    snapshots = _latest_two_snapshots(db)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        'SequenceID', 'PreviousSnapshot', 'CurrentSnapshot', 'PreviousDecision',
        'CurrentDecision', 'DecisionChanged', 'PreviousSetRole', 'CurrentSetRole',
        'PreviousPangenomeGap', 'CurrentPangenomeGap', 'ChangeReason',
    ]
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        if len(snapshots) < 2:
            return str(output)
        (current_id, current), (previous_id, previous) = snapshots
        for sequence_id in sorted(set(previous) | set(current)):
            old = previous.get(sequence_id, {})
            new = current.get(sequence_id, {})
            old_decision = build_selection_decision(old)['decision'] if old else 'NOT_ASSESSED'
            new_decision = build_selection_decision(new)['decision'] if new else 'NOT_ASSESSED'
            changes = []
            for label, key in (
                ('set role', 'sequencing_set_role'), ('pangenome gap', 'pangenome_gap'),
                ('nearest genome', 'nearest_genome_hit'), ('marker QC', 'marker_qc_flag'),
            ):
                if str(old.get(key, 'NA')) != str(new.get(key, 'NA')):
                    changes.append(f'{label}: {old.get(key, "NA")} -> {new.get(key, "NA")}')
            writer.writerow({
                'SequenceID': sequence_id,
                'PreviousSnapshot': previous_id,
                'CurrentSnapshot': current_id,
                'PreviousDecision': old_decision,
                'CurrentDecision': new_decision,
                'DecisionChanged': 'yes' if old_decision != new_decision else 'no',
                'PreviousSetRole': old.get('sequencing_set_role', 'NA'),
                'CurrentSetRole': new.get('sequencing_set_role', 'NA'),
                'PreviousPangenomeGap': old.get('pangenome_gap', 'NA'),
                'CurrentPangenomeGap': new.get('pangenome_gap', 'NA'),
                'ChangeReason': '; '.join(changes) or 'recalculation did not change tracked evidence',
            })
    return str(output)


def write_annual_report(db, outdir: str | Path) -> dict:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    with db.connect() as conn:
        dataset_rows = conn.execute(
            'SELECT s.dataset, COALESCE(r.role, "unregistered"), '
            'CASE WHEN COALESCE(r.genomes_available, 0) = 1 THEN COUNT(s.id) '
            '     ELSE COALESCE(SUM(COALESCE(m.selected_for_wgs, 0)), 0) END, '
            'COUNT(s.id) FROM sequences s '
            'LEFT JOIN dataset_roles r ON r.dataset=s.dataset '
            'LEFT JOIN sequencing_metadata m ON m.id=s.id '
            'GROUP BY s.dataset, r.role, r.genomes_available ORDER BY s.dataset'
        ).fetchall()
        status_rows = conn.execute(
            'SELECT status, COUNT(*) FROM isolate_status GROUP BY status ORDER BY status'
        ).fetchall()
        genome_rows = conn.execute(
            'SELECT genome_id, sequence_id, accession, genome_status, genome_qc_pass, '
            'completeness, contamination, gtdb_taxonomy, ani_cluster FROM genome_records '
            'ORDER BY sequence_id, genome_id'
        ).fetchall()
        run_rows = conn.execute(
            'SELECT run_id, workflow, dataset, status, started_at, completed_at, manifest_path '
            'FROM project_runs ORDER BY started_at, run_id'
        ).fetchall()
        round_rows = conn.execute(
            'SELECT round_id, mode, created_at FROM selection_rounds ORDER BY created_at, round_id'
        ).fetchall()
        isolate_rows = conn.execute(
            'SELECT s.id, s.dataset, COALESCE(m.partner_id, ""), '
            'COALESCE(i.status, "RECEIVED"), COALESCE(p.marker_qc_class, "UNVERIFIED"), '
            'COALESCE(m.selected_for_sequencing, 0), COALESCE(m.selected_for_wgs, 0) FROM sequences s '
            'LEFT JOIN sequencing_metadata m ON m.id=s.id '
            'LEFT JOIN isolate_status i ON i.sequence_id=s.id '
            'LEFT JOIN sequence_provenance p ON p.sequence_id=s.id ORDER BY s.id'
        ).fetchall()
        removal_rows = conn.execute(
            'SELECT sequence_id, original_dataset, partner_id, sequence_length, '
            'taxonomy, reason, removed_at FROM sequence_removals '
            'ORDER BY removed_at, removal_id'
        ).fetchall()

    active = {
        str(row[0]): {
            'dataset': str(row[1] or ''), 'partner_id': str(row[2] or ''),
            'status': str(row[3] or ''), 'marker_qc': str(row[4] or ''),
        }
        for row in isolate_rows
    }
    dataset_roles = db.get_dataset_roles()
    latest_assessments = db.get_latest_assessment_rows()
    candidate_rows = []
    recommendation_counts = {}
    for sequence_id, metadata in active.items():
        if dataset_roles.get(metadata['dataset'], {}).get('role') != 'candidate':
            continue
        assessment = latest_assessments.get(sequence_id, {})
        decision = build_selection_decision(assessment) if assessment else {
            'decision': 'NOT YET ASSESSED', 'evidence_quality': 'NA',
            'decision_reason': 'no stored Performance Review assessment',
        }
        recommendation = decision['decision']
        recommendation_counts[recommendation] = recommendation_counts.get(recommendation, 0) + 1
        candidate_rows.append((
            sequence_id, metadata['dataset'], metadata['partner_id'], metadata['status'],
            metadata['marker_qc'], recommendation, decision.get('evidence_quality', 'NA'),
            assessment.get('sequencing_set_role', 'NA'), assessment.get('taxonomy', 'NA'),
            assessment.get('nearest_hit', 'NA'), assessment.get('nearest_identity', 'NA'),
            assessment.get('project_nearest_hit', 'NA'),
            assessment.get('project_nearest_identity', 'NA'),
            assessment.get('investigation_score', 'NA'),
            decision.get('decision_reason', ''),
        ))

    isolate_path = output / 'project_isolate_ledger.tsv'
    with open(isolate_path, 'w') as handle:
        handle.write('SequenceID\tDataset\tPartnerID\tOperationalStatus\tMarkerQC\tSelectedForGenomeSequencing\tGenomeAlreadySequenced\n')
        for row in isolate_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    genome_path = output / 'project_genome_ledger.tsv'
    with open(genome_path, 'w') as handle:
        handle.write('GenomeID\tSequenceID\tAccession\tGenomeStatus\tGenomeQCPass\tCompleteness\tContamination\tGTDBTaxonomy\tANICluster\n')
        for row in genome_rows:
            handle.write('\t'.join('NA' if value is None else str(value) for value in row) + '\n')
    candidate_path = output / 'project_candidate_overview.tsv'
    with open(candidate_path, 'w') as handle:
        handle.write(
            'SequenceID\tDataset\tPartnerID\tOperationalStatus\tMarkerQC\tRecommendation\t'
            'EvidenceQuality\tSequencingSetRole\tGTDBTaxonomy\tBaselineNearestHit\t'
            'BaselineNearestIdentity\tProjectNearestHit\tProjectNearestIdentity\t'
            'InvestigationScore\tRecommendationReason\n'
        )
        for row in candidate_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    removals_path = output / 'sequence_removal_ledger.tsv'
    with open(removals_path, 'w') as handle:
        handle.write('SequenceID\tOriginalDataset\tPartnerID\tSequenceLength\tTaxonomy\tReason\tRemovedAt\n')
        for row in removal_rows:
            handle.write('\t'.join('NA' if value is None else str(value) for value in row) + '\n')
    changes_path = write_decision_changes(db, output / 'decision_changes.tsv')

    baseline_count = sum(row[3] for row in dataset_rows if row[1] == 'baseline')
    candidate_count = sum(row[3] for row in dataset_rows if row[1] == 'candidate')
    cards = [
        ('Cultured baseline markers', baseline_count),
        ('Active partner candidates', candidate_count),
        ('Withdrawn sequences', len(removal_rows)),
        ('QC-passed genomes', sum(bool(row[4]) for row in genome_rows)),
        ('Genome QC failures', sum(row[3] == 'GENOME_QC_FAILED' for row in genome_rows)),
        ('Selection rounds', len(round_rows)),
    ]
    def table(headers, rows):
        head = ''.join(f'<th>{html.escape(str(item))}</th>' for item in headers)
        body = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(value or ""))}</td>' for value in row) + '</tr>' for row in rows)
        return f'<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'

    report = output / 'annual_report.html'
    report.write_text(
        '<!doctype html><html><head><meta charset="utf-8"><title>BranchManager Project Overview</title>'
        '<style>body{font:15px system-ui;margin:32px;color:#17202a;max-width:1200px}'
        'h1,h2{color:#123b2d}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}'
        '.card{border:1px solid #cad5cf;padding:14px;border-radius:6px}.card b{font-size:26px;display:block}'
        'table{border-collapse:collapse;width:100%;margin:10px 0 28px}th,td{border-bottom:1px solid #d9e1dd;text-align:left;padding:7px}th{background:#eef4f1}</style>'
        '</head><body><h1>Cumulative Project Overview</h1><p>BranchManager point-in-time marker-to-genome project report.</p>'
        '<div class="cards">' + ''.join(f'<div class="card"><b>{value}</b>{html.escape(label)}</div>' for label, value in cards) + '</div>'
        '<h2>Datasets</h2>' + table(('Dataset', 'Role', 'Genomes available', 'Sequences'), dataset_rows) +
        '<h2>Isolate lifecycle</h2>' + table(('Status', 'Count'), status_rows) +
        '<h2>Current recommendations</h2>' + table(('Recommendation', 'Candidates'), sorted(recommendation_counts.items())) +
        '<h2>Workflow runs</h2>' + table(('Run', 'Workflow', 'Dataset', 'Status', 'Started', 'Completed', 'Manifest'), run_rows) +
        '<h2>Selection rounds</h2>' + table(('Round', 'Mode', 'Created'), round_rows) +
        '<h2>Exit Interviews</h2>' + table(('Sequence', 'Dataset', 'Partner', 'Length', 'Taxonomy', 'Reason', 'Removed'), removal_rows) +
        '<p>Detailed ledgers: <a href="project_isolate_ledger.tsv">isolates</a>, '
        '<a href="project_candidate_overview.tsv">current candidate assessments</a>, '
        '<a href="project_genome_ledger.tsv">genomes</a>, '
        '<a href="sequence_removal_ledger.tsv">Exit Interviews</a>, '
        '<a href="decision_changes.tsv">latest decision changes</a>.</p></body></html>\n'
    )
    return {
        'html': str(report), 'isolates': str(isolate_path), 'genomes': str(genome_path),
        'candidates': str(candidate_path), 'removals': str(removals_path),
        'decision_changes': changes_path,
    }
