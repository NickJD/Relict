"""Cumulative project overviews and decision deltas."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import sqlite3
from collections import Counter, defaultdict
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
        'SECONDARY - STRAIN DIVERSITY': 2, 'STRONG CANDIDATE': 3,
        'SECONDARY CANDIDATE': 4, 'REVIEW BEFORE SELECTION': 5,
        'REVIEW - PANGENOME BOUNDARY': 6,
        'EXCLUDE - BASELINE REDUNDANT': 7, 'ALREADY SEQUENCED': 10,
    }
    rows.sort(key=lambda row: (order.get(row.get('Recommendation', ''), 8), row.get('SequenceID', '')))
    columns = [
        'SequenceID', 'PartnerID', 'Recommendation', 'SequencingSetRole',
        'SequencingSetRank', 'SelectionGroupType', 'SelectionDiversityDistance', 'EvidenceQuality', 'MarkerQC',
        'GTDBTaxonomy', 'CulturedGap', 'BaselineExtensionStatus', 'BaselineRedundancyStatus', 'ProjectCoverage',
        'MWLMatchedRank', 'GenomeCoverage', 'SpeciesContext', 'LocalTreeFigure', 'RecommendationReason',
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


def _project_root_from_db(db) -> Path | None:
    db_path = str(getattr(db, 'path', '') or '')
    if not db_path or db_path == ':memory:':
        return None
    try:
        return Path(db_path).expanduser().resolve().parent
    except Exception:
        return None


def _query_dataset_rows(conn, datasets=None) -> list[tuple[str, str, str, int, int, int]]:
    params = ()
    where = ''
    if datasets is not None:
        names = [str(name) for name in datasets if str(name or '').strip()]
        if not names:
            return []
        where = 'WHERE s.dataset IN (' + ','.join('?' for _ in names) + ') '
        params = tuple(names)
    role_columns = {row[1] for row in conn.execute('PRAGMA table_info(dataset_roles)').fetchall()}
    tier_expr = 'r.baseline_tier' if 'baseline_tier' in role_columns else "''"
    rows = conn.execute(
        'SELECT s.dataset, COALESCE(r.role, "unregistered"), COALESCE(' + tier_expr + ', ""), '
        'CASE WHEN COALESCE(r.genomes_available, 0) = 1 THEN COUNT(s.id) '
        '     ELSE COALESCE(SUM(COALESCE(m.selected_for_wgs, 0)), 0) END, '
        "COALESCE(SUM(CASE WHEN UPPER(COALESCE(i.status, '')) = 'PROPOSED' "
        'THEN 1 ELSE 0 END), 0), '
        'COUNT(s.id) FROM sequences s '
        'LEFT JOIN dataset_roles r ON r.dataset=s.dataset '
        'LEFT JOIN sequencing_metadata m ON m.id=s.id '
        'LEFT JOIN isolate_status i ON i.sequence_id=s.id '
        f'{where}'
        'GROUP BY s.dataset, r.role, r.genomes_available, ' + tier_expr + ' ORDER BY s.dataset',
        params,
    ).fetchall()
    return [
        (
            str(dataset or ''),
            str(role or 'unregistered'),
            str(baseline_tier or ''),
            int(genomes or 0),
            int(proposed or 0),
            int(sequences or 0),
        )
        for dataset, role, baseline_tier, genomes, proposed, sequences in rows
        if str(dataset or '').strip()
    ]


def _query_status_rows_for_datasets(conn, datasets) -> list[tuple[str, int]]:
    names = [str(name) for name in datasets if str(name or '').strip()]
    if not names:
        return []
    rows = conn.execute(
        "SELECT COALESCE(i.status, 'RECEIVED'), COUNT(*) FROM sequences s "
        'LEFT JOIN isolate_status i ON i.sequence_id=s.id '
        'WHERE s.dataset IN (' + ','.join('?' for _ in names) + ') '
        "GROUP BY COALESCE(i.status, 'RECEIVED') "
        "ORDER BY COALESCE(i.status, 'RECEIVED')",
        tuple(names),
    ).fetchall()
    return [(str(status or 'RECEIVED'), int(count or 0)) for status, count in rows]


def _query_isolate_rows(conn, datasets=None):
    params = ()
    where = ''
    if datasets is not None:
        names = [str(name) for name in datasets if str(name or '').strip()]
        if not names:
            return []
        where = 'WHERE s.dataset IN (' + ','.join('?' for _ in names) + ') '
        params = tuple(names)
    return conn.execute(
        'SELECT s.id, s.dataset, COALESCE(m.partner_id, ""), '
        'COALESCE(i.status, "RECEIVED"), COALESCE(p.marker_qc_class, "UNVERIFIED"), '
        'COALESCE(m.selected_for_sequencing, 0), COALESCE(m.selected_for_wgs, 0) FROM sequences s '
        'LEFT JOIN sequencing_metadata m ON m.id=s.id '
        'LEFT JOIN isolate_status i ON i.sequence_id=s.id '
        'LEFT JOIN sequence_provenance p ON p.sequence_id=s.id '
        f'{where}'
        'ORDER BY s.id',
        params,
    ).fetchall()


def _command_value(command, flag: str) -> str:
    if not isinstance(command, list):
        return ''
    prefix = f'{flag}='
    for index, token in enumerate(command):
        value = str(token)
        if value == flag and index + 1 < len(command):
            return str(command[index + 1])
        if value.startswith(prefix):
            return value[len(prefix):]
    return ''


def _dataset_from_path_hint(value: object) -> str:
    raw = str(value or '').strip()
    if not raw:
        return ''
    try:
        parts = [part for part in Path(raw).parts if part not in {'', '/'}]
    except TypeError:
        return ''
    for index, part in enumerate(parts):
        if part == 'mailroom' and index > 0:
            return parts[index - 1]
    for part in reversed(parts):
        if part.startswith('01_interview_') and len(part) > len('01_interview_'):
            return part[len('01_interview_'):]
    stage_names = {
        '00_it_desk', '01_onboarding', '02_paper_trail_merge_meeting',
        '03_performance_review_hiring_panel',
    }
    for index, part in enumerate(parts):
        if part in stage_names and index > 0:
            return parts[index - 1]
    return ''


def _manifest_dataset(data: dict, path: Path | None = None) -> str:
    command_dataset = str(_command_value(data.get('command'), '--dataset') or '').strip()
    if command_dataset:
        return command_dataset
    direct = str(data.get('dataset') or '').strip()
    if direct:
        return direct
    mailroom = _command_value(data.get('command'), '--mailroom')
    dataset = _dataset_from_path_hint(mailroom)
    if dataset:
        return dataset
    for record in data.get('inputs', []):
        if str(record.get('role') or '') in {'sample_map', 'mailroom_summary'}:
            dataset = _dataset_from_path_hint(record.get('path'))
            if dataset:
                return dataset
    if path is not None:
        dataset = _dataset_from_path_hint(path)
        if dataset:
            return dataset
    return ''


def _manifest_run_id(path: Path, data: dict, project_root: Path) -> str:
    run_id = str(data.get('run_id') or '').strip()
    if run_id:
        return run_id
    workflow = str(data.get('workflow') or 'workflow').strip()
    dataset = _manifest_dataset(data, path)
    started = str(data.get('started_at') or '').strip()
    try:
        relative = path.relative_to(project_root).as_posix()
    except ValueError:
        relative = path.as_posix()
    return ':'.join(part for part in (workflow, dataset, started, relative) if part)


def _discover_run_manifests(project_root: Path | None, outdir: Path):
    if project_root is None or not project_root.exists():
        return []
    output = outdir.expanduser().resolve()
    rows = []
    for path in sorted(project_root.rglob('run_manifest.json')):
        resolved = path.resolve()
        try:
            resolved.relative_to(output)
            continue
        except ValueError:
            pass
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError, TypeError):
            continue
        workflow = str(data.get('workflow') or '').strip()
        if not workflow or workflow == 'annual_report':
            continue
        rows.append((
            _manifest_run_id(path, data, project_root),
            workflow,
            _manifest_dataset(data, path),
            str(data.get('status') or '').strip(),
            str(data.get('started_at') or '').strip(),
            str(data.get('completed_at') or '').strip(),
            str(path),
        ))
    return rows


def _resolve_manifest_path(value: object, manifest_path: Path, project_root: Path | None) -> Path:
    raw = str(value or '').strip()
    candidate = Path(raw).expanduser()
    if candidate.exists():
        return candidate
    if project_root is not None:
        parts = list(candidate.parts)
        if project_root.name in parts:
            index = parts.index(project_root.name)
            remapped = project_root.joinpath(*parts[index + 1:])
            if remapped.exists():
                return remapped
    local = manifest_path.parent / candidate.name
    if local.exists():
        return local
    return candidate


def _output_paths(data: dict, role: str, manifest_path: Path, project_root: Path | None) -> list[Path]:
    return [
        _resolve_manifest_path(record.get('path'), manifest_path, project_root)
        for record in data.get('outputs', [])
        if str(record.get('role') or '') == role and str(record.get('path') or '').strip()
    ]


def _first_output_path(data: dict, role: str, manifest_path: Path, project_root: Path | None) -> Path | None:
    paths = _output_paths(data, role, manifest_path, project_root)
    return paths[0] if paths else None


def _read_tsv_rows(path: Path | None) -> list[dict]:
    if path is None or not path.exists() or not path.is_file():
        return []
    try:
        with open(path, newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))
    except (OSError, csv.Error):
        return []


def _manifest_input_signature(data: dict) -> str:
    for role in ('sample_map', 'mailroom_summary'):
        for record in data.get('inputs', []):
            if str(record.get('role') or '') == role and record.get('sha256'):
                return str(record.get('sha256'))
    return ''


def _sequence_signature(rows: list[dict]) -> str:
    ids = sorted({str(row.get('SequenceID') or row.get('sequence_id') or '').strip() for row in rows})
    ids = [sid for sid in ids if sid]
    if not ids:
        return ''
    digest = hashlib.sha1('\n'.join(ids).encode('utf-8')).hexdigest()
    return f'sequence_ids:{digest}'


def _qc_summary_from_manifest(path: Path, data: dict, project_root: Path | None) -> dict | None:
    workflow = str(data.get('workflow') or '').strip()
    if workflow not in {'interview', 'paper_trail'}:
        return None
    dataset = _manifest_dataset(data, path)
    assembly_path = _first_output_path(data, 'assembly_tsv', path, project_root)
    read_qc_path = _first_output_path(data, 'read_qc_tsv', path, project_root)
    failed_path = _first_output_path(data, 'failed_manifest_tsv', path, project_root)
    failed_read_path = _first_output_path(data, 'failed_read_manifest_tsv', path, project_root)
    visual_manifest_path = _first_output_path(data, 'visual_manifest_tsv', path, project_root)
    assembly_rows = _read_tsv_rows(assembly_path)
    read_rows = _read_tsv_rows(read_qc_path)
    failed_rows = _read_tsv_rows(failed_path)
    failed_read_rows = _read_tsv_rows(failed_read_path)
    visual_rows = _read_tsv_rows(visual_manifest_path)
    if not any((assembly_rows, read_rows, failed_rows, failed_read_rows, visual_rows)):
        return None
    passed_sequences = sum(
        1 for row in assembly_rows
        if str(row.get('QCClass') or '').strip().upper().startswith('PASS')
    )
    accepted_sequences = sum(
        1 for row in assembly_rows
        if str(row.get('Recommendation') or '').strip().upper() == 'ACCEPT'
    )
    review_sequences = sum(
        1 for row in assembly_rows
        if str(row.get('Recommendation') or '').strip().upper() == 'MANUAL_REVIEW'
        or str(row.get('QCClass') or '').strip().upper() == 'PASS_WITH_WARNINGS'
    )
    failed_sequences = len(failed_rows) or sum(
        1 for row in assembly_rows
        if str(row.get('Recommendation') or '').strip().upper() not in {'', 'ACCEPT'}
    )
    failed_read_files = len(failed_read_rows) or sum(
        1 for row in read_rows
        if str(row.get('Status') or '').strip().lower() not in {'', 'kept'}
    )
    report_counts = Counter(str(row.get('Report') or 'visual_pages') for row in visual_rows)
    sequence_signature = _sequence_signature(assembly_rows)
    signature = sequence_signature or _manifest_input_signature(data) or str(path)
    preference = (
        1 if workflow != 'interview' else 0,
        1 if dataset else 0,
        passed_sequences + failed_sequences,
        len(read_rows),
    )
    return {
        'dataset': dataset,
        'workflow': workflow,
        'started_at': str(data.get('started_at') or ''),
        'manifest_path': str(path),
        'signature': signature,
        'preference': preference,
        'passed_sequence_files': passed_sequences,
        'accepted_sequence_files': accepted_sequences,
        'review_sequence_files': review_sequences,
        'failed_sequence_files': failed_sequences,
        'read_files': len(read_rows),
        'failed_read_files': failed_read_files,
        'visual_pages': len(visual_rows),
        'read_error_pages': report_counts.get('read_error_profiles', len(_output_paths(data, 'read_error_pngs', path, project_root))),
        'chromatogram_pages': report_counts.get('trace_chromatograms', len(_output_paths(data, 'chromatogram_pngs', path, project_root))),
        'assembly_pages': report_counts.get('assembly_overviews', len(_output_paths(data, 'assembly_pngs', path, project_root))),
        'visual_manifest_path': str(visual_manifest_path or ''),
    }


def _discover_qc_summaries(project_root: Path | None, outdir: Path) -> list[dict]:
    if project_root is None or not project_root.exists():
        return []
    output = outdir.expanduser().resolve()
    summaries: dict[tuple[str, str], dict] = {}
    for path in sorted(project_root.rglob('run_manifest.json')):
        resolved = path.resolve()
        try:
            resolved.relative_to(output)
            continue
        except ValueError:
            pass
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError, TypeError):
            continue
        summary = _qc_summary_from_manifest(path, data, project_root)
        if not summary:
            continue
        key = (summary['dataset'], summary['signature'])
        previous = summaries.get(key)
        if previous is None or summary['preference'] > previous['preference']:
            summaries[key] = summary
    return sorted(
        summaries.values(),
        key=lambda row: (str(row.get('dataset') or ''), str(row.get('started_at') or ''), str(row.get('manifest_path') or '')),
    )


def _query_dataset_marker_qc_rows(conn) -> dict[str, dict[str, int]]:
    rows = conn.execute(
        "SELECT s.dataset, "
        "SUM(CASE WHEN UPPER(COALESCE(p.marker_qc_class, '')) LIKE 'PASS%' THEN 1 ELSE 0 END), "
        "SUM(CASE WHEN UPPER(COALESCE(p.marker_qc_class, '')) IN "
        "('PASS_WITH_WARNINGS', 'MARKER_QC_REVIEW_REQUIRED') "
        "OR UPPER(COALESCE(i.status, '')) = 'TRACE_REVIEW' THEN 1 ELSE 0 END), "
        "SUM(CASE WHEN UPPER(COALESCE(p.marker_qc_class, '')) LIKE 'FAIL%' "
        "OR UPPER(COALESCE(i.status, '')) LIKE '%FAILED%' THEN 1 ELSE 0 END) "
        "FROM sequences s "
        "LEFT JOIN sequence_provenance p ON p.sequence_id=s.id "
        "LEFT JOIN isolate_status i ON i.sequence_id=s.id "
        "GROUP BY s.dataset"
    ).fetchall()
    return {
        str(dataset or ''): {
            'active_marker_qc_passed': int(passed or 0),
            'active_marker_qc_review': int(review or 0),
            'active_marker_qc_failed': int(failed or 0),
        }
        for dataset, passed, review, failed in rows
        if str(dataset or '').strip()
    }


def _aggregate_qc_summaries(qc_summaries: list[dict]) -> dict[str, dict[str, int]]:
    totals = defaultdict(lambda: {
        'passed_sequence_files': 0,
        'accepted_sequence_files': 0,
        'review_sequence_files': 0,
        'failed_sequence_files': 0,
        'read_files': 0,
        'failed_read_files': 0,
        'visual_pages': 0,
    })
    for row in qc_summaries:
        dataset = str(row.get('dataset') or '').strip()
        if not dataset:
            continue
        for key in totals[dataset]:
            totals[dataset][key] += int(row.get(key) or 0)
    return dict(totals)


def _html_path_target(path: object, output: Path) -> str:
    raw = str(path or '').strip()
    if not raw:
        return ''
    p = Path(raw).expanduser()
    if p.exists():
        try:
            return Path(os.path.relpath(p.resolve(), output.resolve())).as_posix()
        except (OSError, ValueError):
            return p.resolve().as_posix()
    return raw


def _html_path_link(path: object, output: Path, label: str | None = None) -> str:
    raw = str(path or '').strip()
    if not raw:
        return ''
    target = _html_path_target(raw, output)
    return f'<a href="{html.escape(target)}">{html.escape(label or Path(raw).name or raw)}</a>'


def _resolve_visual_file(manifest_path: Path, file_value: object) -> Path | None:
    raw = str(file_value or '').strip()
    if not raw:
        return None
    candidate = Path(raw).expanduser()
    candidates = [candidate] if candidate.is_absolute() else []
    if not candidate.is_absolute():
        candidates.extend([
            manifest_path.parent.parent / candidate,
            manifest_path.parent / candidate,
            manifest_path.parent / candidate.name,
        ])
    for item in candidates:
        if item.exists() and item.is_file() and item.suffix.lower() == '.png':
            return item
    return None


def _paper_trail_gallery_rows(
    dataset: str, workflow: str, manifest_path: Path, output: Path, *, max_per_report: int = 1,
) -> list[dict]:
    rows = []
    counts: Counter[str] = Counter()
    for row in _read_tsv_rows(manifest_path):
        report = str(row.get('Report') or 'visual_report').strip()
        if counts[report] >= max_per_report:
            continue
        image = _resolve_visual_file(manifest_path, row.get('File'))
        if image is None:
            continue
        page = str(row.get('Page') or '').strip()
        page_count = str(row.get('PageCount') or '').strip()
        first_record = str(row.get('FirstRecord') or '').strip()
        last_record = str(row.get('LastRecord') or '').strip()
        range_text = ''
        if first_record and last_record:
            range_text = f' ({first_record} to {last_record})'
        caption = report.replace('_', ' ').title()
        if page:
            caption += f' page {page}'
            if page_count:
                caption += f' of {page_count}'
        caption += range_text
        rows.append({
            'dataset': dataset,
            'workflow': workflow,
            'report': report.replace('_', ' ').title(),
            'caption': caption,
            'image_path': image.as_posix(),
            'link_path': image.as_posix(),
            'src': _html_path_target(image, output),
            'href': _html_path_target(image, output),
        })
        counts[report] += 1
    return rows


def _neighbourhood_gallery_rows(
    dataset: str, workflow: str, manifest_path: Path, output: Path, *, max_images: int = 2,
) -> list[dict]:
    rows = []
    seen = set()
    for row in _read_tsv_rows(manifest_path):
        if str(row.get('Status') or '').strip().lower() != 'rendered':
            continue
        figure = str(row.get('Figure') or '').strip()
        if not figure or figure in seen:
            continue
        image = _resolve_visual_file(manifest_path, figure)
        if image is None:
            continue
        assessed = str(row.get('AssessedSequencesInFigure') or '').strip()
        leaves = str(row.get('TreeLeavesShown') or '').strip()
        caption = Path(figure).stem.replace('_', ' ').title()
        if assessed:
            caption += f' ({assessed})'
        if leaves:
            caption += f'; {leaves} leaves shown'
        rows.append({
            'dataset': dataset,
            'workflow': workflow,
            'report': 'Local Neighbourhood',
            'caption': caption,
            'image_path': image.as_posix(),
            'link_path': image.as_posix(),
            'src': _html_path_target(image, output),
            'href': _html_path_target(image, output),
        })
        seen.add(figure)
        if len(rows) >= max_images:
            break
    return rows


def _build_visual_gallery_rows(visual_report_rows: list[tuple], output: Path, *, max_images: int = 12) -> list[dict]:
    gallery_rows = []
    for dataset, workflow, report_name, _items, path_value in visual_report_rows:
        manifest_path = Path(str(path_value or '')).expanduser()
        if not manifest_path.exists():
            continue
        if str(report_name) == 'Paper Trail visual pages':
            gallery_rows.extend(_paper_trail_gallery_rows(str(dataset), str(workflow), manifest_path, output))
        elif str(report_name) == 'Local neighbourhood figures':
            gallery_rows.extend(_neighbourhood_gallery_rows(str(dataset), str(workflow), manifest_path, output))
        if len(gallery_rows) >= max_images:
            return gallery_rows[:max_images]
    return gallery_rows[:max_images]


def _visual_gallery_html(rows: list[dict]) -> str:
    if not rows:
        return '<p class="muted">No embeddable PNG previews were discovered for this report.</p>'
    cards = []
    for row in rows:
        title = ' / '.join(part for part in (row.get('dataset'), row.get('report')) if part)
        cards.append(
            '<figure class="visual-card">'
            f'<a href="{html.escape(str(row.get("href") or ""))}">'
            f'<img src="{html.escape(str(row.get("src") or ""))}" alt="{html.escape(str(row.get("caption") or title))}" loading="lazy"></a>'
            f'<figcaption><b>{html.escape(title)}</b><span>{html.escape(str(row.get("caption") or ""))}</span></figcaption>'
            '</figure>'
        )
    return '<div class="gallery">' + ''.join(cards) + '</div>'


def _annual_report_readme_markdown() -> str:
    return """# BranchManager Annual Report Guide

This folder is a point-in-time project overview built from the BranchManager database, workflow manifests, and run outputs.

## How to Read the HTML

- Summary cards give project-wide counts for baseline markers, active candidate markers, proposed sequencing commitments, marker-QC outcomes, and visual pages.
- The Datasets table combines database state with Paper Trail / Interview file-level QC summaries. `QC-passed sequence files` includes accepted files and files that passed with review warnings; `Accepted sequence files` are ready for downstream review without manual marker review; `Review sequence files` passed with warnings; `QC-failed sequence files` should be resequenced or excluded.
- Current recommendations come from the latest Performance Review assessment snapshot. Baseline-redundant candidates are near-identical to cultured baseline markers at the configured identity and coverage thresholds.
- Embedded visual previews are a curated subset of the available PNG outputs. Click any preview to open the full-size PNG, and use `project_visual_reports.tsv` for the full manifest of visual artifacts.
- Workflow runs show reproducibility context: workflow name, inferred dataset, status, timestamps, and manifest path.

## Key TSV Ledgers

- `project_dataset_summary.tsv`: dataset-level counts, marker-QC summaries, and visual page counts.
- `project_candidate_overview.tsv`: current candidate-level recommendations and nearest-hit context.
- `project_visual_reports.tsv`: visual report manifests, dashboards, and local-neighbourhood figure manifests.
- `project_visual_gallery.tsv`: the exact PNG previews embedded in the annual HTML.
- `project_isolate_ledger.tsv`: active isolate/project state from the database.
- `project_genome_ledger.tsv`: genome records and genome QC state.
- `decision_changes.tsv`: latest tracked decision deltas between assessment snapshots.

## Interpretation Notes

BranchManager treats marker quality as decision evidence. Low marker evidence is held for review before sequencing selection; high versus moderate marker evidence is also used when otherwise comparable candidates are ranked. Visuals are supporting evidence, not replacement evidence: the TSV ledgers remain the auditable source of the numerical decisions.
"""


def _annual_report_readme_html() -> str:
    return (
        '<!doctype html><html><head><meta charset="utf-8"><title>BranchManager Annual Report Guide</title>'
        '<style>body{font:15px system-ui;margin:32px;color:#17202a;max-width:920px;line-height:1.45}'
        'h1,h2{color:#123b2d}li{margin:6px 0}code{background:#eef4f1;padding:1px 4px;border-radius:4px}</style></head><body>'
        '<h1>BranchManager Annual Report Guide</h1>'
        '<p>This folder is a point-in-time project overview built from the BranchManager database, workflow manifests, and run outputs.</p>'
        '<h2>How To Read The HTML</h2>'
        '<ul>'
        '<li>Summary cards give project-wide counts for baseline markers, active candidate markers, proposed sequencing commitments, marker-QC outcomes, and visual pages.</li>'
        '<li>The Datasets table combines database state with Paper Trail / Interview file-level QC summaries. QC-passed sequence files include accepted files and files that passed with review warnings.</li>'
        '<li>Current recommendations come from the latest Performance Review assessment snapshot. Baseline-redundant candidates are near-identical to cultured baseline markers at the configured identity and coverage thresholds.</li>'
        '<li>Embedded visual previews are a curated subset of the available PNG outputs. Click any preview to open the full-size PNG.</li>'
        '<li>Workflow runs show reproducibility context: workflow name, inferred dataset, status, timestamps, and manifest path.</li>'
        '</ul>'
        '<h2>Key TSV Ledgers</h2>'
        '<ul>'
        '<li><code>project_dataset_summary.tsv</code>: dataset-level counts, marker-QC summaries, and visual page counts.</li>'
        '<li><code>project_candidate_overview.tsv</code>: current candidate-level recommendations and nearest-hit context.</li>'
        '<li><code>project_visual_reports.tsv</code>: visual report manifests, dashboards, and local-neighbourhood figure manifests.</li>'
        '<li><code>project_visual_gallery.tsv</code>: the exact PNG previews embedded in the annual HTML.</li>'
        '<li><code>project_isolate_ledger.tsv</code>, <code>project_genome_ledger.tsv</code>, and <code>decision_changes.tsv</code>: audit ledgers for active isolates, genome records, and decision changes.</li>'
        '</ul>'
        '<h2>Interpretation Notes</h2>'
        '<p>BranchManager treats marker quality as decision evidence. Low marker evidence is held for review before sequencing selection; high versus moderate marker evidence is also used when otherwise comparable candidates are ranked. Visuals are supporting evidence, not replacement evidence: the TSV ledgers remain the auditable source of the numerical decisions.</p>'
        '</body></html>\n'
    )


def _discover_visual_report_rows(project_root: Path | None, outdir: Path, qc_summaries: list[dict]) -> list[tuple[str, str, str, int, str]]:
    rows: list[tuple[str, str, str, int, str]] = []
    seen: set[tuple[str, str]] = set()
    for summary in qc_summaries:
        visual_manifest = str(summary.get('visual_manifest_path') or '').strip()
        if visual_manifest:
            rows.append((
                str(summary.get('dataset') or ''),
                str(summary.get('workflow') or ''),
                'Paper Trail visual pages',
                int(summary.get('visual_pages') or 0),
                visual_manifest,
            ))
            seen.add(('paper_trail_visuals', visual_manifest))
    if project_root is None or not project_root.exists():
        return rows
    output = outdir.expanduser().resolve()
    for path in sorted(project_root.rglob('run_manifest.json')):
        resolved = path.resolve()
        try:
            resolved.relative_to(output)
            continue
        except ValueError:
            pass
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError, TypeError):
            continue
        workflow = str(data.get('workflow') or '').strip()
        dataset = _manifest_dataset(data, path)
        for dashboard in _output_paths(data, 'performance_review_dashboard', path, project_root):
            key = ('performance_review_dashboard', dashboard.as_posix())
            if key in seen:
                continue
            rows.append((dataset, workflow, 'Performance Review dashboard', 1, dashboard.as_posix()))
            seen.add(key)
        neighbourhood_manifest = path.parent / 'assessment' / 'neighbourhoods' / 'neighbourhood_manifest.tsv'
        if neighbourhood_manifest.exists():
            key = ('local_neighbourhoods', neighbourhood_manifest.as_posix())
            if key not in seen:
                manifest_rows = _read_tsv_rows(neighbourhood_manifest)
                rendered = sum(1 for row in manifest_rows if str(row.get('Status') or '').lower() == 'rendered')
                rows.append((dataset, workflow, 'Local neighbourhood figures', rendered, neighbourhood_manifest.as_posix()))
                seen.add(key)
    return sorted(rows, key=lambda row: (str(row[0]), str(row[1]), str(row[2]), str(row[4])))


def _infer_dataset_role(dataset: str, workflow: str = '') -> str:
    name = str(dataset or '').lower()
    if workflow == 'filing_cabinet' or 'hungate' in name:
        return 'baseline'
    return 'candidate'


def _infer_baseline_tier(dataset: str, role: str = '') -> str:
    if str(role or '').lower() != 'baseline':
        return ''
    return 'priority' if 'hungate' in str(dataset or '').lower() else 'secondary'


def _merge_run_rows(primary_rows, manifest_rows):
    rows = [tuple(row) for row in primary_rows]
    seen = {
        (str(row[1] or ''), str(row[2] or ''), str(row[4] or ''), str(row[5] or ''))
        for row in rows
    }
    for row in manifest_rows:
        key = (str(row[1] or ''), str(row[2] or ''), str(row[4] or ''), str(row[5] or ''))
        if key not in seen:
            rows.append(tuple(row))
            seen.add(key)
    return sorted(rows, key=lambda row: (str(row[4] or ''), str(row[0] or '')))


def _merge_dataset_rows(primary_rows, staged_rows, manifest_rows):
    rows = {str(row[0]): tuple(row) for row in primary_rows if str(row[0] or '').strip()}
    for row in staged_rows:
        dataset = str(row[0] or '')
        if dataset and dataset not in rows:
            rows[dataset] = tuple(row)
    for run in manifest_rows:
        dataset = str(run[2] or '')
        if dataset and dataset not in rows:
            role = _infer_dataset_role(dataset, str(run[1] or ''))
            rows[dataset] = (dataset, role, _infer_baseline_tier(dataset, role), 0, 0, 0)
    return sorted(rows.values(), key=lambda row: str(row[0]))


def _merge_status_rows(primary_rows, extra_rows):
    counts = {str(status): int(count or 0) for status, count in primary_rows}
    for status, count in extra_rows:
        key = str(status or 'RECEIVED')
        counts[key] = counts.get(key, 0) + int(count or 0)
    return sorted(counts.items())


def _discover_staged_state(project_root: Path | None, known_datasets: set[str]):
    if project_root is None or not project_root.exists():
        return [], [], []
    dataset_rows = []
    isolate_rows = []
    status_counts = {}
    added = set()
    for path in sorted(project_root.rglob('.*_project.sqlite')):
        try:
            if not path.is_file() or path.stat().st_size == 0:
                continue
            with sqlite3.connect(str(path)) as conn:
                tables = {
                    str(row[0]) for row in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    ).fetchall()
                }
                if 'sequences' not in tables:
                    continue
                rows = _query_dataset_rows(conn)
                missing = [
                    str(row[0]) for row in rows
                    if str(row[0]) not in known_datasets and str(row[0]) not in added
                ]
                if not missing:
                    continue
                dataset_rows.extend(row for row in rows if str(row[0]) in missing)
                isolate_rows.extend(_query_isolate_rows(conn, missing))
                for status, count in _query_status_rows_for_datasets(conn, missing):
                    status_counts[status] = status_counts.get(status, 0) + int(count or 0)
                added.update(missing)
        except (OSError, sqlite3.Error):
            continue
    return dataset_rows, sorted(status_counts.items()), isolate_rows


def write_annual_report(db, outdir: str | Path) -> dict:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    with db.connect() as conn:
        dataset_rows = _query_dataset_rows(conn)
        marker_qc_rows = _query_dataset_marker_qc_rows(conn)
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
        isolate_rows = _query_isolate_rows(conn)
        removal_rows = conn.execute(
            'SELECT sequence_id, original_dataset, partner_id, sequence_length, '
            'taxonomy, reason, removed_at FROM sequence_removals '
            'ORDER BY removed_at, removal_id'
        ).fetchall()

    project_root = _project_root_from_db(db)
    manifest_rows = _discover_run_manifests(project_root, output)
    qc_summaries = _discover_qc_summaries(project_root, output)
    qc_rows_by_dataset = _aggregate_qc_summaries(qc_summaries)
    visual_report_rows = _discover_visual_report_rows(project_root, output, qc_summaries)
    visual_gallery_rows = _build_visual_gallery_rows(visual_report_rows, output)
    staged_dataset_rows, staged_status_rows, staged_isolate_rows = _discover_staged_state(
        project_root,
        {str(row[0]) for row in dataset_rows},
    )
    run_rows = _merge_run_rows(run_rows, manifest_rows)
    dataset_rows = _merge_dataset_rows(dataset_rows, staged_dataset_rows, manifest_rows)
    status_rows = _merge_status_rows(status_rows, staged_status_rows)
    if staged_isolate_rows:
        isolate_rows = sorted(
            [*isolate_rows, *staged_isolate_rows], key=lambda row: str(row[0] or '')
        )

    active = {
        str(row[0]): {
            'dataset': str(row[1] or ''), 'partner_id': str(row[2] or ''),
            'status': str(row[3] or ''), 'marker_qc': str(row[4] or ''),
        }
        for row in isolate_rows
    }
    dataset_roles = db.get_dataset_roles()
    for dataset, role, baseline_tier, genomes_available, _proposed, _sequences in dataset_rows:
        dataset = str(dataset or '')
        if dataset and dataset not in dataset_roles:
            inferred_role = role if role in {'baseline', 'candidate'} else _infer_dataset_role(dataset)
            dataset_roles[dataset] = {
                'role': inferred_role,
                'genomes_available': bool(genomes_available),
                'baseline_tier': baseline_tier or _infer_baseline_tier(dataset, inferred_role),
            }
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

    dataset_summary_headers = (
        'Dataset', 'Role', 'Baseline tier', 'Genomes available', 'Proposed sequences', 'Sequences',
        'Active marker-QC passed', 'Active marker-QC review', 'Active marker-QC failed',
        'QC-passed sequence files', 'Accepted sequence files', 'Review sequence files',
        'QC-failed sequence files', 'Read files', 'Failed read files', 'Visual report pages',
    )
    dataset_summary_rows = []
    for dataset, role, baseline_tier, genomes_available, proposed, sequences in dataset_rows:
        marker = marker_qc_rows.get(str(dataset), {})
        qc = qc_rows_by_dataset.get(str(dataset), {})
        dataset_summary_rows.append((
            dataset, role, baseline_tier, genomes_available, proposed, sequences,
            marker.get('active_marker_qc_passed', 0),
            marker.get('active_marker_qc_review', 0),
            marker.get('active_marker_qc_failed', 0),
            qc.get('passed_sequence_files', 0),
            qc.get('accepted_sequence_files', 0),
            qc.get('review_sequence_files', 0),
            qc.get('failed_sequence_files', 0),
            qc.get('read_files', 0),
            qc.get('failed_read_files', 0),
            qc.get('visual_pages', 0),
        ))
    dataset_summary_path = output / 'project_dataset_summary.tsv'
    with open(dataset_summary_path, 'w') as handle:
        handle.write('\t'.join(dataset_summary_headers) + '\n')
        for row in dataset_summary_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')

    visual_report_path = output / 'project_visual_reports.tsv'
    with open(visual_report_path, 'w') as handle:
        handle.write('Dataset\tWorkflow\tReport\tItems\tPath\n')
        for row in visual_report_rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    visual_gallery_path = output / 'project_visual_gallery.tsv'
    with open(visual_gallery_path, 'w') as handle:
        handle.write('Dataset\tWorkflow\tReport\tCaption\tImagePath\tLinkPath\n')
        for row in visual_gallery_rows:
            handle.write('\t'.join(str(row.get(key, '')) for key in (
                'dataset', 'workflow', 'report', 'caption', 'image_path', 'link_path',
            )) + '\n')
    readme_path = output / 'README.md'
    readme_path.write_text(_annual_report_readme_markdown())
    readme_html_path = output / 'annual_report_guide.html'
    readme_html_path.write_text(_annual_report_readme_html())

    baseline_count = sum(row[5] for row in dataset_rows if row[1] == 'baseline')
    hungate_count = sum(row[5] for row in dataset_rows if row[1] == 'baseline' and row[2] == 'priority')
    secondary_count = sum(row[5] for row in dataset_rows if row[1] == 'baseline' and row[2] == 'secondary')
    candidate_count = sum(row[5] for row in dataset_rows if row[1] == 'candidate')
    proposed_count = sum(row[4] for row in dataset_rows if row[1] == 'candidate')
    qc_passed_sequence_files = sum(row.get('passed_sequence_files', 0) for row in qc_summaries)
    qc_review_sequence_files = sum(row.get('review_sequence_files', 0) for row in qc_summaries)
    qc_failed_sequence_files = sum(row.get('failed_sequence_files', 0) for row in qc_summaries)
    visual_page_count = sum(row.get('visual_pages', 0) for row in qc_summaries)
    cards = [
        ('Cultured baseline markers', baseline_count),
        ('Hungate baseline markers', hungate_count),
        ('Secondary baseline markers', secondary_count),
        ('Active partner candidates', candidate_count),
        ('Proposed sequences', proposed_count),
        ('QC-passed sequence files', qc_passed_sequence_files),
        ('QC-review sequence files', qc_review_sequence_files),
        ('QC-failed sequence files', qc_failed_sequence_files),
        ('Visual report pages', visual_page_count),
        ('Withdrawn sequences', len(removal_rows)),
        ('QC-passed genomes', sum(bool(row[4]) for row in genome_rows)),
        ('Genome QC failures', sum(row[3] == 'GENOME_QC_FAILED' for row in genome_rows)),
        ('Selection rounds', len(round_rows)),
    ]
    def table(headers, rows, *, raw_columns=None):
        raw_columns = set(raw_columns or ())
        head = ''.join(f'<th>{html.escape(str(item))}</th>' for item in headers)

        def cell(value):
            return '' if value is None else str(value)

        body_rows = []
        for row in rows:
            cells = []
            for header, value in zip(headers, row):
                rendered = cell(value)
                if header in raw_columns:
                    cells.append(f'<td>{rendered}</td>')
                else:
                    cells.append(f'<td>{html.escape(rendered)}</td>')
            body_rows.append('<tr>' + ''.join(cells) + '</tr>')
        return f'<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>'

    report = output / 'annual_report.html'
    report.write_text(
        '<!doctype html><html><head><meta charset="utf-8"><title>BranchManager Project Overview</title>'
        '<style>body{font:15px system-ui;margin:32px;color:#17202a;max-width:1200px;line-height:1.35}'
        'h1,h2{color:#123b2d}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}'
        '.card{border:1px solid #cad5cf;padding:14px;border-radius:6px}.card b{font-size:26px;display:block}'
        '.guide{border:1px solid #cad5cf;background:#f6faf8;border-radius:6px;padding:14px 16px;margin:18px 0 24px}'
        '.guide p{margin:6px 0}.muted{color:#596861}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;margin:10px 0 28px}'
        '.visual-card{border:1px solid #cad5cf;border-radius:6px;margin:0;background:white;overflow:hidden}.visual-card img{display:block;width:100%;height:230px;object-fit:contain;background:#f8faf9}'
        '.visual-card figcaption{border-top:1px solid #d9e1dd;padding:9px 10px}.visual-card figcaption b{display:block;color:#123b2d}.visual-card figcaption span{display:block;color:#44534d;font-size:13px;margin-top:3px}'
        'table{border-collapse:collapse;width:100%;margin:10px 0 28px}th,td{border-bottom:1px solid #d9e1dd;text-align:left;padding:7px;vertical-align:top}th{background:#eef4f1}</style>'
        '</head><body><h1>Cumulative Project Overview</h1><p>BranchManager point-in-time marker-to-genome project report.</p>'
        '<div class="cards">' + ''.join(f'<div class="card"><b>{value}</b>{html.escape(label)}</div>' for label, value in cards) + '</div>'
        '<section class="guide"><h2>How To Read This Report</h2>'
        '<p>Use the summary cards for project-level counts, the Datasets table for dataset/QC file outcomes, and Current recommendations for the latest selection decisions.</p>'
        '<p>Embedded visual previews are sampled from the PNG outputs; click any preview to open the full-size image. The TSV ledgers remain the auditable source for numerical decisions.</p>'
        '<p><a href="annual_report_guide.html">Open the full guide</a> or read <a href="README.md">README.md</a>.</p></section>'
        '<h2>Datasets</h2>' + table(dataset_summary_headers, dataset_summary_rows) +
        '<h2>Isolate lifecycle</h2>' + table(('Status', 'Count'), status_rows) +
        '<h2>Current recommendations</h2>' + table(('Recommendation', 'Candidates'), sorted(recommendation_counts.items())) +
        '<h2>Visual reports</h2>' + table(
            ('Dataset', 'Workflow', 'Report', 'Items', 'Link'),
            [
                (dataset, workflow, report_name, items, _html_path_link(path_value, output))
                for dataset, workflow, report_name, items, path_value in visual_report_rows
            ],
            raw_columns={'Link'},
        ) +
        '<h2>Embedded visual previews</h2>' + _visual_gallery_html(visual_gallery_rows) +
        '<h2>Workflow runs</h2>' + table(('Run', 'Workflow', 'Dataset', 'Status', 'Started', 'Completed', 'Manifest'), run_rows) +
        '<h2>Selection rounds</h2>' + table(('Round', 'Mode', 'Created'), round_rows) +
        '<h2>Exit Interviews</h2>' + table(('Sequence', 'Dataset', 'Partner', 'Length', 'Taxonomy', 'Reason', 'Removed'), removal_rows) +
        '<p>Detailed ledgers: <a href="project_dataset_summary.tsv">datasets</a>, '
        '<a href="project_isolate_ledger.tsv">isolates</a>, '
        '<a href="project_candidate_overview.tsv">current candidate assessments</a>, '
        '<a href="project_genome_ledger.tsv">genomes</a>, '
        '<a href="project_visual_reports.tsv">visual reports</a>, '
        '<a href="project_visual_gallery.tsv">visual gallery</a>, '
        '<a href="annual_report_guide.html">report guide</a>, '
        '<a href="sequence_removal_ledger.tsv">Exit Interviews</a>, '
        '<a href="decision_changes.tsv">latest decision changes</a>.</p></body></html>\n'
    )
    return {
        'html': str(report), 'readme': str(readme_path), 'readme_html': str(readme_html_path),
        'datasets': str(dataset_summary_path), 'isolates': str(isolate_path),
        'genomes': str(genome_path), 'candidates': str(candidate_path),
        'visual_reports': str(visual_report_path), 'visual_gallery': str(visual_gallery_path),
        'removals': str(removals_path), 'decision_changes': changes_path,
    }
