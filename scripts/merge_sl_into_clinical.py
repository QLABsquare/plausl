"""Merge plausl_sl_scores.csv and ASD diagnosis into the clinical trials workbook.

Writes a new workbook (the original is left untouched):
  - Master Sheet: adds diagnosis_group and completion_status after record_id, language and reading composites, moves
    values entered in the unlabeled columns BK-BM to kbit_*_2, fills blank assessment cells from the REDCap exports in
    beh/ (adding orr_mode_2/orr_ss_2/orr_perc_2), clears cells listed in the config, swaps standard scores and
    percentiles entered the wrong way round, fills blank TILLS core sums and disorder flags (cut score 24), and replaces
    the label-only SL columns (ssl1_score, ssl1_score_2, ssl2_score, ssl2_score_2) with the re-scored SL measures.
  - SL score key: one row per SL column with a description and the plausl_with_t2.csv column name.
  - Change log: every assessment cell that differs from the original workbook.

Usage: python3 merge_sl_into_clinical.py                                       (build the workbook)
       python3 merge_sl_into_clinical.py --from-shared <participants.csv> [out]  (recompute composites from shared data)
Reads and writes files in the data folder (PLAUSL_DATA, default: this script's folder). Participant-specific
corrections (invalid scores, swaps, completion status, cells to clear) are read from study_config.json there.
"""
import csv
import json
import os
import re
import sys
import zipfile
from xml.sax.saxutils import escape

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('PLAUSL_DATA', HERE)  # data folder (final_plausl_data/)
SRC = os.path.join(DATA, 'PLAUSL Data - Clinical Trials.xlsx')
OUT = os.path.join(DATA, 'PLAUSL Data - Clinical Trials with SL scores.xlsx')
SCORES = os.path.join(DATA, 'plausl_sl_scores.csv')
DEMOGRAPHICS = os.path.join(DATA, 'clinical_trials_with_demographics.csv')
# participant-specific corrections live in a private file in the data folder (study team only)
# (needed to build the workbook; not needed with --from-shared)
_CONFIG_PATH = os.path.join(DATA, 'study_config.json')
CONFIG = json.load(open(_CONFIG_PATH))['merge'] if os.path.exists(_CONFIG_PATH) else {
    'clear_cells': {'record_id': None, 'columns': [], 'reason': ''}, 'invalid_scores': [], 'completed_without_sl_file': [],
    'not_completed': {}, 'swaps': [], 'moved_note': ['', '', '', '']}
DROP = ['ssl1_score', 'ssl1_score_2', 'ssl2_score', 'ssl2_score_2']
REDCAP = os.path.join(os.path.dirname(DATA), 'beh')
REDCAP_FILES = ['plausl_time1_assessment_scores.csv', 'plausl_time2_assessment_scores.csv']
NEW_ASSESSMENT_COLUMNS = ['orr_mode_2', 'orr_ss_2', 'orr_perc_2']
# cells in the original sheet that hold another child's scores: {record_id, columns, reason}
CLEAR_CELLS = CONFIG['clear_cells']
# children in REDCap but not in the original workbook to add as rows
NEW_CHILDREN = CONFIG.get('new_children', [])

# composites: mean of standard scores on a mean-100/SD-15 metric, blank with fewer than 2 components.
# TILLS subtests (mean 10, SD 3) are converted with 100 + 5 * (score - 10).
COMPOSITES = {
    'language_composite': {'tills': ['tills_va_ss', 'tills_pa_ss', 'tills_nwr_ss'], 'standard': ['picvocab_ss', 'trog_ss']},
    'reading_composite': {'tills': [], 'standard': ['wrmt_letterid_ss', 'wrmt_id_ss', 'wrmt_wa_ss', 'orr_ss']},
    # TILLS was not given at time 2, so the session-2 language composite has only these two tests
    'language_composite_2': {'tills': [], 'standard': ['picvocab_ss_2', 'trog_ss_2']},
}
COMPOSITE_MIN_COMPONENTS = 2
# (record_id, column) scores excluded from composites as invalid (out-of-range or impossible values)
INVALID_SCORES = {tuple(x[:2]) for x in CONFIG['invalid_scores']}


def composite(name, rid, value_of):
    spec, parts = COMPOSITES[name], []
    for col in spec['tills'] + spec['standard']:
        try:
            v = float(value_of(col))
        except ValueError:
            continue
        if (rid, col) not in INVALID_SCORES:
            parts.append(100 + 5 * (v - 10) if col in spec['tills'] else v)
    return (round(sum(parts) / len(parts), 2) if len(parts) >= COMPOSITE_MIN_COMPONENTS else None), len(parts)


# completion = any SL data (a row in plausl_sl_scores.csv); otherwise the stage the child's records reach.
# Children who did an SL task whose raw file is missing count as completed.
COMPLETED_WITHOUT_SL_FILE = set(CONFIG['completed_without_sl_file'])
NOT_COMPLETED = CONFIG['not_completed']  # record_id -> 'not completed: <furthest stage>'
# (record_id, standard-score column, percentile column) entered the wrong way round (percentile above 100)
SWAPS = [tuple(x) for x in CONFIG['swaps']]
# TILLS Identification Core for ages 6;0-7;11 = sum of these standard scores; a sum below the cut score is
# classified as a language/literacy disorder (TILLS Technical Manual, Tables 2.2-2.3). Blank core sums are
# filled where all three subtests exist, and the disorder flag is set from the cut score.
TILLS_CORE = ['tills_va_ss', 'tills_pa_ss', 'tills_nwr_ss']
TILLS_CUT = 24
MOVED_NOTE = tuple(CONFIG['moved_note'])  # change-log row for the values moved out of columns BK-BM

# (new column, plausl_sl_scores.csv column, description)
SL_COLUMNS = []
for s, ssl_sfx, total0 in ((1, '', 'correct_syllables_total_0'), (2, '_t2', 'correct_syllables_total_0_2')):
    SL_COLUMNS.append(('ssl_2afc_acc_s%d' % s, 'mean_accuracy' + ssl_sfx, 'Syllable 2AFC proportion correct, session %d (trial 7 dropped)' % s))
    for tag, tname in ((1, 'target'), (0, 'foil')):
        SL_COLUMNS += [
            ('ssl_sicr_serial_%s_s%d' % (tname, s), 'correct_serial_syllables_%d%s' % (tag, ssl_sfx), 'Syllable SICR, %s trials: proportion of syllables in correct serial position, session %d' % (tname, s)),
            ('ssl_sicr_bigram_%s_s%d' % (tname, s), 'correct_bigrams_%d%s' % (tag, ssl_sfx), 'Syllable SICR, %s trials: proportion of key bigrams recalled, session %d' % (tname, s)),
            ('ssl_sicr_trigram_%s_s%d' % (tname, s), 'correct_trigrams_%d%s' % (tag, ssl_sfx), 'Syllable SICR, %s trials: proportion of trigrams in correct position, session %d' % (tname, s)),
            ('ssl_sicr_total_%s_s%d' % (tname, s), total0 if tag == 0 else 'correct_syllables_total_1' + ssl_sfx, 'Syllable SICR, %s trials: proportion of syllables recalled in any position, session %d' % (tname, s)),
        ]
    SL_COLUMNS += [
        ('ssl_sicr_diff_serial_s%d' % s, 'diff_syllables' + ssl_sfx, 'Syllable SICR difference score, serial position: (target - foil) / (target + foil), session %d' % s),
        ('ssl_sicr_diff_bigram_s%d' % s, 'diff_bigrams' + ssl_sfx, 'Syllable SICR difference score, bigrams: (target - foil) / (target + foil), session %d' % s),
    ]
for s in (1, 2):
    SL_COLUMNS.append(('lsl_2afc_acc_s%d' % s, 'afc_accuracy_s%d' % s, 'Letter 2AFC proportion correct, session %d (practice trial 1 excluded)' % s))
    for tag, tname in ((1, 'target'), (0, 'foil')):
        SL_COLUMNS += [
            ('lsl_sicr_serial_%s_s%d' % (tname, s), 'correct_serial_letters.lsl_%d_s%d' % (tag, s), 'Letter SICR, %s trials: proportion of letters in correct serial position, session %d' % (tname, s)),
            ('lsl_sicr_bigram_%s_s%d' % (tname, s), 'correct_bigrams.lsl_%d_s%d' % (tag, s), 'Letter SICR, %s trials: proportion of key bigrams recalled, session %d' % (tname, s)),
            ('lsl_sicr_trigram_%s_s%d' % (tname, s), 'correct_trigrams.lsl_%d_s%d' % (tag, s), 'Letter SICR, %s trials: proportion of trigrams in correct position, session %d' % (tname, s)),
            ('lsl_sicr_total_%s_s%d' % (tname, s), 'correct_letters_total.lsl_%d_s%d' % (tag, s), 'Letter SICR, %s trials: proportion of letters recalled in any position, session %d' % (tname, s)),
        ]


def col_letter(i):
    s = ''
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def col_index(ref):
    n = 0
    for ch in re.match(r'[A-Z]+', ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def norm_id(record_id):
    """plausl_NNN / pl_NNN / Plausl_NNN / plausl_NNN_M -> 'NNN'; X_ prefix kept."""
    i = record_id.strip().lower()
    excluded = i.startswith('x_')
    i = re.sub(r'^x_', '', i)
    i = re.sub(r'^(plausl|pl)_', '', i)
    i = re.sub(r'_[mt]$', '', i)
    return ('x_' if excluded else '') + i


def read_redcap():
    """Score fields per child from the REDCap exports: {norm_id: {column: (value, file)}}."""
    skip = re.compile(r'dateadmin|examiner|notes|^ssl|^sl[12]_|^bosa|^mri|^osrt|^psrt|^categorical|complete|'
                      r'^redcap|^record_id|scor(er|ing)|second_score|^age_|recordform')
    out = {}
    for name in REDCAP_FILES:
        with open(os.path.join(REDCAP, name), encoding='utf-8-sig', errors='replace') as f:
            for row in csv.DictReader(f):
                child = out.setdefault(norm_id(row['record_id']), {})
                for col, value in row.items():
                    value = (value or '').strip()
                    if skip.search(col) or not value:
                        continue
                    assert child.get(col, (value,))[0] == value, (row['record_id'], col)  # repeat rows agree
                    child.setdefault(col, (value, name))
    return out


def simple_sheet(rows, widths):
    cols = ''.join('<col width="%s" customWidth="1" min="%d" max="%d"/>' % (w, i + 1, i + 1) for i, w in enumerate(widths))
    data = ''.join('<row r="%d">%s</row>' % (r + 1, ''.join(str_cell('%s%d' % (col_letter(c), r + 1), 19 if r == 0 else 1, str(v))
                                                           for c, v in enumerate(row) if v != ''))
                   for r, row in enumerate(rows))
    return ('<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><cols>%s</cols>'
            '<sheetData>%s</sheetData></worksheet>' % (cols, data))


def str_cell(ref, style, text):
    return '<c r="%s" s="%s" t="inlineStr"><is><t>%s</t></is></c>' % (ref, style, escape(text))


def num_cell(ref, style, value):
    return '<c r="%s" s="%s" t="n"><v>%s</v></c>' % (ref, style, value)


def main():
    if not os.path.exists(_CONFIG_PATH):
        sys.exit('study_config.json not found in %s; it is required to build the workbook '
                 '(use --from-shared to recompute composites from a shared participant file)' % DATA)
    with open(SCORES) as f:
        scores = {r['id']: r for r in csv.DictReader(f)}
    with open(DEMOGRAPHICS, encoding='utf-8-sig') as f:
        demog = {r['record_id']: r for r in csv.DictReader(f)}

    zin = zipfile.ZipFile(SRC)
    sheet = zin.read('xl/worksheets/sheet1.xml').decode('utf-8')
    rows = re.findall(r'<row [^>]*?(?:/>|>.*?</row>)', sheet, re.S)
    parsed = []
    for row in rows:
        head = re.match(r'<row [^>]*?(?=/?>)', row).group()
        cells = {col_index(m.group(1)): m.group(0)
                 for m in re.finditer(r'<c r="([A-Z]+)\d+"[^>]*?(?:/>|>.*?</c>)', row, re.S)}
        parsed.append((int(re.search(r'r="(\d+)"', head).group(1)), head, cells))

    def text(cell):
        m = re.search(r'<t[^>]*>(.*?)</t>|<v>(.*?)</v>', cell or '', re.S)
        return (m.group(1) or m.group(2) or '') if m else ''

    # cell styles with a date number format; filled numbers must not inherit them (116 would show as 4/25/1900)
    styles = zin.read('xl/styles.xml').decode('utf-8')
    custom = dict(re.findall(r'<numFmt numFmtId="(\d+)" formatCode="([^"]*)"', styles))
    xfs = re.findall(r'<xf [^>]*?numFmtId="(\d+)"', re.search(r'<cellXfs[^>]*>(.*?)</cellXfs>', styles, re.S).group(1))
    date_styles = {str(i) for i, n in enumerate(xfs) if 14 <= int(n) <= 22 or int(n) in (45, 46, 47)
                   or re.search(r'[dmy]', re.sub(r'"[^"]*"|\[[^\]]*\]', '', custom.get(n, '')), re.I)}

    def number_style(cells, i):
        style = re.search(r's="(\d+)"', cells[i]).group(1) if i in cells else '1'
        return '1' if style in date_styles else style

    header = {i: text(c) for i, c in parsed[0][2].items()}
    by_name = {h: i for i, h in header.items()}

    # a second KBIT entered in unlabeled columns BK-BM; move it under the kbit_*_2 headers
    for stray, target in (('BK', 'kbit_mode_2'), ('BL', 'kbit_raw_2'), ('BM', 'kbit_ss_2')):
        s_i, t_i = col_index(stray), by_name[target]
        for _, _, cells in parsed[1:]:
            if text(cells.get(s_i)):
                assert not text(cells.get(t_i)), (stray, target)
                style = re.search(r's="(\d+)"', cells.get(t_i) or cells[s_i]).group(1)
                cells[t_i] = num_cell('A1', style, text(cells[s_i]))
    unlabeled = [i for i in range(max(max(c) for _, _, c in parsed if c) + 1) if not header.get(i)]

    drop_idx = [i for i, h in header.items() if h in DROP] + unlabeled
    assert len(drop_idx) - len(unlabeled) == len(DROP), header
    insert_at = min(drop_idx)
    n_cols = max(max(c) for _, _, c in parsed if c) + 1
    record_col = 0

    # session-2 oral reading columns exist only in the REDCap export; add them after orr_perc
    header_style = re.search(r's="(\d+)"', parsed[0][2][by_name['orr_perc']]).group(1)
    for k, name in enumerate(NEW_ASSESSMENT_COLUMNS):
        header[n_cols + k] = name
        by_name[name] = n_cols + k
        parsed[0][2][n_cols + k] = str_cell('A1', header_style, name)
    new_after = by_name['orr_perc']

    changes = []
    rows_by_id = {text(cells.get(record_col)): cells for _, _, cells in parsed[1:] if text(cells.get(record_col))}

    # scores entered on the wrong child's row (they duplicate another child's REDCap record); clear them here,
    # the right child's row is filled from REDCap below
    for name in CLEAR_CELLS['columns']:
        cells = rows_by_id[CLEAR_CELLS['record_id']]
        old = text(cells.get(by_name[name]))
        if old:
            del cells[by_name[name]]
            changes.append((CLEAR_CELLS['record_id'], name, old, '', CLEAR_CELLS['reason']))

    # children in REDCap but not in the workbook go into the first empty rows after the data
    id_style = re.search(r's="(\d+)"', parsed[1][2][record_col]).group(1)
    empty_rows = [cells for rnum, _, cells in parsed[1:] if not text(cells.get(record_col))]
    for rid, cells in zip(NEW_CHILDREN, empty_rows):
        cells.clear()
        cells[record_col] = str_cell('A1', id_style, rid)
        rows_by_id[rid] = cells
        changes.append((rid, 'record_id', '', rid, 'row added: child is in the REDCap export but not in the original sheet'))

    # fill blank cells from the REDCap exports (no conflicting values exist; asserted)
    redcap = read_redcap()
    for rid, cells in rows_by_id.items():
        for name, (value, source) in sorted(redcap.get(norm_id(rid), {}).items()):
            if name not in by_name:
                continue
            old = text(cells.get(by_name[name]))
            if old in ('', 'NA', 'N/A'):
                if value in ('NA', 'N/A') or value == old:
                    continue
                reason = 'filled from beh/' + source
                pct = re.fullmatch(r'(\d+(?:\.\d+)?)%', value)
                if pct:
                    value = '%g' % float(pct.group(1))
                    reason += " (converted from '%s%%')" % pct.group(1)
                style = re.search(r's="(\d+)"', cells[by_name[name]]).group(1) if by_name[name] in cells else '1'
                is_num = re.fullmatch(r'-?\d+(\.\d+)?', value)
                cells[by_name[name]] = (num_cell('A1', number_style(cells, by_name[name]), value) if is_num
                                        else str_cell('A1', style, value))
                changes.append((rid, name, old, value, reason))
            else:
                try:
                    same = abs(float(old) - float(value)) < 1e-9
                except ValueError:
                    same = old.lower() == value.lower()
                assert same, (rid, name, old, value)

    for rid, a, b in SWAPS:
        cells = rows_by_id[rid]
        va, vb = text(cells.get(by_name[a])), text(cells.get(by_name[b]))
        cells[by_name[a]], cells[by_name[b]] = (re.sub(r'<v>.*?</v>', '<v>%s</v>' % v, cells[by_name[n]])
                                               for v, n in ((vb, a), (va, b)))
        changes += [(rid, a, va, vb, 'swapped with %s: REDCap value was a percentile above 100' % b),
                    (rid, b, vb, va, 'swapped with %s' % a)]

    def set_num(cells, name, value):
        cells[by_name[name]] = num_cell('A1', number_style(cells, by_name[name]), value)

    def value(cells, name):
        v = text(cells.get(by_name[name]))
        return '' if v in ('NA', 'N/A') else v

    for rid, cells in rows_by_id.items():
        parts = [value(cells, c) for c in TILLS_CORE]
        core, flag = value(cells, 'tills_coresum'), value(cells, 'tills_languagedisorder')
        if not core and all(parts):
            core = '%g' % sum(float(p) for p in parts)
            old = text(cells.get(by_name['tills_coresum']))
            set_num(cells, 'tills_coresum', core)
            changes.append((rid, 'tills_coresum', old, core, 'computed: sum of %s standard scores' % ', '.join(TILLS_CORE)))
        if not core:
            continue
        status = '%d' % (float(core) < TILLS_CUT)
        if flag:
            assert float(flag) == float(status), (rid, core, flag)  # existing flags already follow the cut of 24
        else:
            old = text(cells.get(by_name['tills_languagedisorder']))
            set_num(cells, 'tills_languagedisorder', status)
            changes.append((rid, 'tills_languagedisorder', old, status,
                            'computed: tills_coresum %s %s TILLS cut score of %d (ages 6;0-7;11)'
                            % (core, '<' if status == '1' else '>=', TILLS_CUT)))

    # new column layout: record_id, diagnosis_group, ..., SL scores where the label columns were, ...
    layout = []
    for i in range(n_cols):
        if i == insert_at:
            layout += [('sl', name) for name, _, _ in SL_COLUMNS]
        if i in drop_idx:
            continue
        layout.append(('orig', i))
        if i == record_col:
            layout += [('diag', None), ('status', None)]
        if i == new_after:
            layout += [('orig', by_name[name]) for name in NEW_ASSESSMENT_COLUMNS]
            for comp in COMPOSITES:
                layout += [('comp', comp), ('comp', comp + '_n')]
    sl_source = {name: src for name, src, _ in SL_COLUMNS}

    header_style_sl = re.search(r's="(\d+)"', parsed[0][2][insert_at]).group(1)
    header_style_diag = re.search(r's="(\d+)"', parsed[0][2][1]).group(1)
    unmatched = set(scores)
    out_rows = []
    for rnum, head, cells in parsed:
        rid = text(cells.get(record_col))
        sl_id = demog[rid]['demog_record_id'] if rid in demog else None
        sl = scores.get(sl_id) if sl_id else None
        if sl:
            unmatched.discard(sl_id)
        out = []
        for new_i, (kind, what) in enumerate(layout):
            ref = '%s%d' % (col_letter(new_i), rnum)
            if kind == 'orig':
                if what in cells:
                    out.append(re.sub(r'r="[A-Z]+\d+"', 'r="%s"' % ref, cells[what], count=1))
            elif rnum == 1:
                out.append(str_cell(ref, header_style_diag if kind in ('diag', 'status') else header_style_sl,
                                    {'diag': 'diagnosis_group', 'status': 'completion_status'}.get(kind, what)))
            elif kind == 'comp' and rid:
                score, n = composite(what[:-2] if what.endswith('_n') else what, rid,
                                     lambda name: text(cells.get(by_name[name])))
                if score is not None:
                    out.append(num_cell(ref, 1, n if what.endswith('_n') else score))
            elif kind == 'status' and rid:
                status = 'completed' if sl or rid in COMPLETED_WITHOUT_SL_FILE else NOT_COMPLETED.get(rid)
                assert status and not (sl and rid in NOT_COMPLETED), rid
                out.append(str_cell(ref, 1, status))
            elif kind == 'diag' and rid in demog:
                dx = {'Yes': 'ASD', 'No': 'non-ASD'}.get(demog[rid]['asd_diagnosis'].strip())
                if dx:
                    out.append(str_cell(ref, 1, dx))
            elif kind == 'sl' and sl and sl.get(sl_source[what], '') != '':
                out.append(num_cell(ref, 1, round(float(sl[sl_source[what]]), 4)))
        if cells or out:
            out_rows.append('%s>%s</row>' % (head, ''.join(out)) if out else head + '/>')
        else:
            out_rows.append(head + '/>')
    assert not unmatched, unmatched

    # column widths/styles follow the original columns; new columns get a uniform width
    old_cols = {}
    for m in re.finditer(r'<col ([^>]*)/>', sheet):
        a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        for i in range(int(a['min']) - 1, min(int(a['max']), n_cols + 1)):
            old_cols[i] = (a['width'], a.get('style', '0'))
    col_xml = []
    for new_i, (kind, what) in enumerate(layout):
        width, style = old_cols.get(what, ('15.1640625', '2')) if kind == 'orig' else (('26' if kind == 'sl' else '16'), '2')
        col_xml.append('<col width="%s" customWidth="1" style="%s" min="%d" max="%d"/>' % (width, style, new_i + 1, new_i + 1))

    new_sheet = re.sub(r'<cols>.*?</cols>', '<cols>%s</cols>' % ''.join(col_xml), sheet, flags=re.S)
    new_sheet = re.sub(r'<sheetData>.*</sheetData>', lambda _: '<sheetData>%s</sheetData>' % ''.join(out_rows), new_sheet, flags=re.S)
    new_sheet = re.sub(r'<dimension ref="[^"]*"/>', '<dimension ref="A1:%s%d"/>' % (col_letter(len(layout) - 1), parsed[-1][0]), new_sheet)

    # key sheet
    key_rows = [('column', 'description', 'column in plausl_sl_scores.csv / plausl_with_t2.csv')]
    key_rows.append(('diagnosis_group', 'ASD / non-ASD, from asd_diagnosis in clinical_trials_with_demographics.csv', ''))
    notes = CONFIG.get('key_notes', {})  # study-specific wording (names the excluded scores)
    key_rows.append(('language_composite', notes.get('language_composite',
                     'mean of the available standard scores among TILLS vocabulary, phonemic awareness and nonword '
                     'repetition (converted from mean 10 / SD 3 to mean 100 / SD 15: 100 + 5 x (score - 10)), NIH Toolbox '
                     'Picture Vocabulary and TROG-2 (session 1); blank if fewer than 2 scores; scores listed as invalid '
                     'in the config are excluded. Higher = better language.'), ''))
    key_rows.append(('language_composite_n', 'number of scores averaged into language_composite', ''))
    key_rows.append(('language_composite_2', 'session 2: mean of NIH Toolbox Picture Vocabulary and TROG-2 standard '
                     'scores (TILLS was not given at session 2); blank unless both are present. Higher = better language.', ''))
    key_rows.append(('language_composite_2_n', 'number of scores averaged into language_composite_2', ''))
    key_rows.append(('reading_composite', notes.get('reading_composite',
                     'mean of the available standard scores among WRMT-III letter identification, word identification '
                     'and word attack and NIH Toolbox Oral Reading Recognition (session 1; all mean 100 / SD 15); blank '
                     'if fewer than 2 scores; scores listed as invalid in the config are excluded. Higher = better '
                     'reading.'), ''))
    key_rows.append(('reading_composite_n', 'number of scores averaged into reading_composite', ''))
    key_rows.append(('tills_coresum / tills_languagedisorder', 'TILLS Identification Core for ages 6;0-7;11: sum of '
                     'vocabulary awareness, phonemic awareness and nonword repetition standard scores; '
                     'languagedisorder = 1 if the sum is below the cut score of %d (TILLS Technical Manual, '
                     'Tables 2.2-2.3). Blank sums were computed where all three subtests exist.' % TILLS_CUT, ''))
    key_rows.append(('completion_status', notes.get('completion_status',
                     'completed = has any SL data; otherwise the furthest stage in the records (REDCap assessments, '
                     'surveys, consent forms, Gorilla). Reasons for stopping are not recorded. Children who did an SL '
                     'task whose raw file is missing count as completed; their SL columns are blank.'), ''))
    key_rows += SL_COLUMNS
    key_rows.append(('', 'Scores from score_plausl_sl.py; sources and flags per child in plausl_sl_sources.csv', ''))
    change_rows = [('record_id', 'column', 'original value', 'new value', 'reason')]
    change_rows += [(MOVED_NOTE[0], MOVED_NOTE[1], '', MOVED_NOTE[2], MOVED_NOTE[3])]
    change_rows += changes
    change_rows += [('', ', '.join(DROP), '', '', 'removed: SL task-order labels only, replaced by SL score columns')]

    workbook = zin.read('xl/workbook.xml').decode()
    rels = zin.read('xl/_rels/workbook.xml.rels').decode()
    types = zin.read('[Content_Types].xml').decode()
    extra = {}
    for title, part, rows_, widths in (('SL score key', 'sheet_slkey.xml', key_rows, (30, 95, 40)),
                                       ('Change log', 'sheet_changes.xml', change_rows, (24, 26, 16, 16, 70))):
        rid = 'rId%d' % (max(int(x) for x in re.findall(r'Id="rId(\d+)"', rels)) + 1)
        sheet_id = max(int(x) for x in re.findall(r'sheetId="(\d+)"', workbook)) + 1
        workbook = workbook.replace('</sheets>', '<sheet xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
                                    'name="%s" sheetId="%d" state="visible" r:id="%s"/></sheets>' % (title, sheet_id, rid))
        rels = rels.replace('</Relationships>', '<Relationship Id="%s" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
                            'relationships/worksheet" Target="worksheets/%s"/></Relationships>' % (rid, part))
        types = types.replace('</Types>', '<Override PartName="/xl/worksheets/%s" ContentType="application/'
                              'vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>' % part)
        extra['xl/worksheets/' + part] = simple_sheet(rows_, widths)

    with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            data = {'xl/worksheets/sheet1.xml': new_sheet.encode(), 'xl/workbook.xml': workbook.encode(),
                    'xl/_rels/workbook.xml.rels': rels.encode(), '[Content_Types].xml': types.encode()}.get(item.filename, data)
            zout.writestr(item, data)
        for part, xml in extra.items():
            zout.writestr(part, xml.encode())
    return changes


# session-2 reading composite; computed for the shared dataset (the workbook does not include it)
READING_2 = {'tills': [], 'standard': ['wrmt_letterid_ss_2', 'wrmt_id_ss_2', 'wrmt_wa_ss_2', 'orr_ss_2']}


def from_shared(path, out):
    """Recompute the derived assessment columns (composites, TILLS core sum and disorder status) from the standard
    scores in a shared participant file (plausl_participants.csv), write them to `out`, and report any column whose
    recomputed values differ from the file's."""
    with open(path, encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    specs = dict(COMPOSITES, reading_composite_2=READING_2)
    derived, diffs = [], {}
    for r in rows:
        value_of = lambda col: r.get(col, '')
        d = {'participant': r['participant']}
        for name, spec in specs.items():
            COMPOSITES[name] = spec  # composite() looks names up in COMPOSITES
            score, n = composite(name, r['participant'], value_of)
            d[name], d[name + '_n'] = ('' if score is None else score), (n if score is not None else '')
        parts = [r.get(c, '') for c in TILLS_CORE]
        core = sum(float(p) for p in parts) if all(parts) else None
        d['tills_coresum'] = '' if core is None else '%g' % core
        d['tills_languagedisorder'] = '' if core is None else int(core < TILLS_CUT)
        for col, v in d.items():
            if col == 'participant' or col not in r:
                continue
            a, b = r[col], v
            same = (a == '' and b == '') or (a != '' and b != '' and abs(float(a) - float(b)) < 0.006)
            if not same and not (col == 'tills_languagedisorder' and b == '' and a != ''):
                diffs.setdefault(col, []).append((r['participant'], a, b))
        derived.append(d)
    COMPOSITES.pop('reading_composite_2', None)
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(derived[0]))
        w.writeheader()
        w.writerows(derived)
    print('recomputed %d derived columns for %d participants; wrote %s' % (len(derived[0]) - 1, len(rows), out))
    for col, ds in diffs.items():
        print('  %s differs for %d participant(s): %s' % (col, len(ds), ds[:5]))
    if not diffs:
        print('  all recomputed values match the shared file')


if __name__ == '__main__':
    if '--from-shared' in sys.argv:
        args = sys.argv[sys.argv.index('--from-shared') + 1:]
        if not args:
            sys.exit('usage: python3 merge_sl_into_clinical.py --from-shared <plausl_participants.csv> [output.csv]')
        from_shared(args[0], args[1] if len(args) > 1 else 'derived_scores.csv')
    else:
        print('%d assessment cells changed; see the Change log sheet' % len(main()))
