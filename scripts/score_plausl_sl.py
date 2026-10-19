"""Re-score PLAUSL statistical learning (SSL + LSL; 2AFC + SICR) from raw Gorilla exports.

Follows the scoring rules in all behavioral/plausl_sl_analysis.R, with these corrections:
  - 2AFC answer keys are derived from the stimuli shown on each trial (the option that is a
    familiarization "word" is correct) instead of the hard-coded R vectors, which are misaligned.
  - LSL SICR keys are the letter sequences actually displayed on each trial (R's default_key list
    is out of order for trials 1-4 and uses GJA where version E displayed JGA).
  - Syllable 2AFC trial 7 is dropped (no word option in the older task version).
  - 2AFC: first click per trial is scored (R counted every distinct click); timeouts excluded.
  - Sessions are assigned by date: a child's runs more than 45 days apart start a new session.
    Within a session the earliest run with data for a measure is used, except for RUN_OVERRIDES.

Usage: python3 score_plausl_sl.py                                  (full pipeline from raw exports)
       python3 score_plausl_sl.py --from-trials <trials.csv> [folder] (rescore a shared trial-level file)
Writes plausl_sl_scores.csv, plausl_sl_summary.csv, plausl_sl_sources.csv and
plausl_sl_scores_by_trial.csv to the data folder (PLAUSL_DATA, default: this script's folder).
Participant-specific corrections are read from study_config.json in the data folder.
"""
import csv
import datetime
import glob
import itertools
import json
import os
import re
import statistics
import sys
import xml.etree.ElementTree as ET
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
# data folder (final_plausl_data/); defaults to this script's folder. Outputs are written there too.
DATA = os.environ.get('PLAUSL_DATA', HERE)
RAW_DIRS = [os.path.join(DATA, 'csv_files'),
            os.path.join(os.path.dirname(DATA), 'beh', 'SL data - Copy')]  # second = backups
TRANSCRIPTIONS = os.path.join(DATA, 'sicr_recordings', 'plausl_sicr_transcriptions.xlsx')
# participant-specific corrections live in a private file in the data folder (study team only)
# (needed for the full pipeline from raw exports; not needed with --from-trials)
_CONFIG_PATH = os.path.join(DATA, 'study_config.json')
CONFIG = json.load(open(_CONFIG_PATH))['score'] if os.path.exists(_CONFIG_PATH) else {'run_overrides': [], 'transcription_sessions': {}}
SESSION_GAP_DAYS = 45
# (child, session, measure) -> Gorilla private ID of the run to use instead of the earliest one
RUN_OVERRIDES = {(o['child'], o['session'], o['measure']): o['gorilla_private_id'] for o in CONFIG['run_overrides']}

LSL_WORDS = {'GJA', 'MDH', 'FKC', 'LBE'}
SSL_WORDS = {'da-ro-pi', 'go-la-tu', 'pa-bi-ku', 'ti-bu-do'}
SSL_AFC_DROP_TRIALS = {7}
LSL_AFC_PRACTICE_TRIALS = {1}  # trial 1 follows the practice/instruction screen; trial 17 repeats it

# SSL SICR target/foil tags, from ssl_unscored() in the R script
SSL_TARGET_TAG = {
    'da/ro/pi/go/la/tu': 1, 'pa/bi/ku/da/ro/pi': 1, 'go/tu/ro/pi/la/da': 0, 'tu/ro/la/da/pi/go': 0,
    'bu/da/pi/ti/ro/do': 0, 'la/bi/go/ku/tu/pa': 0, 'do/pa/ti/bu/ku/bi': 0, 'pa/bi/ku/go/la/tu': 1,
    'ti/bu/do/pa/bi/ku': 1, 'tu/bu/do/pa/bi/ku': 1, 'bu/go/la/do/ti/tu': 0, 'da/ro/pi/ti/bu/do': 1,
    'pa/bi/ku/ti/bu/do': 1, 'go/la/tu/ti/bu/do': 1, 'do/bu/ku/bi/pa/ti': 0, 'go/la/tu/da/ro/pi': 1,
    'da/bi/pa/ku/pi/ro': 0,
}

# Transcription label -> session, decided from recording file names (Gorilla private ID -> run date)
# where the recordings exist, otherwise from which sessions contain an SSL SICR run.
TRANSCRIPTION_SESSIONS = {k: tuple(v) for k, v in CONFIG['transcription_sessions'].items()}


def G(row, key):
    return (row.get(key) or '').strip()


def parse_date(s):
    m = re.match(r'(\d{1,2})/(\d{1,2})/(\d{4})', s.strip())
    return datetime.date(int(m.group(3)), int(m.group(2)), int(m.group(1))) if m else None


def child_from_file(path):
    m = re.match(r'plausl_(spark_\d+|\d+)', os.path.basename(path))
    return 'pl_' + m.group(1)


def sheet_value(row, name):
    return G(row, 'Spreadsheet: ' + name) or G(row, name)


# ---------------------------------------------------------------- R scorer (shared by SSL and LSL)
def permutations(tokens):
    if len(tokens) == 6:
        return [tokens]
    if len(tokens) > 6:
        return [tokens[:6]]
    need = 6 - len(tokens)
    out = []
    for comb in itertools.combinations(range(1, len(tokens) + need + 1), need):
        tmp = list(tokens)
        for pos in comb:
            tmp.insert(pos - 1, 'xx')
        out.append(tmp)
    return out


def score_tokens(resp, key):
    serial = sum(a == b for a, b in zip(resp, key))
    resp_bigrams = {resp[i] + '/' + resp[i + 1] for i in range(len(resp) - 1)}
    key_bigrams = [key[0] + '/' + key[1], key[1] + '/' + key[2], key[3] + '/' + key[4], key[4] + '/' + key[5]]
    bigrams = sum(b in resp_bigrams for b in key_bigrams)
    trigrams = int(resp[0:3] == key[0:3]) + int(resp[3:6] == key[3:6])
    total = sum(t in key for t in resp)
    return serial, bigrams, trigrams, total


def best_score(resp_tokens, key_tokens):
    best = None
    for p in permutations(resp_tokens):
        s = score_tokens(p, key_tokens)
        if best is None or s[0] > best[0]:
            best = s
    return best


def sicr_means(trials, unit):
    """trials: list of dicts with tag, counted, serial, bigrams, trigrams, total. R's >2-responses rule."""
    out = {}
    for tag in (0, 1):
        rows = [t for t in trials if t['tag'] == tag]
        n = sum(t['counted'] for t in rows)
        if n <= 2:
            return None
        out[tag] = {'serial': sum(t['serial'] for t in rows) / (6 * n),
                    'total': sum(t['total'] for t in rows) / (6 * n),
                    'bigrams': sum(t['bigrams'] for t in rows) / (4 * n),
                    'trigrams': sum(t['trigrams'] for t in rows) / (2 * n),
                    'n': n}
    return out


# ---------------------------------------------------------------- raw Gorilla runs
def read_runs():
    runs, primary = {}, set()
    for d_index, d in enumerate(RAW_DIRS):
        for path in sorted(glob.glob(os.path.join(d, '*.csv'))):
            with open(path, encoding='utf-8-sig', errors='replace') as f:
                rows = list(csv.DictReader(f))
            groups = {}
            for x in rows:
                task = G(x, 'Task Name')
                kind = 'LSL' if task.startswith('LSL') else 'SSL' if task.startswith('SSL') and 'SICR' not in task and 'end' not in task else None
                if kind:
                    groups.setdefault((G(x, 'Participant Private ID'), kind, task), []).append(x)
            for (priv, kind, task), xs in groups.items():
                key = (priv, kind, task)
                if key in runs or (d_index == 1 and (priv, kind) in primary):
                    continue  # backups duplicate csv_files runs under different task names
                if d_index == 0:
                    primary.add((priv, kind))
                dates = sorted(filter(None, (parse_date(G(x, 'Local Date and Time')) for x in xs)))
                runs[key] = {'child': child_from_file(path), 'kind': kind, 'task': task, 'priv': priv,
                             'date': dates[0] if dates else None,
                             'file': os.path.relpath(path, os.path.dirname(DATA)), 'rows': xs,
                             'backup_only': d_index == 1}
    return list(runs.values())


def lsl_afc(rows):
    options, clicks = {}, {}
    for x in rows:
        t = G(x, 'Trial Number')
        if not t.isdigit():
            continue
        t = int(t)
        if G(x, 'Screen') == '2AFC':
            letters = ''.join(sheet_value(x, 'sicr%d' % i).split('.')[0] for i in range(1, 7))
            if len(letters) == 6:
                options[t] = (letters[:3], letters[3:])
        elif G(x, 'Screen') == '2afc choices' and G(x, 'Response Type') == 'response' and G(x, 'Response').startswith('click'):
            clicks.setdefault(t, G(x, 'Response'))
    trials = []
    for t, (o1, o2) in sorted(options.items()):
        if t in LSL_AFC_PRACTICE_TRIALS:
            continue
        correct = 1 if o1 in LSL_WORDS and o2 not in LSL_WORDS else 2 if o2 in LSL_WORDS and o1 not in LSL_WORDS else None
        if correct is None or t not in clicks:
            continue
        resp = int(clicks[t][-1])
        trials.append({'trial': t, 'option1': o1, 'option2': o2, 'correct_answer': correct,
                       'response': resp, 'score': int(resp == correct)})
    return trials


def ssl_afc(rows):
    options, clicks = {}, {}
    for x in rows:
        t = G(x, 'Trial Number')
        if G(x, 'Screen') != '2afc' or not t.isdigit():
            continue
        t = int(t)
        a = [sheet_value(x, 'audio' + s).replace('.wav', '') for s in ('', '2', '3', '4', '5', '6')]
        if all(a):
            options.setdefault(t, ('-'.join(a[:3]), '-'.join(a[3:])))
        if G(x, 'Response Type') == 'response' and G(x, 'Response').startswith('click'):
            clicks.setdefault(t, G(x, 'Response'))
    trials = []
    for t, (o1, o2) in sorted(options.items()):
        if t in SSL_AFC_DROP_TRIALS:
            continue
        correct = 1 if o1 in SSL_WORDS and o2 not in SSL_WORDS else 2 if o2 in SSL_WORDS and o1 not in SSL_WORDS else None
        if correct is None or t not in clicks:
            continue
        resp = int(clicks[t][-1])
        trials.append({'trial': t, 'option1': o1, 'option2': o2, 'correct_answer': correct,
                       'response': resp, 'score': int(resp == correct)})
    return trials


def lsl_sicr(rows):
    """Response on trial t answers the letters displayed on trial t-1; trial 1's display is a warm-up."""
    boxes, responses = {}, {}
    for x in rows:
        t = G(x, 'Trial Number')
        if not t.isdigit():
            continue
        t = int(t)
        if G(x, 'Screen') == 'sicr boxes' and G(x, 'Response Type') == 'timedOut':
            letters = ''.join(sheet_value(x, 'sicr%d' % i).split('.')[0] for i in range(1, 7)).upper()
            if len(letters) == 6:
                boxes[t] = letters
        elif G(x, 'Screen') == 'sicr practice' and G(x, 'Response Type') == 'response':
            responses.setdefault(t, re.sub(r'\s', '', G(x, 'Response')).upper())
    trials = []
    for t, key in sorted(boxes.items()):
        if t < 2 or t + 1 not in responses:
            continue
        resp = responses[t + 1]
        w1, w2 = key[:3].replace('JGA', 'GJA'), key[3:].replace('JGA', 'GJA')
        tag = int(w1 in LSL_WORDS and w2 in LSL_WORDS)
        tokens = list(resp) if resp else ['X', 'X']
        serial, bigrams, trigrams, total = best_score(tokens, list(key))
        trials.append({'trial': t - 1, 'key': key, 'response': resp, 'tag': tag, 'counted': 1,
                       'serial': serial, 'bigrams': bigrams, 'trigrams': trigrams, 'total': total})
    return trials


# ---------------------------------------------------------------- SSL SICR transcriptions
def read_transcriptions():
    z = zipfile.ZipFile(TRANSCRIPTIONS)
    ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
    strings = [''.join(t.text or '' for t in si.iter(ns + 't'))
               for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall(ns + 'si')]

    def col(ref):
        n = 0
        for ch in re.match(r'[A-Z]+', ref).group():
            n = n * 26 + ord(ch) - 64
        return n - 1

    keys, responses = {}, {}
    for row in ET.fromstring(z.read('xl/worksheets/sheet2.xml')).iter(ns + 'row'):  # 'Sheet9': key/response rows
        cells = {}
        for c in row.findall(ns + 'c'):
            v = c.find(ns + 'v')
            val = '' if v is None else (strings[int(v.text)] if c.get('t') == 's' else v.text)
            cells[col(c.get('r'))] = val.strip()
        label, kind = cells.get(0, ''), cells.get(1, '')
        if label.startswith('plausl_') and kind in ('key', 'response'):
            (keys if kind == 'key' else responses)[label] = [cells.get(i, '') for i in range(2, 18)]
    return keys, responses


def ssl_sicr(key_list, resp_list):
    trials = []
    for i, (key, resp) in enumerate(zip(key_list, resp_list), 1):
        tokens = resp.split('/')
        serial, bigrams, trigrams, total = best_score(tokens, key.split('/'))
        trials.append({'trial': i, 'key': key, 'response': resp, 'tag': SSL_TARGET_TAG[key],
                       'counted': int(resp not in ('', 'N/A')),
                       'serial': serial, 'bigrams': bigrams, 'trigrams': trigrams, 'total': total})
    return trials


def wide_table(chosen):
    """One row per child from {(child, session, measure): trials}; column names as in plausl_with_t2.csv."""
    def afc(trials):
        return round(sum(t['score'] for t in trials) / len(trials), 6) if trials else ''

    children = sorted({c for c, _, _ in chosen})
    wide = []
    for c in children:
        row = {'id': c}
        for s, sfx_ssl in ((1, ''), (2, '_t2')):
            row['mean_accuracy' + ('_t2' if s == 2 else '')] = afc(chosen.get((c, s, 'ssl_afc')))
            m = sicr_means(chosen[(c, s, 'ssl_sicr')], 'syllables') if (c, s, 'ssl_sicr') in chosen else None
            for tag in (0, 1):
                row['correct_serial_syllables_%d%s' % (tag, sfx_ssl)] = round(m[tag]['serial'], 6) if m else ''
                row['correct_bigrams_%d%s' % (tag, sfx_ssl)] = round(m[tag]['bigrams'], 6) if m else ''
                row['correct_trigrams_%d%s' % (tag, sfx_ssl)] = round(m[tag]['trigrams'], 6) if m else ''
                total_name = 'correct_syllables_total_0_2' if (s, tag) == (2, 0) else 'correct_syllables_total_%d%s' % (tag, sfx_ssl)
                row[total_name] = round(m[tag]['total'], 6) if m else ''
            for name, k in (('diff_syllables', 'serial'), ('diff_bigrams', 'bigrams')):
                denom = m[1][k] + m[0][k] if m else 0
                row[name + sfx_ssl] = round((m[1][k] - m[0][k]) / denom, 6) if m and denom else ''
        for s in (1, 2):
            row['afc_accuracy_s%d' % s] = afc(chosen.get((c, s, 'lsl_afc')))
        for s in (1, 2):
            m = sicr_means(chosen[(c, s, 'lsl_sicr')], 'letters') if (c, s, 'lsl_sicr') in chosen else None
            for tag in (0, 1):
                for name, k in (('correct_bigrams', 'bigrams'), ('correct_letters_total', 'total'),
                                ('correct_serial_letters', 'serial'), ('correct_trigrams', 'trigrams')):
                    row['%s.lsl_%d_s%d' % (name, tag, s)] = round(m[tag][k], 6) if m else ''
        for s in (1, 2):
            row['n_ssl_afc_s%d' % s] = len(chosen.get((c, s, 'ssl_afc'), []))
            row['n_lsl_afc_s%d' % s] = len(chosen.get((c, s, 'lsl_afc'), []))
        wide.append(row)
    return wide


# ---------------------------------------------------------------- main
def write(name, rows, folder=None):
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(os.path.join(folder or DATA, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def assign_sessions(dates):
    """Map each date to a session number; a gap of more than SESSION_GAP_DAYS starts a new session."""
    out, session, prev = {}, 0, None
    for d in sorted(set(dates)):
        if prev is None or (d - prev).days > SESSION_GAP_DAYS:
            session += 1
        out[d] = session
        prev = d
    return out


def main():
    if not os.path.exists(_CONFIG_PATH):
        sys.exit('study_config.json not found in %s; it is required to score raw exports '
                 '(use --from-trials to score a shared trial file)' % DATA)
    runs = read_runs()
    for r in runs:
        r['measures'] = {}
        if r['kind'] == 'LSL':
            r['measures']['lsl_afc'] = lsl_afc(r['rows'])
            r['measures']['lsl_sicr'] = lsl_sicr(r['rows'])
        else:
            r['measures']['ssl_afc'] = ssl_afc(r['rows'])
    # SSL SICR runs (dates only) so transcription sessions line up with the same session windows
    sicr_runs = {}
    for d in RAW_DIRS:
        for path in glob.glob(os.path.join(d, '*.csv')):
            with open(path, encoding='utf-8-sig', errors='replace') as f:
                for x in csv.DictReader(f):
                    if G(x, 'Task Name') == 'SSL SICR' and parse_date(G(x, 'Local Date and Time')):
                        sicr_runs.setdefault((child_from_file(path), G(x, 'Participant Private ID')),
                                             parse_date(G(x, 'Local Date and Time')))

    # sessions by date, per child, over every run (including started-but-empty runs, which mark a visit)
    child_dates = {}
    for r in runs:
        if r['date']:
            child_dates.setdefault(r['child'], []).append(r['date'])
    for (child, _), d in sicr_runs.items():
        child_dates.setdefault(child, []).append(d)
    session_of = {c: assign_sessions(ds) for c, ds in child_dates.items()}

    chosen, sources, trial_rows = {}, [], []
    for r in sorted(runs, key=lambda r: (r['child'], r['date'] or datetime.date.min)):
        for measure, trials in r['measures'].items():
            if not trials:
                continue
            session = session_of[r['child']][r['date']]
            slot = (r['child'], session, measure)
            note = []
            if r['backup_only']:
                note.append('run only in beh/SL data - Copy')
            if 'JGA' in ''.join(t.get('key', '') for t in trials):
                note.append('version E: recall sequences show JGA; familiarization word is GJA')
            override = RUN_OVERRIDES.get(slot)
            if slot in chosen or (override and r['priv'] != override):
                sources.append({'participant': r['child'], 'session': session, 'measure': measure,
                                'used': 'no', 'date': r['date'], 'task_version': r['task'],
                                'gorilla_private_id': r['priv'], 'source': r['file'], 'n_trials': len(trials),
                                'note': 'not used: RUN_OVERRIDES selects run %s' % override if override
                                        else 'later run in the same session; earliest run used'})
                continue
            chosen[slot] = trials
            sources.append({'participant': r['child'], 'session': session, 'measure': measure, 'used': 'yes',
                            'date': r['date'], 'task_version': r['task'], 'gorilla_private_id': r['priv'],
                            'source': r['file'], 'n_trials': len(trials), 'note': '; '.join(note)})
            for t in trials:
                trial_rows.append({'participant': r['child'], 'session': session, 'measure': measure, **t})

    keys, responses = read_transcriptions()
    for label, resp in responses.items():
        session, note = TRANSCRIPTION_SESSIONS[label]
        child = 'pl_' + re.match(r'plausl_(.+)_[12]$', label).group(1)
        if session is None:
            sources.append({'participant': child, 'session': '', 'measure': 'ssl_sicr', 'used': 'no',
                            'source': 'sicr_recordings/plausl_sicr_transcriptions.xlsx (%s)' % label, 'note': note})
            continue
        trials = ssl_sicr(keys[label], resp)
        used = 'yes'
        if all(t['response'] in ('xx', '') for t in trials):
            note = '; '.join(filter(None, [note, 'all 16 responses are xx (no response); left blank']))
            used = 'no'
        else:
            chosen[(child, session, 'ssl_sicr')] = trials
        sources.append({'participant': child, 'session': session, 'measure': 'ssl_sicr', 'used': used,
                        'source': 'final_plausl_data/sicr_recordings/plausl_sicr_transcriptions.xlsx (%s)' % label,
                        'n_trials': len(trials), 'note': note})
        for t in trials:
            trial_rows.append({'participant': child, 'session': session, 'measure': 'ssl_sicr', **t})

    wide = wide_table(chosen)

    write('plausl_sl_scores.csv', wide)
    write('plausl_sl_sources.csv', sorted(sources, key=lambda r: (r['participant'], str(r.get('session')), r['measure'])))
    write('plausl_sl_scores_by_trial.csv', trial_rows)
    write('plausl_sl_summary.csv', summarize(wide))


def describe_column(col):
    """(task, measure, session, trial_type) for a plausl_sl_scores.csv column, or None."""
    session = 2 if re.search(r'(_t2|_s2|_0_2)$', col) else 1
    if col.startswith('n_'):
        return None
    if col in ('mean_accuracy', 'mean_accuracy_t2'):
        return 'syllable', '2AFC accuracy', session, 'all'
    if col.startswith('afc_accuracy'):
        return 'letter', '2AFC accuracy', session, 'all'
    m = re.match(r'(correct_\w+?|diff_\w+?)(\.lsl)?(?:_([01]))?(?:_t2|_s[12]|_2)?$', col)
    if not m:
        return None
    name, lsl, tag = m.groups()
    task = 'letter' if lsl else 'syllable'
    labels = {'correct_serial_syllables': 'SICR correct in serial position', 'correct_serial_letters': 'SICR correct in serial position',
              'correct_bigrams': 'SICR correct bigrams', 'correct_trigrams': 'SICR correct trigrams',
              'correct_syllables_total': 'SICR correct items (any position)', 'correct_letters_total': 'SICR correct items (any position)',
              'diff_syllables': 'SICR difference score: serial (target - foil)/(target + foil)',
              'diff_bigrams': 'SICR difference score: bigrams (target - foil)/(target + foil)'}
    return task, labels[name], session, {'0': 'foil', '1': 'target', None: 'difference'}[tag]


def summarize(wide):
    rows = []
    for col in wide[0]:
        d = describe_column(col)
        if not d:
            continue
        vals = [float(r[col]) for r in wide if r.get(col) not in ('', None)]
        task, measure, session, trial_type = d
        rows.append({'task': task, 'measure': measure, 'session': session, 'trial_type': trial_type, 'column': col,
                     'n': len(vals),
                     'mean': round(statistics.mean(vals), 4) if vals else '',
                     'sd': round(statistics.stdev(vals), 4) if len(vals) > 1 else '',
                     'median': round(statistics.median(vals), 4) if vals else '',
                     'min': round(min(vals), 4) if vals else '',
                     'max': round(max(vals), 4) if vals else ''})
    order = {'2AFC accuracy': 0}
    return sorted(rows, key=lambda r: (r['task'] != 'syllable', order.get(r['measure'], 1), r['measure'],
                                       r['session'], ['target', 'foil', 'difference', 'all'].index(r['trial_type'])))


def from_trials(path, folder):
    """Recompute plausl_sl_scores.csv and plausl_sl_summary.csv from a trial-level file (one row per scored trial,
    columns as in plausl_sl_scores_by_trial.csv or the shared plausl_sl_trials.csv)."""
    as_int = lambda v: int(float(v)) if v not in ('', None) else None
    chosen = {}
    with open(path, encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            t = {'trial': as_int(r['trial']), 'response': r.get('response', '')}
            if r['measure'].endswith('_afc'):
                t.update(score=as_int(r['score']))
            else:
                t.update(key=r.get('key', ''), tag=as_int(r['tag']), counted=as_int(r['counted']),
                         **{k: as_int(r[k]) for k in ('serial', 'bigrams', 'trigrams', 'total')})
            chosen.setdefault((r['participant'], as_int(r['session']), r['measure']), []).append(t)
    for slot, trials in list(chosen.items()):  # same rule as main(): an all-"xx" transcription is left blank
        if slot[2] == 'ssl_sicr' and all(t['response'] in ('xx', '') for t in trials):
            del chosen[slot]
    wide = wide_table(chosen)
    os.makedirs(folder, exist_ok=True)
    write('plausl_sl_scores.csv', wide, folder)
    write('plausl_sl_summary.csv', summarize(wide), folder)
    print('scored %d children from %s; wrote plausl_sl_scores.csv and plausl_sl_summary.csv to %s'
          % (len(wide), path, folder))


if __name__ == '__main__':
    if '--from-trials' in sys.argv:
        args = sys.argv[sys.argv.index('--from-trials') + 1:]
        if not args:
            sys.exit('usage: python3 score_plausl_sl.py --from-trials <trials.csv> [output folder]')
        from_trials(args[0], args[1] if len(args) > 1 else '.')
    else:
        main()
