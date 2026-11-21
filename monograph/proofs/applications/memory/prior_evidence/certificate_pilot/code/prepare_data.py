"""Deterministic corpus preparation; no model output or test inspection.
Run after placing the preserved Project Gutenberg download in data/pg15058_raw.txt.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, unicodedata

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
URL = 'https://www.gutenberg.org/cache/epub/15058/pg15058.txt'
EXPECTED_SHA256 = '9af720d0cee3df47311b7dfcec9986c758fa08a2758bb0f9bd58e97de72d5c92'
START = 'Ma bonne, ma bonne, dit un jour Sophie'

def sha(b): return hashlib.sha256(b).hexdigest()
def write_json(path, value): path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')

raw = (DATA / 'pg15058_raw.txt').read_bytes()
assert sha(raw) == EXPECTED_SHA256
s = raw.decode('utf-8-sig').replace('\r\n', '\n').replace('\r', '\n')
assert s.count(START) == 1
start = s.index(START)
end = s.index('*** END OF THE PROJECT GUTENBERG EBOOK', start)
body = unicodedata.normalize('NFC', s[start:end])
paragraphs = [re.sub(r'\s+', ' ', p).strip() for p in re.split(r'\n\s*\n', body) if p.strip()]
normalized = unicodedata.normalize('NFC', '\n\n'.join(paragraphs).lower())
selected = normalized[:180000]
if len(selected) == 180000:
    cuts = [0,140000,160000,180000]
else:
    n = len(selected); cuts = [0, int(.8*n), int(.9*n), n]
splits = {name: selected[cuts[i]:cuts[i+1]] for i,name in enumerate(['train','validation','test'])}
manifest = {}
for i,(name,text) in enumerate(splits.items()):
    b = text.encode('utf-8'); (DATA/f'{name}.txt').write_bytes(b)
    manifest[name] = {'characters':len(text), 'bytes_utf8':len(b), 'sha256':sha(b), 'start_character_in_selected':cuts[i], 'end_character_exclusive':cuts[i+1]}
train64 = {splits['train'][i:i+64] for i in range(max(0,len(splits['train'])-63))}
overlap = {}
for name in ['validation','test']:
    vals = [splits[name][i:i+64] for i in range(max(0,len(splits[name])-63))]
    matches = [s for s in vals if s in train64]
    overlap[name] = {'positions_with_exact_64_character_match_in_train':len(matches), 'distinct_matching_64_character_strings':len(set(matches)), 'positions_checked':len(vals)}
old = ROOT.parent/'multiletter_pilot'/'data'
comparison = {}
for name in splits:
    f = old/f'{name}.txt'
    if f.exists():
        previous = sha(f.read_bytes())
        comparison[name] = {'previous_jules_verne_sha256':previous, 'new_sha256':manifest[name]['sha256'], 'same_bytes':previous==manifest[name]['sha256']}
source = {
 'title':'Les malheurs de Sophie', 'author':'comtesse de Sophie Ségur', 'publication_year':1858, 'language':'French',
 'source_url':URL, 'source_landing_page':'https://www.gutenberg.org/ebooks/15058',
 'source_web_references':['turn93view0','turn94view0'],
 'retrieved_at_utc':'2026-11-01T15:06:00+00:00',
 'retrieved_at_note':'Approximate UTC minute; filesystem mtime of raw file is authoritative below.',
 'raw_file_mtime_utc':datetime.fromtimestamp((DATA/'pg15058_raw.txt').stat().st_mtime, timezone.utc).isoformat(),
 'prepared_at_utc':datetime.now(timezone.utc).isoformat(),
 'source_bytes':len(raw), 'source_sha256':sha(raw), 'raw_file':'pg15058_raw.txt',
 'copyright_note':'Project Gutenberg catalog states public domain in the USA. Original French work published in 1858; author died in 1874. Raw source retains the complete distribution wrapper and license.',
 'transformations':[
  'UTF-8-SIG decode; CRLF and CR normalized to LF',
  'Start at the unique literal first narrative phrase: '+START,
  'Stop before the literal Project Gutenberg END OF EBOOK marker',
  'Unicode NFC normalization',
  'Split paragraphs on blank lines, collapse internal whitespace to one ASCII space, separate paragraphs by two newlines',
  'Lowercase Unicode text; apply NFC again',
  'Select first 180000 normalized Unicode characters, or whole body if shorter',
  'Contiguous deterministic train/validation/test cuts: 140000/20000/20000 if 180000 available, otherwise floor80%/floor90% boundaries',
  'No boundary moving, deduplication, sampling by content, or model-based selection'
 ],
 'full_normalized_body_characters':len(normalized), 'selected_characters':len(selected), 'selected_sha256':sha(selected.encode('utf-8')),
 'is_complete_work':len(selected)==len(normalized),
 'freshness_scope':'New acquisition and different author/work from the preceding Jules Verne pilot. No claim about all historical experiments or pretraining corpora.',
 'test_handling':'Only deterministic preparation, hashing and exact duplicate checks; no model evaluation or content-based inspection of the held-out split.'
}
# Correct retrieval time from the actual preserved download, not a guessed time.
source['retrieved_at_utc'] = source['raw_file_mtime_utc']; del source['retrieved_at_note']
write_json(DATA/'source.json',source)
write_json(DATA/'manifest.json',{'splits':manifest,'overlap_check_64_characters':overlap,'comparison_with_previous_jules_verne':comparison,'selected_characters':len(selected)})
print(json.dumps({'characters':{k:v['characters'] for k,v in manifest.items()},'normalized_body_characters':len(normalized),'overlap_check':overlap,'previous_all_different':all(not v['same_bytes'] for v in comparison.values()),'source_url':URL},indent=2))
