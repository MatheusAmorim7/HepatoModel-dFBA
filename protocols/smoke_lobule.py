"""Run a four-minute numerical smoke test without replacing the normal output."""
import csv
import json
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
out=Path(tempfile.mkdtemp(prefix='hepatomodel-smoke-'))
tree=ET.parse(ROOT/'config/PhysiCell_settings.xml')
root=tree.getroot()
root.find('overall/max_time').text='4'
root.find('save/folder').text=str(out)
settings=out/'settings.xml'
tree.write(settings,encoding='utf-8',xml_declaration=True)
with (out/'console.log').open('w') as log:
    subprocess.run([str(ROOT/'project'),str(settings)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,
                   check=True,timeout=300)
with (out/'fba_summary.csv').open() as f:
    rows=list(csv.DictReader(f))
assert rows and all(r['dfba_status']=='optimal' for r in rows)
assert any(float(r['glucose_flux'])<0 for r in rows)
assert any(float(r['oxygen_flux'])<0 for r in rows)
memory=list(csv.DictReader((out/'memory_usage.csv').open()))
assert memory and max(float(r['rss_MB']) for r in memory)>0
checks=dict(status='passed',scope='numerical smoke, not physiological validation',
            logged_cell_steps=len(rows),all_optimal=True,glucose_and_oxygen_uptake=True,
            memory_report=True,output=str(out))
(out/'smoke_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
