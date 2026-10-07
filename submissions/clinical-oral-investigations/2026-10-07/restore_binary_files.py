from pathlib import Path
import base64,hashlib,json,sys
root=Path(__file__).resolve().parents[3]
folder=root/'submissions/clinical-oral-investigations/2026-10-07'
records=json.loads((folder/'provenance/binary_transport.json').read_text())
for record in records:
    target=(root/record['path']).resolve();encoded=(root/record['encoded_path']).resolve()
    assert target.is_relative_to(folder) and encoded.is_relative_to(folder)
    if encoded.exists():
        raw=base64.b64decode(encoded.read_text(),validate=False)
        assert len(raw)==record['bytes'] and hashlib.sha256(raw).hexdigest()==record['sha256']
        target.write_bytes(raw)
        if '--cleanup' in sys.argv:encoded.unlink()
    else:
        assert hashlib.sha256(target.read_bytes()).hexdigest()==record['sha256']
print('Verified and restored',len(records),'binary files')
