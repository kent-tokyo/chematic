from pathlib import Path
O=Path(__file__).resolve().parent
source=(O/'formal.py').read_text()
marker='\nsource={};external={};negative=[];failures=[];primary={}'
assert source.count(marker)==1
prefix=source.split(marker,1)[0]
exec(compile(prefix,str(O/'formal.py'),'exec'),{'__file__':str(O/'formal.py'),'__name__':'__main__'})
assert not list(O.glob('paired-*.json'))
print('All preflight conditions passed without any formal timing samples',flush=True)
