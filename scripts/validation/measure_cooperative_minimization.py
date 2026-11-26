import hashlib, json, os, statistics, sys, time
import chematic
p='validation/results/2026-11-27-forcefield-hetero/'
a=json.load(open(p+'input.json'))['results']
b={x['name']:x for x in json.load(open(p+'unperturbed.json'))['results']}
rows=[x for x in a if x['category']=='phosphorus_sulfur'][:20]
mode=sys.argv[1]
results=[]
for row in rows:
    mol=chematic.from_smiles(row['smiles']).add_hydrogens()
    coords=b[row['name']]['coordinates']
    for kind in ('MMFF94','UFF'):
        try:
            field=mol.uff_force_field(coords) if kind=='UFF' else mol.mmff_force_field(coords)
            modes=('normal','controlled') if mode=='new' else ('normal',)
            for method in modes:
                for rep in range(3):
                    field.set_positions_(coords)
                    t=time.perf_counter()
                    if method=='normal':
                        status,energy=field.minimize_(200)
                    else:
                        result=field.minimize_with_timeout_(60,max_iters=200)
                        status,energy=result['force_field_status'],result['energy']
                    elapsed=time.perf_counter()-t
                    digest=hashlib.sha256(json.dumps([status,energy,field.positions],separators=(',',':')).encode()).hexdigest()
                    results.append({'name':row['name'],'kind':kind,'method':method,'rep':rep,'elapsed':elapsed,'digest':digest})
            if mode=='new':
                for budget in (0.0001,0.001):
                    field.set_positions_(coords)
                    result=field.minimize_with_timeout_(budget,max_iters=200)
                    actual=field.energy()
                    results.append({'name':row['name'],'kind':kind,'method':'deadline','budget':budget,'result':result,'energy_diff':abs(actual-result['energy']) if result['status']=='timeout' and result['energy'] is not None else None})
        except Exception as e:
            results.append({'name':row['name'],'kind':kind,'error':repr(e)})
json.dump(results,open(sys.argv[2],'w'),indent=2)
