from pathlib import Path
import hashlib,json,shutil,datetime
O=Path(__file__).resolve().parent
original=O/'formal.py';source=original.read_text();plan=json.loads((O/'measurement-plan.json').read_text())
assert not (O/'continuation-seal.json').exists(),'No continuation reruns'
prior={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.glob('paired-*.json')}
assert len(prior)==6
for p in plan['profiles'][:5]:
 data=json.loads((O/('paired-'+p['kind']+'-'+p['profile']+'.json')).read_text());assert data['complete'] and len(data['blocks'])==24
partial=O/'paired-v3000-continued.json';data=json.loads(partial.read_text());assert not data['complete'] and len(data['blocks'])==5
remaining=plan['profiles'][6:];assert len(remaining)==4
for p in remaining:assert not (O/('paired-'+p['kind']+'-'+p['profile']+'.json')).exists()
start=source.index("\nfor profile in plan['profiles']:\n integrity();kind")
end=source.index('\npassed=len(primary)',start)
body=source[start:end];header=source.split("assert not list(O.glob('paired-*.json'))",1)[0]
record={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'First observations only for the four profiles never reached after a transient ENOSPC. Preserve all five completed profiles and the five-block incomplete profile byte-exact. Never rerun the interrupted profile.','original_formal_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'original_timing_loop_sha256':hashlib.sha256(body.encode()).hexdigest(),'continuation_driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'measurement_plan_sha256':hashlib.sha256((O/'measurement-plan.json').read_bytes()).hexdigest(),'existing_raw_files_sha256':prior,'scheduled_keys':[p['kind']+'-'+p['profile'] for p in remaining],'unpersisted_native_rows_in_failed_block':6,'unpersisted_values_not_reconstructed':True,'experiment_remains_incomplete_and_ineligible_for_adoption':True,'failure_log':'formal.log','source_binaries_thresholds_orders_and_timing_loop_unchanged':True,'disk_free_bytes_before':shutil.disk_usage(O).free}
(O/'continuation-seal.json').write_text(json.dumps(record,indent=2)+'\n')
scope={'__file__':str(original),'__name__':'__main__'};exec(compile(header,str(original),'exec'),scope)
scope['integrity']();scope['preflight']=json.loads((O/'preflight-source-equality.json').read_text())
scope.update(source={},external={},negative=[],failures=[],primary={})
for profile in remaining:
 assert shutil.disk_usage(O).free>=512*1024*1024,'Stop before new observations if disk headroom is low'
 for name,h in prior.items():assert hashlib.sha256((O/name).read_bytes()).hexdigest()==h
 scope['plan']={**plan,'profiles':[profile]}
 exec(compile('\n'*source[:start].count('\n')+body,str(original),'exec'),scope)
for name,h in prior.items():assert hashlib.sha256((O/name).read_bytes()).hexdigest()==h
(O/'continuation-finished.json').write_text(json.dumps({'completed_new_profile_keys':record['scheduled_keys'],'all_preexisting_raw_files_unchanged':True,'interrupted_profile_not_rerun':True,'adoption_eligible':False,'reason':'One planned profile remains incomplete; completed SDF source write floors also fail'},indent=2)+'\n')
print('Four previously unmeasured profiles completed; experiment remains incomplete, not adopted',flush=True)
