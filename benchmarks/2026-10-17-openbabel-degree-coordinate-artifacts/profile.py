from pathlib import Path
import subprocess,json,time
O=Path(__file__).resolve().parent
command=[str(O/'baseline-v3000_3d'),str(O/'v3000_3d-branched-spiral.molblocks'),'parse','100','64']
with (O/'diagnostic-profile-worker.json').open('w') as output,(O/'diagnostic-profile-worker.stderr').open('w') as errors:
 process=subprocess.Popen(command,stdout=output,stderr=errors)
 try:
  time.sleep(.3)
  sampled=subprocess.run(['/usr/bin/sample',str(process.pid),'2','-file',str(O/'diagnostic-profile.txt')],capture_output=True,text=True)
  (O/'profile-sampler-output.txt').write_text(sampled.stdout+sampled.stderr)
  assert sampled.returncode==0,sampled.stderr
  assert process.wait(timeout=180)==0
 finally:
  if process.poll() is None:process.terminate();process.wait()
(O/'diagnostic-profile-scope.json').write_text(json.dumps({'scope':'Separate baseline branched-3D parse worker CPU diagnostic, including preparation and repeated parsing. Use parse/classifier frame ancestry; no diagnostic elapsed value enters formal matrices','command':command,'sample_duration_seconds':2},indent=2)+'\n')
