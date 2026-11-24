"""Resumable production entry point; one writer process per run directory."""
import portalocker
import json
import os
import sys
import subprocess
import traceback
from datetime import datetime,timezone
from global_common import *

def main():
    setup()
    lock=(OUT/'progress/run.lock').open('a+')
    try:portalocker.lock(lock,portalocker.LOCK_EX|portalocker.LOCK_NB)
    except portalocker.exceptions.LockException:raise RuntimeError('A global run is already active; do not launch a duplicate')
    atomic_json(OUT/'progress/process.json',{'pid':os.getpid(),'started_utc':datetime.now(timezone.utc).isoformat(),
        'python':sys.executable,'entrypoint':str(Path(__file__).resolve()),'status':'running'})
    try:
        subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/'scripts/tests'),'-v'],check=True)
        from global_prepare import prepare
        from validate_global_model import validate
        from global_rasters import run_years
        prepare();validate();run_years()
        from global_background import build_background
        build_background()
        from global_species import extract_species
        extract_species()
        from global_statistics import calculate
        calculate()
        from finalize_global_tables import main as finalize_tables
        from compare_china_rosters import main as compare_rosters
        finalize_tables()
        compare_rosters()
        atomic_json(OUT/'progress/current.json',{'status':'numerical_analysis_complete','figures':'pending scientific and visual review',
            'completed_utc':datetime.now(timezone.utc).isoformat(),'years':CFG['years'],'global_includes_china':True})
        atomic_json(OUT/'progress/process.json',{'pid':os.getpid(),'status':'complete','completed_utc':datetime.now(timezone.utc).isoformat()})
        event('numerical_analysis_complete',tables=str(OUT/'tables'),figures='pending scientific and visual review')
    except BaseException as error:
        detail={'status':'failed' if not isinstance(error,KeyboardInterrupt) else 'interrupted',
                'time_utc':datetime.now(timezone.utc).isoformat(),'error':str(error),'traceback':traceback.format_exc()}
        atomic_json(OUT/'progress/last_error.json',detail)
        atomic_json(OUT/'progress/process.json',{'pid':os.getpid(),**detail})
        event('run_stopped',error=str(error));raise
    finally:
        portalocker.unlock(lock);lock.close()

if __name__=='__main__':main()
