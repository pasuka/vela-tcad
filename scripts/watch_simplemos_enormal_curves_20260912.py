"""Stream changed status for the current run; no scheduling or solver control."""
import contextlib,io,json,time
import status_simplemos_enormal_curves_20260912 as status

def main():
    output=status.ROOT/'reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912'
    previous=None
    while not (status.LOCAL/'stop_status_watch').exists():
        with contextlib.redirect_stdout(io.StringIO()):snapshot=status.main()
        key=(snapshot['finished_attempts'],snapshot['qualified_states'],snapshot['failed_attempts'],snapshot['paired_currents'])
        if key!=previous:
            print(json.dumps(snapshot,ensure_ascii=False),flush=True);previous=key
        if all((output/('vela_'+arm+'_evidence.json')).exists() for arm in ('continuation','native')):
            print('Both run evidence files exist; independent verification and analysis are still required.',flush=True)
            return
        time.sleep(10)
    print('Status watch stopped by task-local stop file. Solver jobs were not changed.',flush=True)

if __name__=='__main__':main()
