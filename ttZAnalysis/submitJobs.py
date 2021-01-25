import datetime, os, subprocess, sys, time

def get_timestamp():
    return '{:%y%m%d_%H%M%S}'.format(datetime.datetime.now())

def submit(year, process, timestamp=None):
    if timestamp is None:
        timestamp = get_timestamp()
    with open('submit/run_{}_{}_{}.sh'.format(timestamp, year, process), 'w') as f:
        f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        f.write('cd {}\n'.format(os.path.abspath(os.getcwd())))
        f.write('eval `scram runtime -sh`\n')
        f.write('./ttZAnalysis {} {} {}\n'.format(year, process, timestamp))
        f.write('echo "Job done!"\n')
    sub = 'qsub submit/run_{}_{}_{}.sh -l walltime=24:00:00'.format(timestamp, year, process)
    nam = ' -N {}_{}'.format(year, process)
    log = ' -o submit/log_{}_{}_{}.txt'.format(timestamp, year, process)
    err = ' -e submit/err_{}_{}_{}.txt'.format(timestamp, year, process)
    command = sub+nam+log+err
    while True:
        try:
            qsub_output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as error:
            print('Caught error : "{}". Attempting resubmission.'.format(error.output.split('\n')[0]))
            time.sleep(1)
        else:
            print timestamp, year, process, qsub_output.split('\n')[0]
            return

def submit_year(year, timestamp=None):
    if timestamp is None:
        timestamp = get_timestamp()
    processes = (
        "data",
        "ttZ",
        "ttH", "tZq", "tWZ", "ttW", "tHQ", "tHW", "ttZlight", "tttt", "ttWW", "ttWZ", "ttZZ",
        "WZ3L", "WZ2L",
        "DY", "tG", "ttG", "WG", "tt",
        "ZZ4L", "ZZ2E2M", "ZZ2E2T", "ZZ2M2T", "ZZ4E", "ZZ4M", "ggHZZ", "VBFHZZ", "WpHZZ", "WmHZZ", "ZHZZ",
        "WZG", "ZZZ", "WZZ", "WWZ", "WWW", "WWDS", "WW",
    )
    for process in processes:
        submit(year, process, timestamp)

def submit_all():
    timestamp = get_timestamp()
    years = (2016, 2017, 2018)
    for year in years:
        submit_year(year, timestamp)

if __name__ == '__main__' :
    if len(sys.argv)>2:
        submit(sys.argv[1], sys.argv[2])
    elif len(sys.argv)>1:
        submit_year(sys.argv[1])
    else:
        submit_all()
