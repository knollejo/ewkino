import argparse, datetime, getpass, json, os, subprocess, time

def get_timestamp():
    return '{:%y%m%d_%H%M%S}'.format(datetime.datetime.now())

def get_user():
    return getpass.getuser()

def do_qsub(command):
    while True:
        try:
            qsub_output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as error:
            print('Caught error: "{}". Attempting resubmission'.format(error.output.split('\n')[0]))
            time.sleep(1)
        else:
            return qsub_output

def init_submit():
    timestamp = get_timestamp()
    if not os.path.exists('output'):
        os.mkdir('output')
    if not os.path.exists('output/{}'.format(timestamp)):
        os.mkdir('output/{}'.format(timestamp))
    if not os.path.exists('output/{}/submit'.format(timestamp)):
        os.mkdir('output/{}/submit'.format(timestamp))
    print 'submitting to', 'output/{}'.format(timestamp)
    return timestamp

def finalize_submit(timestamp, jobids):
    with open('output/{}/submit/jobids.json'.format(timestamp), 'w') as f:
        json.dump(jobids, f)

def submit(year, process, timestamp=None):
    do_init = bool(timestamp is None)
    if do_init:
        timestamp = init_submit()
    with open('output/{0}/submit/run_{0}_{1}_{2}.sh'.format(timestamp, year, process), 'w') as f:
        f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        f.write('cd {}\n'.format(os.path.abspath(os.getcwd())))
        f.write('eval `scram runtime -sh`\n')
        f.write('./ttZAnalysis {} {} {}\n'.format(year, process, timestamp))
        f.write('echo "Job done!"\n')
    sub = 'qsub output/{0}/submit/run_{0}_{1}_{2}.sh -l walltime=24:00:00'.format(timestamp, year, process)
    nam = ' -N {}_{}'.format(year, process)
    log = ' -o output/{0}/submit/log_{0}_{1}_{2}.txt'.format(timestamp, year, process)
    err = ' -e output/{0}/submit/err_{0}_{1}_{2}.txt'.format(timestamp, year, process)
    command = sub+nam+log+err
    qsub_output = do_qsub(command)
    jobid = qsub_output.split('\n')[0]
    print timestamp, year, process, jobid
    jobids = {'{}_{}'.format(year, process): jobid}
    if do_init:
        finalize_submit(timestamp, jobids)
    else:
        return jobids

def submit_year(year, timestamp=None):
    do_init = bool(timestamp is None)
    if do_init:
        timestamp = init_submit()
    processes = (
        "data",
        "ttZ",
        "ttH", "tZq", "tWZ", "ttW", "tHQ", "tHW", "ttZlight", "tttt", "ttWW", "ttWZ", "ttZZ",
        "WZ3L", "WZ2L",
        "DY", "tG", "ttG", "WG", "tt",
        "ZZ4L", "ZZ2E2M", "ZZ2E2T", "ZZ2M2T", "ZZ4E", "ZZ4M", "ggHZZ", "VBFHZZ", "WpHZZ", "WmHZZ", "ZHZZ",
        "WZG", "ZZZ", "WZZ", "WWZ", "WWW", "WWDS", "WW",
    )
    jobids = {}
    for process in processes:
        jobid = submit(year, process, timestamp)
        jobids.update(jobid)
    if do_init:
        finalize_submit(timestamp, jobids)
    else:
        return jobids

def submit_all():
    timestamp = init_submit()
    years = (2016, 2017, 2018)
    jobids = {}
    for year in years:
        jobids_year = submit_year(year, timestamp)
        jobids.update(jobids_year)
    finalize_submit(timestamp, jobids)

def check(dirname):
    timestamp = dirname.split('/')[-1]
    with open('{}/submit/jobids.json'.format(dirname)) as f:
        jobids = json.load(f)
    user = get_user()
    qstat = do_qsub('qstat -u {}'.format(user)).split('\n')
    running = []
    for line in qstat:
        if not user in line:
            continue
        running.append(line.split('.')[0])
    newids = {}
    for jobname, jobid in jobids.iteritems():
        if jobid.split('.')[0] in running:
            newids[jobname] = jobid
            running.remove(jobid.split('.')[0])
            continue
        year, process = jobname.split('_')
        if os.path.exists('{}/{}_{}.root'.format(dirname, process, year)):
            continue
        newid = submit(year, process, dirname.split('/')[-1])
        newids.update(newid)
    finalize_submit(timestamp, newids)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', '-y', nargs=1, default=None)
    parser.add_argument('--process', '-p', nargs=1, default=None)
    parser.add_argument('--check', '-c', nargs=1, default=None)
    args = parser.parse_args()
    if args.check:
        dirname = args.check[0]
        if dirname[-1] == '/':
            dirname = dirname[:-1]
        check(dirname)
    elif args.process and args.year:
        submit(args.year[0], args.process[0])
    elif args.year:
        submit_year(args.year[0])
    else:
        submit_all()
