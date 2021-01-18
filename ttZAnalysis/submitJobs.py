#import python library classes
import os
import sys

#import other code from framework
from jobSubmission import submitQsubJob, initializeJobScript

if __name__ == '__main__' :

    #WARNING : it is assumed this script is run from the 'skimmer' directory
    current_directory = os.path.dirname( os.path.abspath( __file__ ) )
    print current_directory

    year    = sys.argv[1]
    process = sys.argv[2]
    wall_time = '24:00:00'
    if len( sys.argv ) != 3:
        print( 'Error: submitJobs.py requires additional command-line arguments.' )
        print( 'Usage: python submitJobs.py < year > < process >' )
        sys.exit()

    #make a list of samples (the main directories) and the subdirectories with the latest version of the ntuples ( one for each main directory )
    #make a job script
    script_name = 'submit.sh'
    with open( script_name, 'w') as script:
        initializeJobScript( script )
        script.write('cd {}\n'.format( current_directory ) )
        skim_command = './ttZAnalysis {} {}\n'.format( year, process )
        script.write( skim_command )

    #submit job and catch errors
    submitQsubJob( script_name, wall_time )
