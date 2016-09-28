import drmaa
import os
import time

def main():
    """
    Submit a job and wait for it to finish.
    Note, need file called sleeper.sh in home directory.
    """
    with drmaa.Session() as s:
        print('Creating job template')
        jt = s.createJobTemplate()
        jt.remoteCommand = os.path.join(os.getcwd(), '/mnt/scratch/install/qsub_test.sh')
        jt.args = ['42', 'Simon says:']
        jt.joinFiles = True

        jobid = s.runJob(jt)
        print('Your job has been submitted with ID %s' % jobid)

        retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        print('Job: {0} finished with status {1}'.format(retval.jobId, retval.hasExited))

        print('Cleaning up')
        s.deleteJobTemplate(jt)

def get_drmaa_imformation():
    """ Query the system. """
    with drmaa.Session() as s:
        print('A DRMAA object was created')
        print('Supported contact strings: %s' % s.contact)
        print('Supported DRM systems: %s' % s.drmsInfo)
        print('Supported DRMAA implementations: %s' % s.drmaaImplementation)

        print('Exiting')

if __name__=='__main__':
    main()
    get_drmaa_imformation()

