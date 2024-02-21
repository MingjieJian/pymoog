#!/usr/bin/python
from . import private

class rundir_num(object):

    def __init__(self, pymoog_path, run_type, prefix=''):
        '''

        '''
        self.pymoog_path = pymoog_path

        if prefix == '':
            self.rundir_path = '{}/{}-{}-{}/'.format(self.pymoog_path, run_type, private.datetime.now().strftime("%H:%M:%S.%f"), ''.join(private.secrets.choice(private.string.ascii_uppercase + private.string.ascii_lowercase) for i in range(9)))
        else:
            self.rundir_path = '{}/{}-{}-{}-{}/'.format(self.pymoog_path, prefix, run_type, private.datetime.now().strftime("%H:%M:%S.%f"), ''.join(private.secrets.choice(private.string.ascii_uppercase + private.string.ascii_lowercase) for i in range(9)))

        private.subprocess.run(['mkdir', '-p', self.rundir_path])
        
    def remove(self):
        private.subprocess.run(['rm', '-r', self.rundir_path])
