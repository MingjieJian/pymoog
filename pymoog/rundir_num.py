#!/usr/bin/python
from . import private

class rundir_num(object):

    def __init__(self, pymoog_path, run_type, prefix=''):
        '''

        '''
        # rundir_max_num = 10
        self.pymoog_path = pymoog_path
        # file_list = private.subprocess.run(['ls', pymoog_path], stdout=private.subprocess.PIPE)
        # file_list = str(file_list.stdout, encoding = "utf-8").split('\n')
        # file_list = [i for i in file_list if '.lock' in i]
        if prefix == '':
            self.rundir_path = '{}/{}-{}/'.format(self.pymoog_path, run_type, private.datetime.now().strftime("%H:%M:%S.%f"))
        else:
            self.rundir_path = '{}/{}-{}-{}/'.format(self.pymoog_path, prefix, run_type, private.datetime.now().strftime("%H:%M:%S.%f"))

        private.subprocess.run(['mkdir', '-p', self.rundir_path])
                  
        # dir_exist = [int(private.re.search('rundir(.+)\.lock', i).group(1)) for i in file_list]

        # rundir_list = list(range(1, rundir_max_num+1))
        # for ele in dir_exist:
        #     try:
        #         rundir_list.remove(ele)
        #     except ValueError:
        #         continue

        # if len(rundir_list) == 0:
        #     raise ValueError('No available rundir; please delete the available rundir lock or raise maximum number of rundir.')
        # rundir_num = min(rundir_list)
        # self.rundir_path = '{}/rundir{}/'.format(self.pymoog_path, self.rundir_num)
        
    # def lock(self):
    #     private.subprocess.run(['touch', '{}/rundir{}.lock'.format(self.pymoog_path, self.rundir_num)])
    #     # Create rundir folder
    #     
        
    def remove(self):
        private.subprocess.run(['rm', '-r', self.rundir_path])
        
    # def clear_lock(self, num):
    #     if num != 'all':
    #         private.subprocess.run(['rm', '{}/rundir{}.lock'.format(self.pymoog_path, num)])
    #     elif num == 'all':
    #         private.subprocess.run(['rm', '{}/rundir*.lock'.format(self.pymoog_path)])