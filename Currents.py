#!/usr/bin/env python
# --------------------------------------------------------
#       Script to show the recorded current for the 3D beam test runs
# created on September 4th 2016 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from glob import glob
import time
import ipdb
import json as js
from collections import OrderedDict
import os, sys, shutil

from ConfigParser import ConfigParser
from optparse import OptionParser

import numpy as np
import ROOT as ro
from Utils import *

# ====================================
# CONSTANTS
# ====================================

# ====================================
# CLASS FOR THE DATA
# ====================================
class Currents:
    """reads in information from the keithley log file"""

    def __init__(self, settings_f, diam='', pola='both'):

        # settings
        self.settings_file = Correct_Path(settings_f)
        self.diamond = diam
        self.polarity = pola

        # config
        self.currents_logs_dir = ''
        self.testbeam_log_file = ''
        self.raw_files_dir = ''
        self.ReadSettingsFile()

        self.runs_info = None

        self.Histos = {}
        self.Stuff = []

    def ReadSettingsFile(self):
        if os.path.isfile(self.settings_file):
            pars = ConfigParser()
            pars.read(self.settings_file)
            print 'Loading settings file:', self.settings_file

            if pars.has_section('SETTINGS'):
                if pars.has_option('SETTINGS', 'currents_logs_dir'):
                    self.currents_logs_dir = pars.get('SETTINGS', 'currents_logs_dir')
                else:
                    print 'The option "currents_logs_dir" under the section "SETTINGS" is required'
                if pars.has_option('SETTINGS', 'testbeam_log_file'):
                    self.testbeam_log_file = pars.get('SETTINGS', 'testbeam_log_file')
                else:
                    print 'The option "testbeam_log_file" under the section "SETTINGS" is required'
                if pars.has_option('SETTINGS', 'raw_files_dir'):
                    self.raw_files_dir = pars.get('SETTINGS', 'raw_files_dir')
                else:
                    print 'The option "raw_files_dir" under the section "SETTINGS" is required'
            else:
                print 'The section "SETTINGS" is needed in the settings file'

    # ==========================================================================
    # region INIT

    def LoadRunsInfo(self, force_recreate=False):
        self.runs_info = self.GetRunInfoJson()
        if (not self.runs_info) or (force_recreate):
            if self.testbeam_log_file != '':
                temp_runs_info = OrderedDict()
                with open(self.testbeam_log_file, 'r') as ftb:
                    lines = ftb.readlines()
                    for line in lines:
                        if len(line.split()) > 1 and (IsInt(line.split()[0]) or ('/' in line.split()[0] and IsInt(line.split()[1]))):
                            run_i = line.split()[0] if IsInt(line.split()[0]) else line.split()[1]
                            diam_i = line.split()[1] if IsInt(line.split()[0]) else line.split()[2]
                            volt_i = line.split()[2].strip('V') if IsInt(line.split()[0]) else line.split()[3].strip('V')
                            if '-' in volt_i:
                                volt_i = volt_i.strip('-')
                                volt_i = -1000 * float(volt_i.strip('K').strip('k')) if 'K' in volt_i or 'k' in volt_i else -1 * float(volt_i)
                            else:
                                volt_i = volt_i.strip('+')
                                volt_i = 1000 * float(volt_i.strip('K').strip('k')) if 'K' in volt_i or 'k' in volt_i else float(volt_i)
                            numev_i = line.split()[8] if IsInt(line.split()[0]) else line.split()[9]
                            numev_i = 1000 * int(numev_i.strip('K').strip('k')) if 'K' in numev_i or 'k' in numev_i else int(numev_i)
                            temp_runs_info[run_i] = {'dia': diam_i, 'bias': volt_i, 'num_ev': numev_i}
                            print run_i, diam_i, volt_i, numev_i
                run_raw_paths = glob()

    def GetRunInfoJson(self):
        if os.path.isfile('{d}/{f}.json'.format(d=self.currents_logs_dir, f=self.testbeam_log_file)):
            temp_runs_info = js.load('{d}/{f}.json'.format(d=self.currents_logs_dir, f=self.testbeam_log_file))
            return OrderedDict(sorted(temp_runs_info.iteritems()))

    def RecreateRunsInfo(self):
        self.LoadRunsInfo(True)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-s', '--settings', dest='settings_f', type='string', help='Settings file containing the information of the currents and hv runlogs')
    parser.add_option('-d', '--diamond', dest='diamond', type='string', default='', help='diamond to be analysed. If not given, the default is empty string')
    parser.add_option('-p', '--polarity', dest='polarity', type='string', default='both', help='selects the polarity to analyse. The valid options are: negative, positive, both. The default is "both"')
    parser.add_option('-a', '--automatic', dest='automatic', default=False, action='store_true', help='enables automatic analysis with a given diamond')

    (options, args) = parser.parse_args()
    settings_f = str(options.settings_f)
    diam = str(options.diamond)
    pola = str(options.polarity)
    autom = bool(options.automatic)

    cu = Currents(settings_f, diam, pola)

    if autom:
        pass

