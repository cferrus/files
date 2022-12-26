"""
UpdateH5DataVersionUtils

This file contains the version_log, which holds the version information for
each H5 data file. This file also contains a function for each H5 data file
to apply the update to the data from one version to the next.
"""

import numpy as np
import h5py

#-------------------------------------------------------------------------------
# VERSION LOG
#-------------------------------------------------------------------------------

# UpdateH5DataVersion.py will only recognize updates and files that are recorded
# in this version_log.

version_log = {}

version_log['CCE'] = [
    # For CCE files with filename CceR####.h5
    ['d54f09c4380a4e63aea3b5e6f6803167218a4cbe',
     'Bugfix in CCE radial derivatives (ticket 1096).']
]

version_log['rh_FiniteRadii_CodeUnits.h5'] = [
    ['ef51849550a1d8a5bbdd810c7debf0dd839e86dd',
     'Overall sign change in complex strain h.']
]

version_log['PhiPlus_FiniteRadii_CodeUnits.h5'] = [
    ['75c9eaf42783ccc4ee4eea6b15e8ac5f457c2180',
     'Overall sign change in PhiPlus.']
]

version_log['rhOverM_Asymptotic_GeometricUnits.h5'] = [
    ['ef51849550a1d8a5bbdd810c7debf0dd839e86dd',
     'Overall sign change in complex strain h.']
]

version_log['rhOverM_Asymptotic_GeometricUnits_CoM.h5'] = [
    ['ef51849550a1d8a5bbdd810c7debf0dd839e86dd',
     'Overall sign change in complex strain h.']
]


#-------------------------------------------------------------------------------
# VERSION UPDATE FUNCTIONS
#-------------------------------------------------------------------------------

# The function update_to_most_recent_version() defined in UpdateH5DataVersion.py
# will update a file, e.g. rh_FiniteRadii_CodeUnits.h5, by updating the data
# incrementally (version 0 to 1, then 1 to 2, etc.) until it is on the most
# recent version. The updated data will be stored in a new file, e.g.
# new_rh_FiniteRadii_CodeUnits.h5, and the original file will be unchanged.
#
# In this utils file, we create a dictionary called apply_update that holds
# incremental version update functions which provide instructions on how to take
# the data in a file from one version to the next. There should be one of these
# functions for each new version in each file. To call the appropriate version
# update function, you give the dictionary as keys the filename and the version
# that the file is currently on. For example, if you call:
#
#   > apply_update['rh_FiniteRadii_CodeUnits.h5',0](f_in, f_out)
#
# then this will assume the data in f_in is rh data and apply the update to take
# it from version 0 to 1. This updated data will be stored in a new file called
# f_out. The file f_in will have no change.
#
# To make it simpler to write version update functions, we take care of all the
# file handling with a function called update_wrapper() defined below, which
# wraps the version update functions before they are stored in the dictionary.
# The reason for this wrapper file is that the new file for storing updated data
# is created only when the first incremental update is applied. However, if the
# overall update involves many incremental updates then we don't want to keep
# making a new file for each incremental update. Instead, we keep updating and
# overwriting the data in the new file.
#
# So when writing these version update functions, expect f_in to be an H5 file
# opened with read-only access and write all changes to f_out. All complications
# of creating new files or overwriting is handled elsewhere.
# NOTE: The new file is a copy of the original file, so the version update
# functions have to worry only about the data that must be changed.
#
# The following is an example for adding more updates/files:
#
# def Example_Datafile__0_to_1(f_in, f_out):
#   [MAKE CHANGES]
#   return
# apply_update['Example_Datafile.h5',0] = \
#   lambda A,B : update_wrapper(A,B,
#       Example_Datafile__0_to_1)
#
# def Example_Datafile__1_to_2(f_in, f_out):
#   [MAKE CHANGES]
#   return
# apply_update['Example_Datafile.h5',1] = \
#   lambda A,B : update_wrapper(A,B,
#       Example_Datafile__1_to_2)


apply_update = {}


def PhiPlus_FiniteRadii_CodeUnits__0_to_1(f_in, f_out):
    for group in f_in:
        for mode in f_in[group]:
            if mode.startswith('Y'):
                idx = group + '/' + mode
                f_out[idx][:, 1:] = np.negative(f_in[idx][:, 1:])
    return
apply_update['PhiPlus_FiniteRadii_CodeUnits.h5', 0] = \
    lambda A, B: update_wrapper(A, B,
                                PhiPlus_FiniteRadii_CodeUnits__0_to_1)


def rh_FiniteRadii_CodeUnits__0_to_1(f_in, f_out):
    for group in f_in:
        for mode in f_in[group]:
            if mode.startswith('Y'):
                idx = group + '/' + mode
                f_out[idx][:, 1:] = np.negative(f_in[idx][:, 1:])
    return
apply_update['rh_FiniteRadii_CodeUnits.h5', 0] = \
    lambda A, B: update_wrapper(A, B,
                                rh_FiniteRadii_CodeUnits__0_to_1)


def rhOverM_Asymptotic_GeometricUnits__0_to_1(f_in, f_out):
    for group in f_in:
        for mode in f_in[group]:
            if mode.startswith('Y'):
                idx = group + '/' + mode
                f_out[idx][:, 1:] = np.negative(f_in[idx][:, 1:])
    return
apply_update['rhOverM_Asymptotic_GeometricUnits.h5', 0] = \
    lambda A, B: update_wrapper(A, B,
                                rhOverM_Asymptotic_GeometricUnits__0_to_1)


def rhOverM_Asymptotic_GeometricUnits_CoM__0_to_1(f_in, f_out):
    for group in f_in:
        for mode in f_in[group]:
            if mode.startswith('Y'):
                idx = group + '/' + mode
                f_out[idx][:, 1:] = np.negative(f_in[idx][:, 1:])
    return
apply_update['rhOverM_Asymptotic_GeometricUnits_CoM.h5', 0] = \
    lambda A, B: update_wrapper(A, B,
                                rhOverM_Asymptotic_GeometricUnits_CoM__0_to_1)


#-------------------------------------------------------------------------------
# Needed for TestUpdateH5DataVersion
#-------------------------------------------------------------------------------

version_log['TestUpdateH5DataVersionScript.h5'] = [
# This data file is used in TestUpdateH5DataVersion
    ['0000000000000000000000000000000000000001',
     'Change all signs.'],
    ['0000000000000000000000000000000000000002',
     'Add 10 to each data element.'],
    ['0000000000000000000000000000000000000003',
     'Multiply each data element by 2.']
]

def TestUpdateH5DataVersionScript__0_to_1(f_in, f_out):
    for group in f_in:
        for mode in f_in[group]:
            if mode.startswith('Y'):
                idx = group + '/' + mode
                f_out[idx][:, 1:] = np.negative(f_in[idx][:, 1:])
    return
apply_update['TestUpdateH5DataVersionScript.h5', 0] = \
    lambda A, B: update_wrapper(A, B,
                                TestUpdateH5DataVersionScript__0_to_1)


def TestUpdateH5DataVersionScript__1_to_2(f_in, f_out):
    for group in f_in:
        for mode in f_in[group]:
            if mode.startswith('Y'):
                idx = group + '/' + mode
                f_out[idx][:, 1:] = 10 + f_in[idx][:, 1:]
    return
apply_update['TestUpdateH5DataVersionScript.h5', 1] = \
    lambda A, B: update_wrapper(A, B,
                                TestUpdateH5DataVersionScript__1_to_2)


def TestUpdateH5DataVersionScript__2_to_3(f_in, f_out):
    for group in f_in:
        for mode in f_in[group]:
            if mode.startswith('Y'):
                idx = group + '/' + mode
                f_out[idx][:, 1:] = 2 * f_in[idx][:, 1:]
    return
apply_update['TestUpdateH5DataVersionScript.h5', 2] = \
    lambda A, B: update_wrapper(A, B,
                                TestUpdateH5DataVersionScript__2_to_3)


#-------------------------------------------------------------------------------
# Helper functions to be used in the version update functions
#-------------------------------------------------------------------------------

def update_wrapper(input_file, output_file, update_version_work):
    """
    A wrapper around the update functions for file handling
    """
    if input_file != output_file:
        # For writing updated data to a new file
        copy_H5file(input_file, output_file)
        with h5py.File(input_file, 'r') as f_in, \
             h5py.File(output_file, 'r+') as f_out:

            update_version_work(f_in, f_out)
    else:
        # For overwriting data
        with h5py.File(input_file, 'r+') as f_in:
            f_out = f_in
            update_version_work(f_in, f_out)
    return


def copy_H5file(input_file, output_file):
    """
    Copies a data file to a new location, creating directories as needed.
    """
    from shutil import copyfile
    from os import makedirs
    from os.path import exists, dirname
    if not exists(dirname(output_file)) and not dirname(output_file) == '':
        makedirs(dirname(output_file))
    # We want copyfile, not copy2, because we don't want
    # to preserve permissions (permissions of annex files are read-only).
    copyfile(input_file, output_file)
    return
