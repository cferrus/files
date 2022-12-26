#!/usr/bin/env python

"""
UpdateH5DataVersion

A set of useful tools to handle versioning of waveform files.
"""

from __future__ import print_function

import argparse
import h5py
from os import listdir, walk
from os.path import basename, isdir, isfile, join
from Utils import error

# This imports the version_log and functions needed for data updates:
from UpdateH5DataVersionUtils import *


def get_version(filename):
    """
    Returns the current version number of a data file.
    """
    with h5py.File(filename, 'r') as h5file:
        if 'VersionHist.ver' in h5file:
            version = h5file['VersionHist.ver'].shape[0]
        else:
            version = 0
    return version


def get_version_history(filename):
    """
    Returns version history of a data file as a list of pairs for each
    version change. The first element in the pair is the git commit ID
    before the change and the second element is a description of the change.
    """
    with h5py.File(filename, 'r') as h5file:
        if 'VersionHist.ver' in h5file:
            version_hist = h5file['VersionHist.ver'][:]
            # C++ versioning encodes the VersionHist.ver elements as byte strings
            # so for py2/py3 compatibility we need to convert them to strings as
            # understood by python.
            version_hist = version_hist.astype(str)
        else:
            version_hist = None
    return version_hist


def is_CCE_file_format(filename):
    """
    Checks if the file is a CCE file.
    """
    # For UpdateH5DataVersion to handle more general CCE filenames, update this
    # script to handle names of updated CCE files which can't be identical to
    # the names of not-updated CCE files.
    from re import match
    if match("^CceR[0-9]{4}\.h5$", filename):
        return True
    else:
        return False


def find_newest_available_version(filename):
    """
    Returns the newest available version number of a data file
    from the version_log in UpdateH5DataVersionUtils.
    """
    newest_version = 0
    try:
        if is_CCE_file_format(filename):
            newest_version = len(version_log['CCE'])
        else:
            newest_version = len(version_log[basename(filename)])
    except KeyError:
        newest_version = 0
    return newest_version


def set_COM_params_as_H5_attributes(filename):
    """
    In a COM-corrected data file, the boost velocity and space translation
    parameters used for the COM correction can be found in the History.txt 
    datasets within the data file. 

    This function copies these two COM parameters and creates two HDF5 
    attributes for each HDF5 group that has a History.txt dataset:
        Key               : Value
        boost_velocity    : np.array([v_x, v_y, v_z])
        space_translation : np.array([ dx,  dy,  dz])
    """
    from re import search
    from ast import literal_eval
    tol = 1e-8 # Precision of COM parameters in History.txt 
    with h5py.File(filename, 'r+') as h5file:
        for group in h5file:
            if not isinstance(h5file[group],h5py.Group):
                continue
            try:
                hist = h5file[group+'/History.txt'][()]
            except KeyError:
                error("Attempted to find COM parameters, but no 'History.txt' "
                      "dataset was found in "+filename)

            for param_name in ['boost_velocity', 'space_translation']:
                pattern = r"'{0}': array\((.*?)\)".format(param_name)
                matches = search(pattern, hist)
                if matches:
                    COM_params = literal_eval(matches.group(1))
                else:
                    error("The COM "+param_name+" parameters were not "
                          "found in "+filename)

                preexisting_attr = h5file[group].attrs.get(param_name)
                if preexisting_attr is None:
                    h5file[group].attrs.create(param_name, COM_params)
                elif (np.abs(preexisting_attr-COM_params) >= tol*np.abs(COM_params)).any():
                    error("The attribute '"+param_name+"' of group '"+group+"' "
                          "in the file "+filename+" does not match the values "
                          "reported in the History.txt dataset.")
    return


def set_version(filename, version_hist):
    """
    Creates or modifies the VersionHist.ver dataset. It is up to the
    caller to make sure that VersionHist.ver is set consistently with
    the data in the file.
    """
    with h5py.File(filename, 'r+') as h5file:
        # encoding the version_hist as ascii ensures that the strings are
        # encoded as bytes in the file for both py2 and py3 and this is 
        # what the C++ versioning writes and expects in the file.
        reencoded_version_hist = [
            [n[0].encode('ascii','ignore'), n[1].encode('ascii','ignore')] 
            for n in version_hist]
        dset = h5file.create_dataset('VersionHist.ver',
                                     (len(version_hist), 2),
                                     maxshape=(None,2),
                                     data=reencoded_version_hist,
                                     dtype=h5py.special_dtype(vlen=bytes))
        dset.attrs.create("Version",len(version_hist),shape=(1,),
                          dtype=np.uint64)
    return


def ensure_same_version(list_of_h5filenames):
    """
    For use within CombineSegments.py

    Returns True if all files in list_of_h5filenames have the same
    version history
    """
    reference_hist = get_version_history(list_of_h5filenames[0])
    for h5filename in list_of_h5filenames[1:]:
        hist = get_version_history(h5filename)
        try:
            if hist != reference_hist:
                return False
        except ValueError:
            if (hist != reference_hist).any():
                return False
    return True

def copy_version_history(source_h5filename, dest_h5filename):
    """
    Copies the version history from a source file to a destination file
    """
    version_hist = get_version_history(source_h5filename)
    set_version(dest_h5filename, version_hist)
    return


def get_list_of_h5files(args):
    """
    Given the user's specifications, returns a list of all available h5 files.
    """
    from glob import glob
    if isfile(args.input) and args.input.endswith('.h5'):
        return [args.input]
    elif isdir(args.input):
        list_of_h5files = []
        list_of_h5files += glob('*.h5')
        if args.recursive:
            # Search subdirectories
            dirs = [result[0] for result in walk('.', followlinks=False)]
            dirs.remove('.')
            for dir in dirs:
                list_of_h5files += glob(dir+'/*.h5')
        if not list_of_h5files:
            print('No H5 files found in directory '+args.input)
            exit(0)
        else:
            return sorted(list_of_h5files)
    else:
        error("'"+args.input+"' is neither an H5 file nor a directory "
              "containing H5 files.")


def print_history(filename):
    """
    Prints the version history of an H5 data file.
    """
    if isfile(filename) and filename.endswith('.h5'):
        current_version = get_version(filename)
        print('Currently at version '+str(current_version))
        print('Version change log:')
        print('ver 0: Original')
        history = get_version_history(filename)
        if history is not None:
            for i in range(len(history)):
                # don't print commit ID
                print("ver "+str(i+1)+": "+history[i][1])
        return
    else:
        error('Error: '+filename+' is not an H5 file.')


def print_list_of_file_versions(filename, list_outdated):
    """
    Prints the current version of an H5 data file and checks if an update
    is available. If 'list_outdated' is 'True', then this will only list files
    with an available update.
    """
    current_version = get_version(filename)
    newest_version = find_newest_available_version(filename)
    if newest_version - current_version <= 0:
        new_version = ''
        version_highlight = ''
    else:
        new_version = '  (Outdated by version '+str(newest_version)+')'
        version_highlight = '*'
    if list_outdated and not new_version:
        return
    print('version '+str(current_version)+version_highlight +
          '\t'+filename+new_version)
    return


def update_to_most_recent_version(filename_in, filename_out):
    """
    Updates an H5 data file to the most recent version. The H5 data file
    will be copied and the update made to the output file.
    """
    # See the comment under the section Version Update Functions in
    # UpdateH5DataVersionUtils.py for full documentation on the how
    # updates are performed.
    current_version = get_version(filename_in)
    initial_version = current_version
    newest_version = find_newest_available_version(filename_in)
    if newest_version <= current_version:
        return

    print(filename_in+' is currently on version '+str(current_version))
    # Update versions incrementally
    for ver in range(current_version, newest_version):
        print('updating from version '+str(ver)+' to '+str(ver+1), end="...")
        try:
            if ver == initial_version:
                # The first incremental update should create a new file where the
                # updated data will be stored.
                apply_update[basename(filename_in), ver](
                    filename_in, filename_out)
            else:
                # All further incremental updates will update and overwrite the
                # data in the file created by the first incremental update.
                apply_update[basename(filename_in), ver](
                    filename_out, filename_out)
            current_version = ver+1
        except KeyError:
            set_version(filename_out,
                        version_log[basename(filename_in)][:ver])
            print('failed')
            error('No instructions to update '+basename(filename_in) +
                  ' from version '+str(ver)+' to '+str(ver+1)+'. '+filename_out +
                  ' is currently on version '+str(ver)+'.')
        print('done')
    print(filename_out+' is now on version '+str(current_version)+'.\n')
    if current_version > initial_version:
        if is_CCE_file_format(basename(filename_in)):
            # This branch will not be taken unless the update instructions for
            # CCE data are added to UpdateH5DataVersionUtils.py
            # This is fine though because the code that reads CCE files already
            # handles different versions.
            error('UpdateH5DataVersion does not currently handle updating CCE files.')

            # For UpdateH5DataVersion to handle updating CCE files, add appropriate
            # version update functions in UpdateH5DataVersionUtils.py, remove the
            # error message from the line above, and uncomment the following:
            #
            # set_version(filename_out,
            #             version_log['CCE'][:current_version])
        else:
            set_version(filename_out,
                        version_log[basename(filename_in)][:current_version])
    return

#-------------------------------------------------------------------------------
#  Main routine
#-------------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p2 = p.add_mutually_exclusive_group(required=True)
    p2.add_argument(
        '--history <filename>',
        dest='history',
        action='store_true',
        help='Print a description of the different versions of an H5 data file.')
    p2.add_argument(
        '--list-all <filename|dir>',
        dest='list_all',
        action='store_true',
        help='Print version of <filename> or of all H5 data files in <dir>, and '
             'mark files with available updates. ')
    p2.add_argument(
        '--list-outdated <filename|dir>',
        dest='list_outdated',
        action='store_true',
        help='Same as --list-all but ONLY list files with available updates.')
    p2.add_argument(
        '--update <filename|dir>',
        dest='update',
        action='store_true',
        help='Update files to latest version if possible.')
    p.add_argument(
        'input',
        metavar='filename/dir',
        type=str,
        help=argparse.SUPPRESS)
    p.add_argument(
        '--output',
        '-o',
        metavar='output',
        type=str,
        required=False,
        help='Specify output filename or the directory to output new files into.')
    p.add_argument(
        '--recursive',
        '-r',
        dest='recursive',
        action='store_true',
        default=False,
        help='Perform for all H5 files in the directory and subdirectories.')
    args = p.parse_args()
    if args.update:
        if not args.output:
            error('Error: must specify --output/-o with option --update')
        if args.output == args.input:
            error('Error: output cannot be the same as input')

    list_of_h5files = get_list_of_h5files(args)

    if args.history:
        print_history(args.input)

    else:
        for h5file in list_of_h5files:
            if args.list_all or args.list_outdated:
                print_list_of_file_versions(h5file, args.list_outdated)
            elif args.update:
                if len(list_of_h5files) == 1:
                    update_to_most_recent_version(args.input, args.output)
                else:
                    update_to_most_recent_version(h5file,
                                                  join(args.output, basename(h5file)))


if __name__ == "__main__":
    main()
