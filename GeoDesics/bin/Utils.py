# Python module containing useful wrappers of system utilities

from __future__ import print_function
import os
import sys
import multiprocessing


def ReadH5(file):
    """Opens the file 'file' as read-only H5 file, and returns the
h5py.File object.
"""
    import h5py
    try:
        f = h5py.File(file, 'r')
    except IOError:
        os.sys.stderr.write("\n############# ERROR #############\n")
        os.sys.stderr.write('Attempt to open h5-file "' + file + '" failed.\n')
        raise
    return f


################################################################


# Recursively remove empty subdirectories
def RemoveEmptySubdirs(path):
    if not os.path.isdir(path):
        return

    files = os.listdir(path)
    if len(files):
        for f in files:
            fullpath = os.path.join(path, f)
            if os.path.isdir(fullpath):
                RemoveEmptySubdirs(fullpath)

    files = os.listdir(path)
    if len(files) == 0:
        os.rmdir(path)


################################################################


def ReadDat(file, cols=None):
    """Read an ASCII .dat file, and return the specified columns.
file -- filename of the file to read
cols -- columns to read, in the format expected from numpy.loadtxt's usecols
        option, namely an iterable listing the columns, like (1,4,6,5)
        If ==None, return all columns
RETURNS:
  a numpy array
"""
    from numpy import loadtxt
    try:
        return loadtxt(file, usecols=cols, ndmin=2)  # always return 2D
    except:
        os.sys.stderr.write("\n############# ERROR #############\n")
        os.sys.stderr.write('Attempt to read "' + file +
                            '" as dat file failed.\n')
        raise


# Subroutine to safely run a system bash command
def System(cmd, RequireZeroExitStatus=True, verbose=False):
    from subprocess import call
    if verbose:
        print("System running: {}".format(cmd))
    status = call(cmd, shell=True, executable='/bin/bash')
    if (status != 0) and RequireZeroExitStatus:
        error("The command `" + cmd + "` failed with return status " +
              str(status) + " in " + os.getcwd())
    return status


# Subroutine to safely run a system bash command and capture output
def SystemOutput(cmd, RequireZeroExitStatus=True, verbose=False):
    from subprocess import Popen, PIPE
    if verbose:
        print("SystemOutput running: {}".format(cmd))
    p = Popen(cmd, shell=True, stdout=PIPE)
    stdout, stderr = p.communicate()
    status = p.returncode
    if (status != 0) and RequireZeroExitStatus:
        error("The command `"+cmd+"` failed with return status "+str(status)+\
              " in "+os.getcwd()+"\nStdout=`"+\
              (stdout.rstrip() if stdout else '')+"`\nStderr=`"+\
              (stderr.rstrip() if stderr else '')+"`")
    return stdout.decode("utf-8").rstrip()


def error(msg):
    e_msg = "\n############# ERROR #############\n{}\n".format(msg)
    raise Exception(e_msg)


# Print a warning message to stderr
def warning(msg):
    os.sys.stderr.write("\n############# WARNING #############\n")
    os.sys.stderr.write("{}\n".format(msg))


################################################################


def FirstIndexGEQValue(target, list):
    """Return index of the first element in list with value >= 'target'
"""
    return next(idx for idx, val in enumerate(list) if val >= target)


################################################################


def call_perl(module, func_call, want_array):
    """Call a function with 'func_call' in perl module 'module'.
If want_array=False, calls 'func_call' in scalar context, returning a string.
If want_array=True, calls 'func_call' in array context, returning a list.
WARNING: all perl functions can be silently forced into scalar context.
"""
    import subprocess
    perldir = os.path.realpath(__file__ + "/../../Perl")
    bindir = os.path.realpath(__file__ + "/../../bin")

    # Run the command, using the null character to join elements
    if want_array:
        cmd = '@x=%s::%s; print join("\\0",@x);' % (module, func_call)
    else:
        cmd = '$x=%s::%s; print $x;' % (module, func_call)
    p = subprocess.Popen(
        ['perl', '-I', perldir, '-I', bindir, '-M' + module, '-e', cmd],
        stdout=subprocess.PIPE)
    output, _ = p.communicate()
    status = p.returncode
    if (status != 0):
        error("Perl exited with error " + str(status) + ".")

    # Return the output
    if want_array:
        return output.decode("utf-8").split('\0')
    else:
        return output.decode("utf-8")


################################################################


def SetMaxPrecision(val, p):
    """Sets the precision of a numerical string.
Return a float string with precision p, or just an int string if the
value is indistinguishable from an integer within precision p.
"""
    # TODO: remove ALL trailing zeros, e.g. SetMaxPrecision(0.5005,3)=0.5
    from math import log10
    val = float(val)
    Frac = abs(val - round(val))
    if Frac == 0 or log10(Frac) < -p * 1.00001:  #fudge factor
        val = int(round(val))
        p = 0
    return "%0.*f" % (p, val)


################################################################


def norm(data, axis=None):
    """Returns the L2-norm over a specific axis. The axis option
behaves exactly like in numpy.sum: if axis=None (default), then
sum over all axes. If negative, count from the last axis.
"""
    import numpy as np
    return np.sqrt(np.sum(np.square(data), axis=axis))


################################################################


def SmoothData(t,
               data,
               width,
               Deriv=0,
               conv_factor=4.,
               cut_off=None,
               Tstart=None,
               DT=None,
               dt_tol=1e-12):
    """Smooth the data (t,data) using Gaussian convolution.
Derivatives do not seem to work well when timesteps are much larger than 1.
It also does not work (though it should) for unequally spaced timesteps.
The smooth data will be truncated on both sides based on the smoothing width.
"""
    import math
    import numpy as np

    # Convert to expected objects
    width = float(width)
    try:
        data.shape
    except AttributeError:
        data = np.array(data)
    if len(data.shape) == 1:
        data = data[:, np.newaxis]  #promote the array with singleton dim
    if len(data.shape) != 2:
        error("data is expected to be 2d array, not {}d".format(len(
            data.shape)))

    try:
        t.shape
    except AttributeError:
        t = np.array(t)

    convolution_width = width * conv_factor
    if cut_off == None:
        cut_width = width * conv_factor
    else:
        cut_width = width * cut_off

    # compute average time-spacing
    def ComputeAvgDt(t, dt_tol):
        dt = np.diff(t)
        aveDt = np.mean(dt)
        MaxDtDiff = max(abs(dt - aveDt))
        if MaxDtDiff > dt_tol:
            error(
                "Times are not equally spaced within dt_tol={}.\n".format(
                    dt_tol) +
                "Maximum deviation from the average is {}.".format(MaxDtDiff))
        return aveDt

    aveDt = ComputeAvgDt(t, dt_tol)

    # set up moving window
    if (t[-1] - t[0]) <= 6 * width:
        error("Total time-interval must be at least 6x larger than width")

    # setup output times
    if DT == None:
        t_out = t
    else:
        Tstart = t[0] if (Tstart is None) else Tstart
        if Tstart < t[0]:
            warn("--Tstart smaller than first time-step %g" % t[0])
        t_out = np.arange(Tstart, t[-1], step=DT)

    def K(x, width):
        return 1. / np.sqrt(math.pi) / width * np.exp(-x * x / (width * width))

    def K1(x, width):
        return -2. * x / np.sqrt(np.pi) / width**3 * np.exp(-x * x /
                                                            (width * width))

    def K2(x, width):
        return (4*x*x-2*width**2)/np.sqrt(np.pi)/width**5 \
               * np.exp(-x*x/(width**2))

    idxmin = 0
    idxmax = 0

    # limit t_out to range so that the kernel fits
    t_out = t_out[(t_out >= t[0] + cut_width) * (t_out <= t[-1] - cut_width)]

    output_array = np.empty((len(t_out), 1 + len(data[0])))
    for curr_idx in range(len(t_out)):
        curr_t = t_out[curr_idx]
        # update bounds for convolution
        while (t[idxmin] < curr_t - convolution_width):
            idxmin = idxmin + 1
        while idxmax < len(t) - 1 and t[idxmax] < curr_t + convolution_width:
            idxmax = idxmax + 1

        if (t[idxmax] > curr_t + convolution_width):
            idxmax = idxmax - 1

        if (abs(curr_t - t[0]) < cut_width or abs(curr_t - t[-1]) < cut_width):
            continue

        if (idxmax == idxmin):
            error("For t=%f, no data-points in interval.\n" % curr_t)

        # compute sums
        # s = Sum Kf,  w=Sum K
        f12 = K(curr_t - t[idxmin:idxmax + 1], width)
        w = np.sum(aveDt * (f12[0:-1] + f12[1:]))
        s = np.sum(aveDt * (data[idxmin:idxmax] * f12[0:-1].reshape(
            (-1, 1)) + data[idxmin + 1:idxmax + 1] * f12[1:].reshape((-1, 1))),
                   axis=0)

        # s1= Sum K1f, w1=Sum K1
        if (Deriv >= 1):
            f12 = K1(curr_t - t[idxmin:idxmax + 1], width)
            w1 = np.sum(aveDt * (f12[0:-1] + f12[1:]))
            s1 = np.sum(aveDt * (data[idxmin:idxmax] * f12[0:-1].reshape(
                (-1, 1)) + data[idxmin + 1:idxmax + 1] * f12[1:].reshape(
                    (-1, 1))),
                        axis=0)

        # s2= Sum K2f, w2=Sum K2
        if (Deriv >= 2):
            f12 = K2(curr_t - t[idxmin:idxmax + 1], width)
            w2 = np.sum(aveDt * (f12[0:-1] + f12[1:]))
            s2 = np.sum(aveDt * (data[idxmin:idxmax] * f12[0:-1].reshape(
                (-1, 1)) + data[idxmin + 1:idxmax + 1] * f12[1:].reshape(
                    (-1, 1))),
                        axis=0)

        ### combine for answer ###
        if (Deriv == 0):
            result = np.concatenate(((curr_t, ), s / w))
        elif (Deriv == 1):
            result = np.concatenate(((curr_t, ), s1 / w - s / w**2 * w1))
        elif (Deriv == 2):
            # Note: This formula produces visible noise
            # for non-uniformly spaced data even when the
            # 'w1' terms are included, which sould take
            # correct for non-uniform spacing.
            # (Harald, Nov 29, 2006)
            result = np.concatenate(((curr_t, ), s2 / w - 2. * s1 * w1 / w**2 -
                                     s * w2 / w**2 + 2. * s * w1 * w1 / w**3))
        output_array[curr_idx, :] = result[:]
    return output_array


################################################################


def run_parallel(func, coll, num_workers=None, **kwargs):
    """Run several workers to accomplish a task in parallel.

    Arguments:
    func -- The function to call. This must be a top-level function in a
            module.
    coll -- The collection of arguments to operate over. Must be iterable.
            These are passed as the first argument to func.
    num_workers -- The number of child processes to start. Defaults to the
                   number of concurrent threads a machine can execute. Setting
                   to 1 completely disables parallel infrastructure.

    All remaining keyword arguments get passed to each function call as-is.

    The return value will be equal to [func(x, **kwargs) for x in coll].
    """

    try:
        iter(coll)
    except TypeError:
        raise TypeError("Argument coll must be iterable.")

    if num_workers == 1:
        results = [func(x, **kwargs) for x in coll]
    else:
        pool = multiprocessing.Pool(num_workers)
        mpresults = [pool.apply_async(func, (x, ), kwargs) for x in coll]
        pool.close()

        try:
            # Needs to be a large delay because we don't actually want it to
            # timeout, but Python won't receive KeyboardInterrupt unless
            # we put in the timeout.
            # But don't put in a timeout that is so long that python
            # overflows!
            results = [mp.get(100) for mp in mpresults]
        except KeyboardInterrupt:
            pool.terminate()
            print("Received ^C, exiting...")
            sys.exit(1)

    return results
