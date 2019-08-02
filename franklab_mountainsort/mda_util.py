'''run this from the preprocessing directory containing all dates' data'''
import errno
import os


def make_mda_ntrodeEpoch_links():
    # for each date directory
    for datedir in os.listdir(os.path.curdir):
        date = datedir.split('_')[0]
        # for each ep.mda directory
        for epdirmda in os.listdir(os.path.join(os.path.curdir, date)):
            if '.mda' in epdirmda:
                # for each nt.mda file
                for eptetmda in os.listdir(
                        os.path.join(os.path.curdir, date, epdirmda)):
                    if '.nt' in eptetmda:
                        an = eptetmda.split('_')[1]
                        endf = eptetmda.split('_')[-1]
                        ntr = endf.split('.')[1]
                        cwd = os.getcwd()
                        srclink = os.path.join(
                            cwd, datedir, epdirmda, eptetmda)
                        mntdir = f'{date}_{an}.mnt'
                        ntdir = f'{date}_{an}.{ntr}.mnt'
                        destlink = os.path.join(
                            cwd, datedir, mntdir, ntdir, eptetmda)
                        make_sure_path_exists(
                            os.path.join(cwd, datedir, mntdir))
                        make_sure_path_exists(
                            os.path.join(cwd, datedir, mntdir, ntdir))
                        # to overwrite. remove ntlink if it already exists
                        removeNTfile(destlink)
                        # create directory of sym links to original mda
                        # os.symlink(srclink, destlink)
                        os.symlink(os.path.relpath(
                            srclink, destlink), destlink)


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def removeNTfile(ntrode_filename):
    import os
    import errno
    try:
        os.remove(ntrode_filename)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred


# def createSymlink():
if __name__ == "__main__":
    make_mda_ntrodeEpoch_links()
