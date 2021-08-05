
import os
import os.path
import subprocess
import glob

from ..build import _Build

class _BuildKepler(_Build):

    package = __package__
    path = os.path.dirname(__file__)

    sources = ('kepler.f90', 'data.f')
    objects = ('kepler.o', 'data.o')
    libraries = ('-luuid',)

    # my own stuff
    kepler_source_path  = os.path.join(os.environ['KEPLER_PATH'], 'source')
    kepler_library_path = os.path.join(path, 'lib')
    kepler_library_file = os.path.join(
        kepler_library_path, 'kepler.a')

    project_libraries = (kepler_library_file,)
    include_paths = (
        kepler_source_path,
        kepler_library_path,
        )
    signature_file = '_kepler.pyf'
    module = '_kepler'
    executable_file = 'kepler.exe'

    executable_link_flags = ('-fconvert=big-endian',)
    compile_flags = _Build.compile_flags + executable_link_flags

    def __init__(self):
        """
        Note - all initilisation needed is done in class definition for now.
        Maybe that should go here instead ...
        """
        super().__init__()

    def build_library_check(self, debug = True):
        """
        check whether KEPLER library is up to date
        """
        try:
            library_time = os.path.getctime(self.kepler_library_file)
        except FileNotFoundError:
            library_time = 0

        exclude = ('uuidcom', 'gitcom', 'nburncom', 'gridcom', )
        patterns = ('*com', '*.f', '*.f90', '*.c', 'Makefile*', )
        f = os.path.join(self.kepler_library_path, 'Makefile')
        last_time = os.path.getctime(f)
        for p in patterns:
            for f in glob.glob(os.path.join(self.kepler_source_path, p)):
                if os.path.basename(f) in exclude:
                    continue
                last_time = max(last_time, os.path.getctime(f))
        if last_time > library_time:
            cwd = os.getcwd()
            os.chdir(self.kepler_library_path)
            subprocess.run(['make', '-j'], shell = False, check = True)
            os.chdir(cwd)
            return True
        return False
