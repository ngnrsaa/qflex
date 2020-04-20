import tempfile
import os


class DataStorageInterface:

    def __init__(self, use_temp_files=True):

        self.use_files = use_temp_files

        if self.use_files:
            # Behind the scene, this class creates a temporary file for each object
            # Need to keep a reference in the object to the underlying temp directory
            self._tdir = tempfile.TemporaryDirectory()
            self.storage = tempfile.NamedTemporaryFile(mode='w',
                                                       dir=self._tdir.name,
                                                       delete=False)
            self.fullpath = self.storage.name
        else:
            print("not implemented")
