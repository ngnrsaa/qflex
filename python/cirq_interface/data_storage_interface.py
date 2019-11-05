import tempfile
import os


class DataStorageInterface:

    def __init__(self, use_temp_files=True):

        self.use_files = use_temp_files

        if self.use_files:
            # Behind the scene, this class creates a temporary file for each object
            self._file_handle = tempfile.mkstemp()
        else:
            print("not implemented")

    def __del__(self):

        # The destructor removes the temporary file
        if not self.use_files:
            return

        try:
            # https://stackoverflow.com/questions/36069102/pytest-fixture-finalizer-typeerror-nonetype-object-is-not-callable
            # yield
            os.close(self._file_handle[0])
        except OSError as e:
            if e.errno == 9:
                # if it was closed before
                pass
            else:
                raise e
        except:
            raise ValueError("Something is very wrong!")

        # remove the temporary file from disk
        os.remove(self._file_handle[1])
