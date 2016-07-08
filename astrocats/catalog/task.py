"""Object for each 'task' needed to be carried out by the 'importer'.
"""
import os


class Task():
    """General class for individual catalog operations.

    Currently only used in the import methods.
    Each task instance is loaded with values from the `input/tasks.json` file.

    Attributes
    ----------
    name : str
    nice_name : str or None
    update : bool
    archived : bool
    active : bool
    module : str or None
    groups : list of strings
    repo : str or None
    function : str
    priority : int

    """

    def __init__(self, **kwargs):
        """Class initializer.

        Only existing class attributes are able to be passed and set via the
        constructor `kwargs` dictionary.  Otherwise a `ValueError` is raised.
        """
        # Proper name for this task - used when calling it from command-line
        self.name = None
        # Name for pretty printing
        self.nice_name = None
        # Perform task during update
        self.update = False
        # Use archived data ???
        self.archived = False
        self.active = True    # Whether this task should be performed or not
        # Module in which to find `function` for carrying out this task
        self.module = None
        self.groups = None
        # Repository in which the data comes from
        self.repo = None
        self.function = ''    # Function to execute when carrying out this task
        self.priority = None  # Order in which tasks should be executed

        for key, val in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise ValueError("No attribute '{}'".format(key))

        if self.groups is not None:
            self.groups = [group.lower().strip() for group in self.groups]

        return

    def __repr__(self):
        retval = ("Task(name='{}', nice_name='{}', active='{}', update='{}', "
                  "archived='{}', module='{}', function='{}', repo='{}', "
                  "priority='{}'")
        retval = retval.format(self.name, self.nice_name, self.active,
                               self.update, self.archived, self.module,
                               self.function, self.repo, self.priority)
        return retval

    def current_task(self, args):
        """Name of current action for progress-bar output.

        The specific task string is depends on the configuration via `args`.

        Returns
        -------
        ctask : str
            String representation of this task.
        """
        ctask = self.nice_name if self.nice_name is not None else self.name
        if args is not None:
            if args.update:
                ctask = ctask.replace('%pre', 'Updating')
            else:
                ctask = ctask.replace('%pre', 'Loading')
        return ctask

    def load_archive(self, args):
        """Whether previously archived data should be loaded.
        """
        # If we're running in 'archived' mode, and only loading 'archived'
        # things, then True
        if (args.archived and self.name not in args.refresh_list and not
                args.full_refresh):
            return True
        # For normal running, if we are not sepcifically refreshing this task,
        # then True
        if self.name not in args.refresh_list and not args.full_refresh:
            return True

        return False

    def _get_repo_path(self, base_path):
        """
        """
        return os.path.join(base_path, self.repo, '')
