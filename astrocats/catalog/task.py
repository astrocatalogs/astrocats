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
        Proper name for this task - used when calling it from command-line
    nice_name : str or None
        Name for pretty printing
    update : bool
        Whether this task should be performed during update.
    archived : bool
        Use existing, archived (meta-)data for this task.
    active : bool
        Whether this task should be performed or not by default.
    module : str or None
        Module in which to find `function` for carrying out this task
    groups : list of strings
        Which task groupings this task belongs to.  Allows numerous tasks to
        be run together.
    repo : str or None
        Repository in which the input-/meta- data is stored.  This is *not*
        the final output directory of the entry json file.
    function : str
        Function to execute when carrying out this task.
    priority : int
        Order in which tasks should be executed

    """

    def __init__(self, **kwargs):
        """Class initializer.

        Only existing class attributes are able to be passed and set via the
        constructor `kwargs` dictionary.  Otherwise a `ValueError` is raised.
        """
        self.name = None
        self.nice_name = None
        self.update = False
        self.archived = False
        self.active = True
        self.module = None
        self.groups = None
        self.repo = None
        self.function = ''
        self.priority = None
        self.always_journal = False

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
                  "priority='{}', always_journal='{}'")
        retval = retval.format(self.name, self.nice_name, self.active,
                               self.update, self.archived, self.module,
                               self.function, self.repo, self.priority,
                               self.always_journal)
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
        import warnings
        warnings.warn("`Task.load_archive()` is deprecated!  "
                      "`Catalog.load_url` handles the same functionality.")
        return self.archived or args.archived

    def _get_repo_path(self, base_path):
        """
        """
        return os.path.join(base_path, self.repo, '')
