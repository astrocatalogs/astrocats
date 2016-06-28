class Task:
    name = ''
    nice_name = ''   # Name for pretty printing
    # Perform task during update
    update = False
    # Use archived data ???
    archived = False
    active = True    # Whether this task should be performed or not
    # Module in which to find `function` for carrying out this task
    module = None
    function = ''    # Function to execute when carrying out this task
    priority = None  # Order in which tasks should be executed

    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise ValueError("No attribute '{}'".format(key))

    def __repr__(self):
        retval = ("Task(name='{}', nice_name='{}', active='{}', update='{}', "
                  "archived='{}', module='{}', function='{}', priority='{}'")
        retval = retval.format(self.name, self.nice_name, self.active,
                               self.update, self.archived, self.module,
                               self.function, self.priority)
        return retval

    def current_task(self, args):
        """Name of current action for progress-bar output, depends on run
        configuration.
        """
        ctask = self.nice_name if self.nice_name else self.name
        if args is not None:
            if args.update:
                ctask = ctask.replace('%pre', 'Updating')
            else:
                ctask = ctask.replace('%pre', 'Loading')
        return ctask

    def load_archive(self, args):
        """Depending on run configuration, whether previously archived data
        should be loaded.
        """
        if not self.archived:
            return False
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
