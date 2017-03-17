"""Utilities for working with git.
"""

import os
import subprocess
from glob import glob

import git


class MyProgressPrinter(git.RemoteProgress):
    def update(self, op_code, cur_count, max_count=None, message=''):
        # print(op_code, cur_count, max_count, message, end="\r")
        msg = "{} - {}".format(op_code, cur_count)
        if max_count is not None:
            msg += "/{}".format(max_count)
        if len(msg):
            msg += " - {}".format(message)
        print(msg + "\r")


def fetch(repo_name, progress=True, log=None):
    repo = git.Repo(repo_name)
    _progress = MyProgressPrinter() if progress else None
    for fetch_info in repo.remote().fetch(progress=_progress):
        if log is not None:
            log.debug("Updated {}:{} to {}".format(repo_name, fetch_info.ref,
                                                   fetch_info.commit))


def get_sha(path=None, log=None, short=False, timeout=None):
    """Use `git rev-parse HEAD <REPO>` to get current SHA.
    """
    # git_command = "git rev-parse HEAD {}".format(repo_name).split()
    # git_command = "git rev-parse HEAD".split()
    git_command = ["git", "rev-parse"]
    if short:
        git_command.append("--short")
    git_command.append("HEAD")

    kwargs = {}
    if path is not None:
        kwargs['cwd'] = path
    if timeout is not None:
        kwargs['timeout'] = timeout

    if log is not None:
        log.debug("{} {}".format(git_command, str(kwargs)))

    sha = subprocess.check_output(git_command, **kwargs)
    try:
        sha = sha.decode('ascii').strip()
    except:
        if log is not None:
            log.debug("decode of '{}' failed".format(sha))

    return sha


def git_add_commit_push_all_repos(cat):
    """Add all files in each data repository tree, commit, push.

    Creates a commit message based on the current catalog version info.

    If either the `git add` or `git push` commands fail, an error will be
    raised.  Currently, if `commit` fails an error *WILL NOT* be raised
    because the `commit` command will return a nonzero exit status if
    there are no files to add... which we dont want to raise an error.
    FIX: improve the error checking on this.
    """
    log = cat.log
    log.debug("gitter.git_add_commit_push_all_repos()")

    # Do not commit/push private repos
    all_repos = cat.PATHS.get_all_repo_folders(private=False)
    for repo in all_repos:
        log.info("Repo in: '{}'".format(repo))
        # Get the initial git SHA
        sha_beg = get_sha(repo)
        log.debug("Current SHA: '{}'".format(sha_beg))

        # Get files that should be added, compress and check sizes
        add_files = cat._prep_git_add_file_list(repo,
                                                cat.COMPRESS_ABOVE_FILESIZE)
        log.info("Found {} Files to add.".format(len(add_files)))
        if len(add_files) == 0:
            continue

        try:
            # Add all files in the repository directory tree
            git_comm = ["git", "add"]
            if cat.args.travis:
                git_comm.append("-f")
            git_comm.extend(add_files)
            _call_command_in_repo(
                git_comm, repo, cat.log, fail=True, log_flag=False)

            # Commit these files
            commit_msg = "'push' - adding all files."
            commit_msg = "{} : {}".format(cat._version_long, commit_msg)
            log.info(commit_msg)
            git_comm = ["git", "commit", "-am", commit_msg]
            _call_command_in_repo(git_comm, repo, cat.log)

            # Add all files in the repository directory tree
            git_comm = ["git", "push"]
            if not cat.args.travis:
                _call_command_in_repo(git_comm, repo, cat.log, fail=True)
        except Exception as err:
            try:
                git_comm = ["git", "reset", "HEAD"]
                _call_command_in_repo(git_comm, repo, cat.log, fail=True)
            except:
                pass

            raise err

    return


def git_pull_all_repos(cat, strategy_recursive=True, strategy='theirs'):
    """Perform a 'git pull' in each data repository.

    > `git pull -s recursive -X theirs`
    """
    # raise RuntimeError("THIS DOESNT WORK YET!")
    log = cat.log
    log.debug("gitter.git_pull_all_repos()")
    log.warning("WARNING: using experimental `git_pull_all_repos()`!")

    all_repos = cat.PATHS.get_all_repo_folders()
    for repo_name in all_repos:
        log.info("Repo in: '{}'".format(repo_name))
        # Get the initial git SHA
        sha_beg = get_sha(repo_name)
        log.debug("Current SHA: '{}'".format(sha_beg))

        # Initialize the git repository
        repo = git.Repo(repo_name)
        # Construct the command to call
        git_comm = "git pull --verbose"
        if strategy_recursive:
            git_comm += " -s recursive"
        if strategy is not None:
            git_comm += " -X {:s}".format(strategy)
        log.debug("Calling '{}'".format(git_comm))
        # Call git command (do this manually to use desired options)
        #    Set `with_exceptions=False` to handle errors ourselves (below)
        code, out, err = repo.git.execute(
            git_comm.split(),
            with_stdout=True,
            with_extended_output=True,
            with_exceptions=False)
        # Handle output of git command
        if len(out):
            log.info(out)
        if len(err):
            log.info(err)
        # Hangle error-codes
        if code != 0:
            err_str = "Command '{}' returned exit code '{}'!".format(git_comm,
                                                                     code)
            err_str += "\n\tout: '{}'\n\terr: '{}'".format(out, err)
            log.error(err_str)
            raise RuntimeError(err_str)

        sha_end = get_sha(repo_name)
        if sha_end != sha_beg:
            log.info("Updated SHA: '{}'".format(sha_end))

    return


def git_clone_all_repos(cat):
    """Perform a 'git clone' for each data repository that doesnt exist.
    """
    log = cat.log
    log.debug("gitter.git_clone_all_repos()")

    all_repos = cat.PATHS.get_all_repo_folders()
    out_repos = cat.PATHS.get_repo_output_folders()
    for repo in all_repos:
        log.info("Repo in: '{}'".format(repo))

        if os.path.isdir(repo):
            log.info("Directory exists.")
        else:
            log.debug("Cloning directory...")
            clone(repo, cat.log, depth=max(cat.args.clone_depth, 1))

        if cat.args.purge_outputs and repo in out_repos:
            for fil in glob(os.path.join(repo, '*.json')):
                os.remove(fil)

        grepo = git.cmd.Git(repo)
        try:
            grepo.status()
        except git.GitCommandError:
            log.error("Repository does not exist!")
            raise

        # Get the initial git SHA
        sha_beg = get_sha(repo)
        log.debug("Current SHA: '{}'".format(sha_beg))

    return


def git_reset_all_repos(cat, hard=True, origin=False, clean=True):
    """Perform a 'git reset' in each data repository.
    """
    log = cat.log
    log.debug("gitter.git_reset_all_repos()")

    all_repos = cat.PATHS.get_all_repo_folders()
    for repo in all_repos:
        log.warning("Repo in: '{}'".format(repo))
        # Get the initial git SHA
        sha_beg = get_sha(repo)
        log.debug("Current SHA: '{}'".format(sha_beg))

        grepo = git.cmd.Git(repo)
        # Fetch first
        log.info("fetching")
        grepo.fetch()

        args = []
        if hard:
            args.append('--hard')
        if origin:
            args.append('origin/master')
        log.info("resetting")
        retval = grepo.reset(*args)
        if len(retval):
            log.warning("Git says: '{}'".format(retval))

        # Clean
        if clean:
            log.info("cleaning")
            # [q]uiet, [f]orce, [d]irectories
            retval = grepo.clean('-qdf')
            if len(retval):
                log.warning("Git says: '{}'".format(retval))

        sha_end = get_sha(repo)
        if sha_end != sha_beg:
            log.debug("Updated SHA: '{}'".format(sha_end))

    return


def git_status_all_repos(cat, hard=True, origin=False, clean=True):
    """Perform a 'git status' in each data repository.
    """
    log = cat.log
    log.debug("gitter.git_status_all_repos()")

    all_repos = cat.PATHS.get_all_repo_folders()
    for repo_name in all_repos:
        log.info("Repo in: '{}'".format(repo_name))
        # Get the initial git SHA
        sha_beg = get_sha(repo_name)
        log.debug("Current SHA: '{}'".format(sha_beg))

        log.info("Fetching")
        fetch(repo_name, log=cat.log)

        git_comm = ["git", "status"]
        _call_command_in_repo(
            git_comm, repo_name, cat.log, fail=True, log_flag=True)

        sha_end = get_sha(repo_name)
        if sha_end != sha_beg:
            log.info("Updated SHA: '{}'".format(sha_end))

    return


def clone(repo, log, depth=1):
    """Given a list of repositories, make sure they're all cloned.

    Should be called from the subclassed `Catalog` objects, passed a list
    of specific repository names.

    Arguments
    ---------
    all_repos : list of str
        *Absolute* path specification of each target repository.

    """
    kwargs = {}
    if depth > 0:
        kwargs['depth'] = depth

    try:
        repo_name = os.path.split(repo)[-1]
        repo_name = "https://github.com/astrocatalogs/" + repo_name + ".git"
        log.warning("Cloning '{}' (only needs to be done ".format(repo) +
                    "once, may take few minutes per repo).")
        grepo = git.Repo.clone_from(repo_name, repo, **kwargs)
    except:
        log.error("CLONING '{}' INTERRUPTED".format(repo))
        raise

    return grepo


def _call_command_in_repo(comm, repo, log, fail=False, log_flag=True):
    """Use `subprocess` to call a command in a certain (repo) directory.

    Logs the output (both `stderr` and `stdout`) to the log, and checks the
    return codes to make sure they're valid.  Raises error if not.

    Raises
    ------
    exception `subprocess.CalledProcessError`: if the command fails

    """
    if log_flag:
        log.debug("Running '{}'.".format(" ".join(comm)))
    process = subprocess.Popen(
        comm, cwd=repo, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = process.communicate()
    if stderr is not None:
        err_msg = stderr.decode('ascii').strip().splitlines()
        for em in err_msg:
            log.error(em)
    if stdout is not None:
        out_msg = stdout.decode('ascii').strip().splitlines()
        for om in out_msg:
            log.warning(om)
    # Raises an error if the command failed.
    if fail:
        if process.returncode:
            raise subprocess.CalledProcessError
    return
