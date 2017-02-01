"""Utilities for working with git.
"""

import git

class MyProgressPrinter(git.RemoteProgress):
    def update(self, op_code, cur_count, max_count=None, message=''):
        # print(op_code, cur_count, max_count, message, end="\r")
        msg = "{} - {}".format(op_code, cur_count)
        if max_count is not None:
            msg += "/{}".format(max_count)
        if len(msg):
            msg += " - {}".format(message)
        print(msg)


def fetch(repo_name, progress=True, log=None):
    repo = git.Repo(repo_name)
    _progress = MyProgressPrinter() if progress else None
    for fetch_info in repo.remote().fetch(progress=_progress):
        if log is not None:
            log.debug("Updated {}:{} to {}".format(repo_name, fetch_info.ref, fetch_info.commit))


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
