"""Shared pytest fixtures for Astrocats tests."""
from __future__ import annotations

import json
import os
from argparse import Namespace
from collections import OrderedDict

import pytest

from astrocats.catalog.catalog import Catalog
from astrocats.catalog.task import Task
from astrocats.catalog.utils.logger import get_logger


@pytest.fixture
def catalog(tmp_path):
    """Catalog instance with temporary input/output repositories."""
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    test_input = input_dir / "catalog-test-input"
    test_output = output_dir / "catalog-test-output"
    boneyard = output_dir / "catalog-test-boneyard"
    for path in (test_input, test_output, boneyard):
        path.mkdir(parents=True)

    repos = OrderedDict([
        ("output", ["catalog-test-output"]),
        ("boneyard", ["catalog-test-boneyard"]),
        ("external", ["catalog-test-input"]),
        ("internal", []),
        ("private", []),
    ])
    tasks = {
        "test": {
            "nice_name": "Running tests",
            "active": True,
            "update": False,
            "module": "catalog.tasks.test",
            "function": "do_test",
            "repo": "input/catalog-test-input",
            "priority": 0,
        },
        "merge_duplicates": {
            "nice_name": "Merging duplicates",
            "active": False,
            "update": False,
            "module": "catalog.tasks.merge_duplicates",
            "function": "merge_duplicates",
            "groups": ["meta"],
            "priority": -100,
        },
        "set_preferred_names": {
            "nice_name": "Setting preferred names",
            "active": True,
            "update": False,
            "module": "catalog.tasks.set_preferred_names",
            "function": "set_preferred_names",
            "groups": ["meta"],
            "priority": -10,
        },
        "sanitize": {
            "nice_name": "Cleaning up entries",
            "active": True,
            "update": False,
            "module": "catalog.tasks.sanitize",
            "function": "sanitize",
            "groups": ["meta"],
            "priority": -1,
        },
    }
    (input_dir / "repos.json").write_text(json.dumps(repos), encoding="utf-8")
    (input_dir / "tasks.json").write_text(json.dumps(tasks), encoding="utf-8")

    args = Namespace(
        base_path=str(tmp_path),
        write_entries=True,
        delete_old=False,
        update=False,
        archived=False,
        load_stubs=False,
        travis=True,
        clone_depth=0,
        purge_outputs=False,
        args_task_list=None,
        yes_task_list=None,
        no_task_list=None,
        min_task_priority="test",
        max_task_priority=None,
        task_groups=None,
        verbose=True,
        debug=False,
        log_filename=None,
        private=False,
        count=True,
    )
    log = get_logger(name="astrocats-test", stream_level=40, tostr=True)
    cat = Catalog(args, log, git_clone=False)

    sep = os.sep
    cat.PATHS.PATH_BASE = str(tmp_path) + sep
    cat.PATHS.PATH_INPUT = str(input_dir) + sep
    cat.PATHS.PATH_OUTPUT = str(output_dir) + sep
    cat.PATHS.REPOS_LIST = str(input_dir / "repos.json")
    cat.PATHS.TASK_LIST = str(input_dir / "tasks.json")
    cat.PATHS.repos_dict = repos
    cat.repos_dict = repos
    cat.current_task = Task(
        name="test",
        nice_name="Running tests",
        module="catalog.tasks.test",
        function="do_test",
        repo="input/catalog-test-input",
        priority=0,
        archived=False,
    )
    return cat


@pytest.fixture
def catalog_offline(catalog, monkeypatch):
    """Catalog whose URL downloads never hit the network."""

    def fake_download(self, url, timeout, fail=False, post=None, verify=True):
        if "BAD" in url:
            return None
        return "<html>offline-test</html>"

    monkeypatch.setattr(Catalog, "download_url", fake_download)
    return catalog
