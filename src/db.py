"""SQLite database access layer.

Single-file canonical store, committed to the repo. All schema changes
go through schema.sql. No migrations framework — we just re-run the
schema (it's all CREATE IF NOT EXISTS).
"""

from __future__ import annotations

import logging
import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

from . import config

log = logging.getLogger(__name__)


def init_db(db_path: Path | None = None, schema_path: Path | None = None) -> None:
    """Create the database file (if needed) and apply schema.sql."""
    db_path = db_path or config.DB_PATH
    schema_path = schema_path or config.SCHEMA_PATH
    db_path.parent.mkdir(parents=True, exist_ok=True)
    schema_sql = schema_path.read_text()
    with sqlite3.connect(db_path) as conn:
        conn.executescript(schema_sql)
        conn.commit()
    log.info("Initialized database at %s", db_path)


@contextmanager
def connect(db_path: Path | None = None) -> Iterator[sqlite3.Connection]:
    """Context manager for a database connection with row-as-dict factory."""
    db_path = db_path or config.DB_PATH
    if not db_path.exists():
        init_db(db_path)
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    # Enable foreign keys (we don't have any yet but cheap insurance for later)
    conn.execute("PRAGMA foreign_keys = ON")
    # Write-ahead logging for better concurrency / crash safety
    conn.execute("PRAGMA journal_mode = WAL")
    try:
        yield conn
    finally:
        conn.close()
