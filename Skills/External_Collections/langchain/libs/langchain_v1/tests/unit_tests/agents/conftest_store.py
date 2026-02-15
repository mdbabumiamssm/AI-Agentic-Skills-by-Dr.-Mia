# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from contextlib import asynccontextmanager, contextmanager

from langgraph.store.memory import InMemoryStore


@contextmanager
def _store_memory():
    store = InMemoryStore()
    yield store


@asynccontextmanager
async def _store_memory_aio():
    store = InMemoryStore()
    yield store


# Placeholder functions for other store types that aren't available
@contextmanager
def _store_postgres():
    # Fallback to memory for now
    store = InMemoryStore()
    yield store


@contextmanager
def _store_postgres_pipe():
    # Fallback to memory for now
    store = InMemoryStore()
    yield store


@contextmanager
def _store_postgres_pool():
    # Fallback to memory for now
    store = InMemoryStore()
    yield store


@asynccontextmanager
async def _store_postgres_aio():
    # Fallback to memory for now
    store = InMemoryStore()
    yield store


@asynccontextmanager
async def _store_postgres_aio_pipe():
    # Fallback to memory for now
    store = InMemoryStore()
    yield store


@asynccontextmanager
async def _store_postgres_aio_pool():
    # Fallback to memory for now
    store = InMemoryStore()
    yield store

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
