# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest

from langchain_core.caches import RETURN_VAL_TYPE, InMemoryCache
from langchain_core.outputs import Generation


@pytest.fixture
def cache() -> InMemoryCache:
    """Fixture to provide an instance of InMemoryCache."""
    return InMemoryCache()


def cache_item(item_id: int) -> tuple[str, str, RETURN_VAL_TYPE]:
    """Generate a valid cache item."""
    prompt = f"prompt{item_id}"
    llm_string = f"llm_string{item_id}"
    generations = [Generation(text=f"text{item_id}")]
    return prompt, llm_string, generations


def test_initialization() -> None:
    """Test the initialization of InMemoryCache."""
    cache = InMemoryCache()
    assert cache._cache == {}
    assert cache._maxsize is None

    cache_with_maxsize = InMemoryCache(maxsize=2)
    assert cache_with_maxsize._cache == {}
    assert cache_with_maxsize._maxsize == 2

    with pytest.raises(ValueError, match="maxsize must be greater than 0"):
        InMemoryCache(maxsize=0)


def test_lookup(
    cache: InMemoryCache,
) -> None:
    """Test the lookup method of InMemoryCache."""
    prompt, llm_string, generations = cache_item(1)
    cache.update(prompt, llm_string, generations)
    assert cache.lookup(prompt, llm_string) == generations
    assert cache.lookup("prompt2", "llm_string2") is None


def test_update_with_no_maxsize(cache: InMemoryCache) -> None:
    """Test the update method of InMemoryCache with no maximum size."""
    prompt, llm_string, generations = cache_item(1)
    cache.update(prompt, llm_string, generations)
    assert cache.lookup(prompt, llm_string) == generations


def test_update_with_maxsize() -> None:
    """Test the update method of InMemoryCache with a maximum size."""
    cache = InMemoryCache(maxsize=2)

    prompt1, llm_string1, generations1 = cache_item(1)
    cache.update(prompt1, llm_string1, generations1)
    assert cache.lookup(prompt1, llm_string1) == generations1

    prompt2, llm_string2, generations2 = cache_item(2)
    cache.update(prompt2, llm_string2, generations2)
    assert cache.lookup(prompt2, llm_string2) == generations2

    prompt3, llm_string3, generations3 = cache_item(3)
    cache.update(prompt3, llm_string3, generations3)

    assert cache.lookup(prompt1, llm_string1) is None  # 'prompt1' should be evicted
    assert cache.lookup(prompt2, llm_string2) == generations2
    assert cache.lookup(prompt3, llm_string3) == generations3


def test_clear(cache: InMemoryCache) -> None:
    """Test the clear method of InMemoryCache."""
    prompt, llm_string, generations = cache_item(1)
    cache.update(prompt, llm_string, generations)
    cache.clear()
    assert cache.lookup(prompt, llm_string) is None


async def test_alookup(cache: InMemoryCache) -> None:
    """Test the asynchronous lookup method of InMemoryCache."""
    prompt, llm_string, generations = cache_item(1)
    await cache.aupdate(prompt, llm_string, generations)
    assert await cache.alookup(prompt, llm_string) == generations
    assert await cache.alookup("prompt2", "llm_string2") is None


async def test_aupdate_with_no_maxsize(cache: InMemoryCache) -> None:
    """Test the asynchronous update method of InMemoryCache with no maximum size."""
    prompt, llm_string, generations = cache_item(1)
    await cache.aupdate(prompt, llm_string, generations)
    assert await cache.alookup(prompt, llm_string) == generations


async def test_aupdate_with_maxsize() -> None:
    """Test the asynchronous update method of InMemoryCache with a maximum size."""
    cache = InMemoryCache(maxsize=2)
    prompt, llm_string, generations = cache_item(1)
    await cache.aupdate(prompt, llm_string, generations)
    assert await cache.alookup(prompt, llm_string) == generations

    prompt2, llm_string2, generations2 = cache_item(2)
    await cache.aupdate(prompt2, llm_string2, generations2)
    assert await cache.alookup(prompt2, llm_string2) == generations2

    prompt3, llm_string3, generations3 = cache_item(3)
    await cache.aupdate(prompt3, llm_string3, generations3)

    assert await cache.alookup(prompt, llm_string) is None
    assert await cache.alookup(prompt2, llm_string2) == generations2
    assert await cache.alookup(prompt3, llm_string3) == generations3


async def test_aclear(cache: InMemoryCache) -> None:
    """Test the asynchronous clear method of InMemoryCache."""
    prompt, llm_string, generations = cache_item(1)
    await cache.aupdate(prompt, llm_string, generations)
    await cache.aclear()
    assert await cache.alookup(prompt, llm_string) is None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
