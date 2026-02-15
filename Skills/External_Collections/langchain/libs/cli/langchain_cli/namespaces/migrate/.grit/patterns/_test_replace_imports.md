<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Testing the replace_imports migration

This runs the v0.2 migration with a desired set of rules.

```grit
language python

langchain_all_migrations()
```

## Single import

Before:

```python
from langchain.chat_models import ChatOpenAI
```

After:

```python
from langchain_community.chat_models import ChatOpenAI
```

## Community to partner

```python
from langchain_community.chat_models import ChatOpenAI
```

```python
from langchain_openai import ChatOpenAI
```

## Noop

This file should not match at all.

```python
from foo import ChatOpenAI
```

## Mixed imports

```python
from langchain_community.chat_models import ChatOpenAI, ChatAnthropic, foo
```

```python
from langchain_community.chat_models import foo

from langchain_openai import ChatOpenAI

from langchain_anthropic import ChatAnthropic

```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->