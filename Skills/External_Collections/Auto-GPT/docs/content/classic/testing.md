<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Running tests

To run all tests, use the following command:

```shell
pytest
```

If `pytest` is not found:

```shell
python -m pytest
```

### Running specific test suites

- To run without integration tests:

```shell
pytest --without-integration
```

- To run without *slow* integration tests:

```shell
pytest --without-slow-integration
```

- To run tests and see coverage:

```shell
pytest --cov=autogpt --without-integration --without-slow-integration
```

## Running the linter

This project uses [flake8](https://flake8.pycqa.org/en/latest/) for linting.
We currently use the following rules: `E303,W293,W291,W292,E305,E231,E302`.
See the [flake8 rules](https://www.flake8rules.com/) for more information.

To run the linter:

```shell
flake8 .
```

Or:

```shell
python -m flake8 .
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->