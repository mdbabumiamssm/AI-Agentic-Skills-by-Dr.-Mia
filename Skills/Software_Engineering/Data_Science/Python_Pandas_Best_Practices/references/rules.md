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

# Pandas Best Practices Rules

## 1. Vectorization over Loops

*   **NEVER** iterate over rows using `for index, row in df.iterrows():`. It is extremely slow.
*   **Use Vectorized Operations**:
    ```python
    # BAD
    df['c'] = df.apply(lambda row: row['a'] + row['b'], axis=1)

    # GOOD
    df['c'] = df['a'] + df['b']
    ```

## 2. Method Chaining

*   **Fluent Interface**: Chain methods to create a clear recipe of transformations.
*   **`assign`**: Use `.assign()` for creating new columns within a chain.
    ```python
    df = (
        raw_data
        .query("status == 'active'")
        .assign(total_cost = lambda x: x.price * x.quantity)
        .rename(columns={"total_cost": "Revenue"})
    )
    ```

## 3. Explicit Copies

*   **SettingWithCopyWarning**: Avoid chained indexing for assignment `df[mask]['col'] = val`.
*   **Use `.copy()`**: When creating a subset that you intend to modify, use `.copy()` explicitly.
    ```python
    active_users = df[df['active']].copy()
    active_users.loc[:, 'score'] = 100
    ```

## 4. Dtypes and Memory

*   **Categoricals**: Use `category` dtype for string columns with low cardinality (e.g., 'gender', 'country'). It saves massive amounts of RAM.
*   **Nullable Integers**: Use `Int64` (capital I) to allow `NaN` in integer columns without forcing them to float.

## 5. Input/Output (I/O)

*   **Parquet over CSV**: Prefer `.to_parquet()` and `.read_parquet()` for intermediate storage. It preserves dtypes and is much faster/smaller than CSV.
*   **Date Parsing**: Parse dates on load: `pd.read_csv(..., parse_dates=['date_col'])`.

## 6. Aggregation

*   **Named Aggregation**: Use the explicit syntax for clarity.
    ```python
    df.groupby('group').agg(
        mean_age=('age', 'mean'),
        max_score=('score', 'max')
    )
    ```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->