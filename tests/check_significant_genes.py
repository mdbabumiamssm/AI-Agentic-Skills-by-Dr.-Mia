import pandas as pd

df = pd.read_csv('tests/results_BM_vs_BT/BM_vs_BT_DEGs.csv')
sig = df[df['pval_adj'] < 0.05]
if sig.empty:
    print("No significant DEGs found (adj p < 0.05).")
else:
    print(f"Found {len(sig)} significant DEGs.")
    print("Top 20 by abs(log2fc):")
    # Add abs log2fc column for sorting
    sig['abs_log2fc'] = sig['log2fc'].abs()
    print(sig.sort_values('abs_log2fc', ascending=False).head(20)[['gene', 'log2fc', 'pval_adj', 'comparison_group']])
