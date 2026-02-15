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


# google.generativeai.types.CitationSourceDict

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/google/generative-ai-python/blob/master/google/generativeai/types/citation_types.py#L30-L36">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



A citation to a source for a portion of a specific response.

<!-- Placeholder for "Used in" -->




<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`start_index`<a id="start_index"></a>

</td>
<td>

`int`

Optional. Start of segment of the response
that is attributed to this source.

Index indicates the start of the segment,
measured in bytes.

</td>
</tr><tr>
<td>

`end_index`<a id="end_index"></a>

</td>
<td>

`int`

Optional. End of the attributed segment,
exclusive.

</td>
</tr><tr>
<td>

`uri`<a id="uri"></a>

</td>
<td>

`str`

Optional. URI that is attributed as a source
for a portion of the text.

</td>
</tr><tr>
<td>

`license_`<a id="license_"></a>

</td>
<td>

`str`

Optional. License for the GitHub project that
is attributed as a source for segment.

License info is required for code citations.

</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->