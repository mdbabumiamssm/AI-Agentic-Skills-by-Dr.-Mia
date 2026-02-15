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


# google.generativeai.protos.Segment

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/googleapis/google-cloud-python/tree/main/packages/google-ai-generativelanguage/google/ai/generativelanguage_v1beta/types/generative_service.py#L1074-L1109">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



Segment of the content.

<!-- Placeholder for "Used in" -->




<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`part_index`<a id="part_index"></a>

</td>
<td>

`int`

Output only. The index of a Part object
within its parent Content object.

</td>
</tr><tr>
<td>

`start_index`<a id="start_index"></a>

</td>
<td>

`int`

Output only. Start index in the given Part,
measured in bytes. Offset from the start of the
Part, inclusive, starting at zero.

</td>
</tr><tr>
<td>

`end_index`<a id="end_index"></a>

</td>
<td>

`int`

Output only. End index in the given Part,
measured in bytes. Offset from the start of the
Part, exclusive, starting at zero.

</td>
</tr><tr>
<td>

`text`<a id="text"></a>

</td>
<td>

`str`

Output only. The text corresponding to the
segment from the response.

</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->