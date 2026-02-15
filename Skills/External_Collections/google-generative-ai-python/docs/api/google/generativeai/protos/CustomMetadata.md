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


# google.generativeai.protos.CustomMetadata

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/googleapis/google-cloud-python/tree/main/packages/google-ai-generativelanguage/google/ai/generativelanguage_v1beta/types/retriever.py#L155-L202">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



User provided metadata stored as key-value pairs.

<!-- Placeholder for "Used in" -->

This message has `oneof`_ fields (mutually exclusive fields).
For each oneof, at most one member field can be set at the same time.
Setting any member of the oneof automatically clears all other
members.




<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`string_value`<a id="string_value"></a>

</td>
<td>

`str`

The string value of the metadata to store.

This field is a member of `oneof`_ ``value``.

</td>
</tr><tr>
<td>

`string_list_value`<a id="string_list_value"></a>

</td>
<td>

`google.ai.generativelanguage.StringList`

The StringList value of the metadata to
store.

This field is a member of `oneof`_ ``value``.

</td>
</tr><tr>
<td>

`numeric_value`<a id="numeric_value"></a>

</td>
<td>

`float`

The numeric value of the metadata to store.

This field is a member of `oneof`_ ``value``.

</td>
</tr><tr>
<td>

`key`<a id="key"></a>

</td>
<td>

`str`

Required. The key of the metadata to store.

</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->