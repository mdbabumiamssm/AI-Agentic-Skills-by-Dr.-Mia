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


# google.generativeai.protos.SemanticRetrieverConfig

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/googleapis/google-cloud-python/tree/main/packages/google-ai-generativelanguage/google/ai/generativelanguage_v1beta/types/generative_service.py#L428-L480">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



Configuration for retrieving grounding content from a ``Corpus`` or ``Document`` created using the Semantic Retriever API.

<!-- Placeholder for "Used in" -->





<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`source`<a id="source"></a>

</td>
<td>

`str`

Required. Name of the resource for retrieval. Example:
``corpora/123`` or ``corpora/123/documents/abc``.

</td>
</tr><tr>
<td>

`query`<a id="query"></a>

</td>
<td>

`google.ai.generativelanguage.Content`

Required. Query to use for matching ``Chunk``\ s in the
given resource by similarity.

</td>
</tr><tr>
<td>

`metadata_filters`<a id="metadata_filters"></a>

</td>
<td>

`MutableSequence[google.ai.generativelanguage.MetadataFilter]`

Optional. Filters for selecting ``Document``\ s and/or
``Chunk``\ s from the resource.

</td>
</tr><tr>
<td>

`max_chunks_count`<a id="max_chunks_count"></a>

</td>
<td>

`int`

Optional. Maximum number of relevant ``Chunk``\ s to
retrieve.


</td>
</tr><tr>
<td>

`minimum_relevance_score`<a id="minimum_relevance_score"></a>

</td>
<td>

`float`

Optional. Minimum relevance score for retrieved relevant
``Chunk``\ s.


</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->