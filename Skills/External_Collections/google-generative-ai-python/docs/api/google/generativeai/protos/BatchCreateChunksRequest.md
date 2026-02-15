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


# google.generativeai.protos.BatchCreateChunksRequest

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/googleapis/google-cloud-python/tree/main/packages/google-ai-generativelanguage/google/ai/generativelanguage_v1beta/types/retriever_service.py#L561-L584">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



Request to batch create ``Chunk``\ s.

<!-- Placeholder for "Used in" -->




<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`parent`<a id="parent"></a>

</td>
<td>

`str`

Optional. The name of the ``Document`` where this batch of
``Chunk``\ s will be created. The parent field in every
``CreateChunkRequest`` must match this value. Example:
``corpora/my-corpus-123/documents/the-doc-abc``

</td>
</tr><tr>
<td>

`requests`<a id="requests"></a>

</td>
<td>

`MutableSequence[google.ai.generativelanguage.CreateChunkRequest]`

Required. The request messages specifying the ``Chunk``\ s
to create. A maximum of 100 ``Chunk``\ s can be created in a
batch.

</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->