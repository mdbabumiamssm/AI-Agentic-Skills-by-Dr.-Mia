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


# google.generativeai.protos.ListCorporaResponse

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/googleapis/google-cloud-python/tree/main/packages/google-ai-generativelanguage/google/ai/generativelanguage_v1beta/types/retriever_service.py#L170-L196">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



Response from ``ListCorpora`` containing a paginated list of ``Corpora``.

<!-- Placeholder for "Used in" -->
 The results are sorted by ascending
``corpus.create_time``.



<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`corpora`<a id="corpora"></a>

</td>
<td>

`MutableSequence[google.ai.generativelanguage.Corpus]`

The returned corpora.

</td>
</tr><tr>
<td>

`next_page_token`<a id="next_page_token"></a>

</td>
<td>

`str`

A token, which can be sent as ``page_token`` to retrieve the
next page. If this field is omitted, there are no more
pages.

</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->