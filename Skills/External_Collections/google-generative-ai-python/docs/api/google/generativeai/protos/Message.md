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


# google.generativeai.protos.Message

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/googleapis/google-cloud-python/tree/main/packages/google-ai-generativelanguage/google/ai/generativelanguage_v1beta/types/discuss_service.py#L163-L211">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



The base unit of structured text.

<!-- Placeholder for "Used in" -->

A ``Message`` includes an ``author`` and the ``content`` of the
``Message``.

The ``author`` is used to tag messages when they are fed to the
model as text.





<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`author`<a id="author"></a>

</td>
<td>

`str`

Optional. The author of this Message.

This serves as a key for tagging
the content of this Message when it is fed to
the model as text.

The author can be any alphanumeric string.

</td>
</tr><tr>
<td>

`content`<a id="content"></a>

</td>
<td>

`str`

Required. The text content of the structured ``Message``.

</td>
</tr><tr>
<td>

`citation_metadata`<a id="citation_metadata"></a>

</td>
<td>

`google.ai.generativelanguage.CitationMetadata`

Output only. Citation information for model-generated
``content`` in this ``Message``.

If this ``Message`` was generated as output from the model,
this field may be populated with attribution information for
any text included in the ``content``. This field is used
only on output.


</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->