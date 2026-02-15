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


# google.generativeai.protos.GenerateMessageResponse

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/googleapis/google-cloud-python/tree/main/packages/google-ai-generativelanguage/google/ai/generativelanguage_v1beta/types/discuss_service.py#L124-L160">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



The response from the model.

<!-- Placeholder for "Used in" -->

This includes candidate messages and
conversation history in the form of chronologically-ordered
messages.



<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`candidates`<a id="candidates"></a>

</td>
<td>

`MutableSequence[google.ai.generativelanguage.Message]`

Candidate response messages from the model.

</td>
</tr><tr>
<td>

`messages`<a id="messages"></a>

</td>
<td>

`MutableSequence[google.ai.generativelanguage.Message]`

The conversation history used by the model.

</td>
</tr><tr>
<td>

`filters`<a id="filters"></a>

</td>
<td>

`MutableSequence[google.ai.generativelanguage.ContentFilter]`

A set of content filtering metadata for the prompt and
response text.

This indicates which ``SafetyCategory``\ (s) blocked a
candidate from this response, the lowest ``HarmProbability``
that triggered a block, and the HarmThreshold setting for
that category.

</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->