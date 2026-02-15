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


# google.generativeai.types.SafetyFeedbackDict

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/google/generative-ai-python/blob/master/google/generativeai/types/safety_types.py#L276-L280">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



Safety feedback for an entire request.

<!-- Placeholder for "Used in" -->

This field is populated if content in the input and/or response
is blocked due to safety settings. SafetyFeedback may not exist
for every HarmCategory. Each SafetyFeedback will return the
safety settings used by the request as well as the lowest
HarmProbability that should be allowed in order to return a
result.



<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Attributes</h2></th></tr>

<tr>
<td>

`rating`<a id="rating"></a>

</td>
<td>

`google.ai.generativelanguage.SafetyRating`

Safety rating evaluated from content.

</td>
</tr><tr>
<td>

`setting`<a id="setting"></a>

</td>
<td>

`google.ai.generativelanguage.SafetySetting`

Safety settings applied to the request.

</td>
</tr>
</table>





<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->