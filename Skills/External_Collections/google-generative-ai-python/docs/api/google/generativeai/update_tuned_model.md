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


# google.generativeai.update_tuned_model

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/google/generative-ai-python/blob/master/google/generativeai/models.py#L393-L443">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



Calls the API to push updates to a specified tuned model where only certain attributes are updatable.


<pre class="devsite-click-to-copy prettyprint lang-py tfo-signature-link">
<code>google.generativeai.update_tuned_model(
    tuned_model: (str | protos.TunedModel),
    updates: (dict[str, Any] | None) = None,
    *,
    client: (glm.ModelServiceClient | None) = None,
    request_options: (helper_types.RequestOptionsType | None) = None
) -> model_types.TunedModel
</code></pre>



<!-- Placeholder for "Used in" -->


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->