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


# google.generativeai.list_tuned_models

<!-- Insert buttons and diff -->

<table class="tfo-notebook-buttons tfo-api nocontent">
<td>
  <a target="_blank" href="https://github.com/google/generative-ai-python/blob/master/google/generativeai/models.py#L209-L242">
    <img src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" />
    View source on GitHub
  </a>
</td>
</table>



Calls the API to list all tuned models.


<pre class="devsite-click-to-copy prettyprint lang-py tfo-signature-link">
<code>google.generativeai.list_tuned_models(
    *,
    page_size: (int | None) = 50,
    client: (glm.ModelServiceClient | None) = None,
    request_options: (helper_types.RequestOptionsType | None) = None
) -> model_types.TunedModelsIterable
</code></pre>



<!-- Placeholder for "Used in" -->

```
import pprint
for model in genai.list_tuned_models():
    pprint.pprint(model)
```

<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Args</h2></th></tr>

<tr>
<td>

`page_size`<a id="page_size"></a>

</td>
<td>

How many `types.Models` to fetch per page (api call).

</td>
</tr><tr>
<td>

`client`<a id="client"></a>

</td>
<td>

You may pass a `glm.ModelServiceClient` instead of using the default client.

</td>
</tr><tr>
<td>

`request_options`<a id="request_options"></a>

</td>
<td>

Options for the request.

</td>
</tr>
</table>



<!-- Tabular view -->
 <table class="responsive fixed orange">
<colgroup><col width="214px"><col></colgroup>
<tr><th colspan="2"><h2 class="add-link">Yields</h2></th></tr>
<tr class="alt">
<td colspan="2">

<a href="../../google/generativeai/types/TunedModel.md"><code>types.TunedModel</code></a> objects.

</td>
</tr>

</table>



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->