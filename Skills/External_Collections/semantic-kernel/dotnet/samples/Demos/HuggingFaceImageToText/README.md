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

## HuggingFace ImageToText Service Example

This demonstration is simple WindowsForm Sample application that go thru an **images folder provided at the initialization**, searching for all image files. These images are then displayed in the initial window as soon as the application launches.

The application provides an interactive feature where you can click on each image. Upon clicking, the application employs the Semantic Kernel's HuggingFace ImageToText Service to fetch a descriptive analysis of the clicked image.

A critical aspect of the implementation is how the application captures the binary content of the image and sends a request to the Service, awaiting the descriptive text. This process is a key highlight, showcasing the seamless integration and powerful capabilities of our latest software enhancement.

Required packages to use ImageToText HuggingFace Service:

- Microsoft.SemanticKernel
- Microsoft.SemanticKernel.Connectors.HuggingFace

The following code snippet below shows the most important pieces of code on how to use the ImageToText Service (Hugging Face implementation) to retrieve the descriptive text of an image:

```csharp
// Initializes the Kernel
var kernel = Kernel.CreateBuilder()
	.AddHuggingFaceImageToText("Salesforce/blip-image-captioning-base")
    .Build();

// Gets the ImageToText Service
var service = this._kernel.GetRequiredService<IImageToTextService>();
```

Once one of the images is selected, the binary data of the image is retrieved and sent to the ImageToText Service. The service then returns the descriptive text of the image. The following code snippet demonstrates how to use the ImageToText Service to retrieve the descriptive text of an image:

```csharp
// Get the binary content of a JPEG image:
var imageBinary = File.ReadAllBytes("path/to/file.jpg");

// Prepare the image to be sent to the LLM
var imageContent = new ImageContent(imageBinary) { MimeType = "image/jpeg" };

// Retrieves the image description
var textContent = await service.GetTextContentAsync(imageContent);
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->