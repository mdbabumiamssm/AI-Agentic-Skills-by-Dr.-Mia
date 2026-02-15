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

# Demos

Contains a helper that adds an override `System.Environment` class to a project.
This override version has an enhanced `GetEnvironmentVariable` method that prompts the user
to enter a value if the environment variable is not set.

The code is still fully copyable to another project. These sample projects just allow for a simplified user experience
for users who are new and just getting started.

To use this in your project, add the following to your `.csproj` file:

```xml
<ItemGroup>
  <Using Include="SampleHelpers.SampleEnvironment" Alias="Environment" />
</ItemGroup>

<ItemGroup>
  <Compile Include="$(MSBuildThisFileDirectory)\..\src\Shared\Demos\*.cs" LinkBase="" Visible="false" />
</ItemGroup>
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->