# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

from azure.ai.inference.tracing import AIInferenceInstrumentor
from azure.core.settings import settings

from semantic_kernel.kernel_pydantic import KernelBaseModel
from semantic_kernel.utils.telemetry.model_diagnostics.model_diagnostics_settings import ModelDiagnosticSettings


class AzureAIInferenceTracing(KernelBaseModel):
    """Enable tracing for Azure AI Inference.

    This class is intended to be used as a context manager.
    The instrument() call effect should be scoped to the context manager.
    """

    diagnostics_settings: ModelDiagnosticSettings

    def __init__(self, diagnostics_settings: ModelDiagnosticSettings | None = None) -> None:
        """Initialize the Azure AI Inference Tracing.

        Args:
            diagnostics_settings (ModelDiagnosticSettings, optional): Model diagnostics settings. Defaults to None.
        """
        settings.tracing_implementation = "opentelemetry"
        super().__init__(diagnostics_settings=diagnostics_settings or ModelDiagnosticSettings())

    def __enter__(self) -> None:
        """Enable tracing.

        Both enable_otel_diagnostics and enable_otel_diagnostics_sensitive will enable tracing.
        enable_otel_diagnostics_sensitive will also enable content recording.
        """
        if (
            self.diagnostics_settings.enable_otel_diagnostics
            or self.diagnostics_settings.enable_otel_diagnostics_sensitive
        ):
            AIInferenceInstrumentor().instrument(  # type: ignore
                enable_content_recording=self.diagnostics_settings.enable_otel_diagnostics_sensitive
            )

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Disable tracing."""
        if (
            self.diagnostics_settings.enable_otel_diagnostics
            or self.diagnostics_settings.enable_otel_diagnostics_sensitive
        ):
            AIInferenceInstrumentor().uninstrument()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
