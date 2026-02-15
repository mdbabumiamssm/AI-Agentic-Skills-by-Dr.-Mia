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

from unittest.mock import AsyncMock, Mock

from semantic_kernel.connectors.ai.chat_completion_client_base import ChatCompletionClientBase
from semantic_kernel.connectors.ai.prompt_execution_settings import PromptExecutionSettings
from semantic_kernel.contents.chat_message_content import ChatMessageContent
from semantic_kernel.core_plugins.conversation_summary_plugin import ConversationSummaryPlugin
from semantic_kernel.functions.kernel_arguments import KernelArguments
from semantic_kernel.kernel import Kernel
from semantic_kernel.prompt_template.prompt_template_config import PromptTemplateConfig


def test_conversation_summary_plugin():
    config = PromptTemplateConfig(name="test", description="test")
    plugin = ConversationSummaryPlugin(config)
    assert plugin._summarizeConversationFunction is not None
    assert plugin.return_key == "summary"


def test_conversation_summary_plugin_with_deprecated_value(kernel):
    config = PromptTemplateConfig(name="test", description="test")
    plugin = ConversationSummaryPlugin(config, kernel=kernel)
    assert plugin._summarizeConversationFunction is not None
    assert plugin.return_key == "summary"


async def test_summarize_conversation(kernel: Kernel):
    service = AsyncMock(spec=ChatCompletionClientBase)
    service.service_id = "default"
    service.get_chat_message_contents = AsyncMock(
        return_value=[ChatMessageContent(role="assistant", content="Hello World!")]
    )
    service.get_prompt_execution_settings_class = Mock(return_value=PromptExecutionSettings)
    kernel.add_service(service)
    config = PromptTemplateConfig(
        name="test", description="test", execution_settings={"default": PromptExecutionSettings()}
    )
    kernel.add_plugin(ConversationSummaryPlugin(config), "summarizer")
    args = KernelArguments(input="Hello World!")

    await kernel.invoke(plugin_name="summarizer", function_name="SummarizeConversation", arguments=args)
    args["summary"] == "Hello world"
    service.get_chat_message_contents.assert_called_once()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
